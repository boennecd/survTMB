.min_Y <- sqrt(.Machine$double.eps)

#' @importFrom stats update
.get_Y_names <- function(formula){
  formula <- update(formula, . ~ 1)
  all.vars(formula)
}

.trunc_counting <- function(Y){
  stopifnot(inherits(Y, "Surv") & attr(Y, "type") == "counting")
  tstart <- pmax(Y[, 1], .min_Y)
  keep <- tstart < Y[, 2]

  list(Y = Y[keep, ], keep = keep)
}

.get_joint_knots <- function(n_knots, times)
  stop(".get_joint_knots not implemented")

#' @importFrom stats model.response model.frame terms model.matrix logLik
#' @importFrom reshape2 melt
#' @importFrom splines ns
#' @importFrom lme4 lmer lmerControl fixef VarCorr .makeCC
get_marker_start_params <- function(
  formula, data, mknots, gknots, time_var, id_var, need_start_vals = TRUE){
  mf <- model.frame(formula, data)
  X <- model.matrix(terms(mf), mf)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  out <- list(Y = t(model.response(mf)), X = t(X),
              id = eval(id_var, data),
              time = eval(time_var, data))

  if(is.integer(mknots) && length(mknots) == 1L)
    mknots <- .get_joint_knots(mknots, out$time)
  out$mknots <- mknots
  if(is.integer(gknots) && length(gknots) == 1L)
    gknots <- .get_joint_knots(gknots, out$time)
  out$gknots <- gknots

  if(!need_start_vals){
    out$Psi <- out$Sigma <- out$B <- out$ll <- out$gamma <- numeric()
    return(out)
  }

  # get data.frame to work with
  Y_names <- .get_Y_names(formula)
  X_names <- all.vars(update(formula, 1 ~ .))
  data <- melt(
    data, id.vars = c(deparse(id_var), deparse(time_var), X_names),
    measure.vars = Y_names,
    variable.name = "XXTHEVARIABLEXX", value.name = "XXTHEVALUEXX")

  # fit mixed model
  m_bk <- mknots[ c(1L, length(mknots))]
  m_ik <- mknots[-c(1L, length(mknots))]

  g_bk <- gknots[ c(1L, length(gknots))]
  g_ik <- gknots[-c(1L, length(gknots))]

  frm <- substitute(
    XXTHEVALUEXX ~ -1 +
      XXTHEVARIABLEXX : ns(ti, knots = g_ik, Boundary.knots = g_bk,
                           intercept = TRUE) +
      (XXTHEVARIABLEXX : ns(ti, knots = m_ik, Boundary.knots = m_bk,
                            intercept = TRUE) - 1L | i),
    list(ti = time_var, i = id_var, g_ik = g_ik, g_bk = g_bk,
         m_ik = m_ik, m_bk = m_bk))
  frm <- eval(frm)

  if(length(X_names) > 0)
    for(x_nam in rev(X_names)){
      frm_call <- substitute(
        update(frm, . ~ XXTHEVARIABLEXX : x_var + .),
        list(x_var = as.name(x_nam)))
      frm <- eval(frm_call)
    }

  fit <- lmer(frm, data, control = lmerControl(
    check.conv.grad = .makeCC("ignore", tol = 1e-3, relTol = NULL)))

  n_y <- NROW(out$Y)
  d_x <- length(X_names)

  gamma <- t(matrix(fixef(fit)[seq_len(d_x * n_y)], nrow = n_y))

  B <- matrix(fixef(fit)[-seq_along(gamma)], ncol = n_y)

  vc <- VarCorr(fit)
  Psi <- vc$id
  attr(Psi, "correlation") <- attr(Psi, "stddev") <- NULL
  dimnames(Psi) <- NULL
  d_m <- NROW(Psi) / n_y
  K <- get_commutation(n_y, d_m)
  Psi <- tcrossprod(K %*% Psi, K)

  Sigma <- diag(attr(vc, "sc")^2, n_y)

  out[c("gamma", "B", "Psi", "Sigma", "ll")] <-
    list(gamma, B, Psi, Sigma, ll = c(logLik(fit)))
  out
}

#' @importFrom stats update model.response model.frame model.matrix terms na.omit optim
#' @importFrom survival tmerge coxph Surv
#' @importFrom utils head tail
get_surv_start_params <- function(
  formula, data, mformula, mdata, id_var, time_var, knots, n_nodes,
  need_start_vals = TRUE){
  # get the outcome
  out <- local({
    mf <- model.frame(formula, data)
    outcome <- model.response(mf)
    stopifnot(inherits(outcome, "Surv") &
                attr(outcome, "type") == "counting")
    Z <- model.matrix(terms(mf), mf)
    Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
    id <- eval(id_var, data)

    outcome <- .trunc_counting(outcome)
    keep <- outcome$keep
    list(Y = t(outcome$Y[keep, ]), Z = t(Z[keep, , drop = FALSE]),
         id = id[keep], keep = keep)
  })

  if(is.integer(knots) && length(knots) == 1L)
    knots <- .get_joint_knots(knots, out$Y[out$Y[, 3] == 1, 2])
  out$knots <- knots

  if(!need_start_vals){
    out$delta <- out$alpha <- out$omega <- out$ll <- numeric()
    return(out)
  }

  # make joint data set
  Y_names <- .get_Y_names(mformula)
  n_y <- length(Y_names)
  dY <- sapply(Y_names, as.name)

  tdat <- local({
    dtstart <- formula[[2L]][[2L]]
    dtstop  <- formula[[2L]][[3L]]
    devent  <- formula[[2L]][[4L]]
    tcall <- substitute(tmerge(
      data, data, id = id, tstart = tstart, tstop = dtstop,
      ev = event(dtstop, devent)),
      list(tstart = dtstart, dtstop = dtstop, devent = devent,
           id = id_var))
    out <- eval(tcall, environment())

    add_Y_call <- substitute(tmerge(out, mdata, id = id),
                             list(id = id_var))
    add_Y_call <- as.list(add_Y_call)
    for(i in 1:n_y){
      add_Y_call <- c(add_Y_call, substitute(
        tdc(dtstop, Y), list(dtstop = time_var, Y = dY[[i]])))
      names(add_Y_call)[length(add_Y_call)] <- Y_names[i]
    }
    add_Y_call <- as.call(add_Y_call)

    out <- eval(add_Y_call)
    na.omit(out)
  })

  # fit coxph model
  cformula <- update(formula, Surv(tstart, tstop, ev) ~ . - 1)
  for(dYi in dY)
    cformula <- eval(substitute(update(cformula, . ~ . + XX),
                                list(XX = dYi)))

  cfit <- coxph(cformula, tdat)
  delta <- head(cfit$coefficients, -n_y)
  alpha <- tail(cfit$coefficients,  n_y)
  stopifnot(!anyNA(delta), !anyNA(alpha))

  # fit parameteric baseline
  mf <- model.frame(cformula, tdat)
  Sy <- .trunc_counting(model.response(mf))
  SZ <- t(model.matrix(terms(mf), mf)[Sy$keep, , drop = FALSE])
  bk <- knots[ c(1L, length(knots))]
  ik <- knots[-c(1L, length(knots))]
  omega <- rep(0, length(knots))

  tstart <- Sy$Y[, 1]
  tstop  <- Sy$Y[, 2]
  Y      <- Sy$Y[, 3]

  func <- function(par, ..., grad){
    o <- par[ seq_along(omega)]
    d <- par[-seq_along(omega)]
    -drop(joint_start_ll(
      Y = Y, tstart = tstart, tstop = tstop, omega = o,
      delta = d, Z = SZ, n_nodes = n_nodes, bound_knots = bk,
      inter_knots = ik, grad = grad))
  }
  fn <- func
  formals(fn)$grad <- FALSE
  gr <- func
  formals(gr)$grad <- TRUE

  opt_out <- optim(c(omega, delta, alpha), fn, gr, method = "BFGS")
  omega <- opt_out$par[seq_along(omega)]
  delta <- opt_out$par[seq_along(delta) + length(omega)]
  alpha <- opt_out$par[seq_along(alpha) + length(omega) + length(delta)]

  stopifnot(opt_out$convergence == 0L,
            all(is.finite(omega)), all(is.finite(delta)),
            all(is.finite(alpha)))
  out[c("delta", "alpha", "omega", "ll")] <- list(
    delta, alpha, omega, -opt_out$value)
  out
}

.get_joint_va_start <- function(skew_start, Psi, n_groups) {
  va_pars <- local({
    K <- NCOL(Psi)
    if(length(skew_start) == 1L)
      skew_start <- rep(skew_start, K)
    stopifnot(length(skew_start) == K)

    out <- .cp_to_dp(mu = numeric(K), Sigma = Psi, gamma = skew_start)

    xi <- out$xi
    Omega <- out$Psi
    rho <- out$rho

    c(xi = xi, .cov_to_theta(Omega), rho = rho)
  })

  va_pars <- structure(
    rep(va_pars, n_groups),
    names = outer(names(va_pars), seq_len(n_groups),
                  function(x, y) paste0("g", y, ":", x)))
}

.opt_sub <- function(par, which_par, extract_name, replacement, fn, gr,
                     opt_func, control){
  is_sub_par <- which(grepl(which_par, names(par), perl = TRUE))
  nam_par <- gsub(extract_name, replacement, names(par)[is_sub_par],
                  perl = TRUE)

  nams <- unique(nam_par)
  va_map <- vapply(nams, function(z) is_sub_par[which(z == nam_par)],
                   integer(length(nam_par) / length(nams)))
  if(!is.matrix(va_map))
    va_map <- as.matrix(va_map)

  get_par <- function(x){
    out <- par
    for(i in 1:NCOL(va_map))
      out[va_map[, i]] <- x[i]
    out
  }
  fn_sub <- function(x, ...){
    fn(get_par(x))
  }
  gr_sub <- function(x, ...){
    grad <- gr(get_par(x))
    apply(va_map, 2L, function(indices) sum(grad[indices]))
  }
  par_sub <- apply(
    va_map, 2L, function(indices) mean(par[indices]))

  opt_out <- opt_func(
    par_sub, fn_sub, gr_sub, control = control)
  out <- par
  for(i in seq_len(NCOL(va_map)))
    out[va_map[, i]] <- opt_out$par[i]
  out
}

#' Construct Objective Functions with Derivatives for a Joint Survival and
#' Marker Model
#'
#' @param sformula \code{\link{formula}} for the time to event. The
#'                 left-hand side must be \code{\link{Surv}} object with
#'                 type "counting".
#' @param mformula \code{\link{formula}} for the marker process. The
#'                 left-hand side must \code{\link{cbind}}ed if there are
#'                 more than one marker.
#' @param sdata \code{\link{data.frame}} with the time to event data.
#' @param mdata \code{\link{data.frame}} with the marker data.
#' @param id_var name of the individual identifier. Must be present in both
#'               \code{sdata} and \code{mdata}.
#' @param time_var name of the time variable in \code{mdata}.
#' @param mknots numeric vector with knots for the processes' random mean
#'               term or number of knots to use.
#' @param gknots numeric vector with knots for the processes' fixed mean
#'               term or number of knots to use.
#' @param sknots numeric vector with knots for the baseline hazard
#'               or number of knots to use.
#' @param n_nodes number of nodes to use with Gauss-Legendre quadrature for
#'                the cumulative hazard and nodes to use with Gauss-Hermite
#'                quadrature for the entropy term.
#' @param skew_start starting value for the Pearson's moment coefficient of
#'                   skewness parameter when a SNVA is used. Currently, a
#'                   somewhat arbitrary value.
#' @param opt_func general optimization function to use. It
#'                 needs to have an interface like \code{\link{optim}}.
#' @param n_threads integer with number of threads to use.
#' @param sparse_hess logical for whether to make sparse Hessian computation
#'                    available. Memory and computation time is saved if it is
#'                    \code{FALSE}.
#' @param B staring value for B.
#' @param Psi staring value for Psi.
#' @param Sigma staring value for Sigma.
#' @param omega staring value for omega.
#' @param alpha staring value for alpha.
#' @param delta staring value for delta.
#' @param gamma staring value for gamma.
#' @param va_par staring value for variational parameters.
#'
#' @export
make_joint_ADFun <- function(
  sformula, mformula, sdata, mdata, id_var, time_var, mknots, sknots,
  gknots, n_nodes = 20L, skew_start = -0.02,
  opt_func = .opt_default, n_threads = 1L, sparse_hess = FALSE, B = NULL,
  Psi = NULL, Sigma = NULL, omega = NULL, alpha = NULL, delta = NULL,
  gamma = NULL, va_par = NULL){
  # checks
  check_knots_num <- function(x){
    if(length(x) == 0L)
      stop("Model without terms is not implemented")
    length(x) == 0L ||
      (is.numeric(x) && all(is.finite(x)) && length(x) > 1L)
  }
  check_knots_int <- function(x){
    if(is.integer(x) && x == 0L)
      stop("Model without terms is not implemented")
    is.integer(x) && length(x) == 1L && x != 1L
  }
  stopifnot(
    inherits(sformula, "formula"),
    inherits(mformula, "formula"),
    is.data.frame(sdata),
    is.data.frame(mdata),
    !missing(id_var),
    !missing(time_var),
    check_knots_num(mknots) || check_knots_int(mknots),
    check_knots_num(gknots) || check_knots_int(gknots),
    check_knots_num(sknots) || check_knots_int(sknots),
    is.integer(n_nodes), length(n_nodes) == 1L && n_nodes > 0L,
    is.integer(n_threads), length(n_threads) == 1L, n_threads > 0L,
    is.logical(sparse_hess), length(sparse_hess) == 1L)

  id_var <- substitute(id_var)
  time_var <- substitute(time_var)

  # get data needed for the estimation
  mark <- get_marker_start_params(
    mformula, data = mdata, time_var = time_var, id_var = id_var,
    mknots = mknots, gknots = gknots,
    need_start_vals =
      is.null(B) || is.null(Psi) || is.null(Sigma) || is.null(gamma))
  sr_dat <- get_surv_start_params(
    formula = sformula, data = sdata, mformula = mformula, mdata = mdata,
    id_var = id_var, time_var = time_var, knots = sknots, n_nodes = n_nodes,
    need_start_vals = is.null(omega) || is.null(alpha) || is.null(delta))

  # assign the variables we need
  if(is.null(gamma))
    gamma <- mark$gamma
  if(is.null(B))
    B <- mark$B
  if(is.null(Psi))
    Psi <- mark$Psi
  if(is.null(Sigma))
    Sigma <- mark$Sigma
  if(is.null(omega))
    omega <- sr_dat$omega
  if(is.null(alpha))
    alpha <- sr_dat$alpha
  if(is.null(delta))
    delta <- sr_dat$delta
  mknots <- mark$mknots
  gknots <- mark$gknots
  sknots <- sr_dat$knots

  markers <- mark$Y
  X <- mark$X
  m_id <- mark$id
  m_time <- mark$time

  m_order <- order(m_id)
  markers <- markers[, m_order, drop = FALSE]
  X       <- X[, m_order, drop = FALSE]
  m_id    <- m_id   [m_order]
  m_time  <- m_time [m_order]

  outcomes <- sr_dat$Y
  s_id <- sr_dat$id
  Z <- sr_dat$Z

  s_order  <- order(s_id)
  outcomes <- outcomes[, s_order]
  s_id     <- s_id    [s_order]
  Z        <- Z       [, s_order, drop = FALSE]

  # checks
  n_y <- NROW(mark$Y)
  dim_m <- length(mknots)
  dim_g <- length(gknots)
  K <- n_y * dim_m
  n_Z <- NROW(Z)
  n_groups <- length(unique(s_id))

  c_num <- function(x)
    is.numeric(x) && all(is.finite(x))
  stopifnot(
    c_num(gamma), is.matrix(gamma), NCOL(gamma) == n_y,
    c_num(B), is.matrix(B), NROW(B) == dim_g, NCOL(B) == n_y,
    c_num(Psi), is.matrix(Psi), NROW(Psi) == K, NCOL(Psi) == K,
    c_num(Sigma), is.matrix(Sigma), NROW(Sigma) == n_y, NCOL(Sigma) == n_y,
    c_num(alpha), length(alpha) == n_y,
    c_num(delta), length(delta) == n_Z,
    c_num(omega), length(omega) == length(sknots),
    c_num(markers), is.matrix(markers),
    c_num(X), is.matrix(X), NCOL(X) == NCOL(markers),
    setequal(unique(s_id), unique(m_id)),
    n_groups == length(s_id),
    check_knots_num(mknots),
    check_knots_num(sknots),
    check_knots_num(gknots))

  # prepare VA parameters
  if(missed_va <- is.null(va_par))
    va_par <- .get_joint_va_start(skew_start = skew_start, Psi = Psi,
                                  n_groups = n_groups)
  stopifnot(length(va_par) == n_groups * (2L * K + (K * (K + 1L)) / 2L),
            all(is.finite(va_par)))

  # make AD func
  n_markers <- table(m_id)

  data <- list(
    markers = markers, n_markers = n_markers, m_time = m_time, X = X,
    mknots = mknots, gknots = gknots,
    tstart = outcomes[1, ], tstop = outcomes[2, ], outcomes = outcomes[3, ],
    sknots = sknots, Z = Z,
    n_threads = n_threads, sparse_hess = sparse_hess, n_nodes = n_nodes)
  parameters <- list(
    gamma = gamma, B = B, Psi = .cov_to_theta(Psi),
    Sigma = .cov_to_theta(Sigma), delta = delta, omega = omega,
    alpha = alpha, va_par = va_par)

  func <- get_joint_funcs(data = data, parameters = parameters)

  # create parameter vector and return
  par <- c(gamma, B, .cov_to_theta(Psi), .cov_to_theta(Sigma), delta,
           omega, alpha, va_par)
  names(par) <- c(
    outer(paste0("gamma:", rownames(X)), rownames(mark$Y),
          function(x, y) paste0(x, ".", y)),
    outer(paste0("B:g", seq_len(NROW(B))), rownames(mark$Y),
          function(x, y) paste0(x, ".", y)),
    paste0("Psi:", names(.cov_to_theta(Psi))),
    paste0("Sigma:", names(.cov_to_theta(Sigma))),
    paste0("delta:", rownames(Z)),
    paste0("omega:b", seq_along(omega)),
    paste0("alpha:", rownames(mark$Y)),
    names(va_par))

  out <- list(
    fn = function(x, ...){
      -joint_funcs_eval_lb(p = func, x)
    },
    gr = function(x, ...){
      -joint_funcs_eval_grad(p = func, x)
    },
    he = function(x, ...){
      stop("he not implemented")
    },
    get_params = function(x)
      stop("get_params not implemented"))

  if(missed_va){
    # make a quick search where we set all VA parameters to the same value
    is_va_regex <- "^g\\d+"
    par <- .opt_sub(
      par = par, which_par = is_va_regex, extract_name = "^(g\\d+:)(.+)",
      replacement = "\\2", fn = out$fn, gr = out$gr,
      opt_func = opt_func, control = .opt_func_quick_control)

    # refine by optimizing each individually
    is_va <- which(grepl(is_va_regex, names(par)))
    fn_va <- function(x, ...){
      par[is_va] <- x
      out$fn(par)
    }
    gr_va <- function(x, ...){
      par[is_va] <- x
      out$gr(par)[is_va]
    }

    opt_va <- opt_func(
      par[is_va], fn_va, gr_va, control = .opt_func_quick_control)
    par[is_va] <- opt_va$par
  }

  out$par <- par

  # clean up enviroment
  rm(list = ls()[!ls() %in% c("func", "out")])
  out
}

