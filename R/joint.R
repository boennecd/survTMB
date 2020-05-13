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

#' @importFrom stats model.response model.frame
#' @importFrom reshape2 melt
#' @importFrom splines ns
#' @importFrom lme4 lmer lmerControl fixef VarCorr .makeCC
get_marker_start_params <- function(formula, data, knots, time_var, id_var,
                                    need_start_vals = TRUE){
  out <- list(Y = model.response(model.frame(formula, data)),
              id = eval(id_var, data),
              time = eval(time_var, data))

  if(is.integer(knots) && length(knots) == 1L)
    knots <- .get_joint_knots(knots, out$time)
  out$knots <- knots

  if(!need_start_vals){
    out$Psi <- out$Sigma <- out$B <- out$ll <- numeric()
    return(out)
  }

  # get data.frame to work with
  Y_names <- .get_Y_names(formula)
  data <- melt(
    data, id.vars = c(deparse(id_var), deparse(time_var)),
    measure.vars = Y_names)

  # fit mixed model
  bk <- knots[ c(1L, length(knots))]
  ik <- knots[-c(1L, length(knots))]
  frm <- substitute(
    value ~
      ns(ti, knots = ik, Boundary.knots = bk, intercept = TRUE) : variable - 1L +
      (ns(ti, knots = ik, Boundary.knots = bk, intercept = TRUE) : variable - 1L | i),
    list(ti = time_var, i = id_var, ik = ik, bk = bk))

  fit <- lmer(frm, data, control = lmerControl(
    check.conv.grad = .makeCC("ignore", tol = 1e-3, relTol = NULL)))

  n_y <- NCOL(out$Y)
  B <- matrix(fixef(fit), nc = n_y)
  vc <- VarCorr(fit)
  Psi <- vc$id
  attr(Psi, "correlation") <- attr(Psi, "stddev") <- NULL
  dimnames(Psi) <- NULL

  Sigma <- diag(attr(vc, "sc")^2, n_y)

  out[c("B", "Psi", "Sigma", "ll")] <-
    list(B, Psi, Sigma, ll = c(logLik(fit)))
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
    list(Y = outcome$Y[keep, ], Z = Z[keep, , drop = FALSE], id = id[keep],
         keep = keep)
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
      sdata, sdata, id = id, tstart = tstart, tstop = dtstop,
      ev = event(dtstop, devent)),
      list(tstart = dtstart, dtstop = dtstop, devent = devent,
           id = id_var))
    out <- eval(tcall)

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
  cformula <- update(sformula, Surv(tstart, tstop, ev) ~ . - 1)
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
  SZ <- model.matrix(terms(mf), mf)[Sy$keep, ]
  offsets <- drop(SZ %*% c(delta, alpha))
  bk <- knots[ c(1L, length(knots))]
  ik <- knots[-c(1L, length(knots))]
  omega <- rep(0, length(knots))

  tstart <- Sy$Y[, 1]
  tstop  <- Sy$Y[, 2]
  Y      <- Sy$Y[, 3]

  func <- function(omega, ..., grad)
    -drop(joint_start_baseline(
      Y = Y, tstart = tstart, tstop = tstop, omega = omega,
      offsets = offsets, n_nodes = n_nodes, bound_knots = bk,
      inter_knots = ik, grad = grad))
  fn <- func
  formals(fn)$grad <- FALSE
  gr <- func
  formals(gr)$grad <- TRUE

  opt_out <- optim(omega, fn, gr, method = "BFGS")
  stopifnot(opt_out$convergence == 0L)

  omega <- opt_out$par
  out[c("delta", "alpha", "omega", "ll")] <- list(
    delta, alpha, omega, -opt_out$value)
  out
}

.get_joint_va_start <- function(skew_start, K, Psi, n_groups) {
  va_pars <- local({
    alpha <- .gamma_to_alpha(skew_start)
    if(length(alpha) == 1L)
      alpha <- rep(alpha, K)
    stopifnot(length(alpha) == K)

    nu <- sqrt(2 / pi) * alpha / sqrt(1 + alpha^2)
    omega_denom <- sqrt(1 - nu^2)

    sds <- sqrt(diag(Psi))

    omega <- sds / omega_denom
    dnu <- omega * nu
    Omega <- Psi + outer(dnu, dnu)

    xi <- - dnu
    c(xi = xi, .cov_to_theta(Omega), alpha = alpha)
  })

  va_pars <- structure(
    rep(va_pars, n_groups),
    names = outer(names(va_pars), seq_len(n_groups),
                  function(x, y) paste0("g", y, ":", x)))
}

#' @importFrom stats optim
.opt_default_CG <- function(par, fn, gr, ...){
  cl <- match.call()
  cl[[1L]] <- quote(stats::optim)
  if(is.null(cl$method))
    cl$method <- "CG"

  out <- eval(cl, parent.frame())
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
#' @param mknots vector with knots for the process mean or number of knots
#'               to use.
#' @param mknots vector with knots for the baseline hazard or number of
#'               knots to use.
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
#' @param va_par staring value for variational parameters.
#'
#' @export
make_joint_ADFun <- function(
  sformula, mformula, sdata, mdata, id_var, time_var, mknots, sknots,
  n_nodes = 20L, skew_start = .alpha_to_gamma(-1),
  opt_func = .opt_default_CG, n_threads = 1L, sparse_hess = FALSE, B = NULL,
  Psi = NULL, Sigma = NULL, omega = NULL, alpha = NULL, delta = NULL,
  va_par = NULL){
  # checks
  check_knots_num <- function(x)
    is.numeric(x) && all(is.finite(x)) && length(x) > 1L
  check_knots_int <- function(x)
    is.integer(x) && length(x) == 1L && x > 1L
  stopifnot(
    inherits(sformula, "formula"),
    inherits(mformula, "formula"),
    is.data.frame(sdata),
    is.data.frame(mdata),
    !missing(id_var),
    !missing(time_var),
    check_knots_num(mknots) || check_knots_int(mknots),
    check_knots_num(sknots) || check_knots_int(sknots),
    is.integer(n_nodes), length(n_nodes) == 1L && n_nodes > 0L,
    is.integer(n_threads), length(n_threads) == 1L, n_threads > 0L,
    is.logical(sparse_hess), length(sparse_hess) == 1L)

  id_var <- substitute(id_var)
  time_var <- substitute(time_var)

  # get data needed for the estimation
  mark <- get_marker_start_params(
    mformula, data = mdata, time_var = time_var, id_var = id_var,
    knots = mknots,
    need_start_vals = is.null(B) || is.null(Psi) || is.null(Sigma))
  sr_dat <- get_surv_start_params(
    formula = sformula, data = sdata, mformula = mformula, mdata = mdata,
    id_var = id_var, time_var = time_var, knots = sknots, n_nodes = n_nodes,
    need_start_vals = is.null(omega) || is.null(alpha) || is.null(delta))

  # assign the variables we need
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
  mknots <- mark$knots
  sknots <- sr_dat$knots

  markers <- mark$Y
  m_id <- mark$id
  m_time <- mark$time

  m_order <- order(m_id)
  markers <- markers[m_order, , drop = FALSE]
  m_id    <- m_id   [m_order]
  m_time  <- m_time [m_order]

  outcomes <- sr_dat$Y
  s_id <- sr_dat$id
  Z <- sr_dat$Z

  s_order  <- order(s_id)
  outcomes <- outcomes[s_order, ]
  s_id     <- s_id    [s_order]
  Z        <- Z       [s_order, , drop = FALSE]

  # checks
  n_y <- NCOL(mark$Y)
  dim_m <- length(mknots)
  K <- n_y * dim_m
  n_Z <- NCOL(Z)
  n_groups <- length(unique(s_id))

  c_num <- function(x)
    is.numeric(x) && all(is.finite(x))
  stopifnot(
    c_num(B), is.matrix(B), NROW(B) == dim_m, NCOL(B) == n_y,
    c_num(Psi), is.matrix(Psi), NROW(Psi) == K, NCOL(Psi) == K,
    c_num(Sigma), is.matrix(Sigma), NROW(Sigma) == n_y, NCOL(Sigma) == n_y,
    c_num(alpha), length(alpha) == n_y,
    c_num(delta), length(delta) == n_Z,
    c_num(omega), length(omega) == length(sknots),
    setequal(unique(s_id), unique(m_id)),
    n_groups == length(s_id),
    check_knots_num(mknots),
    check_knots_num(sknots))

  # prepare VA parameters
  if(missed_va <- is.null(va_par))
    va_par <- .get_joint_va_start(skew_start = skew_start, K = K, Psi = Psi,
                                  n_groups = n_groups)
  stopifnot(length(va_par) == n_groups * (2L * K + (K * (K + 1L)) / 2L),
            all(is.finite(va_par)))

  # make AD func
  n_markers <- table(m_id)
  markers <- t(markers)
  Z <- t(Z)

  data <- list(
    markers = markers, n_markers = n_markers, m_time = m_time,
    mknots = mknots,
    tstart = outcomes[, 1], tstop = outcomes[, 2], outcomes = outcomes[, 3],
    sknots = sknots, Z = Z,
    n_threads = n_threads, sparse_hess = sparse_hess, n_nodes = n_nodes)
  parameters <- list(
    B = B, Psi = .cov_to_theta(Psi), Sigma = .cov_to_theta(Sigma),
    delta = delta, omega = omega, alpha = alpha,
    va_par = va_par)

  func <- get_joint_funcs(data = data, parameters = parameters)

  # create parameter vector and return
  par <- c(B, .cov_to_theta(Psi), .cov_to_theta(Sigma), delta, omega, alpha,
           va_par)
  names(par) <- c(
    outer(paste0("B:m", seq_len(NROW(B))), colnames(mark$Y),
          function(x, y) paste0(x, ".", y)),
    paste0("Psi:", names(.cov_to_theta(Psi))),
    paste0("Sigma:", names(.cov_to_theta(Sigma))),
    paste0("delta:", rownames(Z)),
    paste0("omega:b", seq_along(omega)),
    paste0("alpha:", colnames(mark$Y)),
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
    })

  if(missed_va){
    ctrl_use <- list(
      reltol = .Machine$double.eps^(1/5), maxit = 100L, abstol = 0)

    # find better starting values for the VA parameters
    is_va    <- which(grepl("^g\\d+:", names(par)))
    is_alpha <- which(grepl("^g\\d+:alpha", names(par)))
    is_xi    <- which(grepl("^g\\d+:xi", names(par)))
    is_psi   <- which(grepl("^g\\d+:(log_sd|L)", names(par), perl = TRUE))
    stopifnot(length(is_va) > 0,
              length(is_alpha) > 0,
              length(is_xi) > 0,
              length(is_psi) > 0)

    # make a quick search where we set all VA parameters to the same value
    va_par_start <- local({
      nam_par <- gsub("^(g\\d+:)(.+)", "\\2", names(par), perl = TRUE)
      . <- function(x_idx){
        nams <- unique(nam_par[x_idx])
        sapply(nams, function(z) x_idx[z == nam_par[x_idx]])
      }

      va_map <- .(is_va)
      par_sub <- apply(va_map, 2L, function(indices) mean(par[indices]))

      get_par <- function(x){
        out <- par
        for(i in 1:NCOL(va_map))
          out[va_map[, i]] <- x[i]
        out
      }

      fn_va <- function(x, ...){
        out$fn(get_par(x))
      }
      gr_va <- function(x, ...){
        grad <- out$gr(get_par(x))
        apply(va_map, 2L, function(indices) sum(grad[indices]))
      }

      opt_va <- opt_func(
        par_sub, fn_va, gr_va, control = ctrl_use)

      opt_va$par
    })

    # refine by optimizing each individually
    par[is_va] <- va_par_start

    fn_va <- function(x, ...){
      par[is_va] <- x
      out$fn(par)
    }
    gr_va <- function(x, ...){
      par[is_va] <- x
      out$gr(par)[is_va]
    }

    opt_va <- opt_func(par[is_va], fn_va, gr_va, control = ctrl_use)
    par[is_va] <- opt_va$par
  }

  out$par <- par
  out
}

