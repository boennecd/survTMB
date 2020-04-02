.laplace_char <- "Laplace"
.gva_char     <- "GVA"
.snva_char    <- "SNVA"

.alpha_to_gamma <- function(alpha){
  mu <- sqrt(2 / pi) * alpha / sqrt(1 + alpha^2)
  (4 - pi) / 2 * mu^3 / (1 - mu^2)^(3/2)
}
.gamma_to_alpha <- function(gamma){
  cf <- 2 * gamma / (4 - pi)
  v1 <- local({
    cf <- ( cf)^(1/3)
    cf / sqrt(1 + cf^(2))
  })
  v2 <- local({
    cf <- (-cf)^(1/3)
    -cf / sqrt(1 + cf^(2))
  })
  mu <- ifelse(gamma > 0, v1, v2)
  o <- sqrt(pi) * mu
  o / sqrt(2 - pi * mu^2)
}

#' @importFrom stats optim
.opt_default <- function(par, fn, gr, ...){
  cl <- match.call()
  cl[[1L]] <- quote(stats::optim)
  if(is.null(cl$method))
    cl$method <- "BFGS"
  if(is.null(cl$control))
    cl$control <- list(reltol = .Machine$double.eps^(1/4), maxit = 1000)

  out <- eval(cl, parent.frame())
  stopifnot(out$convergence == 0L)
  out
}

#' @importFrom survival survfit
Shat <- function(obj){
  ## predicted survival for individuals (adjusted for covariates)
  newobj <- survfit(obj, se.fit = FALSE)
  surv <- newobj$surv
  rr <- try(predict(obj, type = "risk"), silent = TRUE)
  if (inherits(rr, "try-error"))
    rr <- 1
  surv2 <- surv[match(obj$y[, ncol(obj$y) - 1], newobj$time)]

  surv2^rr
}

.cov_to_theta <- function(theta){
  n_rng <- NCOL(theta)
  ch <- t(chol(theta))
  log_sd <- structure(log(diag(ch)), names = paste0("log_sd", 1:n_rng))
  keep <- lower.tri(ch)
  lower_tri <-(diag(diag(ch)^(-1), n_rng) %*% ch)[keep]
  if(n_rng > 1L)
    names(lower_tri) <-
      outer(1:n_rng, 1:n_rng, function(x, y) paste0("L", x, ".", y))[keep]
  c(log_sd, lower_tri)
}

.theta_to_cov <- function(theta){
  dim <- .5 * (sqrt(8 * length(theta) + 1) - 1)
  if(dim < 2L)
    return(exp(2 * theta))

  L <- diag(dim)
  L[lower.tri(L)] <- tail(theta, -dim)
  L <- diag(exp(head(theta, dim))) %*% L
  tcrossprod(L)
}


.set_use_own_VA_method <- with(new.env(), {
  .use_own_VA_method <- TRUE
  function(use_own_VA_method){
    stopifnot(
      is.logical(use_own_VA_method),
      length(use_own_VA_method) == 1L,
      !is.na(use_own_VA_method))
    .use_own_VA_method <<- use_own_VA_method
    invisible(use_own_VA_method)
  }
})

.get_use_own_VA_method <- with(
  environment(.set_use_own_VA_method),
  function()
    .use_own_VA_method)

#####
# Define function to get objective function.
#
# Args:
#   formula: two-sided formula where the left-handside is a `Surv` object
#            and the right-handside is the fixed effects.
#   data: `data.frame` with variables used in the model.
#   df: integer scalar with the degrees of freedom used for the baseline
#       spline.
#   Z: one-sided formula where the right-handside are the random effects.
#   cluster: vector with integers or factors for group identifiers.
#   do_setup: character vector indicating which approximation to setup.
#             It is included to test computation time.
#   n_nodes: number of nodes to use in (adaptive) Gauss–-Hermite quadrature.
#   param_type: characters for the parameterization used with the SNVA.
#   skew_start: starting value for the Pearson's moment coefficient of
#               skewness parameter when a SNVA is used.
#   link: character specifying the link function.
#   theta: starting values for covariance matrix.
#   beta: starting values for fixed effect coefficients.
#   opt_func: general optimization function to get starting values.
#
# Returns:
#   TODO: what?
#
#' @importFrom TMB MakeADFun
#' @importFrom stats model.frame model.response terms model.matrix lm lm.fit predict qnorm
#' @importFrom rstpm2 nsx
#' @importFrom survival coxph frailty Surv
#' @export
make_gsm_ADFun <- function(
  formula, data, df, Z, cluster, do_setup = c("Laplace", "GVA", "SNVA"),
  n_nodes = 20L, param_type = c("DP", "CP_trans", "CP"),
  link = c("PH", "PO", "probit"), theta = NULL, beta = NULL,
  opt_func = .opt_default, n_threads = 1L,
  skew_start = .alpha_to_gamma(-1)){
  link <- link[1]
  param_type <- param_type[1]
  skew_boundary <- 0.99527
  stopifnot(
    is.integer(df), df > 0L, inherits(formula, "formula"),
    inherits(Z, "formula"), is.data.frame(data),
    !missing(cluster),
    all(do_setup %in% c(.laplace_char, .gva_char, .snva_char)),
    is.integer(n_nodes), length(n_nodes) == 1L, n_nodes > 0L,
    link %in% c("PH", "PO", "probit"),
    param_type %in% c("DP", "CP_trans", "CP"),
    is.integer(n_threads) && n_threads > 0L && length(n_threads) == 1L,
    is.numeric(skew_start), length(skew_start) == 1L)
  eval(bquote(stopifnot(
    .(-skew_boundary) < skew_start && skew_start < .(skew_boundary))))

  #####
  # get the cluster variable
  cluster <- substitute(cluster)
  grp <- eval(cluster, data)
  stopifnot(is.factor(grp) || is.integer(grp) || is.character(grp),
            isTRUE(length(grp) == NROW(data)))
  if(!is.factor(grp))
    grp <- as.factor(grp)
  grp <- as.integer(grp)
  n_grp <- length(unique(grp))

  #####
  # change order of data set and assign group size variables
  ord <- order(grp)
  data <- data[ord, ]
  grp <- grp[ord]
  grp_size <- table(grp)

  #####
  # get design matrices and outcome
  mf_X <- model.frame(formula, data = data)
  mt_X <- terms(mf_X)
  X <- model.matrix(mt_X, mf_X)
  n_x_fix <- NCOL(X)

  y <- model.response(mf_X)
  stopifnot(inherits(y, "Surv"), isTRUE(attr(y, "type") == "right"))
  event <- y[, 2]
  tobs  <- y[, 1]

  formula_Z <- Z
  mf_Z <- model.frame(formula_Z, data = data)
  mt_Z <- terms(mf_Z)
  Z <- model.matrix(mt_Z, mf_Z)
  n_rng <- NCOL(Z)
  stopifnot(NCOL(X) > 0, NCOL(Z) > 0, NROW(X) == NROW(Z),
            length(y) == NROW(X))

  #####
  # add time-varying baseline
  time_var <- formula[[2L]][[2L]]
  formula_b <- eval(
    bquote(~ nsx(log(.(time_var)), df = .(df), intercept = FALSE) - 1))
  mt_b <- terms(model.frame(formula_b, data = data[event > 0, ]))
  X <- cbind(X, model.matrix(mt_b, data))

  #####
  # approximate derivatives w/ finite difference
  is_fix <- 1:n_x_fix
  XD <- X
  XD[,  is_fix] <- 0.
  XD[, -is_fix] <- local({
    dt <- .Machine$double.eps^(1/3)
    dat_m1 <- dat_p1 <- eval(bquote(data.frame(.(time_var))), data)
    x <- dat_m1[[1L]]
    dat_p1[[1L]] <- x * exp( dt)
    dat_m1[[1L]] <- x * exp(-dt)
    (model.matrix(mt_b, dat_p1) - model.matrix(mt_b, dat_m1)) / 2 / dt / x
  })

  #####
  # get starting values w/ coxph and then a lm fit
  need_theta <- is.null(theta)
  need_beta  <- is.null(beta)

  stopifnot(isTRUE(attr(mt_Z, "intercept") == 1L),
            colnames(Z)[1] == "(Intercept)")

  inits <- list()

  if(need_beta){
    inits$link_Hat <- local({
        cox_fit <- coxph(formula, data, model = TRUE)
        S_hat <- Shat(cox_fit)
        if(link == "PH")
          pmax(log(.Machine$double.eps) / 4, log(-log(S_hat)))
        else if(link == "PO")
          log((1 - S_hat) / S_hat)
        else if(link == "probit")
          -qnorm(S_hat)
        else
          stop(sprintf("%s not implemented", sQuote(link)))
      })

    inits$coef <- local({
      keep <- event > 0
      lm.fit(x = X[keep, , drop = FALSE],
             y = inits$link_Hat[keep])$coefficients
    })
  }

  if(need_theta)
    inits$theta <- local({
      cox_frm <- eval(bquote(
        update(formula, . ~ . + frailty(.(cluster),
                                        distribution = "gaussian"))))

      cox_fit <- coxph(cox_frm, data)

      cox_fit$history[[1L]]$theta
    })

  #####
  # setup ADFun object for the Laplace approximation
  data_ad_func <- list(
    tobs = tobs, event = event, X = X, XD = XD, Z = Z, grp = grp - 1L,
    link = link, grp_size = grp_size, n_threads = n_threads)

  # the user may have provided values
  theta <- if(!need_theta){
    stopifnot(all(dim(theta) == n_rng))
    theta

  } else local({
    out <- matrix(0., n_rng, n_rng)
    out[1L] <- inits$theta
    if(n_rng > 1L)
      diag(out)[-1L] <- rep(1, n_rng - 1L)
    out
  })

  theta <- .cov_to_theta(theta)

  beta <- if(!is.null(beta)){
    stopifnot(length(beta) == NCOL(X))
    beta
  } else
    inits$coef
  names(beta) <- colnames(X)

  # assign parameter list
  params = list(
    eps = .Machine$double.eps^(1/2), kappa = 1e8, b = beta,
    theta = theta)

  if(.laplace_char %in% do_setup){
  # get Laplace AD function
  adfunc_laplace <- local({
    data_ad_func <- c(
      list(app_type = .laplace_char), data_ad_func)
    # TODO: initialize random effects in a smarter way...
    params$u <- matrix(0., n_rng, n_grp)
    MakeADFun(
      data = data_ad_func, parameters = params, DLL = "survTMB",
      silent = TRUE, random = "u")
  })

  # we make a wrapper object to account for the eps and kappa and allow the
  # user to change these
  laplace_out <- with(new.env(), {
    eps <- adfunc_laplace$par["eps"]
    kappa <- adfunc_laplace$par["kappa"]
    fn <- adfunc_laplace$fn
    gr <- adfunc_laplace$gr
    he <- adfunc_laplace$he
    get_x <- function(x)
      c(eps = eps, kappa = kappa, x)

    out <- adfunc_laplace[
      !names(adfunc_laplace) %in% c("par", "fn", "gr", "he")]

    par <- adfunc_laplace$par[-(1:2)]
    names(par)[seq_along(inits$coef)] <- names(inits$coef)

    c(list(
      par = par,
      fn = function(x, ...){ fn(get_x(x))                               },
      gr = function(x, ...){ gr(get_x(x))[-(1:2)]                       },
      he = function(x, ...){ he(get_x(x))[-(1:2), -(1:2), drop = FALSE] },
      # function to set penalty parameters
      update_pen = function(eps, kappa){
        p_env <- parent.env(environment())
        if(!missing(eps))
          assign("eps"  , eps  , p_env)
        if(!missing(kappa))
          assign("kappa", kappa, p_env)
        invisible(with(p_env, c(eps = eps, kappa = kappa)))
      }
    ), out)
  })
  } else
    laplace_out <- NULL

  get_par_va <- function(params)
    with(params, {
      names(b)        <- rep("b", length(b))
      names(theta)    <- rep("theta", length(theta))
      names(theta_VA) <- rep("theta_VA", length(theta_VA))
      c(eps = eps, kappa = kappa, b, theta, theta_VA)
    })

  #####
  # setup ADFun object for the GVA
  if(.gva_char %in% do_setup || .snva_char %in% do_setup){
  # set the initial values
  n_mu     <- n_rng
  n_Lambda <- (n_rng * (n_rng + 1L)) / 2L
  n_p_grp  <- n_mu + n_Lambda
  theta_VA <- rep(NA_real_, n_p_grp * n_grp)

  # set means to zero
  idx_mean <- sapply(1:n_grp - 1L, function(i) i * n_p_grp + 1:n_mu)
  theta_VA[idx_mean] <- 0.
  # set parameters for Lambda
  idx_Lambda <- sapply(1:n_grp - 1L, function(i)
    i * n_p_grp + n_mu + 1:n_Lambda)
  theta_VA[idx_Lambda] <- params$theta
  # set names
  theta_VA_names <- c(paste0("mu", 1:n_rng), names(params$theta))
  theta_VA_names <- c(outer(theta_VA_names, 1:n_grp, paste, sep = ":"))
  names(theta_VA) <- theta_VA_names

  get_gva_out <- function(theta_VA){
    adfunc_VA <- local({
      data_ad_func <- c(
        list(app_type = .gva_char), data_ad_func,
        list(n_nodes = n_nodes))
      params$theta_VA <- theta_VA

      if(.get_use_own_VA_method()){
        data_ad_func$param_type <- "DP"
        within(list(), {
          ptr <- get_VA_funcs(data = data_ad_func, parameters = params)
          fn <- function(par)
            VA_funcs_eval_lb(ptr, par)
          gr <- function(par)
            drop(VA_funcs_eval_grad(ptr, par))
          he <- function(par)
            stop("he not implemented")
          par <- get_par_va(params)
        })

      } else
        MakeADFun(
          data = data_ad_func, parameters = params, DLL = "survTMB",
          silent = TRUE)
    })

    # we make a wrapper object to account for the eps and kappa and allow the
    # user to change these
    gva_out <- with(new.env(), {
      eps <- adfunc_VA$par["eps"]
      kappa <- adfunc_VA$par["kappa"]
      fn <- adfunc_VA$fn
      gr <- adfunc_VA$gr
      he <- adfunc_VA$he
      get_x <- function(x)
        c(eps = eps, kappa = kappa, x)

      out <- adfunc_VA[
        !names(adfunc_VA) %in% c("par", "fn", "gr", "he")]

      par <- adfunc_VA$par[-(1:2)]
      names(par)[seq_along(inits$coef)] <- names(inits$coef)
      idx_va <- (length(par) - length(theta_VA_names) + 1):length(par)
      names(par)[idx_va] <-
        theta_VA_names

      c(list(
        par = par,
        fn = function(x, ...){ fn(get_x(x))                               },
        gr = function(x, ...){ gr(get_x(x))[-(1:2)]                       },
        he = function(x, ...){ he(get_x(x))[-(1:2), -(1:2), drop = FALSE] },
        # function to set penalty parameters
        update_pen = function(eps, kappa){
          p_env <- parent.env(environment())
          if(!missing(eps))
            assign("eps"  , eps  , p_env)
          if(!missing(kappa))
            assign("kappa", kappa, p_env)
          invisible(with(p_env, c(eps = eps, kappa = kappa)))
        },
        # function to get parameters
        get_params = function(x)
          x[-idx_va],
        control = list(maxit = 1000L)
      ), out)
    })
  }

  # find initial VA params
  func <- get_gva_out(theta_VA)
  coefs_start <- c(params$b, params$theta)
  do_drop <- seq_along(coefs_start)
  get_x <- function(x)
    c(coefs_start, x)

  fn <- function(x)
    func$fn(get_x(x))
  gr <- function(x)
    func$gr(get_x(x))[-do_drop]

  opt_out <- opt_func(theta_VA, fn = fn, gr = gr)
  gva_out <- get_gva_out(opt_out$par)

  } else
    gva_out <- NULL

  #####
  # setup ADFun object for the SNVA
  if(.snva_char %in% do_setup){
  # set the initial values
  n_mu     <- n_rng
  n_rho    <- n_rng
  n_Lambda <- (n_rng * (n_rng + 1L)) / 2L
  n_p_grp  <- n_mu + n_Lambda + n_rho
  theta_VA <- rep(NA_real_, n_p_grp * n_grp)

  #####
  # setup starting values for VA parameters
  # find GVA solution
  gva_opt <- with(gva_out, opt_func(par, fn, gr))
  gva_va_vals <- gva_opt$par[-seq_len(length(beta) + length(theta))]

  params$b     <- gva_opt$par[1:length(beta)]
  params$theta <- gva_opt$par[1:length(theta) + length(beta)]

  if(length(skew_start) == 1L && n_rng > 1L)
    skew_start <- rep(skew_start, n_rng)
  n_p_grp_gva <-  n_p_grp - n_rho

  theta_VA <- if(param_type == "DP"){
    alpha <- .gamma_to_alpha(skew_start)
    nu <- sqrt(2 / pi) * alpha / sqrt(1 + alpha^2)
    omega_denom <- sqrt(1 - nu^2)
    n_lower_tri <- (n_rng * (n_rng - 1L)) / 2L

    vapply(1:n_grp, function(i){
      # setup mean and VA variance
      gva_par <- gva_va_vals[(i - 1L) * n_p_grp_gva + 1:n_p_grp_gva]
      sds <- exp(gva_par[1:n_mu + n_mu])

      Sig <- diag(n_rng)
      if(n_rng > 1L)
        Sig[lower.tri(Sig)] <- tail(gva_par, n_lower_tri)
      Sig <- diag(sds, n_rng) %*% Sig
      Sig <- tcrossprod(Sig)

      # compute VA parameters and return
      omega <- sds / omega_denom
      dnu <- omega * nu
      Omega <- Sig + outer(dnu, dnu)

      xi <- gva_par[1:n_mu] - dnu
      c(xi, .cov_to_theta(Omega), alpha)
    }, numeric(n_p_grp), USE.NAMES = FALSE)
  } else {
    stopifnot(param_type == "CP_trans")
    skew_trans <-
      log((skew_boundary + skew_start) / (skew_boundary - skew_start))
    vapply(1:n_grp, function(i){
      # setup mean and VA variance
      gva_par <- gva_va_vals[(i - 1L) * n_p_grp_gva + 1:n_p_grp_gva]
      c(gva_par, skew_trans)
    }, numeric(n_p_grp), USE.NAMES = FALSE)
  }

  # set names
  theta_VA_names <- local({
    keep <- which(lower.tri(diag(n_rng)))
    L <-
      outer(1:n_rng, 1:n_rng, function(x, y) paste0("L", x, ".", y))[keep]

    proto <- c(paste0("mu", 1:n_rng), paste0("log_sd", 1:n_rng), L,
               paste0("skew", 1:n_rng))

    c(outer(proto, 1:n_grp, paste, sep = ":"))
  })
  theta_VA <- structure(c(theta_VA), names = theta_VA_names)

  adfunc_VA <- local({
    data_ad_func <- c(
      list(app_type = .snva_char), data_ad_func,
      list(n_nodes = n_nodes, param_type = param_type))
    params$theta_VA <- theta_VA

    if(.get_use_own_VA_method())
      within(list(), {
        ptr <- get_VA_funcs(data = data_ad_func, parameters = params)
        fn <- function(par)
          VA_funcs_eval_lb(ptr, par)
        gr <- function(par)
          drop(VA_funcs_eval_grad(ptr, par))
        he <- function(par)
          stop("he not implemented")
        par <- get_par_va(params)
      })
    else
      MakeADFun(
        data = data_ad_func, parameters = params, DLL = "survTMB",
        silent = TRUE)
  })

  # we make a wrapper object to account for the eps and kappa and allow the
  # user to change these
  get_snva_out <- function(theta_VA){
    with(new.env(), {
      eps <- adfunc_VA$par["eps"]
      kappa <- adfunc_VA$par["kappa"]
      fn <- adfunc_VA$fn
      gr <- adfunc_VA$gr
      he <- adfunc_VA$he
      get_x <- function(x)
        c(eps = eps, kappa = kappa, x)

      out <- adfunc_VA[
        !names(adfunc_VA) %in% c("par", "fn", "gr", "he")]

      par <- adfunc_VA$par[-(1:2)]
      names(par)[seq_along(beta)] <- names(beta)
      idx_va <- (length(par) - length(theta_VA_names) + 1):length(par)
      names(par)[idx_va] <- theta_VA_names

      c(list(
        par = par,
        fn = function(x, ...){ fn(get_x(x))                               },
        gr = function(x, ...){ gr(get_x(x))[-(1:2)]                       },
        he = function(x, ...){ he(get_x(x))[-(1:2), -(1:2), drop = FALSE] },
        # function to set penalty parameters
        update_pen = function(eps, kappa){
          p_env <- parent.env(environment())
          if(!missing(eps))
            assign("eps"  , eps  , p_env)
          if(!missing(kappa))
            assign("kappa", kappa, p_env)
          invisible(with(p_env, c(eps = eps, kappa = kappa)))
        },
        # function to get parameters
        get_params = function(x)
          x[-idx_va],
        control = list(maxit = 1000L)
      ), out)
    })
  }

  # find initial VA params
  func <- get_snva_out(theta_VA)
  coefs_start <- c(params$b, params$theta)
  do_drop <- seq_along(coefs_start)
  get_x <- function(x)
    c(coefs_start, x)

  fn <- function(x)
    func$fn(get_x(x))
  gr <- function(x)
    func$gr(get_x(x))[-do_drop]

  opt_out <- opt_func(theta_VA, fn = fn, gr = gr)
  snva_out <- get_snva_out(opt_out$par)

  } else
    snva_out <- NULL

  structure(
    list(laplace = laplace_out, gva = gva_out, snva = snva_out, y = y,
         event = event, X = X, XD = XD, Z = Z, grp = grp, terms = list(
           X = mt_X, Z = mt_Z, baseline = mt_b), cl = match.call(),
         link = link, opt_func = opt_func),
    class = "GSM_ADFun")
}
