.laplace_char <- "Laplace"
.gva_char     <- "GVA"
.snva_char    <- "SNVA"

#' @importFrom stats optim
.opt_default <- function(par, fn, gr, ...){
  out <- optim(par, fn = fn, gr = gr, method = "BFGS",
               control = list(reltol = .Machine$double.eps^(1/4),
                              maxit = 1000))
  stopifnot(out$convergence == 0L)
  out
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
  link = c("PH", "PO", "probit"),
  theta = NULL, beta = NULL, opt_func = .opt_default){
  link <- link[1]
  param_type <- param_type[1]
  stopifnot(
    is.integer(df), df > 0L, inherits(formula, "formula"),
    inherits(Z, "formula"), is.data.frame(data),
    !missing(cluster),
    all(do_setup %in% c(.laplace_char, .gva_char, .snva_char)),
    is.integer(n_nodes), length(n_nodes) == 1L, n_nodes > 0L,
    link %in% c("PH", "PO", "probit"),
    param_type %in% c("DP", "CP_trans", "CP"))

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
    link = link, grp_size = grp_size)

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

  list(laplace = laplace_out, gva = NULL, snva = NULL, y = y,
       event = event, X = X, XD = XD, Z = Z, grp = grp, terms = list(
         X = mt_X, Z = mt_Z, baseline = mt_b))
}
