#' Construct Objective Functions with Derivatives for a Mixed Generalized
#' Survival Model for Pedigree Data
#'
#' @param formula two-sided \code{\link{formula}} where the left-hand side is a
#'                \code{\link{Surv}} object and the right-hand side is the
#'                fixed effects. The formula is passed on each cluster so
#'                the \code{\link{model.matrix}} cannot depend on the
#'                data set (e.g. you cannot use \code{\link{ns}} without
#'                specifying the knots).
#' @param tformula \code{\link{formula}} with baseline survival function.
#'                 The time variable must be the same
#'                 symbol as used in the left-hand-side of \code{formula}.
#'                 The formula is passed on each cluster so
#'                 the \code{\link{model.matrix}} cannot depend on the
#'                 data set (e.g. you cannot use \code{\link{ns}} without
#'                 specifying the knots).
#' @param c_data \code{\link{list}} with cluster data.
#' @param n_nodes integer with the number of nodes to use in (adaptive)
#'                Gauss-Hermite quadrature and Gauss–Legendre quadrature.
#' @param omega starting value for baseline survival function's parameters.
#' @param beta starting values for fixed effects coefficients.
#' @param sds starting values for scale matrix scales.
#' @param trace logical for whether to print tracing information.
#' @param kappa numeric scalar with the penalty in the relaxed problem
#' ensuring the monotonicity of the survival curve.
#' @param method character vector indicating which approximation to use.
#' @inheritParams make_mgsm_ADFun
#'
#' @export
make_pedigree_ADFun <- function(
  c_data, formula, tformula, n_nodes = 20L, n_threads = 1L,
  sparse_hess = FALSE, link = c("PH", "PO", "probit"),
  opt_func = .opt_default, skew_start = -.0001, omega = NULL,
  beta = NULL, sds = NULL, trace = FALSE, kappa = .MGSM_default_kappa,
  dtformula = NULL, method = c("SNVA", "GVA")){
  # checks
  link <- link[1L]
  method <- method[1L]
  stopifnot(
    is.list(c_data), length(c_data) > 0L,
    is.integer(n_nodes), length(n_nodes) == 1L, n_nodes > 0L,
    is.integer(n_threads), length(n_threads) == 1L, n_threads > 0L,
    is.logical(sparse_hess), length(sparse_hess) == 1L,
    is.character(link), link %in% c(c("PH", "PO", "probit")),
    is.null(omega) || (
      is.numeric(omega) && is.vector(omega) && all(is.finite(omega))),
    is.null(beta) || (
      is.numeric(beta) && is.vector(beta) && all(is.finite(beta))),
    is.null(sds) || (
      is.numeric(sds) && is.vector(sds) && all(sds > 0)),
    is.logical(trace), length(trace) == 1L,
    is.double(kappa), length(kappa) == 1L, kappa >= 0,
    method %in% c("SNVA", "GVA"))
  eval(bquote(stopifnot(
    .(-.skew_boundary) < skew_start && skew_start < .(.skew_boundary))))

  # get design matrices etc.
  c_data <- lapply(c_data, function(x){
    gsm_dat <- gsm(formula = formula, data = x$data, tformula = tformula,
                   link = link, n_threads = n_threads, do_fit = FALSE,
                   opt_func = opt_func, dtformula = dtformula)
    y <- gsm_dat$y
    stopifnot(inherits(y, "Surv"), isTRUE(attr(y, "type") == "right"))

    list(cor_mats = x$cor_mats, X = t(gsm_dat$X), XD = t(gsm_dat$XD),
         Z = t(gsm_dat$Z), event = y[, 2], y_surv = y)
  })

  #####
  # get starting values
  miss_model_params <- is.null(beta) || is.null(omega) || is.null(sds)
  if(is.null(beta) || is.null(omega)){
    if(trace)
      cat("Finding starting values for the fixed effects...\n")

    fix_start <- local({
      X  <- t(do.call(cbind, lapply(c_data, `[[`, "X")))
      XD <- t(do.call(cbind, lapply(c_data, `[[`, "XD")))
      Z  <- t(do.call(cbind, lapply(c_data, `[[`, "Z")))
      y  <- do.call(rbind, lapply(c_data, `[[`, "y_surv"))

      class(y) <- "Surv"
      attr(y, "type") <- "right"

      fit <- gsm_fit(X = X, XD = XD, Z = Z, y = y, link = link,
                     n_threads = n_threads, opt_func = opt_func,
                     offset_eta = numeric(), offset_etaD = numeric())
      if(trace)
        cat(sprintf(
          "Maximum log-likelihood without random effects is: %.3f\n",
          -fit$optim$value))

      list(omega = fit$beta, beta = fit$gamma)
    })

    omega <- fix_start$omega
    beta <- fix_start$beta
  }

  if(is.null(sds))
    sds <- rep(sqrt(.5 / length(c_data[[1L]]$cor_mats)),
               length(c_data[[1L]]$cor_mats))

  #####
  # checks
  for(i in seq_along(c_data)){
    n_obs <- c_data[[i]]$n_obs
    expr <- substitute({
      stopifnot(
        is.numeric(x$event), length(x$event) == n_obs,
        all(is.finite(x$event)),
        is.matrix(x$X), NCOL(x$X) == n_obs, NROW(x$X) == length(omega),
        all(is.finite(x$X)),
        is.matrix(x$XD), NCOL(x$XD) == n_obs, NROW(x$XD) == length(omega),
        all(is.finite(x$XD)),
        is.matrix(x$Z), NCOL(x$Z) == n_obs, NROW(x$Z) == length(beta),
        all(is.finite(x$Z)),
        is.list(x$cor_mats),
        length(x$cor_mats) == length(sds))
      for(j in seq_along(sds))
        stopifnot(is.matrix(x$cor_mats[[j]]),
                  NROW(x$cor_mats[[j]]) == n_obs,
                  NCOL(x$cor_mats[[j]]) == n_obs)
    }, list(x = bquote(c_data[[.(i)]])))
    eval(expr, environment())
  }

  #####
  # setup variational parameters
  g <- 0L
  va_par <- sapply(c_data, function(c_dat){
    g <<- g + 1L

    n_obs <- length(c_dat$event)
    sigma <- matrix(0., n_obs, n_obs)
    for(i in seq_along(sds))
      sigma <- sigma + sds[i]^2 * c_dat$cor_mats[[i]]

    offset_eta  <- drop(omega %*% c_dat$X + beta %*% c_dat$Z)
    offset_etaD <- drop(omega %*% c_dat$XD)

    X <- matrix(nrow = 0, ncol = n_obs)
    # TODO: this is very inefficient
    opt_obj <- get_gsm_pointer(
      X = X, XD = X, Z = diag(n_obs), y = c_dat$event, eps = 1e-16,
      kappa = 1e8, link = link, n_threads = 1L, offset_eta = offset_eta,
      offset_etaD = offset_etaD)

    par <- numeric(n_obs)
    sig_inv <- solve(sigma)
    chol_sig_inv <- chol(sig_inv)

    fn <- function(x, ...)
      -gsm_eval_ll(ptr = opt_obj, beta = numeric(), gamma = x) +
      sum((chol_sig_inv %*% x)^2) / 2
    gr <- function(x, ...)
      -gsm_eval_grad(ptr = opt_obj, beta = numeric(), gamma = x) +
      drop(sig_inv %*% x)
    he <- function(x, ...)
      -gsm_eval_hess(ptr = opt_obj, beta = numeric(), gamma = x) + sig_inv

    opt_ret <- opt_func(par, fn, gr)
    if(opt_ret$ok){
      mu <- opt_ret$par
      sig_use <- solve(he(mu))

    } else {
      mu <- par
      sig_use <- sigma

    }

    sig_use <- .rescale_cov(sig_use)

    if(method == "GVA"){
      sig_use <- cov_to_theta(sig_use)
      out <- c(setNames(mu, paste0("mu", seq_along(mu))),
               sig_use)
      return(setNames(out, paste0(sprintf("g%d:", g), names(out))))
    }

    out <- cp_to_dp(
      mu = mu, Sigma = sig_use, gamma = rep(skew_start, n_obs))

    xi <- out$xi
    names(xi) <- paste0("xi", seq_along(xi))
    Psi <- cov_to_theta(out$Psi)
    alpha <- out$alpha
    names(xi) <- paste0("xi", seq_along(xi))
    names(alpha) <- paste0("alpha", seq_along(alpha))

    out <- c(xi, Psi, alpha)
    names(out) <- paste0(sprintf("g%d:", g), names(out))
    out
  }, simplify = FALSE)
  va_par <- unlist(va_par)

  #####
  # get AD funcs
  if(length(omega) > 0){
    rnames <- rownames(c_data[[1L]]$X)
    if(!is.null(rnames))
      names(omega) <- paste0("omega:", rnames)
    else
      names(omega) <- paste0("omega:", seq_along(omega))
  }
  if(length(beta) > 0){
    rnames <- rownames(c_data[[1L]]$Z)
    if(!is.null(rnames))
      names(beta) <- paste0("beta:", rnames)
    else
      names(beta) <- paste0("beta:", seq_along(beta))
  }
  if(length(sds) > 0)
    names(sds) = paste0("log_sds", seq_along(sds))

  if(trace)
    cat("Creating ADFun...\n")

  # setup cache
  setup_atomic_cache(n_nodes = n_nodes, type = .snva_char, link = link)
  adfun <- get_pedigree_funcs(
    data = c_data, n_nodes = n_nodes, link = link, omega = omega,
    beta = beta, log_sds = log(sds), va_par = va_par, method = method,
    eps = .MGSM_defaul_eps, kappa = kappa, n_threads = n_threads)
  par <- c(omega, beta, log(sds), va_par)

  #####
  # find variational parameters
  if(trace)
    cat("Finding starting values for the variational parameters...\n")
  par <- psqn_optim_pedigree_private(
    val = par, ptr = adfun, rel_eps = sqrt(.Machine$double.eps),
    max_it = 1000L, n_threads = n_threads, c1 = 1e-4, c2 = .9,
    method = method)
  if(trace)
    cat(sprintf(
      "The lower bound at the starting values is: %.3f\n",
      -attr(par, "value")))
  par <- c(par)

  #####
  # create output list
  out <- list(
    par = par,
    fn = function(x, ...){
      eval_psqn_pedigree(x, ptr = adfun, n_threads = n_threads,
                         method = method)
    },
    gr = function(x, ...){
      grad_psqn_pedigree(x, ptr = adfun, n_threads = n_threads,
                         method = method)
    },
    he = function(x, ...){
      stop("he not implemented")
    },
    get_params = function(x)
      stop("get_params not implemented"),
    opt_func = opt_func,
    sparse_hess = sparse_hess,
    cl = match.call()
    # TODO: save terms
  )
  rm(list = setdiff(ls(), c("out", "adfun", "n_threads", "method")))
  out
}
