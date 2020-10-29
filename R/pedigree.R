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
#' @param args_gva_opt arguments to pass to \code{\link{optim_pedigree_psqn}}
#' when optimizing a GVA to find initial values for the SNVA.
#' @inheritParams make_mgsm_ADFun
#'
#' @export
make_pedigree_ADFun <- function(
  c_data, formula, tformula, n_nodes = 20L, n_threads = 1L,
  sparse_hess = FALSE, link = c("PH", "PO", "probit"),
  opt_func = .opt_default, skew_start = -.0001, omega = NULL,
  beta = NULL, sds = NULL, trace = FALSE, kappa = .MGSM_default_kappa,
  dtformula = NULL, method = c("SNVA", "GVA"), args_gva_opt = list()){
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
    method %in% c("SNVA", "GVA"),
    is.list(args_gva_opt))
  eval(bquote(stopifnot(
    .(-.skew_boundary) < skew_start && skew_start < .(.skew_boundary))))

  # first we get the GVA. Then we either return it or use it to find the
  # starting values for the SNVA

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
    sds <- rep(sqrt(.2 / length(c_data[[1L]]$cor_mats)),
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

    sig_use <- cov_to_theta(sig_use)
    out <- c(setNames(mu, paste0("mu", seq_along(mu))),
             sig_use)
    setNames(out, paste0(sprintf("g%d:", g), names(out)))
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

  # setup cache
  if(trace)
    cat("Creating ADFun for the GVA...\n")
  setup_atomic_cache(n_nodes = n_nodes, type = .gva_char, link = link)
  adfun <- get_pedigree_funcs(
    data = c_data, n_nodes = n_nodes, link = link, omega = omega,
    beta = beta, log_sds = log(sds), va_par = va_par, method = "GVA",
    eps = .MGSM_defaul_eps, kappa = kappa, n_threads = n_threads)
  par <- c(omega, beta, log(sds), va_par)

  #####
  # find variational parameters
  if(trace)
    cat("Finding starting values for the variational parameters in the GVA...\n")
  par <- psqn_optim_pedigree_private(
    val = par, ptr = adfun, rel_eps = sqrt(.Machine$double.eps),
    max_it = 1000L, n_threads = n_threads, c1 = 1e-4, c2 = .9,
    method = "GVA")
  if(trace)
    cat(sprintf(
      "The lower bound at the starting values with the GVA is: %.3f\n",
      -attr(par, "value")))
  par <- c(par)

  out <- list(
    par = par,
    fn = function(x, ...){
      eval_psqn_pedigree(x, ptr = adfun, n_threads = n_threads,
                         method = "GVA")
    },
    gr = function(x, ...){
      grad_psqn_pedigree(x, ptr = adfun, n_threads = n_threads,
                         method = "GVA")
    },
    he = function(x, ...){
      stop("he not implemented")
    },
    get_params = function(x)
      stop("get_params not implemented"),
    opt_func = opt_func,
    sparse_hess = sparse_hess,
    cl = match.call(),
    method = method
    # TODO: save terms
  )

  if(method == "GVA"){
    #####
    # return the GVA output
    rm(list = setdiff(ls(), c("out", "adfun", "n_threads", "method")))
    return(out)
  }

  #####
  # optimize the GVA and setup the starting values
  if(trace)
    cat("Optimizing the GVA to get starting values...\n")
  gva_defaults <- list(
    max_cg = max(min(100L, length(par)), as.integer(length(par) * .001)))
  opt_args <- c(args_gva_opt,
                gva_defaults[setdiff(names(gva_defaults),
                                     names(args_gva_opt))])
  opt_args$object <- out
  gva_res <- do.call(optim_pedigree_psqn, opt_args)

  if(trace)
    cat(sprintf("Maximum lower bound with the GVA is %.2f...\n",
                -gva_res$value))

  # get the model parameters and the VA parameters. First the global
  # parameters
  m_par    <- gva_res$params
  va_start <- gva_res$va_params

  omega <- head(m_par, length(omega))
  m_par <- m_par[-seq_along(omega)]
  beta <- head(m_par, length(beta))
  sds <- exp(tail(m_par, -length(beta)))

  idx_start <- 1L
  g <- 0L
  va_start <- sapply(c_data, function(c_dat){
    n_obs <- length(c_dat$event)
    na_va <- n_obs + (n_obs * (n_obs + 1L)) / 2L
    va_par <- va_start[idx_start - 1L + 1:na_va]
    idx_start <<- idx_start + na_va
    g         <<- g + 1L

    mu <-               head(va_par,  n_obs)
    vcv <- theta_to_cov(tail(va_par, -n_obs))
    out <- cp_to_dp(mu, vcv, rep(skew_start, n_obs))

    nams <- mgsm_get_snva_names(n_rng = n_obs, n_grp = 1L,
                                param_type = "DP")
    nams <- gsub("(^g\\d+)(.+)", paste0("g", g, "\\2"), nams)
    setNames(c(out$xi, cov_to_theta(out$Psi), out$alpha),
             nams)
  })
  va_start <- unlist(va_start)

  # setup AD function
  if(trace)
    cat("Creating ADFun for the SNVA...\n")
  setup_atomic_cache(n_nodes = n_nodes, type = .snva_char, link = link)
  adfun <- get_pedigree_funcs(
    data = c_data, n_nodes = n_nodes, link = link, omega = omega,
    beta = beta, log_sds = log(sds), va_par = va_start, method = "SNVA",
    eps = .MGSM_defaul_eps, kappa = kappa, n_threads = n_threads)
  par <- c(omega, beta, log(sds), va_start)

  # improve starting values
  if(trace)
    cat("Finding starting values for the variational parameters in the SNVA...\n")
  par <- psqn_optim_pedigree_private(
    val = par, ptr = adfun, rel_eps = sqrt(.Machine$double.eps),
    max_it = 1000L, n_threads = n_threads, c1 = 1e-4, c2 = .9,
    method = "SNVA")
  if(trace)
    cat(sprintf(
      "The lower bound at the starting values for the SNVA is: %.3f\n",
      -attr(par, "value")))
  par <- c(par)

  # return
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
    cl = match.call(),
    method = method
    # TODO: save terms
  )
  rm(list = setdiff(ls(), c("out", "adfun", "n_threads", "method")))
  out
}

#' Optimize Mixed Generalized Survival Model for Pedigree Data Using the
#' psqn Method
#'
#' @param object object from \code{\link{make_pedigree_ADFun}}.
#' @inheritParams optim_mgsm_psqn
#'
#' @seealso
#' \code{\link{make_pedigree_ADFun}}
#'
#' @importFrom stats setNames
#' @export
optim_pedigree_psqn <- function(object, par = NULL,
                                rel_eps = sqrt(.Machine$double.eps),
                                max_it = 1000L, c1 = 1e-4, c2 = .1,
                                trace = 0L, use_bfgs = TRUE,
                                cg_tol = .5, strong_wolfe = TRUE,
                                max_cg = NULL, pre_method = 1L){
  #####
  # setup and checks
  if(is.null(par))
    par <- object$par
  if(is.null(max_cg))
    max_cg <- length(par)
  env <- environment(object$fn)
  stopifnot(
    is.numeric(par), all(is.finite(par)),
    is.numeric(rel_eps), is.finite(rel_eps), length(rel_eps) == 1L,
    is.numeric(c1), is.finite(c1), length(c1) == 1L,
    is.numeric(c2), is.finite(c2), length(c2) == 1L,
    is.numeric(cg_tol), is.finite(cg_tol), length(cg_tol) == 1L,
    is.integer(max_it), is.finite(max_it), length(max_it) == 1L,
    is.integer(trace), is.finite(trace), length(trace) == 1L,
    is.logical(strong_wolfe), is.finite(strong_wolfe),
    length(strong_wolfe) == 1L,
    is.logical(use_bfgs), is.finite(use_bfgs),
    length(use_bfgs) == 1L,
    is.integer(max_cg), is.finite(max_cg), length(max_cg) == 1L,
    is.integer(pre_method), length(pre_method) == 1L, pre_method %in% 0:1)

  #####
  # optimize
  res <- psqn_optim_pedigree(
    par, ptr = env$adfun, rel_eps = rel_eps, max_it = max_it,
    n_threads = env$n_threads, c1 = c1, c2 = c2, trace = trace,
    use_bfgs = use_bfgs, cg_tol = cg_tol, strong_wolfe = strong_wolfe,
    method = object$method, max_cg = max_cg, pre_method)

  # sort out variational parameters from model parameters
  par <- setNames(res$par, names(object$par))
  res$par <- par

  is_val <- grepl("^g\\d+:", names(par))
  res$params    <- res$par[!is_val]
  res$va_params <- res$par[ is_val]

  res
}
