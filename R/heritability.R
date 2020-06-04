#' @export
make_heritability_ADFun <- function(
  c_data, n_nodes = 20L, n_threads = 1L, sparse_hess = FALSE,
  link = c("PH", "PO", "probit"), opt_func = .opt_default,
  skew_start = -.0001, omega, beta, sds){
  # checks
  link <- link[1L]
  stopifnot(
    is.list(c_data), length(c_data) > 0L,
    is.integer(n_nodes), length(n_nodes) == 1L, n_nodes > 0L,
    is.integer(n_threads), length(n_threads) == 1L, n_threads > 0L,
    is.logical(sparse_hess), length(sparse_hess) == 1L,
    is.character(link), link %in% c(c("PH", "PO", "probit")),
    is.numeric(omega), is.vector(omega), all(is.finite(omega)),
    is.numeric(beta), is.vector(beta), all(is.finite(beta)),
    is.numeric(sds), is.vector(sds), all(sds > 0),
    is.numeric(skew_start), all(is.finite(skew_start)))

  for(i in seq_along(c_data)){
    n_obs <- c_data[[i]]$n_obs
    expr <- substitute({
      stopifnot(
        is.numeric(x$y), length(x$y) == n_obs, all(is.finite(x$y)),
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
        stopifnot(is.matrix(x$cor_mats[[j]]), NROW(x$cor_mats[[j]]) == n_obs,
                  NCOL(x$cor_mats[[j]]) == n_obs)
    }, list(x = bquote(c_data[[.(i)]])))
    eval(expr, environment())
  }

  #####
  # setup variational parameters
  g <- 0L
  va_par <- sapply(c_data, function(c_dat){
    g <<- g + 1L

    n_obs <- length(c_dat$y)
    sigma <- matrix(0., n_obs, n_obs)
    for(i in seq_along(sds))
      sigma <- sigma + sds[i]^2 * c_dat$cor_mats[[i]]

    out <- cp_to_dp(mu = numeric(n_obs), Sigma = sigma,
                    gamma = rep(skew_start, n_obs))

    xi <- out$xi
    names(xi) <- paste0("xi", seq_along(xi))
    Psi <- cov_to_theta(out$Psi)
    alpha <- out$alpha
    names(alpha) <- paste0("alpha", seq_along(alpha))

    out <- c(out$xi, Psi, alpha)
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

  data <- list(n_threads = n_threads, sparse_hess = sparse_hess,
               n_nodes = n_nodes, link = link, c_data = c_data)
  parameters <- list(omega = omega, beta = beta, log_sds = log(sds),
                     va_par = va_par, eps = .MGSM_defaul_eps,
                     kappa = .MGSM_default_kappa)

  adfun <- get_herita_funcs(data = data, parameters = parameters)
  par <- c(parameters$omega, parameters$beta, parameters$log_sds,
           parameters$va_par)

  #####
  # find variational parameters
  is_va <- -seq_len(length(omega) + length(beta) + length(sds))
  new_vas <- local({
    fn_va <- function(x, ...){
      par[is_va] <- x
      herita_funcs_eval_lb(p = adfun, par)
    }
    gr_va <- function(x, ...){
      par[is_va] <- x
      herita_funcs_eval_grad(p = adfun, par)[is_va]
    }

    va_opt <-
      opt_func(par[is_va], fn_va, gr_va, control = list(maxit = 1000L))

    va_opt$par
  })
  par[is_va] <- new_vas

  #####
  # create output list
  out <- list(
    par = par,
    fn = function(x, ...){
      herita_funcs_eval_lb(p = adfun, x)
    },
    gr = function(x, ...){
      herita_funcs_eval_grad(p = adfun, x)
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

  out
}
