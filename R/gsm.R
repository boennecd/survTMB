gsm_get_XD <- function(time_var, mt_X, data){
  stopifnot(is.symbol(time_var))

  #####
  # approximate derivatives w/ finite difference
  dt <- .Machine$double.eps^(1/3)
  dat_m1 <- dat_p1 <- data
  x <- data[[deparse(time_var)]]
  dat_p1[[time_var]] <- x * exp( dt)
  dat_m1[[time_var]] <- x * exp(-dt)
  (model.matrix(mt_X, dat_p1) - model.matrix(mt_X, dat_m1)) / 2 / dt / x
}

# fits a GSM
#' @importFrom stats ecdf
gsm <- function(formula, data, df, tformula = NULL, link, n_threads,
                do_fit, opt_func = .opt_default){
  # checks
  stopifnot(
    inherits(formula, "formula"),
    is.data.frame(data),
    missing(df) || (is.integer(df) && df > 0 && length(df) == 1L),
    is.null(tformula) || inherits(tformula, "formula"),
    !(missing(df) && is.null(tformula)),
    is.integer(n_threads) && n_threads > 0 && length(n_threads) == 1L,
    link %in% c("PH", "PO", "probit"),
    is.logical(do_fit), length(do_fit) == 1L)

  #####
  # get design matrices and outcome
  mf_Z <- model.frame(formula, data = data)
  mt_Z <- terms(mf_Z)
  Z <- model.matrix(mt_Z, mf_Z)

  y <- model.response(mf_Z)
  stopifnot(inherits(y, "Surv"), isTRUE(attr(y, "type") == "right"))
  event <- y[, 2]

  #####
  # get time-varying baseline
  time_var <- formula[[2L]][[2L]]
  if(is.null(tformula))
    tformula <- eval(
      bquote(~ nsx(log(.(time_var)), df = .(df), intercept = FALSE) - 1))
  mt_X <- terms(model.frame(tformula, data = data[event > 0, ]))
  X <- model.matrix(mt_X, data)

  XD <- gsm_get_XD(time_var = time_var, mt_X = mt_X, data = data)

  out <- list(Z = Z, X = X, XD = XD, y = y, mt_Z = mt_Z, mt_X = mt_X)
  if(!do_fit)
    return(out)

  out$fit <- gsm_fit(X, XD, Z, y, link, n_threads, opt_func)
  out
}

# fits a GSM
gsm_fit <- function(X, XD, Z, y, link, n_threads, opt_func = .opt_default){
  # checks
  n <- NROW(y)
  stopifnot(NROW(X) == n, is.matrix(X),
            NROW(XD) == n, is.matrix(XD),
            NROW(Z) == n, is.matrix(Z),
            inherits(y, "Surv"), isTRUE(attr(y, "type") == "right"),
            is.integer(n_threads), length(n_threads) == 1L, n_threads > 0L)
  event <- y[, 2]

  # get starting values
  start_coef <- local({
    keep <- event > 0
    y_pass <- y[, 1][keep]
    fit <- ecdf(y_pass)
    S_hat <- 1 - fit(y_pass)
    n <- length(S_hat)
    S_hat <- pmax(.25 / n, pmin(S_hat, 1 - .25 / n))

    link_hat <- if(link == "PH")
      pmax(log(.Machine$double.eps) / 4, log(-log(S_hat)))
    else if(link == "PO")
      log((1 - S_hat) / S_hat)
    else if(link == "probit")
      -qnorm(S_hat)
    else
      stop(sprintf("%s not implemented", sQuote(link)))

    sfit <-
      lm.fit(x = cbind(X[keep, , drop = FALSE], Z[keep, , drop = FALSE]),
             y = link_hat)

    list(beta  = sfit$coefficients[ seq_len(NCOL(X))],
         gamma = sfit$coefficients[-seq_len(NCOL(X))])
  })

  # get object for ML estimation and perform the estimation
  opt_obj <- get_gsm_pointer(X = t(X), XD = t(XD), Z = t(Z),
                             y = event, eps = 1e-16, kappa = 1e8,
                             link = link, n_threads = n_threads)

  is_beta <- seq_along(start_coef$beta)
  is_gamma <- with(start_coef, seq_along(gamma) + length(beta))
  fn <- function(x){
    b <- x[is_beta]
    g <- x[is_gamma]
    -gsm_eval_ll(ptr = opt_obj, beta = b, gamma = g)
  }
  gr <- function(x){
    b <- x[is_beta]
    g <- x[is_gamma]
    -gsm_eval_grad(ptr = opt_obj, beta = b, gamma = g)
  }

  par <- with(start_coef, c(beta, gamma))
  opt_out <- opt_func(par, fn = fn, gr = gr)

  list(beta  = opt_out$par[is_beta],
       gamma = opt_out$par[is_gamma],
       optim = opt_out, mlogli = fn, grad = gr)
}
