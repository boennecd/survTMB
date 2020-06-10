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
    cox_fit <- coxph(y ~ Z - 1, model = TRUE)
    S_hat <- Shat(cox_fit)
    link_hat <- if(link == "PH")
      pmax(log(.Machine$double.eps) / 4, log(-log(S_hat)))
    else if(link == "PO")
      log((1 - S_hat) / S_hat)
    else if(link == "probit")
      -qnorm(S_hat)
    else
      stop(sprintf("%s not implemented", sQuote(link)))

    keep <- event > 0
    sfit <-
      lm.fit(x = cbind(X[keep, , drop = FALSE], Z[keep, , drop = FALSE]),
             y = link_hat[keep])

    list(beta  = sfit$coefficients[ seq_len(NCOL(X))],
         gamma = sfit$coefficients[-seq_len(NCOL(X))])
  })

  # get object for ML estimation and perform the estimation
  opt_obj <- get_gsm_pointer(X = t(X), XD = t(XD), Z = t(Z),
                             y = event, eps = 1e-16, kappa = 1e8,
                             link = link, n_threads = n_threads)

  is_beta <- seq_along(start_coef$beta)
  fn <- function(x){
    b <- x[ is_beta]
    g <- x[-is_beta]
    -gsm_eval_ll(ptr = opt_obj, beta = b, gamma = g)
  }
  gr <- function(x){
    b <- x[ is_beta]
    g <- x[-is_beta]
    -gsm_eval_grad(ptr = opt_obj, beta = b, gamma = g)
  }

  par <- with(start_coef, c(beta, gamma))
  opt_out <- opt_func(par, fn = fn, gr = gr)

  list(beta  = opt_out$par[is_beta],
       gamma = opt_out$par[-is_beta],
       optim = par)
}
