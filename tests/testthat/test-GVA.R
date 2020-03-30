context("testing GVA")

get_par_val_eortc <- function(x){
  x$par <- head(x$par, 6)
  x[c("par", "value")]
}

get_func_eortc <- function(link, n_threads)
  make_gsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
    df = 3L, data = eortc, link = link, do_setup = "GVA",
    n_threads = n_threads)

for(link in c("PH", "PO", "probit"))
  for(n_threads in 1:2)
    test_that(sprintf("GVA gives previous results (%s, %d)", sQuote(link),
                      n_threads), {
      func <- get_func_eortc(link, n_threads)

      eps <- .Machine$double.eps^(3/5)
      func$gva$control <- eps
      res <- do.call(optim, func$gva)

      expect_known_value(
        get_par_val_eortc(res), sprintf("test-res/GVA-%s.RDS", link),
        tolerance = sqrt(eps))
    })
