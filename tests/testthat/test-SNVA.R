context("testing SNVA")

get_par_val_eortc <- function(x){
  x <- x$optim
  x$par <- head(x$par, 6)
  x[c("par", "value")]
}

get_func_eortc <- function(link, n_threads, param_type)
  make_gsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
    df = 3L, data = eortc, link = link, do_setup = "SNVA",
    n_threads = n_threads, param_type = param_type)

for(link in c("PH", "PO", "probit"))
  for(n_threads in 1:2)
    for(param_type in c("DP", "CP_trans"))
      test_that(sprintf("SNVA gives previous results (%s, %s, %d)",
                        sQuote(link), sQuote(param_type), n_threads), {
        func <- get_func_eortc(link, n_threads, param_type)
        expect_s3_class(func, "GSM_ADFun")

        eps <- .Machine$double.eps^(3/5)
        res <- fit_mgsm(func, "SNVA", control = list(reltol = eps))
        expect_s3_class(res, "GSM_ADFit")

        expect_known_value(
          get_par_val_eortc(res), sprintf("test-res/SNVA-%s-%s.RDS",
                                          link, param_type),
          tolerance = sqrt(eps))
        expect_known_output(res, sprintf("test-res/SNVA-%s-%s.txt",
                                         link, param_type),
                            print = TRUE)
      })
