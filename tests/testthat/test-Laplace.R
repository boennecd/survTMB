context("testing Laplace approximation")

get_par_val_eortc <- function(x){
  x <- x$optim
  x$par <- head(x$par, 6)
  x[c("par", "value")]
}

get_func_eortc <- function(link, n_threads)
  make_mgsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
    df = 3L, data = eortc, link = link, do_setup = "Laplace",
    n_threads = n_threads)

for(link in c("PH", "PO", "probit"))
  for(n_threads in 1:2)
    test_that(sprintf("Laplace gives previous results (%s, %d)", sQuote(link),
                      n_threads), {
      func <- get_func_eortc(link, n_threads)
      expect_s3_class(func, "MGSM_ADFun")
      expect_setequal(names(func), .MGSM_ADFun_members)
      expect_known_output(
        func, sprintf(file.path(test_res_dir, "Laplace-func-%s.txt"), link),
        print = TRUE)

      eps <- .Machine$double.eps^(3/5)
      res <- fit_mgsm(func, "Laplace", control = list(reltol = eps))
      expect_s3_class(res, "MGSM_ADFit")
      expect_setequal(names(res), .MGSM_fit_members)

      expect_known_value(
        get_par_val_eortc(res),
        sprintf(file.path(test_res_dir, "Laplace-%s.RDS"), link),
        tolerance = sqrt(eps))
      expect_known_output(
        res, sprintf(file.path(test_res_dir, "Laplace-%s.txt"), link),
        print = TRUE)

      free_laplace(func)
    })
