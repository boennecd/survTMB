context("testing GVA")

get_par_val_eortc <- function(x){
  x <- x$optim
  x$par <- head(x$par, 6)
  x[c("par", "value")]
}

get_func_eortc <- function(link, n_threads, dense_hess = FALSE,
                           sparse_hess = FALSE)
  make_mgsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
    df = 3L, data = eortc, link = link, do_setup = "GVA",
    n_threads = n_threads, dense_hess = dense_hess,
    sparse_hess = sparse_hess)

for(link in c("PH", "PO", "probit"))
  for(n_threads in 1:2)
    for(use_own in c(TRUE, FALSE))
      test_that(sprintf("GVA gives previous results (%s, %d, %d)",
                        sQuote(link), n_threads, use_own), {
        old_val <- survTMB:::.get_use_own_VA_method()
        on.exit(survTMB:::.set_use_own_VA_method(old_val))
        survTMB:::.set_use_own_VA_method(use_own)

        func <- get_func_eortc(link, n_threads)
        expect_s3_class(func, "MGSM_ADFun")
        expect_setequal(names(func), .MGSM_ADFun_members)

        eps <- .Machine$double.eps^(3/5)
        res <- fit_mgsm(func, "GVA", control = list(reltol = eps))
        expect_s3_class(res, "MGSM_ADFit")
        expect_setequal(names(res), .MGSM_fit_members)

        expect_known_value(
          get_par_val_eortc(res),
          sprintf(file.path(test_res_dir, "GVA-%s.RDS"), link),
          tolerance = sqrt(eps))
        expect_known_output(
          res, sprintf(file.path(test_res_dir, "GVA-%s.txt"), link),
          print = TRUE)
      })

for(link in c("PH", "PO", "probit"))
  test_that(
    sprintf("GVA gives correct Hessian estimates (%s)", sQuote(link)), {
      skip_if_not_installed("numDeriv")

      old_val <- survTMB:::.get_use_own_VA_method()
      on.exit(survTMB:::.set_use_own_VA_method(old_val))

      survTMB:::.set_use_own_VA_method(TRUE)
      my_func <- get_func_eortc(link = link, 2L, dense_hess = TRUE)
      expect_known_output(
        my_func, sprintf(file.path(test_res_dir, "GVA-func-dense-%s.txt"),
                         link),
        print = TRUE)
      sp_func <- get_func_eortc(link = link, 2L, sparse_hess = TRUE)
      survTMB:::.set_use_own_VA_method(FALSE)
      tm_func <- get_func_eortc(link = link, 2L)

      my_fit <- fit_mgsm(my_func, "GVA")

      par <- with(my_fit, c(params, va_params))

      # check Hessian
      eps <- .Machine$double.eps^(3/5)
      my_hes <- my_func$gva$he(par)
      sp_hes <- sp_func$gva$he_sp(par)
      expect_s4_class(sp_hes, "dsCMatrix")
      tm_hes <- tm_func$gva$he(par)

      nu_hes <- numDeriv::jacobian(
        my_func$gva$gr, par, method.args = list(eps = eps))

      expect_equal(my_hes, tm_hes)
      expect_equal(my_hes, as.matrix(sp_hes), check.attributes = FALSE)
      expect_equal(my_hes, nu_hes, tolerance = sqrt(eps))
    })
