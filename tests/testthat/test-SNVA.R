context("testing SNVA")

get_par_val_eortc <- function(x){
  x <- x$optim
  x$par <- head(x$par, 6)
  x[c("par", "value")]
}

get_func_eortc <- function(link, n_threads, param_type, dense_hess = FALSE,
                           sparse_hess = FALSE)
  make_mgsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
    df = 3L, data = eortc, link = link, do_setup = "SNVA",
    n_threads = n_threads, param_type = param_type, dense_hess = dense_hess,
    sparse_hess = sparse_hess)

for(link in c("PH", "PO", "probit"))
  for(n_threads in 1:2)
    for(param_type in c("DP", "CP_trans"))
      for(use_own in c(TRUE, FALSE))
        test_that(sprintf("SNVA gives previous results (%s, %s, %d, %d)",
                          sQuote(link), sQuote(param_type), n_threads,
                          use_own), {
          old_val <- survTMB:::.get_use_own_VA_method()
          on.exit(survTMB:::.set_use_own_VA_method(old_val))
          survTMB:::.set_use_own_VA_method(use_own)

          func <- get_func_eortc(link, n_threads, param_type)
          expect_s3_class(func, "MGSM_ADFun")
          expect_setequal(names(func), .MGSM_ADFun_members)

          eps <- .Machine$double.eps^(3/5)
          res <- fit_mgsm(func, "SNVA", control = list(reltol = eps))
          expect_s3_class(res, "MGSM_ADFit")
          expect_setequal(names(res), .MGSM_fit_members)

          expect_known_value(
            get_par_val_eortc(res),
            sprintf(file.path(test_res_dir, "SNVA-%s-%s.RDS"),
                    link, param_type),
            tolerance = sqrt(eps), check.attributes = use_own)
          expect_known_output(
            res, sprintf(file.path(test_res_dir, "SNVA-%s-%s.txt"),
                         link, param_type),
            print = TRUE)

          clear_cppad_mem(2L)
        })

for(link in c("PH", "PO", "probit"))
  for(param_type in c("DP", "CP_trans"))
    test_that(
      sprintf("SNVA gives correct Hessian estimates (%s, %s)",
              sQuote(link), sQuote(param_type)), {
        skip_if_not_installed("numDeriv")

        old_val <- survTMB:::.get_use_own_VA_method()
        on.exit(survTMB:::.set_use_own_VA_method(old_val))

        survTMB:::.set_use_own_VA_method(TRUE)
        my_func <- get_func_eortc(
          link = link, 2L, param_type, dense_hess = TRUE)
        sp_func <- get_func_eortc(
          link = link, 2L, param_type, sparse_hess = TRUE)
        expect_known_output(
          sp_func, sprintf(file.path(test_res_dir,
                                     "SNVA-func-sparse-%s-%s.txt"),
                           link, param_type),
          print = TRUE)
        survTMB:::.set_use_own_VA_method(FALSE)
        tm_func <- get_func_eortc(
          link = link, 2L, param_type)

        my_fit <- fit_mgsm(my_func, "SNVA")

        par <- with(my_fit, c(params, va_params))
        my_hes <- my_func$snva$he(par)
        sp_hes <- sp_func$snva$he_sp(par)
        expect_s4_class(sp_hes, "dsCMatrix")
        tm_hes <- tm_func$snva$he(par)

        eps <- .Machine$double.eps^(3/5)
        nu_hes <- numDeriv::jacobian(
          my_func$snva$gr, par, method.args = list(eps = eps))

        expect_equal(my_hes, tm_hes)
        expect_equal(my_hes, as.matrix(sp_hes), check.attributes = FALSE)
        expect_equal(my_hes, nu_hes, tolerance = sqrt(eps))

        clear_cppad_mem(2L)
      })

