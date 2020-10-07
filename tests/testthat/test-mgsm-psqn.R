context("testing mgsm-psqn")

get_func_eortc_adfun <- function(link, n_threads)
  make_mgsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
    df = 3L, data = eortc, link = link, do_setup = "SNVA",
    n_threads = n_threads, dense_hess = FALSE,
    sparse_hess = FALSE, param_type = "DP")

get_psqn_obj <-  function(link, n_threads)
  make_mgsm_psqn_obj(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
    df = 3L, data = eortc, link = link, do_setup = "SNVA",
    n_threads = n_threads)

for(link in c("PH", "PO", "probit"))
  test_that(sprintf("psqn method gives something similar to adfun version with link '%s'", link), {
    obj <- get_psqn_obj(link, 2L)
    expect_s3_class(obj, "mgsm_psqn")

    org_vers <- get_func_eortc_adfun(link, n_threads = 2L)
    org_res <- fit_mgsm(org_vers, method = "SNVA")
    new_vers <- optim_mgsm_psqn(object = obj)

    expect_equal(names(org_res$optim$par), names(new_vers$par))
    expect_equal(org_res$optim$value, new_vers$value,
                 tolerance = 1e-4)

    # check the gradient and function value
    val <- obj$par
    expect_equal(org_vers$snva$fn(val),   obj$fn(val) )
    expect_equal(org_vers$snva$gr(val), c(obj$gr(val)))
  })
