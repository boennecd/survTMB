context("testing mgsm-psqn")

get_func_eortc_adfun <- function(link, n_threads, method)
  make_mgsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
    df = 3L, data = eortc, link = link, do_setup = method,
    n_threads = n_threads, dense_hess = FALSE,
    sparse_hess = FALSE, param_type = "DP")

get_psqn_obj <-  function(link, n_threads, method)
  make_mgsm_psqn_obj(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
    df = 3L, data = eortc, link = link, method = method,
    n_threads = n_threads)

for(m in c("SNVA", "GVA"))
for(link in c("PH", "PO", "probit"))
  test_that(sprintf("psqn method gives something similar to adfun version with link '%s' and method '%s'", link, m), {
    obj <- get_psqn_obj(link, 2L, method = m)
    expect_s3_class(obj, "mgsm_psqn")

    org_vers <- get_func_eortc_adfun(link, n_threads = 2L, method = m)
    org_res <- fit_mgsm(org_vers, method = m)
    new_vers <- optim_mgsm_psqn(object = obj)

    expect_equal(names(org_res$optim$par), names(new_vers$par))
    expect_equal(org_res$optim$value, new_vers$value,
                 tolerance = 1e-4)

    # check the gradient and function value
    val <- obj$par
    expect_equal(org_vers[[tolower(m)]]$fn(val),   obj$fn(val) )
    expect_equal(org_vers[[tolower(m)]]$gr(val), c(obj$gr(val)))
  })

test_that("make_mgsm_psqn_obj gives the same with the dtformula argument", {
  skip_if_not_installed("splines")
  library(splines)
  ks <- sort(c(rep(c(0, 3000), each = 4), 1500))
  f <- function(x)
    splineDesign(ks, x)
  df <- function(x)
    splineDesign(ks, x, derivs = 1)

  r1 <- make_mgsm_psqn_obj(
    Surv(y, uncens) ~ trt - 1, cluster = as.factor(center), Z = ~ 1,
    tformula = ~ f(y) - 1, dtformula = ~ df(y) - 1,
    df = 3L, data = eortc, link = "PH", method = "GVA", n_threads = 1L)
  r2 <- make_mgsm_psqn_obj(
    Surv(y, uncens) ~ trt - 1, cluster = as.factor(center), Z = ~ 1,
    tformula = ~ f(y) - 1,
    df = 3L, data = eortc, link = "PH", method = "GVA", n_threads = 1L)

  expect_equal(r1$par, r2$par)
})
