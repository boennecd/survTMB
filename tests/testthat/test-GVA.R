context("testing Laplace approximation")

get_par_val_eortc <- function(x){
  x$par <- head(x$par, 6)
  x[c("par", "value")]
}

get_func_eortc <- function(link)
  make_gsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
    df = 3L, data = eortc, link = link, do_setup = "GVA")

test_that("GVA gives previous results (PH)", {
  func <- get_func_eortc("PH")
  res <- do.call(optim, func$gva)
  get_par_val_eortc(res)

  expect_equal(
    get_par_val_eortc(res),
    list(par = c(`(Intercept)` = -8.74960475750348,
                 trt = 0.627588644748086, `nsx(log(y), df = 3, intercept = FALSE)1` = 6.17807627728185,
                 `nsx(log(y), df = 3, intercept = FALSE)2` = 13.1246066637928,
                 `nsx(log(y), df = 3, intercept = FALSE)3` = 4.77284091742587,
                 theta = -1.34470293032808), value = 3629.76651939251))
})

test_that("GVA gives previous results (PO)", {
  func <- get_func_eortc("PO")
  res <- do.call(optim, func$gva)
  get_par_val_eortc(res)

  expect_equal(
    get_par_val_eortc(res),
    list(par = c(`(Intercept)` = -8.91480725104131,
                 trt = 0.877105340100762, `nsx(log(y), df = 3, intercept = FALSE)1` = 6.42859409671986,
                 `nsx(log(y), df = 3, intercept = FALSE)2` = 13.4490486270861,
                 `nsx(log(y), df = 3, intercept = FALSE)3` = 5.41894268592975,
                 theta = -1.02233305567567), value = 3629.28331674209))
})

test_that("GVA gives previous results (probit)", {
  func <- get_func_eortc("probit")
  res <- do.call(optim, func$gva)
  get_par_val_eortc(res)

  expect_equal(
    get_par_val_eortc(res),
    list(par = c(`(Intercept)` = -4.05421750814269,
                 trt = 0.524883805380364, `nsx(log(y), df = 3, intercept = FALSE)1` = 2.9017701774533,
                 `nsx(log(y), df = 3, intercept = FALSE)2` = 5.55262972082832,
                 `nsx(log(y), df = 3, intercept = FALSE)3` = 2.83841845310392,
                 theta = -1.55358670744936), value = 3629.13051119782))
})
