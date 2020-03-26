context("testing Laplace approximation")

get_par_val_eortc <- function(x){
  x$par <- head(x$par, 6)
  x[c("par", "value")]
}

get_func_eortc <- function(link)
  make_gsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
    df = 3L, data = eortc, link = link, do_setup = "Laplace")

test_that("Laplace gives previous results (PH)", {
  func <- get_func_eortc("PH")
  res <- do.call(optim, func$laplace)
  get_par_val_eortc(res)

  expect_equal(
    get_par_val_eortc(res),
    list(par = c(`(Intercept)` = -8.74959213009747,
                 trt = 0.62759066229636, `nsx(log(y), df = 3, intercept = FALSE)1` = 6.17807221835789,
                 `nsx(log(y), df = 3, intercept = FALSE)2` = 13.1245390254641,
                 `nsx(log(y), df = 3, intercept = FALSE)3` = 4.77285859376998,
                 theta = -1.34453684515036), value = 3629.76101490319))
})

test_that("Laplace gives previous results (PO)", {
  func <- get_func_eortc("PO")
  res <- do.call(optim, func$laplace)
  get_par_val_eortc(res)

  expect_equal(
    get_par_val_eortc(res),
    list(par = c(`(Intercept)` = -8.91505349480837,
                 trt = 0.877004876875972, `nsx(log(y), df = 3, intercept = FALSE)1` = 6.42872669657199,
                 `nsx(log(y), df = 3, intercept = FALSE)2` = 13.4497398630251,
                 `nsx(log(y), df = 3, intercept = FALSE)3` = 5.41887041787906,
                 theta = -1.02363322587053), value = 3629.29314256304))
})

test_that("Laplace gives previous results (probit)", {
  func <- get_func_eortc("probit")
  res <- do.call(optim, func$laplace)
  get_par_val_eortc(res)

  expect_equal(
    get_par_val_eortc(res),
    list(par = c(`(Intercept)` = -4.05448057657354,
                 trt = 0.524894660123169, `nsx(log(y), df = 3, intercept = FALSE)1` = 2.90168535341615,
                 `nsx(log(y), df = 3, intercept = FALSE)2` = 5.55253403403667,
                 `nsx(log(y), df = 3, intercept = FALSE)3` = 2.83834140526638,
                 theta = -1.55451789664845), value = 3629.13112243673))
})
