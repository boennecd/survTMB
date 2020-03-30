context("testing SNVA")

get_par_val_eortc <- function(x){
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

        eps <- .Machine$double.eps^(3/5)
        func$snva$control <- eps
        res <- do.call(optim, func$snva)

        expect_known_value(
          get_par_val_eortc(res), sprintf("test-res/SNVA-%s-%s.RDS",
                                          link, param_type),
          tolerance = sqrt(eps))
      })

test_that("PH", {
  func <- get_func_eortc("PH", 1L, "CP_trans")

  eps <- .Machine$double.eps^(3/5)
  func$snva$control <- eps
  res <- do.call(optim, func$snva)

  expect_equal(
    get_par_val_eortc(res),
    list(par = c(`(Intercept)` = -8.74958889264355,
                 trt = 0.627618058968419, `nsx(log(y), df = 3, intercept = FALSE)1` = 6.17805616140279,
                 `nsx(log(y), df = 3, intercept = FALSE)2` = 13.1244767484505,
                 `nsx(log(y), df = 3, intercept = FALSE)3` = 4.77286713362206,
                 theta = -1.34336688466204), value = 3629.75909884551),
    tolerance = sqrt(eps))
})

test_that("PO", {
  func <- get_func_eortc("PO", 1L, "CP_trans")

  eps <- .Machine$double.eps^(3/5)
  func$snva$control <- eps
  res <- do.call(optim, func$snva)

  expect_equal(
    get_par_val_eortc(res),
    list(par = c(`(Intercept)` = -8.9148514143319,
                 trt = 0.87706216683176, `nsx(log(y), df = 3, intercept = FALSE)1` = 6.42861488790485,
                 `nsx(log(y), df = 3, intercept = FALSE)2` = 13.4492463731085,
                 `nsx(log(y), df = 3, intercept = FALSE)3` = 5.41890519741321,
                 theta = -1.02243720013472), value = 3629.28307310362),
    tolerance = sqrt(eps))
})

test_that("probit", {
  func <- get_func_eortc("probit", 1L, "CP_trans")

  eps <- .Machine$double.eps^(3/5)
  func$snva$control <- eps
  res <- do.call(optim, func$snva)

  expect_equal(
    get_par_val_eortc(res),
    list(par = c(`(Intercept)` = -4.05421750814269,
                 trt = 0.524883805380364, `nsx(log(y), df = 3, intercept = FALSE)1` = 2.9017701774533,
                 `nsx(log(y), df = 3, intercept = FALSE)2` = 5.55262972082832,
                 `nsx(log(y), df = 3, intercept = FALSE)3` = 2.83841845310392,
                 theta = -1.55358670744936), value = 3629.13051119782),
    tolerance = sqrt(eps))
})
