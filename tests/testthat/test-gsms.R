context("testing gsm functions")

test_that("Gradients are correct with a high epsilon", {
  f <- file.path("local-tests", "fail-mgsm.RDS")
  if(!file.exists(f))
    f <- file.path("tests", "testthat", f)
  skip_if_not(file.exists(f))

  # this data set revealed a bug
  dat <- readRDS(f)
  get_fit <- function(link, eps = 1e-16){
    obj <- with(dat, survTMB:::gsm(
      formula = Surv(y, event) ~ x + treatment, data = dframe, df = 3L,
      link = link, n_threads = 1L, do_fit = FALSE))

    with(obj, survTMB:::gsm_fit(
      X = X, XD = XD, Z = Z, y = y, link = link, n_threads = 1L,
      offset_eta = numeric(), offset_etaD = numeric(),
      eps = eps))
  }

  fit <- get_fit("PH", eps = 1)
  # dput(fit$start_coef, control = c("keepNA", "keepInteger"))
  expect_equal(
    unname(fit$start_coef),
    c(4.86255677159865, 6.86769765659779, 6.33629031801487,
      -5.08757251566875, -0.00538590235340096, -0.0289071515215028))
  # dput(with(fit, numDeriv::grad(mlogli, start_coef)))
  expect_equal(
    drop(with(fit, grad(start_coef))),
    c(-4462.51544524202, -1477.02030303319, 695.465371884567,
      -5634.08636135467, -1501.32718063916, -3658.98920614737),
    tolerance = 1e-5)

  fit <- get_fit("PO", eps = 1)
  # dput(fit$start_coef, control = c("keepNA", "keepInteger"))
  expect_equal(
    unname(fit$start_coef),
    c(5.06327194574772, 9.1761419521409, 9.21245263420721, -5.2133840102223,
      0.00572938602129334, -0.0228969047980255))
  # dput(with(fit, numDeriv::grad(mlogli, start_coef)))
  expect_equal(
    drop(with(fit, grad(start_coef))),
    c(-2987.58250085956, -932.702279672211, 512.066111198917,
      -3635.17501093579, -993.50437332079, -2417.89958934265),
    tolerance = 1e-5)

  fit <- get_fit("probit")
  # dput(fit$start_coef, control = c("keepNA", "keepInteger"))
  expect_equal(
    unname(fit$start_coef),
    c(2.55891961359501, 3.78842886000736, 4.98626484286831,
      -2.33733568340046, -0.00486370289392708, -0.0216345308794209))
  # dput(with(fit, numDeriv::grad(mlogli, start_coef)))
  expect_equal(
    drop(with(fit, grad(start_coef))),
    c(11115.8432497678, -33621.8717162036, 22883.8510911714,
      32976.8315701619, 36842.5634464707, 32730.3964620881))

  clear_cppad_mem(2L)
})
