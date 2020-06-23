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
    c(-267.988148455606, 76.1320534565749, 298.740066381521,
      78.4837230074173, -190.722172637887, -105.395726242187),
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
    c(-200.08884765849, -15.4083728748582, 106.693482823233,
      -123.895063426942, -99.7090137336823, -133.106028981174),
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
    c(1180.59297436751, -3238.7442226338, 2435.97387283985,
      3644.49430227336, 3638.23404178443, 3398.05918852532))
})
