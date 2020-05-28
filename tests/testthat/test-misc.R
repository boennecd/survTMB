context("testing miscellaneous functions")

test_that("cov -> .cov_to_theta -> .theta_to_cov gives the same matrix", {
  xcov <- 1.157
  expect_equal(.theta_to_cov(.cov_to_theta(xcov)), xcov,
               check.attributes = FALSE)

  xcov <- structure(c(0.43, 0.96, -0.27, 1.18, -0.4, 0.96, 3.84, -1.32,
                      2.45, 1.17, -0.27, -1.32, 2.98, -1.02, 2.32, 1.18, 2.45, -1.02,
                      3.86, -0.96, -0.4, 1.17, 2.32, -0.96, 7.76), .Dim = c(5L, 5L))
  expect_equal(.theta_to_cov(.cov_to_theta(xcov)), xcov,
               check.attributes = FALSE)
})

test_that(".cp_to_dp -> .dp_to_cp gives the same parameters", {
  set.seed(1)
  K <- 4L
  xi <- rnorm(K)
  Psi <- drop(rWishart(1, K, diag(K)))
  rho <- runif(K, -1, 1)
  alpha <- rho * sqrt(diag(Psi))

  cps <- .dp_to_cp(xi = xi, Psi = Psi, alpha = alpha)
  dps_back <- .cp_to_dp(mu = cps$mu, Sigma = cps$Sigma, gamma = cps$gamma)

  expect_equal(xi , dps_back$xi)
  expect_equal(Psi, dps_back$Psi)
  expect_equal(alpha, dps_back$alpha)
})
