context("pedigree functions")

test_that("SNVA gives previous results", {
  f <- "pedigree-test-data.RDS"
  if(!file.exists("pedigree-test-data.RDS"))
    f <- file.path("tests", "testthat", f)
  dat <- readRDS(f)
  c_data <- lapply(dat$sim_data, function(x){
    data <- data.frame(Z = x$Z, y = x$y, event = x$event)
    cor_mats <- list(x$rel_mat)
    list(data = data, cor_mats = cor_mats)
  })

  sbase_haz <- function(x){
    x <- log(x)
    cbind(cubed = x^3, squared = x^2, x = x)
  }
  dsbase_haz <- function(x){
    y <- log(x)
    cbind(3 * y^2, 2 * y, 1) / x
  }

  rel_eps <- sqrt(.Machine$double.eps)
  func <- make_pedigree_ADFun(
    formula = Surv(y, event) ~ Z.1 + Z.2 - 1, skew_start = -.001,
    tformula  = ~  sbase_haz(y) - 1, n_nodes = 15L,
    dtformula = ~ dsbase_haz(y) - 1, method = "SNVA",
    c_data = c_data, link = "probit", n_threads = 1L,
    args_gva_opt = list(max_cg = 100L, c2 = .01, rel_eps = sqrt(rel_eps)))

  expect_known_value(func$par, update = FALSE,
                     file.path(test_res_dir, "pedigree-SNVA-start.RDS"),
                     tolerance = rel_eps^(1/4))

  opt <- optim_pedigree_psqn(
    func, max_cg = 100L, rel_eps = rel_eps,
    c2 = .01, cg_tol = .1)

  expect_known_value(opt[c("par", "params", "va_params")],
                     file.path(test_res_dir, "pedigree-SNVA-opt.RDS"),
                     tolerance = sqrt(rel_eps))
})
