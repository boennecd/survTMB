context("testing 'make_joint_ADFun'")

for(n_threads in 1:2){
  test_that(sprintf(
    "fn and gr gives the same with all components (n_threads: %d)",
    n_threads), {
      dat <- readRDS(get_test_file_name("joint-all.RDS"))

      eps <- .Machine$double.eps^(1/4)
      opt_func <- function(par, fn, gr, ...)
        optim(par, fn, gr, control = list(reltol = eps, maxit = 10000L),
              method = "BFGS")

      out <- make_joint_ADFun(
        sformula =  Surv(left_trunc, y, event) ~ Z1 + Z2,
        mformula = cbind(Y1, Y2) ~ X1,
        id_var = id, time_var = obs_time, skew_start = -1e-16,
        sdata = dat$survival_data, mdata = dat$marker_data,
        mknots = dat$params$m_attr$knots, sknots = dat$params$b_attr$knots,
        gknots = dat$params$g_attr$knots, n_nodes = 15L,
        n_threads = n_threads, opt_func = opt_func)

      expect_known_value(
        head(out$par        , 45),
        file.path(test_res_dir, "joint-all-par.RDS"),
        tolerance = sqrt(eps) * 10)
      # dput(out$fn(out$par))
      expect_equal(out$fn(out$par), 342.626333272812,
                   tolerance = sqrt(eps) * 10)
      expect_known_value(
        head(out$gr(out$par), 45),
        file.path(test_res_dir, "joint-all-gr.RDS"),
        tolerance = sqrt(eps) * 10)
  })
}
