# install package
if(FALSE)
  devtools::install_github("boennecd/SimSurvNMarker",
                           ref = "cb2f777e065b229e5739c372460f88b18b5f5d38")

#####
# x-all
library(SimSurvNMarker)
alpha <- c(0.32, -0.31)
omega <- c(-2.6, -1.32)
delta <- c(.2, -0.17)
gamma <- structure(c(0.14, -0.8), .Dim = 1:2)
B <- structure(c(-0.96, 0.33, 0.39, 0.26, -0.76, 0.19), .Dim = 3:2)
sig <- structure(c(0.03, 0, 0, 0.05), .Dim = c(2L, 2L))
Psi <- structure(c(
  1.57, -0.37, -0.08, -0.17, -0.37, 0.98, -0.05, 0.09,
  -0.08, -0.05, 0.87, 0.53, -0.17, 0.09, 0.53, 1.17), .Dim = c(4L, 4L))
n_obs <- 1000L

n_y <- length(alpha)
d_m <- NROW(Psi) / n_y
d_g <- NROW(B)
d_b <- length(omega)
d_z <- length(delta)
d_x <- NROW(gamma)

r_n_marker <- function(id)
  rpois(1, 10) + 1L
r_obs_time <- function(id, n_markes)
  sort(runif(n_markes, 0, 10))
r_z <- function(id)
  as.numeric(runif(d_z) > .5)
r_x <- function(id)
  as.numeric(runif(d_x) > .5)
r_left_trunc <- function(id)
  rbeta(1, 1, 2) * 3
r_right_cens <- function(id)
  rbeta(1, 2, 1) * 6 + 4

b_ks <- seq(log(1e-1), log(10), length.out = d_b)
m_ks <- seq(0     , 10     , length.out = d_m)
g_ks <- seq(0     , 10     , length.out = d_g)
b_func <- get_ns_spline(b_ks, do_log = TRUE)
m_func <- get_ns_spline(m_ks, do_log = FALSE)
g_func <- get_ns_spline(g_ks, do_log = FALSE)

gl_dat <- get_gl_rule(50L)

set.seed(95002011L)
dat <- sim_joint_data_set(
  n_obs = n_obs, B = B, Psi = Psi, omega = omega, delta = delta,
  alpha = alpha, sigma = sig, gamma = gamma, b_func = b_func,
  m_func = m_func, g_func = g_func, gl_dat = gl_dat, r_z = r_z,
  r_left_trunc = r_left_trunc, r_right_cens = r_right_cens,
  r_n_marker = r_n_marker, r_x = r_x, r_obs_time = r_obs_time,
  y_max = 10.)
saveRDS(dat, file.path("inst", "test-data", "large-joint-all.RDS"))

set.seed(95002011L)
dat <- sim_joint_data_set(
  n_obs = 25L, B = B, Psi = Psi, omega = omega, delta = delta,
  alpha = alpha, sigma = sig, gamma = gamma, b_func = b_func,
  m_func = m_func, g_func = g_func, gl_dat = gl_dat, r_z = r_z,
  r_left_trunc = r_left_trunc, r_right_cens = r_right_cens,
  r_n_marker = r_n_marker, r_x = r_x, r_obs_time = r_obs_time,
  y_max = 10.)
saveRDS(dat, file.path("inst", "test-data", "joint-all.RDS"))
