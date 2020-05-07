#####
# example used in test file

# parameters
omega <- c(-0.74, 0.76, 0.20)
Psi <- structure(c(1.08, 0.12, -0.36, -0.48, 0.36, -0.12, 0.12, 0.36,
                   0, 0.12, 0, -0.12, -0.36, 0, 0.84, 0.12, 0.12, 0.12, -0.48, 0.12,
                   0.12, 0.84, -0.12, 0.24, 0.36, 0, 0.12, -0.12, 0.84, -0.12, -0.12,
                   -0.12, 0.12, 0.24, -0.12, 0.6), .Dim = c(6L, 6L))
B <- structure(c(0.97, 0.01, -0.07, -0.78, -0.86, 0.98), .Dim = 3:2)
alpha <- c(0.7, 0.6)

mu     <- c(-0.79, -0.67, 0.05, -0.71, -0.13, -0.33)
rho    <- c(0.13, -0.93, 0.93, -0.11, -0.04, -0.44)
Lambda <- Psi

# spline functions
library(splines)

dput(b_knots <- seq(log(1), log(10), length.out = 4))
b_func <- with(new.env(), {
  ks <- b_knots

  function(x)
    ns(log(x), knots = ks[-c(1L, length(ks))],
       Boundary.knots = ks[ c(1L, length(ks))],
       intercept = FALSE)
})

dput(fix_knots <- seq(0, 10, length.out = 3))
mark_fix_func <- with(new.env(), {
  ks <- fix_knots

  function(x)
    ns(x, knots = ks[-c(1L, length(ks))],
       Boundary.knots = ks[ c(1L, length(ks))],
       intercept = TRUE)
})

# get nodes to use for Gauss–Legendre quadrature
n_nodes <- 30L
wdat <- survTMB:::get_gl_rule(n_nodes)
.time_eps <- sqrt(.Machine$double.eps)

# approximate integral in variational approximation
k <- Lambda %*% rho
dput(k <- drop(k / sqrt(1 + drop(rho %*% k))))
dput(U <- c(B) + mu)

library(compiler)
.integrant_inner <- function(ti, omega, U, Lambda, k, alpha){
  b  <- c(b_func(ti))
  mb <- c(mark_fix_func(ti))
  m  <- rep(mb, length(alpha))
  ma <- rep(alpha, each = length(mb)) * m

  v1 <- drop(omega %*% b + ma %*% U + (ma %*% (Lambda %*% ma)) * .5)
  v2 <- pnorm(drop(ma %*% k), log.p = TRUE)
  2 * exp(v1 + v2)
}
.integrant_inner <- cmpfun(.integrant_inner)

integrant <- function(ti, omega, U, Lambda, k, alpha)
  vapply(ti, .integrant_inner, FUN.VALUE = numeric(1L),
         omega = omega, U = U, Lambda = Lambda, k = k, alpha = alpha)

va_integral <- function(lb, ub, omega, U, Lambda, k, alpha){
  nodes <- (ub - lb) / 2 * wdat$node + (ub + lb) / 2
  f <- integrant(ti = nodes, omega = omega, U = U, Lambda = Lambda, k = k,
                 alpha = alpha)
  (ub - lb) / 2 * drop(wdat$weight %*% f)
}

lb <- 2.4
ub <- 6.7

dput(va_integral(lb = lb, ub = ub, omega = omega, U = U, Lambda = Lambda,
                 k = k, alpha = alpha))

# approximate the gradient
.integrant_inner_grad <- function(ti, omega, U, Lambda, k, alpha){
  b  <- c(b_func(ti))
  mb <- c(mark_fix_func(ti))
  m  <- rep(mb, length(alpha))
  M  <- (diag(length(alpha)) %x% t(mb))
  ma <- rep(alpha, each = length(mb)) * m

  v1 <- drop(omega %*% b + ma %*% U + (ma %*% (Lambda %*% ma)) * .5)
  v2 <- pnorm(drop(ma %*% k), log.p = TRUE)
  v3 <- dnorm(drop(ma %*% k), log = TRUE)
  g <- 2 * exp(v1 + v2)
  h <- 2 * exp(v1 + v3)

  do <- g * b
  da <-
    M %*% (g * (U + Lambda %*% crossprod(M, alpha)) + h * k)
  dU <- g * ma
  dLambda <- .5 * g * c(ma %o% ma)
  dk <- h * ma

  c(omega = do, alpha = da, U = dU, k = dk, Lambda = dLambda)
}
.integrant_inner_grad <- cmpfun(.integrant_inner_grad)

integrant_grad <- function(ti, omega, U, Lambda, k, alpha){
  len <- length(omega) + length(U) + length(Lambda) + length(k) + length(alpha)
  vapply(
    ti, .integrant_inner_grad, FUN.VALUE = numeric(len),
    omega = omega, U = U, Lambda = Lambda, k = k, alpha = alpha)
}

va_integral_grad <- function(lb, ub, omega, U, Lambda, k, alpha){
  nodes <- (ub - lb) / 2 * wdat$node + (ub + lb) / 2
  f <- integrant_grad(ti = nodes, omega = omega, U = U, Lambda = Lambda, k = k,
                      alpha = alpha)

  rowSums(rep((ub - lb) / 2 * wdat$weight, each = NROW(f)) * f)
}

my_grad <- va_integral_grad(lb = lb, ub = ub, omega = omega, U = U, Lambda = Lambda,
                            k = k, alpha = alpha)

# check that the result is correct
library(numDeriv)
wrap <- function(par){
  # get the parameters
  lo <- length(omega)
  o <- head(par, lo)
  par <- par[-(1:lo)]

  la <- length(alpha)
  a <- par[1:la]
  par <- par[-(1:la)]

  lU <- length(U)
  u <- par[1:lU]
  par <- par[-(1:lU)]

  lk <- length(k)
  k <- par[1:lk]
  par <- par[-(1:lk)]

  n <- NCOL(Lambda)
  L <- matrix(nr = n, nc = n)
  L[lower.tri(L, TRUE)] <- par
  L[upper.tri(L)] <- t(L)[upper.tri(L)]

  va_integral(lb = lb, ub = ub, omega = o, alpha = a, U = u,
              Lambda = L, k = k)
}

num_grad <- grad(wrap, c(omega, alpha, U, k, Lambda[lower.tri(Lambda, TRUE)]))

lamb <- tail(my_grad, length(Lambda))
lamb[lower.tri(Lambda)] <- 2 * lamb[lower.tri(Lambda)]
lamb <- lamb[lower.tri(Lambda, TRUE)]

tmp <- head(my_grad, length(num_grad))
tmp[1:length(lamb) + length(tmp) - length(lamb)] <- lamb

rbind(tmp, num_grad)
tmp - num_grad

dput(unname(my_grad))
