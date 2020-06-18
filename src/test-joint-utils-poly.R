#####
# example used in test file

without_b <- FALSE
without_g <- FALSE
without_m <- FALSE

# parameters
bas <- with(new.env(), {
  bas <- poly(1:10, 2)
  function(x){
    out <- predict(bas, x)
    cbind(1, out)
  }
})
omega <- if(without_b) numeric() else c(-0.74, 0.76, 0.10)
Psi <- structure(c(1.08, 0.12, -0.36, -0.48, 0.36, -0.12, 0.12, 0.36,
                   0, 0.12, 0, -0.12, -0.36, 0, 0.84, 0.12, 0.12, 0.12, -0.48, 0.12,
                   0.12, 0.84, -0.12, 0.24, 0.36, 0, 0.12, -0.12, 0.84, -0.12, -0.12,
                   -0.12, 0.12, 0.24, -0.12, 0.6), .Dim = c(6L, 6L))
B <- if(without_g)
  numeric() else
    structure(c(0.97, 0.01, .02, -0.78, -0.86, -0.03),
              .Dim = c(3, 2))
dput(t(B))

alpha <- c(0.7, 0.6)

if(without_m){
  mu <- rho <- Lambda <- numeric()
} else {
  mu     <- c(-0.79, -0.67, 0.05, -0.71, -0.13, -0.33)
  rho    <- c(0.13, -0.93, 0.93, -0.11, -0.04, -0.44)
  Lambda <- Psi
}

# spline functions
library(splines)

b_func <- bas
if(without_b)
  b_func <- NULL

m_func <- bas
if(without_m)
  m_func <- NULL

g_func <- bas
if(without_g)
  g_func <- NULL

# get nodes to use for Gauss–Legendre quadrature
n_nodes <- 30L
wdat <- survTMB:::get_gl_rule(n_nodes)
.time_eps <- sqrt(.Machine$double.eps)

# approximate integral in variational approximation
if(!is.null(m_func)){
  k <- Lambda %*% rho
  dput(k <- drop(k / sqrt(1 + drop(rho %*% k))))
  dput(U <- mu)
} else {
  k <- U <- numeric()
}

library(compiler)
.integrant_inner <- function(ti, omega, B, U, Lambda, k, alpha){
  if(!is.null(b_func))
    b  <- c(b_func(ti))

  if(!is.null(m_func)){
    mb <- c(m_func(ti))
    m  <- rep(mb, length(alpha))
    ma <- rep(alpha, each = length(mb)) * m
  }

  if(!is.null(g_func)){
    gb <- c(g_func(ti))
    g <- rep(gb, length(alpha))
    ga <- rep(alpha, each = length(gb)) * g
  }

  v1 <- 0.
  if(!is.null(b_func))
    v1 <- v1 + drop(omega %*% b)
  if(!is.null(g_func))
    v1 <- v1 + drop(ga %*% c(B))
  if(!is.null(m_func)){
    v1 <- v1 + drop(ma %*% U + (ma %*% (Lambda %*% ma)) * .5)
    v2 <- pnorm(drop(ma %*% k), log.p = TRUE)
  } else
    v2 <- 0.
  2 * exp(v1 + v2)
}
.integrant_inner <- cmpfun(.integrant_inner)

integrant <- function(ti, omega, B, U, Lambda, k, alpha)
  vapply(ti, .integrant_inner, FUN.VALUE = numeric(1L),
         omega = omega, U = U, Lambda = Lambda, k = k, alpha = alpha, B = B)

va_integral <- function(lb, ub, omega, B, U, Lambda, k, alpha){
  nodes <- (ub - lb) / 2 * wdat$node + (ub + lb) / 2
  f <- integrant(ti = nodes, omega = omega, U = U, Lambda = Lambda, k = k,
                 alpha = alpha, B = B)
  (ub - lb) / 2 * drop(wdat$weight %*% f)
}

lb <- 1.2
ub <- 9.5

dput(va_integral(lb = lb, ub = ub, omega = omega, U = U, Lambda = Lambda,
                 k = k, alpha = alpha, B = B))

# approximate the gradient
.integrant_inner_grad <- function(ti, omega, B, U, Lambda, k, alpha){
  if(!is.null(b_func))
    b  <- c(b_func(ti))

  if(!is.null(m_func)){
    mb <- c(m_func(ti))
    m  <- rep(mb, length(alpha))
    M  <- (diag(length(alpha)) %x% t(mb))
    ma <- rep(alpha, each = length(mb)) * m
  }

  if(!is.null(g_func)){
    gb <- c(g_func(ti))
    g  <- rep(gb, length(alpha))
    G  <- (diag(length(alpha)) %x% t(gb))
    ga <- rep(alpha, each = length(gb)) * g
  }

  v1 <- 0.
  if(!is.null(b_func))
    v1 <- v1 + drop(omega %*% b)
  if(!is.null(g_func))
    v1 <- v1 + drop(ga %*% c(B))
  if(!is.null(m_func)){
    v1 <- v1 + drop(ma %*% U + (ma %*% (Lambda %*% ma)) * .5)
    v2 <- pnorm(drop(ma %*% k), log.p = TRUE)
    v3 <- dnorm(drop(ma %*% k), log = TRUE)
  } else
    v2 <- v3 <- 0.

  g <- 2 * exp(v1 + v2)
  h <- 2 * exp(v1 + v3)

  do <- if(!is.null(b_func))
    g * b else numeric()
  da <- numeric(length(alpha))
  if(!is.null(m_func))
    da <- da + M %*% (g * (U + Lambda %*% crossprod(M, alpha)) + h * k)
  if(!is.null(g_func))
    da <- da + G %*% (g  * c(B))
  dB <- if(!is.null(g_func))
    g * ga else numeric()
  dU <- if(!is.null(m_func)) g * ma else numeric()
  dLambda <- if(!is.null(m_func))
    .5 * g * c(ma %o% ma) else numeric()
  dk <- if(!is.null(m_func))
    h * ma else numeric()

  c(omega = do, alpha = da, dB = dB, U = dU, k = dk, Lambda = dLambda)
}
.integrant_inner_grad <- cmpfun(.integrant_inner_grad)

integrant_grad <- function(ti, omega, B, U, Lambda, k, alpha){
  len <- length(omega) + length(B) + length(U) + length(Lambda) +
    length(k) + length(alpha)
  vapply(
    ti, .integrant_inner_grad, FUN.VALUE = numeric(len), B = B,
    omega = omega, U = U, Lambda = Lambda, k = k, alpha = alpha)
}

va_integral_grad <- function(lb, ub, omega, B, U, Lambda, k, alpha){
  nodes <- (ub - lb) / 2 * wdat$node + (ub + lb) / 2
  f <- integrant_grad(ti = nodes, omega = omega, U = U, Lambda = Lambda, k = k,
                      alpha = alpha, B = B)

  rowSums(rep((ub - lb) / 2 * wdat$weight, each = NROW(f)) * f)
}

my_grad <- va_integral_grad(
  lb = lb, ub = ub, omega = omega, B = B, U = U, Lambda = Lambda, k = k,
  alpha = alpha)

# check that the result is correct
library(numDeriv)
wrap <- function(par){
  # get the parameters
  if(!is.null(b_func)){
    lo <- length(omega)
    o <- head(par, lo)
    par <- par[-(1:lo)]
  }

  la <- length(alpha)
  a <- par[1:la]
  par <- par[-(1:la)]

  if(!is.null(g_func)){
    lB <- length(B)
    B_val <- par[1:lB]
    par <- par[-(1:lB)]
  }

  if(!is.null(m_func)){
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
  }

  va_integral(lb = lb, ub = ub, omega = o, alpha = a, U = u,
              Lambda = L, k = k, B = B_val)
}

num_grad <- grad(wrap, c(omega, alpha, B, U, k,
                         Lambda[lower.tri(Lambda, TRUE)]))

if(!is.null(m_func)){
  lamb <- tail(my_grad, length(Lambda))
  lamb[lower.tri(Lambda)] <- 2 * lamb[lower.tri(Lambda)]
  lamb <- lamb[lower.tri(Lambda, TRUE)]

  tmp <- head(my_grad, length(num_grad))
  tmp[1:length(lamb) + length(tmp) - length(lamb)] <- lamb
} else
  tmp <- my_grad

rbind(tmp, num_grad)
tmp - num_grad
max(abs(tmp - num_grad))

length(my_grad)
dput(unname(my_grad))
