
# survTMB

[![Build Status on
Travis](https://travis-ci.org/boennecd/survTMB.svg?branch=master,osx)](https://travis-ci.org/boennecd/survTMB)
<!-- [![](https://www.r-pkg.org/badges/version/survTMB)](https://www.r-pkg.org/badges/version/survTMB) -->
<!-- [![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/survTMB)](https://cran.r-project.org/package=survTMB) -->

This package contains methods to estimated mixed generalized survival
models (Liu, Pawitan, and Clements 2016, 2017). All methods use
automatic differentiation using the CppAD library (Bell 2019) through
[the TMB package](https://github.com/kaskr/adcomp) (Kristensen et al.
2016). The estimation methods are

  - a Laplace approximation using
    [TMB](https://github.com/kaskr/adcomp).
  - Gaussian variational approximation (GVA) similar to the method shown
    by Ormerod and Wand (2012).
  - Skew-normal variational approximation (SNVA) similar to the method
    shown by Ormerod (2011).

The [example](#example) section shows an example of how to use the
package with different methods. The [benchmark](#benchmark) section
shows a comparison of the computation time of the methods.

Joint marker and survival models are also available in the package. We
show an example of estimating a joint model in the [joint
models](#joint-models) section.

## Example

We estimate a GSM a below with the proportional odds (PO) link function
using both a Laplace approximation, a GVA, and a SNVA. First, we define
a function to perform the estimation.

``` r
# assign variable with data 
dat <- coxme::eortc

# assign function to estimate the model
library(survTMB)
#> Loading required package: splines
#> Loading required package: survival
library(survival)
fit_model <- function(link, n_threads = 2L, method = "Laplace", 
                      param_type = "DP", dense_hess = FALSE, 
                      sparse_hess = FALSE, do_fit = TRUE, 
                      do_free = FALSE)
  eval(bquote({
    adfun <- make_mgsm_ADFun(
      Surv(y, uncens) ~ trt, cluster = as.factor(center), 
      Z = ~ trt, df = 3L, data = dat, link = .(link), do_setup = .(method), 
      n_threads = .(n_threads), param_type = .(param_type), n_nodes = 15L, 
      dense_hess = .(dense_hess), sparse_hess = .(sparse_hess))
    fit <- if(.(do_fit))
      fit_mgsm(adfun, method = .(method)) else NULL
    if(.(do_free)){
      free_laplace(adfun)
      clear_cppad_mem(.(n_threads), keep_work_space = TRUE)
    }
    list(fit = fit, fun = adfun)
  }), parent.frame())

# # estimate the model using different methods. Start w/ Laplace
# (lap_ph <- fit_model("PO"))$fit

# w/ GVA
(gva_fit <- fit_model("PO", method = "GVA"))$fit
#> 
#> MGSM estimated with method 'GVA' with link 'PO' from call:
#>   make_mgsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, 
#>       df = 3L, Z = ~trt, cluster = as.factor(center), do_setup = "GVA", 
#>       n_nodes = 15L, param_type = "DP", link = "PO", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "GVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                   -8.05                                    1.03 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                    5.70                                   11.82 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                    5.59 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0446 0.0585             0.211 0.806
#> trt              0.0585 0.1181             0.806 0.344
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13031.10

# w/ SNVA
(snva_fit <- fit_model("PO", method = "SNVA", param_type = "DP"))$fit
#> 
#> MGSM estimated with method 'SNVA' with link 'PO' from call:
#>   make_mgsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, 
#>       df = 3L, Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       n_nodes = 15L, param_type = "DP", link = "PO", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                   -8.05                                    1.03 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                    5.70                                   11.82 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                    5.59 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0446 0.0585             0.211 0.806
#> trt              0.0585 0.1182             0.806 0.344
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13031.10
```

### Computing the Hessian

The Hessian using a variational approximation (VA) can be computed as
both a dense matrix and as a sparse matrix. We show an example below
where we compare the two approaches.

``` r
library(microbenchmark) # needed for benchmarking
```

``` r
# fit model w/ GVA
fit <- fit_model("PO", method = "GVA", dense_hess = TRUE, 
                 sparse_hess = TRUE)

# compute dense Hessian
par <- with(fit$fit, c(params, va_params))
dense_hess <- fit$fun$gva$he(par)
num_hess <- numDeriv::jacobian(fit$fun$gva$gr, par)
all.equal(dense_hess, num_hess, tolerance = 1e-5)
#> [1] TRUE

# has many zeros (i.e. it is sparse)
mean(abs(dense_hess) > .Machine$double.eps) # fraction of non-zeros
#> [1] 0.102

# plot non-zero entries (black block's are non-zero; ignore upper triangle)
par(mar = c(1, 1, 1, 1))
is_non_zero <- t(abs(dense_hess) > .Machine$double.eps)
is_non_zero[upper.tri(is_non_zero)] <- FALSE
image(is_non_zero, xaxt = "n", yaxt = "n", 
      col = gray.colors(2, 1, 0))
```

<img src="man/figures/README-comp_hess-1.png" width="100%" />

``` r

# compute sparse Hessian
sparse_hess <- fit$fun$gva$he_sp(par)

# they are identical 
stopifnot(isTRUE(
  all.equal(as.matrix(sparse_hess), dense_hess, check.attributes = FALSE)))

# compare storage cost
as.numeric(object.size(dense_hess) / object.size(sparse_hess))
#> [1] 11

# we usually want the first part the inverse negative Hessian for the model 
# parameters. This can be computed as follows
library(Matrix)
n_vars <- length(fit$fit$params)
naiv_vcov <- function(hess)
  solve(hess)[1:n_vars, 1:n_vars]
alte_vcov <- function(hess){
  idx <- 1:n_vars
  A <- hess[ idx,  idx]
  C <- hess[-idx,  idx]
  D <- hess[-idx, -idx]
  solve(A - crossprod(C, solve(D, C)))
}

# these are the asymptotic standard deviations
structure(sqrt(diag(alte_vcov(dense_hess))), names = names(fit$fit$params))
#>                             (Intercept)                                     trt 
#>                                   0.420                                   0.109 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   0.282                                   0.816 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                   0.162                                   0.527 
#>                                   theta                                   theta 
#>                                   0.303                                   2.506

# check output is the same
stopifnot(
  isTRUE(all.equal(naiv_vcov(dense_hess), alte_vcov(dense_hess))),
  isTRUE(all.equal(naiv_vcov(dense_hess), as.matrix(alte_vcov(sparse_hess)), 
                   check.attributes = FALSE)),
  isTRUE(all.equal(naiv_vcov(dense_hess), as.matrix(naiv_vcov(sparse_hess)), 
                   check.attributes = FALSE)))

# compare computation time
microbenchmark(
  `Compute dense Hessian`               = fit$fun$gva$he(par), 
  `Compute sparse Hessian`              = fit$fun$gva$he_sp(par), 
  `Invert dense Hessian (naive)`        = naiv_vcov(dense_hess), 
  `Invert sparse Hessian (naive)`       = naiv_vcov(sparse_hess),
  `Invert dense Hessian (alternative)`  = alte_vcov(dense_hess), 
  `Invert sparse Hessian (alternative)` = alte_vcov(sparse_hess),
  times = 10)
#> Unit: milliseconds
#>                                 expr    min     lq   mean median     uq    max
#>                Compute dense Hessian 313.42 317.32 320.51 318.23 324.65 329.31
#>               Compute sparse Hessian  18.37  18.69  19.02  19.01  19.28  19.77
#>         Invert dense Hessian (naive)   5.26   5.30   5.39   5.35   5.42   5.73
#>        Invert sparse Hessian (naive)   1.06   1.07   1.13   1.09   1.12   1.30
#>   Invert dense Hessian (alternative)   1.27   1.33   1.35   1.35   1.37   1.47
#>  Invert sparse Hessian (alternative)   2.71   2.77   2.96   2.87   3.09   3.52
#>  neval
#>     10
#>     10
#>     10
#>     10
#>     10
#>     10
```

The sparse matrix only becomes more favorable for larger data sets (that
is, in terms of the number of clusters). However,
[recording](https://www.coin-or.org/CppAD/Doc/independent.htm) takes
some time and requires additional memory. We illustrate the additional
time below.

``` r
microbenchmark(
  `W/o Hessians     ` = fit_model("PO", method = "GVA", do_fit = FALSE), 
  `W/ dense Hessian ` = fit_model("PO", method = "GVA", do_fit = FALSE, 
                                  dense_hess = TRUE), 
  `W/ sparse Hessian` = fit_model("PO", method = "GVA", do_fit = FALSE, 
                                  sparse_hess = TRUE), 
  times = 10)
#> Unit: milliseconds
#>               expr min    lq  mean median    uq   max neval
#>  W/o Hessians       28  28.5  29.4   28.9  30.9  31.7    10
#>  W/ dense Hessian   78  80.1  80.6   80.6  81.2  84.1    10
#>  W/ sparse Hessian 634 638.0 643.8  645.3 648.8 651.6    10
```

### Approximation of the Conditional Distribution

The variational parameters provide an approximation of the conditional
distribution given the data and parameters or the posterior in a
Bayesian view. As an example, we can look at the multivariate normal
distribution approximation which is made by the GVA for the first group
below.

``` r
va_params <- gva_fit$fit$va_params
is_this_group <- which(grepl("^g1:", names(va_params)))
n_random_effects <- 2L

# conditional mean of random effects
va_params[is_this_group][seq_len(n_random_effects)]
#> g1:mu1 g1:mu2 
#>  0.367  0.614

# conditional covariance matrix of random effects
theta_to_cov(va_params[is_this_group][-seq_len(n_random_effects)])
#>         [,1]    [,2]
#> [1,] 0.01314 0.00623
#> [2,] 0.00623 0.02920
```

We can compare this with the multivariate skew-normal distribution
approximation from the SNVA.

``` r
va_params <- snva_fit$fit$va_params
is_this_group <- which(grepl("^g1:", names(va_params)))
n_random_effects <- 2L

xi <- va_params[is_this_group][seq_len(n_random_effects)]
Psi <- head(tail(va_params[is_this_group], -n_random_effects), 
            -n_random_effects)
Psi <- theta_to_cov(Psi)
alpha <- tail(va_params[is_this_group], n_random_effects)

# conditional mean, covariance matrix, and Pearson's moment coefficient of 
# skewness
dp_to_cp(xi = xi, Psi = Psi, alpha = alpha)
#> $mu
#> g1:mu1 g1:mu2 
#>  0.367  0.614 
#> 
#> $Sigma
#>         [,1]    [,2]
#> [1,] 0.01314 0.00622
#> [2,] 0.00622 0.02920
#> 
#> $gamma
#> [1] -1e-04 -1e-04
```

from the default values possibly because the lower bound is quite flat
in these parameters in this area.

``` r
skews <- sapply(1:37, function(id){
  va_params <- snva_fit$fit$va_params
  is_this_group <- which(grepl(paste0("^g", id, ":"), names(va_params)))
  
  xi <- va_params[is_this_group][seq_len(n_random_effects)]
  Psi <- head(tail(va_params[is_this_group], -n_random_effects), 
              -n_random_effects)
  Psi <- theta_to_cov(Psi)
  alpha <- tail(va_params[is_this_group], n_random_effects)
  dp_to_cp(xi = xi, Psi = Psi, alpha = alpha)$gamma
})

apply(skews, 1L, quantile, probs = seq(0, 1, by = .25))
#>           [,1]      [,2]
#> 0%   -1.00e-04 -1.00e-04
#> 25%  -1.00e-04 -1.00e-04
#> 50%  -1.00e-04 -1.00e-04
#> 75%  -1.00e-04 -1.00e-04
#> 100% -9.99e-05 -9.99e-05
```

Again, the skewness parameter have not moved much from the defaults.

### Other link functions

We estimate the same model below with other link functions.

``` r
# ######
# # w/ Laplace
# fit_model("PH"    , do_free = TRUE)$fit
# fit_model("PO"    , do_free = TRUE)$fit
# fit_model("probit", do_free = TRUE)$fit

######
# w/ GVA
fit_model("PH"    , method = "GVA")$fit
#> 
#> MGSM estimated with method 'GVA' with link 'PH' from call:
#>   make_mgsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, 
#>       df = 3L, Z = ~trt, cluster = as.factor(center), do_setup = "GVA", 
#>       n_nodes = 15L, param_type = "DP", link = "PH", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "GVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -7.827                                   0.724 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   5.394                                  11.390 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   4.800 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0292 0.0256             0.171 0.640
#> trt              0.0256 0.0550             0.640 0.234
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13026.74
fit_model("PO"    , method = "GVA")$fit
#> 
#> MGSM estimated with method 'GVA' with link 'PO' from call:
#>   make_mgsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, 
#>       df = 3L, Z = ~trt, cluster = as.factor(center), do_setup = "GVA", 
#>       n_nodes = 15L, param_type = "DP", link = "PO", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "GVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                   -8.05                                    1.03 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                    5.70                                   11.82 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                    5.59 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0446 0.0585             0.211 0.806
#> trt              0.0585 0.1181             0.806 0.344
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13031.10
fit_model("probit", method = "GVA")$fit
#> 
#> MGSM estimated with method 'GVA' with link 'probit' from call:
#>   make_mgsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, 
#>       df = 3L, Z = ~trt, cluster = as.factor(center), do_setup = "GVA", 
#>       n_nodes = 15L, param_type = "DP", link = "probit", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "GVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -3.742                                   0.602 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   2.660                                   5.004 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   2.972 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0195 0.0149             0.140 0.513
#> trt              0.0149 0.0434             0.513 0.208
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13035.14

######
# w/ SNVA (DP: direct parameterization)
fit_model("PH"    , method = "SNVA", param_type = "DP")$fit
#> 
#> MGSM estimated with method 'SNVA' with link 'PH' from call:
#>   make_mgsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, 
#>       df = 3L, Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       n_nodes = 15L, param_type = "DP", link = "PH", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -7.827                                   0.724 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   5.394                                  11.390 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   4.800 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0292 0.0256             0.171 0.638
#> trt              0.0256 0.0548             0.638 0.234
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13026.74
fit_model("PO"    , method = "SNVA", param_type = "DP")$fit
#> 
#> MGSM estimated with method 'SNVA' with link 'PO' from call:
#>   make_mgsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, 
#>       df = 3L, Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       n_nodes = 15L, param_type = "DP", link = "PO", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                   -8.05                                    1.03 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                    5.70                                   11.82 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                    5.59 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0446 0.0585             0.211 0.806
#> trt              0.0585 0.1182             0.806 0.344
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13031.10
fit_model("probit", method = "SNVA", param_type = "DP")$fit
#> 
#> MGSM estimated with method 'SNVA' with link 'probit' from call:
#>   make_mgsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, 
#>       df = 3L, Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       n_nodes = 15L, param_type = "DP", link = "probit", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -3.742                                   0.602 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   2.660                                   5.004 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   2.973 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0195 0.0149             0.140 0.513
#> trt              0.0149 0.0434             0.513 0.208
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13035.14

######
# w/ SNVA (CP: centralized parameterization)
fit_model("PH"    , method = "SNVA", param_type = "CP_trans")$fit
#> 
#> MGSM estimated with method 'SNVA' with link 'PH' from call:
#>   make_mgsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, 
#>       df = 3L, Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       n_nodes = 15L, param_type = "CP_trans", link = "PH", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -7.827                                   0.724 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   5.394                                  11.390 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   4.800 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0292 0.0256             0.171 0.638
#> trt              0.0256 0.0548             0.638 0.234
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13026.74
fit_model("PO"    , method = "SNVA", param_type = "CP_trans")$fit
#> 
#> MGSM estimated with method 'SNVA' with link 'PO' from call:
#>   make_mgsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, 
#>       df = 3L, Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       n_nodes = 15L, param_type = "CP_trans", link = "PO", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                   -8.05                                    1.03 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                    5.70                                   11.82 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                    5.59 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0446 0.0585             0.211 0.806
#> trt              0.0585 0.1182             0.806 0.344
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13031.10
fit_model("probit", method = "SNVA", param_type = "CP_trans")$fit
#> 
#> MGSM estimated with method 'SNVA' with link 'probit' from call:
#>   make_mgsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, 
#>       df = 3L, Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       n_nodes = 15L, param_type = "CP_trans", link = "probit", 
#>       n_threads = 2L, dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -3.742                                   0.602 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   2.660                                   5.004 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   2.973 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0195 0.0149             0.140 0.513
#> trt              0.0149 0.0434             0.513 0.208
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13035.14
```

## Benchmark

We provide a benchmark of the estimation methods used in section
[example](#example) below.

``` r
clear_cppad_mem(4L)
#> [1] 1
```

``` r
for(mth in c("GVA")){
# for(mth in c("Laplace", "GVA")){
  msg <- sprintf("Method: %s", mth)
  cat(sprintf("\n%s\n%s\n", msg, 
              paste0(rep("-", nchar(msg)), collapse = "")))
  print(microbenchmark(
    `PH         ` = fit_model("PH"    , 1L, mth, do_free = TRUE),
    `PH     (2L)` = fit_model("PH"    , 2L, mth, do_free = TRUE),
    `PH     (4L)` = fit_model("PH"    , 4L, mth, do_free = TRUE),
    
    `PO         ` = fit_model("PO"    , 1L, mth, do_free = TRUE),
    `PO     (2L)` = fit_model("PO"    , 2L, mth, do_free = TRUE),
    `PO     (4L)` = fit_model("PO"    , 4L, mth, do_free = TRUE), 
    
    `probit     ` = fit_model("probit", 1L, mth, do_free = TRUE),
    `probit (2L)` = fit_model("probit", 2L, mth, do_free = TRUE),
    `probit (4L)` = fit_model("probit", 4L, mth, do_free = TRUE),
    times = 5))
}
#> 
#> Method: GVA
#> -----------
#> Unit: milliseconds
#>         expr   min    lq mean median    uq max neval
#>  PH          180.7 181.7  182  182.1 183.1 184     5
#>  PH     (2L) 114.5 114.9  116  116.8 117.0 119     5
#>  PH     (4L)  86.8  86.9  108   91.2  94.2 180     5
#>  PO          573.7 576.2  577  577.1 578.7 580     5
#>  PO     (2L) 337.6 339.5  341  340.7 341.9 344     5
#>  PO     (4L) 230.4 230.5  250  231.4 232.6 324     5
#>  probit      716.9 719.2  720  719.2 721.6 722     5
#>  probit (2L) 417.7 419.0  420  419.1 420.4 424     5
#>  probit (4L) 281.8 285.8  285  285.9 286.6 287     5
```

``` r
for(param_type in c("DP", "CP_trans")){
  mth <- "SNVA"
  msg <- sprintf("Method: %s (%s)", mth, param_type)
  cat(sprintf("\n%s\n%s\n", msg, 
              paste0(rep("-", nchar(msg)), collapse = "")))
  print(suppressWarnings(microbenchmark(
    `PH         ` = fit_model("PH"    , 1L, mth, param_type = param_type),
    `PH     (2L)` = fit_model("PH"    , 2L, mth, param_type = param_type),
    `PH     (4L)` = fit_model("PH"    , 4L, mth, param_type = param_type),
    
    `PO         ` = fit_model("PO"    , 1L, mth, param_type = param_type),
    `PO     (2L)` = fit_model("PO"    , 2L, mth, param_type = param_type),
    `PO     (4L)` = fit_model("PO"    , 4L, mth, param_type = param_type), 
    
    `probit     ` = fit_model("probit", 1L, mth, param_type = param_type),
    `probit (2L)` = fit_model("probit", 2L, mth, param_type = param_type),
    `probit (4L)` = fit_model("probit", 4L, mth, param_type = param_type),
    times = 5)))
}
#> 
#> Method: SNVA (DP)
#> -----------------
#> Unit: milliseconds
#>         expr min  lq mean median  uq max neval
#>  PH          213 215  216    215 218 218     5
#>  PH     (2L) 141 142  143    144 144 146     5
#>  PH     (4L) 105 107  109    108 108 116     5
#>  PO          736 740  744    740 741 765     5
#>  PO     (2L) 453 457  464    461 472 478     5
#>  PO     (4L) 306 309  311    313 313 316     5
#>  probit      910 914  922    918 922 943     5
#>  probit (2L) 554 559  565    567 571 577     5
#>  probit (4L) 366 374  374    374 375 382     5
#> 
#> Method: SNVA (CP_trans)
#> -----------------------
#> Unit: milliseconds
#>         expr min  lq mean median  uq max neval
#>  PH          284 288  291    292 295 296     5
#>  PH     (2L) 181 187  206    188 188 284     5
#>  PH     (4L) 143 144  144    144 144 145     5
#>  PO          741 742  744    742 744 750     5
#>  PO     (2L) 461 467  469    470 470 479     5
#>  PO     (4L) 309 310  317    316 323 325     5
#>  probit      920 922  929    924 936 943     5
#>  probit (2L) 551 553  564    566 567 581     5
#>  probit (4L) 370 373  376    375 375 385     5
```

``` r
rm(list = ls())
gc()
#>           used  (Mb) gc trigger  (Mb) max used  (Mb)
#> Ncells 2105482 112.5    4304134 229.9  2751727 147.0
#> Vcells 3609966  27.6    8412869  64.2  8412869  64.2
clear_cppad_mem(4L)
#> [1] 1
```

## Joint Models

We will use one of the test data sets in the
[inst/test-data](inst/test-data) directory. The data is generated with
the [inst/test-data/gen-test-data.R](inst/test-data/gen-test-data.R)
file which is available on Github. The file uses the
[SimSurvNMarker](https://github.com/boennecd/SimSurvNMarker) package to
simulate a data set. The model is simulated
from

<!-- $$\begin{align*} -->

<!-- \vec Y_{ij} \mid \vec U_i = \vec u_i -->

<!--   &\sim N^{(r)}(\vec \mu_i(s_{ij}, \vec u_i), \Sigma) -->

<!--   \\ -->

<!-- \vec\mu(s, \vec u) &= -->

<!--   \Gamma^\top \vec x_i + B^\top\vec g(s) + U^\top\vec m(s) -->

<!--   \\ -->

<!-- &= \left(I \otimes \vec x_i^\top\right)\text{vec}\Gamma -->

<!--      + \left(I \otimes \vec g(s)^\top\right)\text{vec} B -->

<!--      + \left(I \otimes \vec m(s)^\top\right) \vec u -->

<!--   \\ -->

<!-- \vec U_i &\sim N^{(K)}(\vec 0, \Psi) -->

<!--   \\ -->

<!-- h(t\mid \vec u) &= \exp\left( -->

<!--   \vec\omega^\top\vec b(t) + -->

<!--   \vec z_i^\top\vec\delta + -->

<!--   \vec\alpha^\top\vec\mu(t, \vec u) -->

<!--   \right) -->

<!--   \\ -->

<!-- &= \exp\Bigg( -->

<!--   \vec\omega^\top\vec b(t) + -->

<!--   \vec z_i^\top\vec\delta -->

<!--   + \vec 1^\top\left( -->

<!--   \text{diag}(\vec \alpha) \otimes \vec x_i^\top\right)\text{vec}\Gamma -->

<!--   + \vec 1^\top\left( -->

<!--   \text{diag}(\vec \alpha) \otimes \vec g(t)^\top\right)\text{vec} B \\ -->

<!-- &\hspace{50pt}+ \vec 1^\top\left( -->

<!--   \text{diag}(\vec \alpha) \otimes \vec m(t)^\top\right)\vec u -->

<!--   \Bigg) -->

<!-- \end{align*}$$ -->

  
![\\begin{align\*} \\vec Y\_{ij} \\mid \\vec U\_i = \\vec u\_i &\\sim
N^{(r)}(\\vec \\mu\_i(s\_{ij}, \\vec u\_i), \\Sigma) \\\\ \\vec\\mu(s,
\\vec u) &= \\Gamma^\\top \\vec x\_i + B^\\top\\vec g(s) + U^\\top\\vec
m(s) \\\\ &= \\left(I \\otimes \\vec
x\_i^\\top\\right)\\text{vec}\\Gamma + \\left(I \\otimes \\vec
g(s)^\\top\\right)\\text{vec} B + \\left(I \\otimes \\vec
m(s)^\\top\\right) \\vec u \\\\ \\vec U\_i &\\sim N^{(K)}(\\vec 0,
\\Psi) \\\\ h(t\\mid \\vec u) &= \\exp\\left( \\vec\\omega^\\top\\vec
b(t) + \\vec z\_i^\\top\\vec\\delta + \\vec\\alpha^\\top\\vec\\mu(t,
\\vec u) \\right) \\\\ &= \\exp\\Bigg( \\vec\\omega^\\top\\vec b(t) +
\\vec z\_i^\\top\\vec\\delta + \\vec 1^\\top\\left( \\text{diag}(\\vec
\\alpha) \\otimes \\vec x\_i^\\top\\right)\\text{vec}\\Gamma +
\\vec 1^\\top\\left( \\text{diag}(\\vec \\alpha) \\otimes \\vec
g(t)^\\top\\right)\\text{vec} B \\\\ &\\hspace{50pt}+
\\vec 1^\\top\\left( \\text{diag}(\\vec \\alpha) \\otimes \\vec
m(t)^\\top\\right)\\vec u \\Bigg)
\\end{align\*}](https://latex.codecogs.com/svg.latex?%5Cbegin%7Balign%2A%7D%20%20%5Cvec%20Y_%7Bij%7D%20%5Cmid%20%5Cvec%20U_i%20%3D%20%5Cvec%20u_i%20%20%20%20%26%5Csim%20N%5E%7B%28r%29%7D%28%5Cvec%20%5Cmu_i%28s_%7Bij%7D%2C%20%5Cvec%20u_i%29%2C%20%5CSigma%29%20%20%20%20%5C%5C%20%20%5Cvec%5Cmu%28s%2C%20%5Cvec%20u%29%20%26%3D%20%20%20%20%5CGamma%5E%5Ctop%20%5Cvec%20x_i%20%2B%20B%5E%5Ctop%5Cvec%20g%28s%29%20%2B%20U%5E%5Ctop%5Cvec%20m%28s%29%20%20%20%20%5C%5C%20%20%26%3D%20%5Cleft%28I%20%5Cotimes%20%5Cvec%20x_i%5E%5Ctop%5Cright%29%5Ctext%7Bvec%7D%5CGamma%20%20%20%20%20%20%20%2B%20%5Cleft%28I%20%5Cotimes%20%5Cvec%20g%28s%29%5E%5Ctop%5Cright%29%5Ctext%7Bvec%7D%20B%20%20%20%20%20%20%20%2B%20%5Cleft%28I%20%5Cotimes%20%5Cvec%20m%28s%29%5E%5Ctop%5Cright%29%20%5Cvec%20u%20%20%20%20%5C%5C%20%20%5Cvec%20U_i%20%26%5Csim%20N%5E%7B%28K%29%7D%28%5Cvec%200%2C%20%5CPsi%29%20%20%20%20%5C%5C%20%20h%28t%5Cmid%20%5Cvec%20u%29%20%26%3D%20%5Cexp%5Cleft%28%20%20%20%20%5Cvec%5Comega%5E%5Ctop%5Cvec%20b%28t%29%20%2B%20%20%20%20%5Cvec%20z_i%5E%5Ctop%5Cvec%5Cdelta%20%2B%20%20%20%20%5Cvec%5Calpha%5E%5Ctop%5Cvec%5Cmu%28t%2C%20%5Cvec%20u%29%20%20%20%20%5Cright%29%20%20%20%20%5C%5C%20%20%26%3D%20%5Cexp%5CBigg%28%20%20%20%20%5Cvec%5Comega%5E%5Ctop%5Cvec%20b%28t%29%20%2B%20%20%20%20%5Cvec%20z_i%5E%5Ctop%5Cvec%5Cdelta%20%20%20%20%2B%20%5Cvec%201%5E%5Ctop%5Cleft%28%20%20%20%20%5Ctext%7Bdiag%7D%28%5Cvec%20%5Calpha%29%20%5Cotimes%20%5Cvec%20x_i%5E%5Ctop%5Cright%29%5Ctext%7Bvec%7D%5CGamma%20%20%20%20%2B%20%5Cvec%201%5E%5Ctop%5Cleft%28%20%20%20%20%5Ctext%7Bdiag%7D%28%5Cvec%20%5Calpha%29%20%5Cotimes%20%5Cvec%20g%28t%29%5E%5Ctop%5Cright%29%5Ctext%7Bvec%7D%20B%20%5C%5C%20%20%26%5Chspace%7B50pt%7D%2B%20%5Cvec%201%5E%5Ctop%5Cleft%28%20%20%20%20%5Ctext%7Bdiag%7D%28%5Cvec%20%5Calpha%29%20%5Cotimes%20%5Cvec%20m%28t%29%5E%5Ctop%5Cright%29%5Cvec%20u%20%20%20%20%5CBigg%29%20%20%5Cend%7Balign%2A%7D
"\\begin{align*}  \\vec Y_{ij} \\mid \\vec U_i = \\vec u_i    &\\sim N^{(r)}(\\vec \\mu_i(s_{ij}, \\vec u_i), \\Sigma)    \\\\  \\vec\\mu(s, \\vec u) &=    \\Gamma^\\top \\vec x_i + B^\\top\\vec g(s) + U^\\top\\vec m(s)    \\\\  &= \\left(I \\otimes \\vec x_i^\\top\\right)\\text{vec}\\Gamma       + \\left(I \\otimes \\vec g(s)^\\top\\right)\\text{vec} B       + \\left(I \\otimes \\vec m(s)^\\top\\right) \\vec u    \\\\  \\vec U_i &\\sim N^{(K)}(\\vec 0, \\Psi)    \\\\  h(t\\mid \\vec u) &= \\exp\\left(    \\vec\\omega^\\top\\vec b(t) +    \\vec z_i^\\top\\vec\\delta +    \\vec\\alpha^\\top\\vec\\mu(t, \\vec u)    \\right)    \\\\  &= \\exp\\Bigg(    \\vec\\omega^\\top\\vec b(t) +    \\vec z_i^\\top\\vec\\delta    + \\vec 1^\\top\\left(    \\text{diag}(\\vec \\alpha) \\otimes \\vec x_i^\\top\\right)\\text{vec}\\Gamma    + \\vec 1^\\top\\left(    \\text{diag}(\\vec \\alpha) \\otimes \\vec g(t)^\\top\\right)\\text{vec} B \\\\  &\\hspace{50pt}+ \\vec 1^\\top\\left(    \\text{diag}(\\vec \\alpha) \\otimes \\vec m(t)^\\top\\right)\\vec u    \\Bigg)  \\end{align*}")  

where ![\\vec Y\_{ij}\\in\\mathbb
R^{n\_y}](https://latex.codecogs.com/svg.latex?%5Cvec%20Y_%7Bij%7D%5Cin%5Cmathbb%20R%5E%7Bn_y%7D
"\\vec Y_{ij}\\in\\mathbb R^{n_y}") is individual
![i](https://latex.codecogs.com/svg.latex?i "i")’s
![j](https://latex.codecogs.com/svg.latex?j "j")th observed marker at
time ![s\_{ij}](https://latex.codecogs.com/svg.latex?s_%7Bij%7D
"s_{ij}"), ![U\_i\\in\\mathbb
R^K](https://latex.codecogs.com/svg.latex?U_i%5Cin%5Cmathbb%20R%5EK
"U_i\\in\\mathbb R^K") is individual
![i](https://latex.codecogs.com/svg.latex?i "i")’s random effect, and
![h](https://latex.codecogs.com/svg.latex?h "h") is the instantaneous
hazard rate for the time-to-event outcome.
![\\vec\\alpha](https://latex.codecogs.com/svg.latex?%5Cvec%5Calpha
"\\vec\\alpha") is the so-called association parameter. It shows the
strength of the relation between the latent mean function,
![\\vec\\mu(t,\\vec
u)](https://latex.codecogs.com/svg.latex?%5Cvec%5Cmu%28t%2C%5Cvec%20u%29
"\\vec\\mu(t,\\vec u)"), and the log of the instantaneous rate,
![h(t\\mid \\vec
u)](https://latex.codecogs.com/svg.latex?h%28t%5Cmid%20%5Cvec%20u%29
"h(t\\mid \\vec u)"). ![\\vec
m(t)](https://latex.codecogs.com/svg.latex?%5Cvec%20m%28t%29
"\\vec m(t)"), ![\\vec
g(t)](https://latex.codecogs.com/svg.latex?%5Cvec%20g%28t%29
"\\vec g(t)") and ![\\vec
b(t)](https://latex.codecogs.com/svg.latex?%5Cvec%20b%28t%29
"\\vec b(t)") are basis expansions of time. As an example, these can be
a polynomial, a B-spline, or a natural cubic spline. The expansion for
the baseline hazard, ![\\vec
b(t)](https://latex.codecogs.com/svg.latex?%5Cvec%20b%28t%29
"\\vec b(t)"), is typically made on ![\\log
t](https://latex.codecogs.com/svg.latex?%5Clog%20t "\\log t") instead of
![t](https://latex.codecogs.com/svg.latex?t "t"). One reason is that the
model reduces to a Weibull distribution when a first polynomial is used
and ![\\vec\\alpha =
\\vec 0](https://latex.codecogs.com/svg.latex?%5Cvec%5Calpha%20%3D%20%5Cvec%200
"\\vec\\alpha = \\vec 0"). ![\\vec
x\_i](https://latex.codecogs.com/svg.latex?%5Cvec%20x_i "\\vec x_i") and
![\\vec z\_i](https://latex.codecogs.com/svg.latex?%5Cvec%20z_i
"\\vec z_i") are individual specific known covariates.

We start by loading the simulated data set.

``` r
dat <- readRDS(file.path("inst", "test-data", "large-joint-all.RDS"))

# the marker data
m_data <- dat$marker_data
head(m_data, 10)
#>    obs_time      Y1      Y2 X1 id
#> 1     2.358 -0.0313 -1.4347  1  1
#> 2     1.941 -1.6829  0.0705  0  2
#> 3     2.469 -1.1793 -0.2994  0  2
#> 4     0.679  0.5640 -0.4988  0  3
#> 5     0.633 -1.0614  0.1318  0  4
#> 6     0.704 -0.9367 -0.4324  0  4
#> 7     0.740 -0.6472  0.0842  0  4
#> 8     0.989 -1.0603  0.4740  0  4
#> 9     1.291 -0.5059  0.4103  0  4
#> 10    2.404 -0.6974  0.2074  0  4

# the survival data
s_data <- dat$survival_data
head(s_data, 10)
#>    Z1 Z2 left_trunc     y event id
#> 1   0  1      2.198 2.905  TRUE  1
#> 2   0  0      0.107 3.350  TRUE  2
#> 3   0  1      0.593 0.797  TRUE  3
#> 4   1  1      0.279 2.507  TRUE  4
#> 5   0  0      2.062 4.669  TRUE  5
#> 6   1  1      2.215 8.058 FALSE  6
#> 7   1  0      1.061 9.123 FALSE  7
#> 8   1  1      0.214 8.894 FALSE  8
#> 9   1  1      1.106 7.498 FALSE  9
#> 10  1  0      1.237 1.734  TRUE 10
```

There is

``` r
length(unique(s_data$id))
#> [1] 1000
length(unique(s_data$id)) == NROW(s_data) # one row per id
#> [1] TRUE
```

individuals who each has an average of

``` r
NROW(m_data) / length(unique(s_data$id))
#> [1] 3.74
```

observed markers. The data is simulated. Thus, we know the true
parameters. These are

``` r
dat$params[c("gamma", "B", "Psi", "omega", "delta", "alpha", "sigma")]
#> $gamma
#>      [,1] [,2]
#> [1,] 0.14 -0.8
#> 
#> $B
#>       [,1]  [,2]
#> [1,] -0.96  0.26
#> [2,]  0.33 -0.76
#> [3,]  0.39  0.19
#> 
#> $Psi
#>       [,1]  [,2]  [,3]  [,4]
#> [1,]  1.57 -0.37 -0.08 -0.17
#> [2,] -0.37  0.98 -0.05  0.09
#> [3,] -0.08 -0.05  0.87  0.53
#> [4,] -0.17  0.09  0.53  1.17
#> 
#> $omega
#> [1] -2.60 -1.32
#> 
#> $delta
#> [1]  0.20 -0.17
#> 
#> $alpha
#> [1]  0.32 -0.31
#> 
#> $sigma
#>      [,1] [,2]
#> [1,] 0.03 0.00
#> [2,] 0.00 0.05
```

We start by constructing the objective function in order to estimate the
model.

``` r
system.time(
  out <- make_joint_ADFun(
    sformula =  Surv(left_trunc, y, event) ~ Z1 + Z2, 
    mformula = cbind(Y1, Y2) ~ X1, 
    id_var = id, time_var = obs_time, 
    sdata = s_data, mdata = m_data, m_coefs = dat$params$m_attr$knots,
    s_coefs = dat$params$b_attr$knots, g_coefs = dat$params$g_attr$knots, 
    n_nodes = 30L, n_threads = 6L))
#>    user  system elapsed 
#>  11.942   0.032   4.297
```

Next, we fit the model using the default optimization function.

``` r
system.time(
  opt_out <- out$opt_func(
    out$par, out$fn, out$gr, control = list(maxit = 10000L)))
#>    user  system elapsed 
#>  69.163   0.016  11.535
```

The estimated lower bound of the log marginal likelihood at the optimum
is shown
below.

<!-- with(environment(out$fn), c(mark$ll, sr_dat$ll, mark$ll + sr_dat$ll)) -->

``` r
-opt_out$value
#> [1] -4124
```

Further, we can compare the estimated model parameters with the true
model parameters as follows.

``` r
names(opt_out$par) <- names(out$par)
true_params <- with(dat$params, c(
  gamma, B, cov_to_theta(Psi), cov_to_theta(sigma),
  delta, omega, alpha))
n_params <- length(true_params)
names(true_params) <- names(out$par)[seq_along(true_params)]
rbind(Estimate = opt_out$par[1:n_params], 
      `True value` = true_params)
#>            gamma:X1.Y1 gamma:X1.Y2 B:g1.Y1 B:g2.Y1 B:g3.Y1 B:g1.Y2 B:g2.Y2
#> Estimate         0.133      -0.802  -0.942   0.367   0.403   0.256  -0.843
#> True value       0.140      -0.800  -0.960   0.330   0.390   0.260  -0.760
#>            B:g3.Y2 Psi:L1.1 Psi:L2.1 Psi:L3.1 Psi:L4.1 Psi:L2.2 Psi:L3.2
#> Estimate    0.0686    0.236   -0.303  -0.1158   -0.225  -0.0412  -0.1171
#> True value  0.1900    0.226   -0.295  -0.0638   -0.136  -0.0567  -0.0729
#>            Psi:L4.2 Psi:L3.3 Psi:L4.3 Psi:L4.4 Sigma:L1.1 Sigma:L2.1 Sigma:L2.2
#> Estimate    -0.0226  -0.1241    0.555  -0.0875      -1.76    0.00302      -1.53
#> True value   0.0528  -0.0751    0.566  -0.0942      -1.75    0.00000      -1.50
#>            delta:Z1 delta:Z2 omega:b1 omega:b2 alpha:Y1 alpha:Y2
#> Estimate      0.237   -0.219    -2.51    -1.25    0.338   -0.277
#> True value    0.200   -0.170    -2.60    -1.32    0.320   -0.310
```

Next, we compare the estimated covariance matrix of the random effects
with the true
values.

``` r
# random effect covariance matrix (first estimated and then the true values)
is_psi <- which(grepl("Psi", names(true_params)))
theta_to_cov(opt_out$par[is_psi]) 
#>        [,1]    [,2]    [,3]    [,4]
#> [1,]  1.603 -0.3842 -0.1466 -0.2848
#> [2,] -0.384  1.0129 -0.0772  0.0466
#> [3,] -0.147 -0.0772  0.8074  0.5188
#> [4,] -0.285  0.0466  0.5188  1.1983
dat$params$Psi
#>       [,1]  [,2]  [,3]  [,4]
#> [1,]  1.57 -0.37 -0.08 -0.17
#> [2,] -0.37  0.98 -0.05  0.09
#> [3,] -0.08 -0.05  0.87  0.53
#> [4,] -0.17  0.09  0.53  1.17
cov2cor(theta_to_cov(opt_out$par[is_psi]))
#>        [,1]    [,2]    [,3]    [,4]
#> [1,]  1.000 -0.3015 -0.1289 -0.2055
#> [2,] -0.302  1.0000 -0.0854  0.0423
#> [3,] -0.129 -0.0854  1.0000  0.5274
#> [4,] -0.206  0.0423  0.5274  1.0000
cov2cor(dat$params$Psi)
#>         [,1]    [,2]    [,3]   [,4]
#> [1,]  1.0000 -0.2983 -0.0685 -0.125
#> [2,] -0.2983  1.0000 -0.0541  0.084
#> [3,] -0.0685 -0.0541  1.0000  0.525
#> [4,] -0.1254  0.0840  0.5253  1.000
```

Further, we compare the estimated covariance matrix of the noise with
the true values.

``` r
# noise covariance matrix (first estimated and then the true values)
is_sigma <- which(grepl("Sigma", names(true_params)))
theta_to_cov(opt_out$par[is_sigma])
#>          [,1]     [,2]
#> [1,] 0.029376 0.000517
#> [2,] 0.000517 0.047144
dat$params$sigma
#>      [,1] [,2]
#> [1,] 0.03 0.00
#> [2,] 0.00 0.05
cov2cor(theta_to_cov(opt_out$par[is_sigma]))
#>        [,1]   [,2]
#> [1,] 1.0000 0.0139
#> [2,] 0.0139 1.0000
cov2cor(dat$params$sigma)
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
```

We can look at quantiles of mean, standard deviations, and Pearson’s
moment coefficient of skewness for each individuals estimated
variational distribution as follows.

``` r
va_stats <- lapply(1:1000, function(id){
  is_grp_x <- which(grepl(paste0("^g", id, ":"), names(opt_out$par)))
  x_va_pars <- opt_out$par[is_grp_x]
  xi <- x_va_pars[grepl(":xi", names(x_va_pars))]
  Lambda <- theta_to_cov(
    x_va_pars[grepl(":(log_sd|L)", names(x_va_pars))])
  alpha <- x_va_pars[grepl(":alpha", names(x_va_pars))]
  
  dp_to_cp(xi = xi, Psi = Lambda, alpha = alpha)
})

sum_func <- function(x)
  apply(x, 2L, quantile, probs = seq(0, 1, by = .1))

# mean 
sum_func(do.call(rbind, lapply(va_stats, `[[`, "mu")))
#>       g1:xi1  g1:xi2   g1:xi3  g1:xi4
#> 0%   -3.4738 -3.6373 -2.37207 -2.7812
#> 10%  -1.4551 -0.9531 -0.82561 -0.9635
#> 20%  -1.0003 -0.6354 -0.48732 -0.5692
#> 30%  -0.6377 -0.3933 -0.28593 -0.3441
#> 40%  -0.3140 -0.1956 -0.13987 -0.1612
#> 50%  -0.0261 -0.0138 -0.00653 -0.0174
#> 60%   0.2929  0.1955  0.12421  0.1363
#> 70%   0.6008  0.3913  0.27279  0.3103
#> 80%   0.9946  0.6001  0.46744  0.5549
#> 90%   1.5250  1.0146  0.84856  0.9648
#> 100%  3.9340  2.9911  2.31430  2.6400

# standard deviation
sum_func(do.call(rbind, lapply(va_stats, 
                               function(x) sqrt(diag(x[["Sigma"]])))))
#>        [,1]  [,2]  [,3]  [,4]
#> 0%   0.0851 0.107 0.107 0.135
#> 10%  0.1286 0.236 0.160 0.295
#> 20%  0.1644 0.335 0.204 0.416
#> 30%  0.2265 0.462 0.280 0.569
#> 40%  0.3345 0.619 0.430 0.793
#> 50%  0.4549 0.659 0.607 0.858
#> 60%  0.5380 0.680 0.714 0.900
#> 70%  0.5984 0.704 0.785 0.936
#> 80%  0.6433 0.733 0.827 0.986
#> 90%  0.6816 0.773 0.855 1.034
#> 100% 0.7246 0.998 0.878 1.080

# skewness
skews <-  sum_func(do.call(rbind, lapply(va_stats, `[[`, "gamma")))
skews[] <- sprintf("%8.4f", skews)
print(skews, quote = FALSE)
#>      [,1]     [,2]     [,3]     [,4]    
#> 0%    -0.0122  -0.0126  -0.0041  -0.0041
#> 10%   -0.0042  -0.0041  -0.0003  -0.0003
#> 20%   -0.0031  -0.0030  -0.0002  -0.0002
#> 30%   -0.0025  -0.0024  -0.0002  -0.0002
#> 40%   -0.0022  -0.0022  -0.0002  -0.0002
#> 50%   -0.0021  -0.0021  -0.0001  -0.0001
#> 60%   -0.0020  -0.0020  -0.0001  -0.0001
#> 70%   -0.0018  -0.0018  -0.0001  -0.0001
#> 80%   -0.0012  -0.0011  -0.0001  -0.0000
#> 90%   -0.0004  -0.0003  -0.0000  -0.0000
#> 100%   0.0000   0.0001   0.0001   0.0001
```

We only see a low amount of skewness.

``` r
rm(list = ls())
gc()
#>           used (Mb) gc trigger  (Mb) max used  (Mb)
#> Ncells 2183618  117    4304134 229.9  2827218 151.0
#> Vcells 3791940   29    8412869  64.2  8412869  64.2
clear_cppad_mem(6L)
#> [1] 1
```

## Pedigree Data

<!-- $$\begin{align*}  -->

<!-- g(S(t\mid \vec x_{ij}, \epsilon_{ij})) &=  -->

<!--   \vec\omega^\top\vec f(t) + \vec\beta^\top\vec x_{ij} +\epsilon_{ij} \\ -->

<!-- \vec\epsilon_i &= (\epsilon_{i1}, \dots, \epsilon_{in_i})^\top \sim  -->

<!--   N^{(n_i)}\left(\vec 0, \sum_{l = 1}^K\sigma_l^2 C_{il} -->

<!--   \right) -->

<!-- \end{align*}$$ -->

The package contains an implementation of models which can used to
estimate heritability using pedigree data. These are GSMs of the
following form

  
![\\begin{align\*} g(S(t\\mid \\vec x\_{ij}, \\epsilon\_{ij})) &=
\\vec\\omega^\\top\\vec f(t) + \\vec\\beta^\\top\\vec x\_{ij}
+\\epsilon\_{ij} \\\\\\vec\\epsilon\_i &= (\\epsilon\_{i1}, \\dots,
\\epsilon\_{in\_i})^\\top \\sim N^{(n\_i)}\\left(\\vec 0, \\sum\_{l
= 1}^K\\sigma\_l^2 C\_{il}
\\right)\\end{align\*}](https://latex.codecogs.com/svg.latex?%5Cbegin%7Balign%2A%7D%20g%28S%28t%5Cmid%20%5Cvec%20x_%7Bij%7D%2C%20%5Cepsilon_%7Bij%7D%29%29%20%26%3D%20%20%20%5Cvec%5Comega%5E%5Ctop%5Cvec%20f%28t%29%20%2B%20%5Cvec%5Cbeta%5E%5Ctop%5Cvec%20x_%7Bij%7D%20%2B%5Cepsilon_%7Bij%7D%20%5C%5C%5Cvec%5Cepsilon_i%20%26%3D%20%28%5Cepsilon_%7Bi1%7D%2C%20%5Cdots%2C%20%5Cepsilon_%7Bin_i%7D%29%5E%5Ctop%20%5Csim%20%20%20N%5E%7B%28n_i%29%7D%5Cleft%28%5Cvec%200%2C%20%5Csum_%7Bl%20%3D%201%7D%5EK%5Csigma_l%5E2%20C_%7Bil%7D%20%20%5Cright%29%5Cend%7Balign%2A%7D
"\\begin{align*} g(S(t\\mid \\vec x_{ij}, \\epsilon_{ij})) &=   \\vec\\omega^\\top\\vec f(t) + \\vec\\beta^\\top\\vec x_{ij} +\\epsilon_{ij} \\\\\\vec\\epsilon_i &= (\\epsilon_{i1}, \\dots, \\epsilon_{in_i})^\\top \\sim   N^{(n_i)}\\left(\\vec 0, \\sum_{l = 1}^K\\sigma_l^2 C_{il}  \\right)\\end{align*}")  

where ![g](https://latex.codecogs.com/svg.latex?g "g") is a given link
function, ![\\vec f](https://latex.codecogs.com/svg.latex?%5Cvec%20f
"\\vec f") is a given function, the
![\\epsilon\_{ij}](https://latex.codecogs.com/svg.latex?%5Cepsilon_%7Bij%7D
"\\epsilon_{ij}")s are individual specific random effects, and the
![K](https://latex.codecogs.com/svg.latex?K "K")
![C\_{il}](https://latex.codecogs.com/svg.latex?C_%7Bil%7D "C_{il}")
matrices are known. Various types of
![C\_{il}](https://latex.codecogs.com/svg.latex?C_%7Bil%7D "C_{il}")
matrices can be used. A typical example is to use a kinship matrix to
estimate genetic effect Other examples are to include maternal effects,
paternal effects, shared environment etc.

As an example, we will use the `pedigree.RDS` in the
[inst/test-data](inst/test-data) directory.

``` r
# load the data
library(survTMB)
dat <- readRDS(file.path("inst", "test-data", "pedigree.RDS"))

# prepare the cluster data
c_data <- lapply(dat$sim_data, function(x){
  data <- data.frame(Z = x$Z, y = x$y, event = x$event)
  cor_mats <- list(x$rel_mat)
  list(data = data, cor_mats = cor_mats)
})
length(c_data)    # number of clusters/families
#> [1] 100
str(c_data[[1L]]) # example with the first cluster/family
#> List of 2
#>  $ data    :'data.frame':    32 obs. of  4 variables:
#>   ..$ Z.1  : num [1:32] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ Z.2  : num [1:32] -0.04471 -1.73322 0.00213 -0.6303 -0.34097 ...
#>   ..$ y    : num [1:32] 6.58058 2.74164 4.72234 0.04791 0.00399 ...
#>   ..$ event: logi [1:32] FALSE TRUE FALSE TRUE TRUE FALSE ...
#>  $ cor_mats:List of 1
#>   ..$ : num [1:32, 1:32] 1 0.5 0.5 0.125 0.125 0.125 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:32] "8" "9" "10" "16" ...
#>   .. .. ..$ : chr [1:32] "8" "9" "10" "16" ...
sapply(c_data, function(x) NCOL(x$cor_mats[[1L]]))
#>   [1] 32 42 37 34 36 35 33 38 28 32 32 35  8 35  3 31 31 35 29 39 34 36 33  4  4
#>  [26] 35 38  3 33 34 34 35 36  5 35 29 30 39  1 33  2 34  1 38 12 35 35 33 33 29
#>  [51] 31 33 31 43 38 32 33 40 35 28 32 37 32 34  2 34 33  2  2 35 38 36 34 40 37
#>  [76] 38 30 13 30  2  5 33 35 35 38 37 37  4 35 33 33 32 35 33 35 31 37 36 31  6

# use a third order polynomial as in the true model
sbase_haz <- function(x){
  x <- log(x)
  cbind(x^3, x^2, x)
}

# create ADFun
system.time(
  func <- make_pedigree_ADFun(
    formula = Surv(y, event) ~ Z.1 + Z.2 - 1,
    tformula = ~ sbase_haz(y) - 1, trace = TRUE,
    c_data = c_data, link = "probit", n_threads = 6L))
#> Finding starting values for fixed effects...
#> Maximum log-likelihood without random effects is: -1583.452
#> Creating ADFun...
#> Finding starting values for variational parameters...
#>    user  system elapsed 
#>   3.289   0.056   0.792

-func$fn(func$par) # lower bound of the log-likelihood
#> [1] -1575

# check memory usage
# TODO: export the function
siz <- survTMB:::pedigree_get_size(environment(func$fn)$adfun)

# the amount of work and memory necessary for computing function values 
# and derivatives using f is roughly proportional to ...
siz$size_var 
#> [1] 1509225

# the number of parameters in the operation sequence
siz$size_par 
#> [1] 20429

# optimize and compare the results with the true parameters
library(lbfgs)
system.time(
  opt_out <- lbfgs(func$fn, func$gr, func$par, m = 10, 
                   max_iterations = 5000L, invisible = 1))
#>    user  system elapsed 
#> 1037.43    0.06  173.00
```

We show the estimates below and compare them with the true values.

``` r
-opt_out$value # lower bound on the marginal log-likelihood in the end
#> [1] -1550

rbind(
  `Starting values` = head(func$par, 6), 
  Estimates = head(opt_out$par, 6),
  `True values` = c(dat$omega, dat$beta, log(dat$sds)))
#>                 omega:sbase_haz(y) omega:sbase_haz(y) omega:sbase_haz(y)x
#> Starting values             0.0150             0.0760               0.140
#> Estimates                   0.0197             0.0996               0.184
#> True values                 0.0200             0.1000               0.175
#>                 beta:Z.1 beta:Z.2 log_sds1
#> Starting values   -0.792    0.217   -0.693
#> Estimates         -1.029    0.278   -0.179
#> True values       -1.000    0.250   -0.223

# check the skew parameters
names(opt_out$par) <- names(func$par)
reg_exp <- "(^g\\d+)(:.+$)"
va_ests <- opt_out$par[grepl(reg_exp, names(func$par), perl = TRUE)]
grp <- gsub(reg_exp, "\\1", names(va_ests), perl = TRUE)
cps <- tapply(va_ests, grp, function(x){
  n <- length(x)
  n <- .5 * (sqrt(8 * n + 25) - 5)
  dp_to_cp(xi = head(x, n), 
           Psi = theta_to_cov(tail(head(x, -n), -n)), 
           alpha = tail(x, n))
})

summary(unlist(lapply(cps, `[[`, "gamma")))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -0.0150 -0.0009  0.0000 -0.0002  0.0001  0.0345
```

## References

<div id="refs" class="references">

<div id="ref-Bell19">

Bell, B. 2019. *CppAD: A Package for C++ Algorithmic Differentiation*.
<http://www.coin-or.org/CppAD>.

</div>

<div id="ref-Kristensen16">

Kristensen, Kasper, Anders Nielsen, Casper W. Berg, Hans Skaug, and
Bradley M. Bell. 2016. “TMB: Automatic Differentiation and Laplace
Approximation.” *Journal of Statistical Software* 70 (5): 1–21.
<https://doi.org/10.18637/jss.v070.i05>.

</div>

<div id="ref-Liu16">

Liu, Xing-Rong, Yudi Pawitan, and Mark Clements. 2016. “Parametric and
Penalized Generalized Survival Models.” *Statistical Methods in Medical
Research* 27 (5): 1531–46. <https://doi.org/10.1177/0962280216664760>.

</div>

<div id="ref-Liu17">

Liu, Xing-Rong, Yudi Pawitan, and Mark S. Clements. 2017. “Generalized
Survival Models for Correlated Time-to-Event Data.” *Statistics in
Medicine* 36 (29): 4743–62. <https://doi.org/10.1002/sim.7451>.

</div>

<div id="ref-Ormerod11">

Ormerod, J. T. 2011. “Skew-Normal Variational Approximations for
Bayesian Inference.” *Unpublished Article*.

</div>

<div id="ref-Ormerod12">

Ormerod, J. T., and M. P. Wand. 2012. “Gaussian Variational Approximate
Inference for Generalized Linear Mixed Models.” *Journal of
Computational and Graphical Statistics* 21 (1). Taylor & Francis: 2–17.
<https://doi.org/10.1198/jcgs.2011.09118>.

</div>

</div>
