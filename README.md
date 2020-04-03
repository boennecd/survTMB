
survTMB
=======

This package contains methods to estimated mixed generalized survival models (Liu, Pawitan, and Clements 2016; Liu, Pawitan, and Clements 2017). All methods use automatic differentiation using the CppAD library (B. Bell 2019) through the TMB package (Kristensen et al. 2016). The estimation methods are

-   a Laplace approximation using TMB.
-   Gaussian variational approximation (GVA) similar to the method shown by Ormerod and Wand (2012).
-   Skew-normal variational approximation (SNVA) similar to the method shown by Ormerod (2011).

The [example](#example) section shows an example of how to use the package with different methods. The [benchmark](#benchmark) section shows a comparison of the computation time of the methods.

Example
-------

We estimate a GSM a below with the proportional odds (PO) link function using both a Laplace approximation, a GVA, and a SNVA. First, we define a function to perform the estimation.

``` r
# assign variable with data 
dat <- coxme::eortc

# assign function to estimate the model
library(survTMB)
library(survival)
fit_model <- function(link, n_threads = 2L, method = "Laplace", 
                      param_type = "DP", with_hess = FALSE){
  eval(bquote({
    adfun <- make_gsm_ADFun(
      Surv(y, uncens) ~ trt, cluster = as.factor(center), 
      Z = ~ trt, df = 3L, data = dat, link = .(link), do_setup = .(method), 
      n_threads = .(n_threads), param_type = .(param_type), n_nodes = 15L, 
      dense_hess = .(with_hess), sparse_hess = .(with_hess))
    fit <- fit_mgsm(adfun, method = .(method))
    list(fit = fit, fun = adfun)
  }), parent.frame())
}

# estimate the model using different methods. Start w/ Laplace
(lap_ph <- fit_model("PO"))$fit
#> 
#> GSM estimated with method 'Laplace' with link 'PO' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "Laplace", 
#>       n_nodes = 15L, param_type = "DP", link = "PO", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "Laplace")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                   -7.62                                    1.11 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                    5.37                                   11.00 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                    5.39 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)      trt       (Intercept)     trt
#> (Intercept)     0.09336 -0.00461            0.3055 -0.0386
#> trt            -0.00461  0.15272           -0.0386  0.3908
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated log-likelihood is -13033.04

# w/ GVA
(gva_fit <- fit_model("PO", method = "GVA"))$fit
#> 
#> GSM estimated with method 'GVA' with link 'PO' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "GVA", 
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
#>             (Intercept)   trt       (Intercept)   trt
#> (Intercept)       0.048 0.054             0.219 0.702
#> trt               0.054 0.123             0.702 0.351
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13031.11

# the GVA solution seems to be better. First, we print a log-likelihood 
# approximation using the Laplace approximation at the GVA solution
print(-lap_ph$fun$laplace$fn(gva_fit$fit$params), digits = 7)
#> [1] -13031.17
#> attr(,"logarithm")
#> [1] TRUE
# then we print the log-likelihood approximation at the Laplace solution
-lap_ph$fit$optim$value
#> [1] -13033

# w/ SNVA
fit_model("PO", method = "SNVA", param_type = "DP")$fit
#> 
#> GSM estimated with method 'SNVA' with link 'PO' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
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
#> (Intercept)      0.0474 0.0559             0.218 0.739
#> trt              0.0559 0.1207             0.739 0.347
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13031.11
```

### Computing the Hessian

The Hessian using a variational approximation (VA) can be computed as both a dense and as sparse matrix. We show an example below where we compare the two approaches.

``` r
library(microbenchmark) # needed for benchmarking
```

``` r
# fit model w/ GVA
fit <- fit_model("PO", method = "GVA", with_hess = TRUE)

# compute dense Hessian
par <- with(fit$fit, c(params, va_params))
dense_hess <- fit$fun$gva$he(par)

# has many zeros and is sparse
mean(abs(dense_hess) > 0) # fraction of non-zeros
#> [1] 0.105

# plot non-zero entries (zero block's are non-zero; ignore upper triangle)
par(mar = c(1, 1, 1, 1))
is_non_zero <- t(abs(dense_hess) > 0)
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
#>                                   0.162                                   0.466 
#>                                   theta                                   theta 
#>                                   1.230                                   2.116

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
#>                Compute dense Hessian 144.46 144.77 148.86 147.85 151.54 158.54
#>               Compute sparse Hessian  18.91  19.47  20.00  19.83  20.38  21.65
#>         Invert dense Hessian (naive)   5.16   5.32   5.35   5.35   5.41   5.49
#>        Invert sparse Hessian (naive)   1.12   1.14   1.27   1.26   1.36   1.43
#>   Invert dense Hessian (alternative)   1.32   1.33   1.39   1.37   1.45   1.49
#>  Invert sparse Hessian (alternative)   2.78   2.94   3.03   2.99   3.18   3.27
#>  neval
#>     10
#>     10
#>     10
#>     10
#>     10
#>     10

# compare storage cost
as.numeric(object.size(dense_hess) / object.size(sparse_hess))
#> [1] 5.05
```

The sparse matrix only becomes more favorable for larger data sets (that is, in terms of the number of clusters).

### Other link functions

We estimate the same model below with other link functions.

``` r
######
# w/ Laplace
fit_model("PH"    )$fit
#> 
#> GSM estimated with method 'Laplace' with link 'PH' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "Laplace", 
#>       n_nodes = 15L, param_type = "DP", link = "PH", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "Laplace")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -7.439                                   0.752 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   5.110                                  10.636 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   4.695 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)     trt       (Intercept)    trt
#> (Intercept)      0.0579 -0.0132             0.241 -0.190
#> trt             -0.0132  0.0837            -0.190  0.289
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated log-likelihood is -13028.98
fit_model("PO"    )$fit
#> 
#> GSM estimated with method 'Laplace' with link 'PO' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "Laplace", 
#>       n_nodes = 15L, param_type = "DP", link = "PO", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "Laplace")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                   -7.62                                    1.11 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                    5.37                                   11.00 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                    5.39 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)      trt       (Intercept)     trt
#> (Intercept)     0.09336 -0.00461            0.3055 -0.0386
#> trt            -0.00461  0.15272           -0.0386  0.3908
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated log-likelihood is -13033.04
fit_model("probit")$fit
#> 
#> GSM estimated with method 'Laplace' with link 'probit' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "Laplace", 
#>       n_nodes = 15L, param_type = "DP", link = "probit", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "Laplace")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -3.630                                   0.587 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   2.647                                   4.714 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   2.988 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)     trt       (Intercept)    trt
#> (Intercept)      0.0499 -0.0153             0.223 -0.204
#> trt             -0.0153  0.1115            -0.204  0.334
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated log-likelihood is -13039.10

######
# w/ GVA
fit_model("PH"    , method = "GVA")$fit
#> 
#> GSM estimated with method 'GVA' with link 'PH' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "GVA", 
#>       n_nodes = 15L, param_type = "DP", link = "PH", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "GVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -7.827                                   0.723 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   5.394                                  11.390 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   4.799 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0280 0.0259             0.167 0.660
#> trt              0.0259 0.0550             0.660 0.235
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13026.74
fit_model("PO"    , method = "GVA")$fit
#> 
#> GSM estimated with method 'GVA' with link 'PO' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "GVA", 
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
#>             (Intercept)   trt       (Intercept)   trt
#> (Intercept)       0.048 0.054             0.219 0.702
#> trt               0.054 0.123             0.702 0.351
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13031.11
fit_model("probit", method = "GVA")$fit
#> 
#> GSM estimated with method 'GVA' with link 'probit' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "GVA", 
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
#> (Intercept)      0.0195 0.0144             0.140 0.495
#> trt              0.0144 0.0437             0.495 0.209
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13035.15

######
# w/ SNVA (DP: direct parameterization)
fit_model("PH"    , method = "SNVA", param_type = "DP")$fit
#> 
#> GSM estimated with method 'SNVA' with link 'PH' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       n_nodes = 15L, param_type = "DP", link = "PH", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -7.827                                   0.723 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   5.394                                  11.390 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   4.799 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0271 0.0262             0.165 0.675
#> trt              0.0262 0.0555             0.675 0.236
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13026.74
fit_model("PO"    , method = "SNVA", param_type = "DP")$fit
#> 
#> GSM estimated with method 'SNVA' with link 'PO' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
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
#> (Intercept)      0.0474 0.0559             0.218 0.739
#> trt              0.0559 0.1207             0.739 0.347
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13031.11
fit_model("probit", method = "SNVA", param_type = "DP")$fit
#> 
#> GSM estimated with method 'SNVA' with link 'probit' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       n_nodes = 15L, param_type = "DP", link = "probit", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -3.747                                   0.606 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   2.661                                   5.006 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   2.974 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)     trt       (Intercept)   trt
#> (Intercept)     0.02728 0.00668             0.165 0.179
#> trt             0.00668 0.05107             0.179 0.226
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13035.31

######
# w/ SNVA (CP: centralized parameterization)
fit_model("PH"    , method = "SNVA", param_type = "CP_trans")$fit
#> 
#> GSM estimated with method 'SNVA' with link 'PH' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
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
#> (Intercept)      0.0294 0.0257             0.172 0.639
#> trt              0.0257 0.0549             0.639 0.234
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13026.70
fit_model("PO"    , method = "SNVA", param_type = "CP_trans")$fit
#> 
#> GSM estimated with method 'SNVA' with link 'PO' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       n_nodes = 15L, param_type = "CP_trans", link = "PO", n_threads = 2L, 
#>       dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                   -8.06                                    1.04 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                    5.70                                   11.82 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                    5.60 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0624 0.0373             0.250 0.396
#> trt              0.0373 0.1422             0.396 0.377
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13031.23
fit_model("probit", method = "SNVA", param_type = "CP_trans")$fit
#> 
#> GSM estimated with method 'SNVA' with link 'probit' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       n_nodes = 15L, param_type = "CP_trans", link = "probit", 
#>       n_threads = 2L, dense_hess = FALSE, sparse_hess = FALSE)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -3.745                                   0.605 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   2.660                                   5.005 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   2.973 
#> 
#> Estimated random effect covariance matrix (correlation matrix) is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0236 0.0108             0.154 0.325
#> trt              0.0108 0.0469             0.325 0.217
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13035.20
```

Benchmark
---------

We provide a benchmark of the estimation methods used in section [example](#example) below.

``` r
for(mth in c("Laplace", "GVA")){
  msg <- sprintf("Method: %s", mth)
  cat(sprintf("\n%s\n%s\n", msg, 
              paste0(rep("-", nchar(msg)), collapse = "")))
  print(microbenchmark(
    `PH         ` = fit_model("PH"    , 1L, mth),
    `PH     (2L)` = fit_model("PH"    , 2L, mth),
    `PH     (4L)` = fit_model("PH"    , 4L, mth),
    
    `PO         ` = fit_model("PO"    , 1L, mth),
    `PO     (2L)` = fit_model("PO"    , 2L, mth),
    `PO     (4L)` = fit_model("PO"    , 4L, mth), 
    
    `probit     ` = fit_model("probit", 1L, mth),
    `probit (2L)` = fit_model("probit", 2L, mth),
    `probit (4L)` = fit_model("probit", 4L, mth),
    times = 5))
}
#> 
#> Method: Laplace
#> ---------------
#> Unit: milliseconds
#>         expr  min   lq mean median   uq  max neval
#>  PH          1025 1044 1044   1048 1049 1054     5
#>  PH     (2L)  647  653  658    658  664  669     5
#>  PH     (4L)  484  484  485    485  487  487     5
#>  PO          1606 1617 1619   1621 1626 1626     5
#>  PO     (2L)  995 1008 1010   1012 1015 1017     5
#>  PO     (4L)  701  701  714    718  724  727     5
#>  probit      1003 1006 1042   1010 1074 1117     5
#>  probit (2L)  619  619  622    620  626  626     5
#>  probit (4L)  431  435  452    443  466  484     5
#> 
#> Method: GVA
#> -----------
#> Unit: milliseconds
#>         expr  min   lq mean median   uq  max neval
#>  PH           272  272  274    273  274  277     5
#>  PH     (2L)  181  182  186    182  183  204     5
#>  PH     (4L)  142  143  144    144  145  147     5
#>  PO           848  850  862    855  867  892     5
#>  PO     (2L)  511  516  524    525  533  533     5
#>  PO     (4L)  352  354  359    355  355  375     5
#>  probit      1434 1440 1442   1441 1444 1453     5
#>  probit (2L)  847  849  857    849  850  888     5
#>  probit (4L)  566  566  574    570  571  598     5
```

``` r
for(param_type in c("DP", "CP_trans")){
  mth <- "SNVA"
  msg <- sprintf("Method: %s (%s)", mth, param_type)
  cat(sprintf("\n%s\n%s\n", msg, 
              paste0(rep("-", nchar(msg)), collapse = "")))
  print(microbenchmark(
    `PH         ` = fit_model("PH"    , 1L, mth, param_type = param_type),
    `PH     (2L)` = fit_model("PH"    , 2L, mth, param_type = param_type),
    `PH     (4L)` = fit_model("PH"    , 4L, mth, param_type = param_type),
    
    `PO         ` = fit_model("PO"    , 1L, mth, param_type = param_type),
    `PO     (2L)` = fit_model("PO"    , 2L, mth, param_type = param_type),
    `PO     (4L)` = fit_model("PO"    , 4L, mth, param_type = param_type), 
    
    `probit     ` = fit_model("probit", 1L, mth, param_type = param_type),
    `probit (2L)` = fit_model("probit", 2L, mth, param_type = param_type),
    `probit (4L)` = fit_model("probit", 4L, mth, param_type = param_type),
    times = 5))
}
#> 
#> Method: SNVA (DP)
#> -----------------
#> Unit: milliseconds
#>         expr  min   lq mean median   uq  max neval
#>  PH           454  457  462    461  463  477     5
#>  PH     (2L)  290  293  295    295  298  300     5
#>  PH     (4L)  214  214  223    215  223  250     5
#>  PO          3073 3091 3122   3106 3138 3202     5
#>  PO     (2L) 1808 1838 1916   1853 1990 2089     5
#>  PO     (4L) 1197 1199 1210   1208 1219 1228     5
#>  probit      3247 3249 3274   3267 3286 3322     5
#>  probit (2L) 1881 1932 1986   2019 2032 2065     5
#>  probit (4L) 1255 1262 1279   1280 1293 1307     5
#> 
#> Method: SNVA (CP_trans)
#> -----------------------
#> Unit: milliseconds
#>         expr  min   lq mean median   uq  max neval
#>  PH           534  544  552    548  565  569     5
#>  PH     (2L)  369  372  377    375  378  390     5
#>  PH     (4L)  267  267  271    268  268  285     5
#>  PO          5664 5738 5830   5779 5899 6071     5
#>  PO     (2L) 1587 1621 1655   1639 1671 1760     5
#>  PO     (4L) 1009 1015 1059   1020 1043 1209     5
#>  probit      5780 5802 5875   5852 5856 6086     5
#>  probit (2L) 3171 3433 3392   3447 3452 3457     5
#>  probit (4L) 2224 2231 2340   2329 2428 2488     5
```

References
----------

Bell, B. 2019. *CppAD: A Package for C++ Algorithmic Differentiation*. <http://www.coin-or.org/CppAD>.

Kristensen, Kasper, Anders Nielsen, Casper W. Berg, Hans Skaug, and Bradley M. Bell. 2016. “TMB: Automatic Differentiation and Laplace Approximation.” *Journal of Statistical Software* 70 (5): 1–21. doi:[10.18637/jss.v070.i05](https://doi.org/10.18637/jss.v070.i05).

Liu, Xing-Rong, Yudi Pawitan, and Mark Clements. 2016. “Parametric and Penalized Generalized Survival Models.” *Statistical Methods in Medical Research* 27 (5): 1531–46. doi:[10.1177/0962280216664760](https://doi.org/10.1177/0962280216664760).

Liu, Xing-Rong, Yudi Pawitan, and Mark S. Clements. 2017. “Generalized Survival Models for Correlated Time-to-Event Data.” *Statistics in Medicine* 36 (29): 4743–62. doi:[10.1002/sim.7451](https://doi.org/10.1002/sim.7451).

Ormerod, J. T. 2011. “Skew-Normal Variational Approximations for Bayesian Inference.” *Unpublished Article*.

Ormerod, J. T., and M. P. Wand. 2012. “Gaussian Variational Approximate Inference for Generalized Linear Mixed Models.” *Journal of Computational and Graphical Statistics* 21 (1). Taylor & Francis: 2–17. doi:[10.1198/jcgs.2011.09118](https://doi.org/10.1198/jcgs.2011.09118).
