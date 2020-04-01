
survTMB
=======

This package contains methods to estimated mixed generalized survival models (Liu, Pawitan, and Clements 2016; Liu, Pawitan, and Clements 2017). All methods use automatic differentiation using the CppAD library (B. Bell 2019) through the TMB package (Kristensen et al. 2016). The estimation methods are

-   A Laplace approximation method using TMB.
-   Gaussian variational approximation (GVA) similar to the method shown by Ormerod and Wand (2012).
-   Skew-normal variational approximation (SNVA) similar to the method shown by Ormerod (2011).

The [example](#example) section shows an example of how to use the package with different methods. The [benchmark](#benchmark) section shows a comparison of the computation time of the methods.

Example
-------

We estimate a GSM using different link functions and different methods below. First, we define function to perform the estimation. Then we use the different methods with models using various link functions.

``` r
dat <- coxme::eortc

library(survTMB)
library(survival)
fit_model <- function(link, n_threads = 2L, method = "Laplace", 
                      param_type = "DP"){
  eval(bquote({
    adfun <- make_gsm_ADFun(
      Surv(y, uncens) ~ trt, cluster = as.factor(center), 
      Z = ~ trt, df = 3L, data = dat, link = .(link), do_setup = .(method), 
      n_threads = .(n_threads), param_type = .(param_type))
    fit <- fit_mgsm(adfun, method = .(method))
    list(fit = fit, fun = adfun)
  }), parent.frame())
}

######
# w/ Laplace
(lap_ph <- fit_model("PH"    ))$fit
#> 
#> GSM estimated with method 'Laplace' with link 'PH' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "Laplace", 
#>       param_type = "DP", link = "PH", n_threads = 2L)
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
#> Estimated random effect covariance matrix (correlation matrix is:
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
#>       param_type = "DP", link = "PO", n_threads = 2L)
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
#> Estimated random effect covariance matrix (correlation matrix is:
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
#>       param_type = "DP", link = "probit", n_threads = 2L)
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
#> Estimated random effect covariance matrix (correlation matrix is:
#>             (Intercept)     trt       (Intercept)    trt
#> (Intercept)      0.0499 -0.0153             0.223 -0.204
#> trt             -0.0153  0.1115            -0.204  0.334
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated log-likelihood is -13039.10

######
# w/ GVA
(gva_fit <- fit_model("PH"    , method = "GVA"))$fit
#> 
#> GSM estimated with method 'GVA' with link 'PH' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "GVA", 
#>       param_type = "DP", link = "PH", n_threads = 2L)
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
#> Estimated random effect covariance matrix (correlation matrix is:
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
#>       param_type = "DP", link = "PO", n_threads = 2L)
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
#> Estimated random effect covariance matrix (correlation matrix is:
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
#>       param_type = "DP", link = "probit", n_threads = 2L)
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
#> Estimated random effect covariance matrix (correlation matrix is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0195 0.0144             0.140 0.495
#> trt              0.0144 0.0437             0.495 0.209
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13035.15

# the GVA solutions seems to be better. First, we print a log-likelihood 
# approximation using the Laplace approximation at the GVA solution
print(-lap_ph$fun$laplace$fn(gva_fit$fit$params), digits = 7)
#> [1] -13026.71
#> attr(,"logarithm")
#> [1] TRUE
# then we print the log-likelihood approximation at the Laplace solution
-lap_ph$fit$optim$value
#> [1] -13029

######
# w/ SNVA (DP: direct parameterization)
fit_model("PH"    , method = "SNVA", param_type = "DP")$fit
#> 
#> GSM estimated with method 'SNVA' with link 'PH' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       param_type = "DP", link = "PH", n_threads = 2L)
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
#> Estimated random effect covariance matrix (correlation matrix is:
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
#>       param_type = "DP", link = "PO", n_threads = 2L)
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
#> Estimated random effect covariance matrix (correlation matrix is:
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
#>       param_type = "DP", link = "probit", n_threads = 2L)
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
#> Estimated random effect covariance matrix (correlation matrix is:
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
#>       param_type = "CP_trans", link = "PH", n_threads = 2L)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -7.829                                   0.725 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   5.394                                  11.390 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   4.800 
#> 
#> Estimated random effect covariance matrix (correlation matrix is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0315 0.0236             0.177 0.557
#> trt              0.0236 0.0568             0.557 0.238
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13026.73
fit_model("PO"    , method = "SNVA", param_type = "CP_trans")$fit
#> 
#> GSM estimated with method 'SNVA' with link 'PO' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       param_type = "CP_trans", link = "PO", n_threads = 2L)
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
#> Estimated random effect covariance matrix (correlation matrix is:
#>             (Intercept)    trt       (Intercept)   trt
#> (Intercept)      0.0617 0.0375             0.248 0.402
#> trt              0.0375 0.1409             0.402 0.375
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13031.23
fit_model("probit", method = "SNVA", param_type = "CP_trans")$fit
#> 
#> GSM estimated with method 'SNVA' with link 'probit' from call:
#>   make_gsm_ADFun(formula = Surv(y, uncens) ~ trt, data = dat, df = 3L, 
#>       Z = ~trt, cluster = as.factor(center), do_setup = "SNVA", 
#>       param_type = "CP_trans", link = "probit", n_threads = 2L)
#>   fit_mgsm(object = adfun, method = "SNVA")
#> 
#> Estimated fixed effects:
#>                             (Intercept)                                     trt 
#>                                  -3.745                                   0.606 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   2.660                                   5.003 
#> nsx(log(y), df = 3, intercept = FALSE)3 
#>                                   2.973 
#> 
#> Estimated random effect covariance matrix (correlation matrix is:
#>             (Intercept)     trt       (Intercept)   trt
#> (Intercept)     0.02478 0.00967             0.157 0.282
#> trt             0.00967 0.04758             0.282 0.218
#> (standard deviations are in the diagonal of the correlation matrix)
#> 
#> Estimated lower bound is -13035.31
```

Benchmark
---------

We provide a benchmark of the estimation methods used in section [example](#example) below.

``` r
library(microbenchmark)
```

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
#>  PH           832  838  848    840  858  869     5
#>  PH     (2L)  546  555  558    559  559  573     5
#>  PH     (4L)  414  419  422    421  426  430     5
#>  PO          1346 1362 1366   1363 1373 1386     5
#>  PO     (2L)  860  861  871    867  879  886     5
#>  PO     (4L)  609  633  641    634  658  673     5
#>  probit       902  903  909    912  914  916     5
#>  probit (2L)  562  565  574    572  582  590     5
#>  probit (4L)  403  408  410    411  411  418     5
#> 
#> Method: GVA
#> -----------
#> Unit: milliseconds
#>         expr  min   lq mean median   uq  max neval
#>  PH           245  245  246    246  246  248     5
#>  PH     (2L)  167  167  170    169  170  176     5
#>  PH     (4L)  138  138  139    139  140  141     5
#>  PO           990  993  999    995 1006 1010     5
#>  PO     (2L)  592  592  593    594  594  594     5
#>  PO     (4L)  414  417  420    418  425  428     5
#>  probit      1776 1784 1788   1787 1788 1803     5
#>  probit (2L) 1037 1047 1053   1048 1056 1078     5
#>  probit (4L)  708  715  725    722  732  748     5
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
#>  PH           434  435  438    436  438  446     5
#>  PH     (2L)  302  304  309    308  314  320     5
#>  PH     (4L)  235  236  239    237  242  247     5
#>  PO          3655 3671 3697   3689 3693 3777     5
#>  PO     (2L) 2375 2564 2562   2592 2610 2667     5
#>  PO     (4L) 1596 1649 1706   1650 1743 1892     5
#>  probit      3939 3961 4003   3984 3993 4140     5
#>  probit (2L) 2474 2689 2716   2754 2776 2887     5
#>  probit (4L) 1675 1711 1728   1727 1754 1776     5
#> 
#> Method: SNVA (CP_trans)
#> -----------------------
#> Unit: milliseconds
#>         expr  min   lq mean median   uq  max neval
#>  PH           607  612  615    617  620  620     5
#>  PH     (2L)  319  325  328    331  332  332     5
#>  PH     (4L)  284  286  288    286  288  296     5
#>  PO          5574 5626 5628   5638 5645 5657     5
#>  PO     (2L) 2404 2498 2551   2584 2628 2638     5
#>  PO     (4L) 3126 3135 3242   3158 3248 3540     5
#>  probit      8704 8715 8762   8733 8791 8868     5
#>  probit (2L) 3249 3528 3559   3529 3599 3887     5
#>  probit (4L) 2443 2450 2570   2509 2575 2873     5
```

References
----------

Bell, B. 2019. *CppAD: A Package for C++ Algorithmic Differentiation*. <http://www.coin-or.org/CppAD>.

Kristensen, Kasper, Anders Nielsen, Casper W. Berg, Hans Skaug, and Bradley M. Bell. 2016. “TMB: Automatic Differentiation and Laplace Approximation.” *Journal of Statistical Software* 70 (5): 1–21. doi:[10.18637/jss.v070.i05](https://doi.org/10.18637/jss.v070.i05).

Liu, Xing-Rong, Yudi Pawitan, and Mark Clements. 2016. “Parametric and Penalized Generalized Survival Models.” *Statistical Methods in Medical Research* 27 (5): 1531–46. doi:[10.1177/0962280216664760](https://doi.org/10.1177/0962280216664760).

Liu, Xing-Rong, Yudi Pawitan, and Mark S. Clements. 2017. “Generalized Survival Models for Correlated Time-to-Event Data.” *Statistics in Medicine* 36 (29): 4743–62. doi:[10.1002/sim.7451](https://doi.org/10.1002/sim.7451).

Ormerod, J. T. 2011. “Skew-Normal Variational Approximations for Bayesian Inference.” *Unpublished Article*.

Ormerod, J. T., and M. P. Wand. 2012. “Gaussian Variational Approximate Inference for Generalized Linear Mixed Models.” *Journal of Computational and Graphical Statistics* 21 (1). Taylor & Francis: 2–17. doi:[10.1198/jcgs.2011.09118](https://doi.org/10.1198/jcgs.2011.09118).
