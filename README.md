
survTMB
=======

This package contains methods to estimated mixed generalized survival models (Liu, Pawitan, and Clements 2016; Liu, Pawitan, and Clements 2017). All methods use automatic differentiation using the CppAD library (B. Bell 2019) through the TMB package (Kristensen et al. 2016). The estimation methods are

-   A Laplace approximation method using TMB.
-   Gaussian variational approximation (GVA) similar to the method shown by Ormerod and Wand (2012).
-   Skew-normal variational approximation (SNVA) similar to the method shown by Ormerod (2011).

The [example](#example) section shows an example of how to use the package with different methods. The [benchmark](#benchmark) section shows a comparison of the computation time of the methods.

Example
-------

We estimate a GSM using different link functions and different methods below. First, we define a function to perform the estimation. Then we use the different methods with models using various link functions.

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
      n_threads = .(n_threads), param_type = .(param_type), n_nodes = 15L)
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
#>       n_nodes = 15L, param_type = "DP", link = "PH", n_threads = 2L)
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
#>       n_nodes = 15L, param_type = "DP", link = "PO", n_threads = 2L)
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
#>       n_nodes = 15L, param_type = "DP", link = "probit", n_threads = 2L)
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
#>       n_nodes = 15L, param_type = "DP", link = "PH", n_threads = 2L)
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
#>       n_nodes = 15L, param_type = "DP", link = "PO", n_threads = 2L)
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
#>       n_nodes = 15L, param_type = "DP", link = "probit", n_threads = 2L)
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
#>       n_nodes = 15L, param_type = "DP", link = "PH", n_threads = 2L)
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
#>       n_nodes = 15L, param_type = "DP", link = "PO", n_threads = 2L)
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
#>       n_nodes = 15L, param_type = "DP", link = "probit", n_threads = 2L)
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
#>       n_nodes = 15L, param_type = "CP_trans", link = "PH", n_threads = 2L)
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
#> Estimated random effect covariance matrix (correlation matrix is:
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
#>       n_nodes = 15L, param_type = "CP_trans", link = "PO", n_threads = 2L)
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
#>       n_threads = 2L)
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
#> Estimated random effect covariance matrix (correlation matrix is:
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
#>  PH           872  877  899    892  905  948     5
#>  PH     (2L)  568  575  580    578  583  595     5
#>  PH     (4L)  422  428  433    432  436  448     5
#>  PO          1380 1386 1393   1387 1397 1414     5
#>  PO     (2L)  895  898  917    911  913  967     5
#>  PO     (4L)  641  646  653    652  661  664     5
#>  probit       914  937  951    957  968  979     5
#>  probit (2L)  580  587  609    611  629  637     5
#>  probit (4L)  419  423  428    430  432  438     5
#> 
#> Method: GVA
#> -----------
#> Unit: milliseconds
#>         expr  min   lq mean median   uq  max neval
#>  PH           247  248  251    252  253  255     5
#>  PH     (2L)  169  169  171    171  172  172     5
#>  PH     (4L)  140  141  143    142  142  152     5
#>  PO           812  814  819    820  820  827     5
#>  PO     (2L)  496  497  500    501  502  502     5
#>  PO     (4L)  351  355  358    356  358  368     5
#>  probit      1399 1399 1405   1403 1408 1414     5
#>  probit (2L)  830  834  838    841  842  842     5
#>  probit (4L)  579  580  582    582  584  585     5
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
#>  PH           427  435  435    435  435  444     5
#>  PH     (2L)  279  279  285    283  290  293     5
#>  PH     (4L)  219  219  222    220  224  228     5
#>  PO          2888 2897 2931   2901 2954 3016     5
#>  PO     (2L) 1963 2100 2156   2170 2228 2319     5
#>  PO     (4L) 1259 1311 1312   1318 1320 1354     5
#>  probit      3066 3067 3092   3102 3111 3114     5
#>  probit (2L) 1979 2043 2149   2237 2241 2247     5
#>  probit (4L) 1313 1342 1354   1349 1370 1394     5
#> 
#> Method: SNVA (CP_trans)
#> -----------------------
#> Unit: milliseconds
#>         expr  min   lq mean median   uq  max neval
#>  PH           492  495  498    497  502  505     5
#>  PH     (2L)  349  352  353    354  355  358     5
#>  PH     (4L)  262  267  269    268  270  276     5
#>  PO          5263 5272 5336   5307 5320 5516     5
#>  PO     (2L) 1672 1679 1791   1800 1892 1910     5
#>  PO     (4L) 1421 1439 1467   1472 1479 1525     5
#>  probit      5429 5433 5454   5440 5445 5521     5
#>  probit (2L) 3194 3504 3577   3716 3720 3749     5
#>  probit (4L) 2432 2500 2544   2517 2567 2702     5
```

References
----------

Bell, B. 2019. *CppAD: A Package for C++ Algorithmic Differentiation*. <http://www.coin-or.org/CppAD>.

Kristensen, Kasper, Anders Nielsen, Casper W. Berg, Hans Skaug, and Bradley M. Bell. 2016. “TMB: Automatic Differentiation and Laplace Approximation.” *Journal of Statistical Software* 70 (5): 1–21. doi:[10.18637/jss.v070.i05](https://doi.org/10.18637/jss.v070.i05).

Liu, Xing-Rong, Yudi Pawitan, and Mark Clements. 2016. “Parametric and Penalized Generalized Survival Models.” *Statistical Methods in Medical Research* 27 (5): 1531–46. doi:[10.1177/0962280216664760](https://doi.org/10.1177/0962280216664760).

Liu, Xing-Rong, Yudi Pawitan, and Mark S. Clements. 2017. “Generalized Survival Models for Correlated Time-to-Event Data.” *Statistics in Medicine* 36 (29): 4743–62. doi:[10.1002/sim.7451](https://doi.org/10.1002/sim.7451).

Ormerod, J. T. 2011. “Skew-Normal Variational Approximations for Bayesian Inference.” *Unpublished Article*.

Ormerod, J. T., and M. P. Wand. 2012. “Gaussian Variational Approximate Inference for Generalized Linear Mixed Models.” *Journal of Computational and Graphical Statistics* 21 (1). Taylor & Francis: 2–17. doi:[10.1198/jcgs.2011.09118](https://doi.org/10.1198/jcgs.2011.09118).
