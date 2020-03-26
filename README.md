
survTMB
=======

TODO: Write intro to the package...

Example
-------

TODO: write description

``` r
dat <- coxme::eortc

library(survTMB)
fit_model <- function(link){
  adfun <- make_gsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), 
    Z = ~ trt,
    df = 3L, data = dat, link = link, do_setup = "Laplace")
  fit <- do.call(optim, adfun$laplace)
  list(logLik = -fit$value, par = fit$par)
}

fit_model("PH"    )
#> Loading required package: rstpm2
#> Loading required package: survival
#> Loading required package: splines
#> 
#> Attaching package: 'rstpm2'
#> The following object is masked from 'package:survival':
#> 
#>     colon
#> Loading required package: TMB
#> $logLik
#> [1] -13026.68
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                              -7.8232329                               0.7198517 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                               5.3941547                              11.3890933 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                               4.7992329                              -1.8429549 
#>                                   theta                                   theta 
#>                              -2.2318986                               1.7952915
fit_model("PO"    )
#> $logLik
#> [1] -13031.16
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                               -8.049407                                1.029742 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                5.697061                               11.821001 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                5.593478                               -1.584165 
#>                                   theta                                   theta 
#>                               -1.869971                                1.947962
fit_model("probit")
#> $logLik
#> [1] -13035.14
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                              -3.7399189                               0.5992700 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                               2.6601204                               5.0043914 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                               2.9722065                              -1.9988364 
#>                                   theta                                   theta 
#>                              -1.8502383                               0.7924196
```

``` r
library(microbenchmark)
microbenchmark(
  PH     = fit_model("PH"    ), 
  PO     = fit_model("PO"    ), 
  probit = fit_model("probit"), times = 5)
#> Unit: milliseconds
#>    expr       min        lq      mean    median        uq       max neval
#>      PH  857.3811  858.7811  878.1654  880.1165  888.5134  906.0347     5
#>      PO 1245.5381 1245.9614 1257.7621 1251.1682 1257.1723 1288.9706     5
#>  probit 1598.4720 1604.4196 1614.1876 1613.0919 1613.6504 1641.3039     5
```
