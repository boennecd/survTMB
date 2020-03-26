
survTMB
=======

TODO: Write intro to the package...

Example
-------

TODO: write description

``` r
dat <- coxme::eortc

library(survTMB)
library(survival)
fit_model <- function(link, n_threads = 2L){
  adfun <- make_gsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), 
    Z = ~ trt, df = 3L, data = dat, link = link, do_setup = "Laplace", 
    n_threads = n_threads)
  fit <- do.call(optim, adfun$laplace)
  list(logLik = -fit$value, par = fit$par)
}

fit_model("PH"    )
#> $logLik
#> [1] -13027
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                 -7.8232                                  0.7199 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                  5.3942                                 11.3891 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                  4.7992                                 -1.8430 
#>                                   theta                                   theta 
#>                                 -2.2319                                  1.7953
fit_model("PO"    )
#> $logLik
#> [1] -13031
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                  -8.049                                   1.030 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                   5.697                                  11.821 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                   5.593                                  -1.584 
#>                                   theta                                   theta 
#>                                  -1.870                                   1.948
fit_model("probit")
#> $logLik
#> [1] -13035
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                 -3.7399                                  0.5993 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                  2.6601                                  5.0044 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                  2.9722                                 -1.9988 
#>                                   theta                                   theta 
#>                                 -1.8502                                  0.7924
```

``` r
library(microbenchmark)
microbenchmark(
   PH           = fit_model("PH"    , 1L),
  `PH     (2L)` = fit_model("PH"    , 2L),
  `PH     (4L)` = fit_model("PH"    , 4L),
  
   PO           = fit_model("PO"    , 1L),
  `PO     (2L)` = fit_model("PO"    , 2L),
  `PO     (4L)` = fit_model("PO"    , 4L), 
  
  
   probit       = fit_model("probit", 1L),
  `probit (2L)` = fit_model("probit", 2L),
  `probit (4L)` = fit_model("probit", 4L),
  times = 5)
#> Unit: milliseconds
#>         expr    min     lq   mean median     uq    max neval
#>           PH  884.2  895.4  911.9  909.9  913.2  957.0     5
#>  PH     (2L)  568.5  570.8  609.2  571.4  640.1  695.1     5
#>  PH     (4L)  413.3  419.1  436.1  420.7  445.6  481.5     5
#>           PO 1265.5 1269.5 1294.2 1279.4 1286.2 1370.4     5
#>  PO     (2L)  781.9  783.6  811.8  786.9  840.2  866.4     5
#>  PO     (4L)  565.2  570.1  584.8  575.2  593.8  619.7     5
#>       probit 1627.5 1669.5 1703.8 1681.7 1766.8 1773.2     5
#>  probit (2L)  959.4  960.9  981.7  979.6  993.1 1015.8     5
#>  probit (4L)  689.9  716.0  732.4  724.2  734.7  797.3     5
```
