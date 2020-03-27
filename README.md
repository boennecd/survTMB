
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
fit_model <- function(link, n_threads = 2L, method = "Laplace"){
  adfun <- make_gsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), 
    Z = ~ trt, df = 3L, data = dat, link = link, do_setup = method, 
    n_threads = n_threads)
  fit <- 
    if(method == "Laplace")
      do.call(optim, adfun$laplace) else if(method == "GVA")
        do.call(optim, adfun$gva)
  
  out <- list(-fit$value, par = head(fit$par, 8))
  
  names(out)[1L] <- if(method == "Laplace")
    "logLik" else "lower bound"
  out
}

######
# w/ Laplace
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

######
# w/ GVA
fit_model("PH"    , method = "GVA")
#> $`lower bound`
#> [1] -13027
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                 -7.8269                                  0.7233 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                  5.3940                                 11.3896 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                  4.7994                                 -1.7872 
#>                                   theta                                   theta 
#>                                 -1.7361                                  0.8784
fit_model("PO"    , method = "GVA")
#> $`lower bound`
#> [1] -13031
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                 -8.0542                                  1.0345 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                  5.6972                                 11.8220 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                  5.5943                                 -1.5178 
#>                                   theta                                   theta 
#>                                 -1.3850                                  0.9845
fit_model("probit", method = "GVA")
#> $`lower bound`
#> [1] -13035
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                 -3.7418                                  0.6019 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                  2.6602                                  5.0039 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                  2.9725                                 -1.9681 
#>                                   theta                                   theta 
#>                                 -1.7060                                  0.5691
```

``` r
library(microbenchmark)

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
#>         expr    min     lq   mean median     uq    max neval
#>  PH           897.6  902.3  906.1  905.5  910.0  915.0     5
#>  PH     (2L)  569.9  578.3  579.6  580.8  583.5  585.3     5
#>  PH     (4L)  412.8  419.4  420.6  421.3  423.4  426.1     5
#>  PO          1282.6 1291.1 1297.8 1302.5 1303.6 1309.3     5
#>  PO     (2L)  788.8  797.2  804.1  801.8  805.7  827.0     5
#>  PO     (4L)  572.8  575.8  576.1  576.0  577.8  578.2     5
#>  probit      1644.0 1668.9 1679.2 1673.4 1687.0 1722.4     5
#>  probit (2L)  982.8  983.4  995.5  996.7 1004.1 1010.5     5
#>  probit (4L)  695.8  700.0  706.6  707.5  708.8  721.1     5
#> 
#> Method: GVA
#> -----------
#> Unit: milliseconds
#>         expr    min     lq   mean median     uq    max neval
#>  PH           287.2  287.4  291.4  291.3  294.8  296.4     5
#>  PH     (2L)  191.1  191.8  197.3  192.4  200.0  211.2     5
#>  PH     (4L)  153.8  155.0  155.0  155.2  155.5  155.5     5
#>  PO          1940.5 1943.5 1957.9 1960.9 1968.5 1975.8     5
#>  PO     (2L) 1124.2 1130.4 1143.5 1143.8 1153.5 1165.4     5
#>  PO     (4L)  795.7  798.1  806.4  803.6  805.8  828.6     5
#>  probit      3113.9 3147.0 3196.8 3160.9 3177.6 3384.7     5
#>  probit (2L) 1811.0 1815.1 1833.7 1823.8 1831.1 1887.3     5
#>  probit (4L) 1226.8 1231.8 1239.0 1239.4 1243.3 1253.7     5
```
