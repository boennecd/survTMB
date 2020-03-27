
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
# w/ Laplace approximatio
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
    
    `PO         `  = fit_model("PO"    , 1L, mth),
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
#>  PH           926.1  931.6  932.4  932.0  934.4  938.1     5
#>  PH     (2L)  586.7  591.6  599.7  598.2  598.7  623.5     5
#>  PH     (4L)  421.3  427.7  433.2  433.8  435.4  447.8     5
#>  PO          1298.6 1314.4 1318.8 1323.2 1323.7 1334.2     5
#>  PO     (2L)  818.1  818.1  823.7  821.8  824.1  836.4     5
#>  PO     (4L)  581.1  581.4  588.7  585.8  596.0  598.9     5
#>  probit      1658.5 1672.1 1674.8 1672.7 1678.8 1691.9     5
#>  probit (2L)  984.7  994.8 1004.1  995.7  996.6 1048.8     5
#>  probit (4L)  694.5  698.5  705.3  700.0  708.6  724.9     5
#> 
#> Method: GVA
#> -----------
#> Unit: milliseconds
#>         expr    min     lq   mean median     uq    max neval
#>  PH           288.9  290.4  290.9  290.7  291.6  293.0     5
#>  PH     (2L)  192.8  193.0  196.5  196.9  199.8  200.0     5
#>  PH     (4L)  148.5  149.4  150.4  149.8  151.0  153.4     5
#>  PO          1932.0 1935.2 1963.5 1967.0 1991.3 1991.8     5
#>  PO     (2L) 1111.1 1111.7 1114.5 1113.8 1113.9 1122.0     5
#>  PO     (4L)  735.4  741.7  743.7  741.8  746.4  753.3     5
#>  probit      3128.7 3132.7 3167.9 3151.0 3177.8 3249.1     5
#>  probit (2L) 1746.6 1751.3 1763.6 1753.0 1766.4 1800.9     5
#>  probit (4L) 1072.5 1080.8 1088.9 1087.3 1101.6 1102.2     5
```
