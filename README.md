
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
#>  PH           736.2  736.3  744.7  748.8  750.7  751.8     5
#>  PH     (2L)  482.9  484.5  487.1  488.2  489.7  490.4     5
#>  PH     (4L)  358.4  370.7  369.0  371.4  371.7  373.1     5
#>  PO          1088.3 1098.6 1102.4 1100.1 1106.6 1118.2     5
#>  PO     (2L)  698.4  703.4  706.3  706.1  708.4  715.1     5
#>  PO     (4L)  498.7  500.4  510.5  514.5  518.5  520.1     5
#>  probit      1478.2 1482.8 1484.9 1483.4 1487.0 1493.3     5
#>  probit (2L)  884.6  885.3  889.3  886.9  890.3  899.5     5
#>  probit (4L)  646.1  649.4  652.4  651.3  653.7  661.2     5
#> 
#> Method: GVA
#> -----------
#> Unit: milliseconds
#>         expr    min     lq   mean median     uq    max neval
#>  PH           242.2  245.2  247.3  246.8  247.6  254.5     5
#>  PH     (2L)  167.7  168.3  169.0  169.3  169.5  170.2     5
#>  PH     (4L)  138.5  139.1  140.1  139.4  140.8  143.0     5
#>  PO          1822.3 1823.5 1826.6 1824.8 1830.9 1831.5     5
#>  PO     (2L) 1090.1 1090.2 1097.5 1095.5 1101.9 1109.5     5
#>  PO     (4L)  755.4  760.5  764.5  764.5  765.7  776.4     5
#>  probit      2950.6 2974.0 2979.9 2979.7 2979.8 3015.5     5
#>  probit (2L) 1720.1 1731.6 1733.3 1732.3 1735.5 1747.0     5
#>  probit (4L) 1173.2 1173.9 1176.4 1175.8 1176.4 1182.7     5
```
