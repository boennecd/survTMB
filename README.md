
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
fit_model <- function(link, n_threads = 2L, method = "Laplace", 
                      param_type = "DP"){
  adfun <- make_gsm_ADFun(
    Surv(y, uncens) ~ trt, cluster = as.factor(center), 
    Z = ~ trt, df = 3L, data = dat, link = link, do_setup = method, 
    n_threads = n_threads, param_type = param_type)
  fit <- 
    if(method == "Laplace")
      do.call(optim, adfun$laplace) else if(method == "GVA") 
        do.call(optim, adfun$gva) else if(method == "SNVA")
          do.call(optim, adfun$snva)
  
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
#>                                -7.82323                                 0.71985 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                 5.39415                                11.38909 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                 4.79923                                -1.84295 
#>                                   theta                                   theta 
#>                                -2.23190                                 1.79529
fit_model("PO"    )
#> $logLik
#> [1] -13031
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                 -8.0494                                  1.0297 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                  5.6971                                 11.8210 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                  5.5935                                 -1.5842 
#>                                   theta                                   theta 
#>                                 -1.8700                                  1.9480
fit_model("probit")
#> $logLik
#> [1] -13035
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                -3.73992                                 0.59927 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                 2.66012                                 5.00439 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                 2.97221                                -1.99884 
#>                                   theta                                   theta 
#>                                -1.85024                                 0.79242

######
# w/ GVA
fit_model("PH"    , method = "GVA")
#> $`lower bound`
#> [1] -13027
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                -7.82687                                 0.72329 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                 5.39399                                11.38961 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                 4.79935                                -1.78716 
#>                                   theta                                   theta 
#>                                -1.73611                                 0.87839
fit_model("PO"    , method = "GVA")
#> $`lower bound`
#> [1] -13031
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                -8.05425                                 1.03454 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                 5.69722                                11.82202 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                 5.59432                                -1.51778 
#>                                   theta                                   theta 
#>                                -1.38503                                 0.98453
fit_model("probit", method = "GVA")
#> $`lower bound`
#> [1] -13035
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                -3.74182                                 0.60186 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                 2.66017                                 5.00389 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                 2.97248                                -1.96807 
#>                                   theta                                   theta 
#>                                -1.70599                                 0.56907

######
# w/ SNVA (DP: direct parameterization)
fit_model("PH"    , method = "SNVA", param_type = "DP")
#> $`lower bound`
#> [1] -13027
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                -7.82680                                 0.72339 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                 5.39410                                11.38990 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                 4.79947                                -1.80457 
#>                                   theta                                   theta 
#>                                -1.75024                                 0.91588
fit_model("PO"    , method = "SNVA", param_type = "DP")
#> $`lower bound`
#> [1] -13031
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                 -8.0537                                  1.0329 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                  5.6973                                 11.8220 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                  5.5942                                 -1.5246 
#>                                   theta                                   theta 
#>                                 -1.4517                                  1.0957
fit_model("probit", method = "SNVA", param_type = "DP")
#> $`lower bound`
#> [1] -13035
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                -3.74731                                 0.60649 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                 2.66109                                 5.00553 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                 2.97377                                -1.80078 
#>                                   theta                                   theta 
#>                                -1.50356                                 0.18189

######
# w/ SNVA (CP: centralized parameterization)
fit_model("PH"    , method = "SNVA", param_type = "CP_trans")
#> $`lower bound`
#> [1] -13027
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                -7.82860                                 0.72525 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                 5.39404                                11.39005 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                 4.79951                                -1.72912 
#>                                   theta                                   theta 
#>                                -1.61982                                 0.67150
fit_model("PO"    , method = "SNVA", param_type = "CP_trans")
#> $`lower bound`
#> [1] -13031
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                -8.06200                                 1.04236 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                 5.69761                                11.82252 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                 5.59589                                -1.41066 
#>                                   theta                                   theta 
#>                                -1.09963                                 0.50051
fit_model("probit", method = "SNVA", param_type = "CP_trans")
#> $`lower bound`
#> [1] -13035
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                                -3.73762                                 0.60509 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                                 2.65679                                 4.98911 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                                 2.97249                                -1.85121 
#>                                   theta                                   theta 
#>                                -1.56902                                 0.30720
```

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
#>         expr     min      lq    mean  median      uq     max neval
#>  PH           738.26  742.60  745.22  747.59  748.75  748.91     5
#>  PH     (2L)  481.78  483.67  483.67  483.72  484.02  485.18     5
#>  PH     (4L)  353.13  361.74  360.62  361.91  362.58  363.76     5
#>  PO          1097.56 1099.71 1116.59 1116.12 1122.17 1147.41     5
#>  PO     (2L)  683.62  692.32  694.46  697.47  699.23  699.63     5
#>  PO     (4L)  501.71  502.78  504.71  503.73  506.68  508.62     5
#>  probit      1477.58 1477.86 1482.35 1478.21 1479.61 1498.50     5
#>  probit (2L)  874.49  883.35  887.60  883.51  884.05  912.59     5
#>  probit (4L)  630.24  630.63  632.21  632.70  633.51  633.97     5
#> 
#> Method: GVA
#> -----------
#> Unit: milliseconds
#>         expr     min      lq    mean  median      uq     max neval
#>  PH           240.14  242.23  244.17  244.85  245.91  247.71     5
#>  PH     (2L)  163.87  164.98  166.96  165.25  170.07  170.64     5
#>  PH     (4L)  135.83  136.28  137.66  136.90  139.20  140.12     5
#>  PO          1802.56 1805.13 1809.05 1806.78 1809.41 1821.35     5
#>  PO     (2L) 1074.84 1080.06 1079.92 1080.83 1080.98 1082.90     5
#>  PO     (4L)  746.32  747.97  753.97  751.83  756.25  767.49     5
#>  probit      2925.66 2928.30 2944.67 2929.57 2954.34 2985.47     5
#>  probit (2L) 1705.91 1707.02 1713.25 1710.41 1713.55 1729.37     5
#>  probit (4L) 1160.91 1162.11 1164.53 1163.75 1165.15 1170.72     5
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
#>         expr     min      lq    mean  median      uq     max neval
#>  PH           423.00  423.29  426.72  427.25  429.46  430.59     5
#>  PH     (2L)  290.08  290.20  295.03  296.54  298.23  300.08     5
#>  PH     (4L)  230.61  230.67  231.54  230.89  232.75  232.76     5
#>  PO          4818.53 4836.45 4856.13 4845.70 4849.59 4930.37     5
#>  PO     (2L) 2815.99 2826.59 2832.49 2827.17 2836.92 2855.79     5
#>  PO     (4L) 1940.50 1942.48 1955.05 1949.79 1960.15 1982.32     5
#>  probit      5237.08 5253.59 5257.73 5258.01 5267.07 5272.90     5
#>  probit (2L) 3016.84 3022.58 3034.02 3030.32 3040.88 3059.48     5
#>  probit (4L) 2058.96 2063.00 2071.73 2069.38 2072.80 2094.50     5
#> 
#> Method: SNVA (CP_trans)
#> -----------------------
#> Unit: milliseconds
#>         expr     min      lq    mean  median      uq     max neval
#>  PH           588.39  590.38  592.95  593.94  594.70  597.32     5
#>  PH     (2L)  321.20  325.31  326.59  325.90  327.59  332.93     5
#>  PH     (4L)  277.03  281.59  282.55  281.94  284.96  287.21     5
#>  PO          4529.42 4546.81 4554.70 4547.64 4574.04 4575.58     5
#>  PO     (2L) 2949.73 2950.58 2960.60 2955.59 2970.63 2976.46     5
#>  PO     (4L) 2021.87 2025.01 2026.50 2026.39 2028.22 2031.02     5
#>  probit      7966.02 7981.01 7985.19 7986.18 7994.27 7998.46     5
#>  probit (2L) 5133.06 5139.84 5152.83 5156.47 5158.47 5176.33     5
#>  probit (4L) 2732.39 2748.96 2752.41 2753.71 2754.35 2772.61     5
```
