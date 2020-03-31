
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
#> [1] -13026.6757
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -7.8232328523                            0.7198516887 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            5.3941546610                           11.3890932863 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            4.7992329025                           -1.8429548771 
#>                                   theta                                   theta 
#>                           -2.2318985734                            1.7952915040
fit_model("PO"    )
#> $logLik
#> [1] -13031.16041
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                            -8.049406749                             1.029741911 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                             5.697060912                            11.821001229 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                             5.593477692                            -1.584164581 
#>                                   theta                                   theta 
#>                            -1.869970860                             1.947962080
fit_model("probit")
#> $logLik
#> [1] -13035.13832
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -3.7399188924                            0.5992700212 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            2.6601203807                            5.0043913748 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            2.9722064746                           -1.9988364130 
#>                                   theta                                   theta 
#>                           -1.8502382811                            0.7924195633

######
# w/ GVA
fit_model("PH"    , method = "GVA")
#> $`lower bound`
#> [1] -13026.73806
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -7.8268712791                            0.7232878681 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            5.3939925692                           11.3896128017 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            4.7993519927                           -1.7871605349 
#>                                   theta                                   theta 
#>                           -1.7361101036                            0.8783932750
fit_model("PO"    , method = "GVA")
#> $`lower bound`
#> [1] -13031.10877
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -8.0542478654                            1.0345436621 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            5.6972245460                           11.8220193230 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            5.5943198048                           -1.5177841973 
#>                                   theta                                   theta 
#>                           -1.3850294115                            0.9845349819
fit_model("probit", method = "GVA")
#> $`lower bound`
#> [1] -13035.14892
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -3.7418205373                            0.6018560988 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            2.6601732516                            5.0038905656 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            2.9724793475                           -1.9680664220 
#>                                   theta                                   theta 
#>                           -1.7059861559                            0.5690654814

######
# w/ SNVA (DP: direct parameterization)
fit_model("PH"    , method = "SNVA", param_type = "DP")
#> $`lower bound`
#> [1] -13026.73553
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -7.8267969825                            0.7233878832 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            5.3940985676                           11.3899028519 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            4.7994743473                           -1.8045742722 
#>                                   theta                                   theta 
#>                           -1.7502359919                            0.9158786061
fit_model("PO"    , method = "SNVA", param_type = "DP")
#> $`lower bound`
#> [1] -13031.10591
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                            -8.053660182                             1.032859368 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                             5.697300624                            11.822048489 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                             5.594198518                            -1.524556393 
#>                                   theta                                   theta 
#>                            -1.451694088                             1.095708146
fit_model("probit", method = "SNVA", param_type = "DP")
#> $`lower bound`
#> [1] -13035.3148
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -3.7473128000                            0.6064949783 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            2.6610868775                            5.0055287310 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            2.9737651686                           -1.8007789312 
#>                                   theta                                   theta 
#>                           -1.5035552927                            0.1818880866

######
# w/ SNVA (CP: centralized parameterization)
fit_model("PH"    , method = "SNVA", param_type = "CP_trans")
#> $`lower bound`
#> [1] -13026.72826
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -7.8286041580                            0.7252509681 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            5.3940371090                           11.3900490246 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            4.7995132667                           -1.7291160711 
#>                                   theta                                   theta 
#>                           -1.6198187270                            0.6715015168
fit_model("PO"    , method = "SNVA", param_type = "CP_trans")
#> $`lower bound`
#> [1] -13031.26088
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -8.0671606397                            1.0446280211 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            5.6987652826                           11.8245162174 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            5.5972058647                           -1.3689499765 
#>                                   theta                                   theta 
#>                           -1.0405385420                            0.3930502906
fit_model("probit", method = "SNVA", param_type = "CP_trans")
#> $`lower bound`
#> [1] -13035.21828
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -3.7480724987                            0.6049313457 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            2.6623278032                            5.0096127002 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            2.9737886481                           -1.8604772282 
#>                                   theta                                   theta 
#>                           -1.5753973824                            0.3220013039
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
#>         expr    min     lq   mean median     uq    max neval
#>  PH           774.4  775.5  782.4  783.3  788.2  790.5     5
#>  PH     (2L)  503.1  510.3  510.7  511.2  512.3  516.7     5
#>  PH     (4L)  360.2  372.3  375.7  376.7  383.7  385.5     5
#>  PO          1103.2 1138.8 1143.9 1143.7 1159.7 1174.4     5
#>  PO     (2L)  700.4  704.2  713.1  719.3  720.0  721.3     5
#>  PO     (4L)  518.3  523.8  525.6  526.8  527.8  531.2     5
#>  probit      1509.6 1540.7 1539.8 1546.2 1550.7 1551.9     5
#>  probit (2L)  902.8  911.0  920.8  913.1  931.0  946.1     5
#>  probit (4L)  656.0  658.6  664.7  663.2  672.1  673.7     5
#> 
#> Method: GVA
#> -----------
#> Unit: milliseconds
#>         expr    min     lq   mean median     uq    max neval
#>  PH           249.7  252.1  254.4  252.3  252.8  265.3     5
#>  PH     (2L)  169.4  170.5  171.6  171.2  172.8  174.2     5
#>  PH     (4L)  138.2  139.5  140.2  139.7  139.8  143.7     5
#>  PO          1012.3 1014.8 1015.3 1014.8 1015.1 1019.3     5
#>  PO     (2L)  620.8  622.6  625.1  624.6  625.4  631.9     5
#>  PO     (4L)  429.8  429.9  431.8  430.8  433.9  434.7     5
#>  probit      1780.8 1787.8 1789.3 1790.7 1791.7 1795.5     5
#>  probit (2L) 1045.7 1053.6 1055.6 1056.0 1060.3 1062.6     5
#>  probit (4L)  704.4  712.9  714.7  715.7  717.8  722.9     5
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
#>         expr    min     lq   mean median     uq    max neval
#>  PH           445.5  446.9  451.8  453.0  455.8  457.9     5
#>  PH     (2L)  294.4  295.6  297.8  296.6  300.7  301.6     5
#>  PH     (4L)  232.2  232.3  234.0  232.6  236.2  236.5     5
#>  PO          4637.5 4640.4 4668.0 4644.7 4674.2 4743.1     5
#>  PO     (2L) 2765.5 2773.0 2805.6 2783.4 2823.9 2882.0     5
#>  PO     (4L) 1872.3 1874.7 1902.8 1878.0 1883.1 2006.2     5
#>  probit      4904.0 4913.8 4934.4 4922.7 4932.9 4998.5     5
#>  probit (2L) 2855.5 2862.2 2876.5 2868.8 2888.1 2908.0     5
#>  probit (4L) 1942.8 1943.5 1962.1 1960.3 1979.1 1984.7     5
#> 
#> Method: SNVA (CP_trans)
#> -----------------------
#> Unit: milliseconds
#>         expr    min     lq   mean median     uq    max neval
#>  PH           623.8  632.1  641.1  637.5  655.1  657.3     5
#>  PH     (2L)  324.3  332.0  339.5  335.5  345.4  360.3     5
#>  PH     (4L)  279.9  282.5  287.2  283.9  291.3  298.2     5
#>  PO          5479.2 5494.1 5577.8 5516.1 5693.0 5706.8     5
#>  PO     (2L) 2195.8 2231.8 2246.3 2246.5 2258.5 2299.2     5
#>  PO     (4L) 2377.7 2398.9 2468.4 2444.3 2457.9 2663.4     5
#>  probit      7774.5 7819.4 7823.1 7825.8 7841.0 7855.0     5
#>  probit (2L) 4371.3 4430.6 4459.0 4462.4 4463.8 4567.1     5
#>  probit (4L) 2924.5 2941.9 2997.6 2973.8 3070.9 3076.8     5
```
