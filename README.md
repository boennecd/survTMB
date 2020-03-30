
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
#> [1] -13031.19526
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -8.0620008745                            1.0423618353 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            5.6976091208                           11.8225239481 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            5.5958936252                           -1.4106585136 
#>                                   theta                                   theta 
#>                           -1.0996280311                            0.5005085105
fit_model("probit", method = "SNVA", param_type = "CP_trans")
#> $`lower bound`
#> [1] -13035.23606
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -3.7376240206                            0.6050923633 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            2.6567877023                            4.9891095080 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            2.9724926394                           -1.8512086029 
#>                                   theta                                   theta 
#>                           -1.5690173942                            0.3072032997
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
#>  PH           759.0  760.5  782.1  762.1  790.1  838.9     5
#>  PH     (2L)  470.8  481.7  488.4  488.0  496.3  505.3     5
#>  PH     (4L)  350.0  375.4  375.5  378.3  378.5  395.0     5
#>  PO          1116.8 1161.1 1189.6 1168.7 1196.0 1305.2     5
#>  PO     (2L)  689.6  706.0  722.3  730.2  731.6  754.2     5
#>  PO     (4L)  529.8  541.5  570.5  554.3  562.9  664.0     5
#>  probit      1505.9 1523.4 1560.7 1538.8 1555.2 1679.9     5
#>  probit (2L)  880.9  898.2  935.5  934.1  945.5 1019.0     5
#>  probit (4L)  631.1  641.0  654.4  643.0  645.0  711.8     5
#> 
#> Method: GVA
#> -----------
#> Unit: milliseconds
#>         expr    min     lq   mean median     uq    max neval
#>  PH           250.1  250.2  253.8  252.1  257.1  259.7     5
#>  PH     (2L)  168.5  170.3  170.7  170.8  171.1  172.8     5
#>  PH     (4L)  138.1  139.4  140.3  139.6  142.0  142.6     5
#>  PO          1804.9 1829.0 1842.6 1840.9 1843.2 1895.0     5
#>  PO     (2L) 1085.7 1089.3 1113.8 1106.5 1119.6 1167.8     5
#>  PO     (4L)  756.9  765.5  772.4  770.3  772.3  797.2     5
#>  probit      2948.1 2954.8 3016.5 2979.9 3047.3 3152.3     5
#>  probit (2L) 1727.6 1733.1 1748.1 1738.2 1752.7 1788.8     5
#>  probit (4L) 1159.8 1163.4 1176.9 1171.8 1177.8 1211.7     5
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
#>  PH           432.1  433.1  434.8  434.0  436.4  438.1     5
#>  PH     (2L)  288.3  288.9  290.2  289.2  291.4  293.1     5
#>  PH     (4L)  227.6  227.7  228.7  229.2  229.3  229.7     5
#>  PO          4843.0 4861.0 4879.5 4873.5 4908.3 4911.7     5
#>  PO     (2L) 2861.4 2868.7 2882.4 2876.7 2902.4 2902.9     5
#>  PO     (4L) 1937.2 1972.9 1972.0 1979.6 1980.2 1990.1     5
#>  probit      5275.8 5279.0 5282.5 5283.6 5285.8 5288.2     5
#>  probit (2L) 3060.1 3065.4 3077.1 3073.8 3077.4 3108.7     5
#>  probit (4L) 2075.0 2077.7 2087.3 2081.5 2088.3 2113.8     5
#> 
#> Method: SNVA (CP_trans)
#> -----------------------
#> Unit: milliseconds
#>         expr    min     lq   mean median     uq    max neval
#>  PH           609.0  610.4  621.9  622.5  627.5  639.9     5
#>  PH     (2L)  323.9  324.3  331.7  325.0  332.3  352.8     5
#>  PH     (4L)  275.5  276.7  282.9  278.9  291.5  291.8     5
#>  PO          4610.9 4621.0 4643.3 4633.8 4662.0 4688.8     5
#>  PO     (2L) 3035.9 3087.6 3174.6 3171.2 3256.0 3322.0     5
#>  PO     (4L) 2130.0 2136.6 2177.6 2152.4 2170.5 2298.6     5
#>  probit      8103.7 8136.1 8210.7 8171.0 8302.8 8339.8     5
#>  probit (2L) 5216.9 5251.1 5294.9 5252.7 5305.9 5447.8     5
#>  probit (4L) 2924.7 2925.7 2947.9 2939.2 2968.0 2981.8     5
```
