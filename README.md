
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
#> [1] -13031.22757
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -8.0638727322                            1.0425044459 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            5.6983653823                           11.8233679076 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            5.5964995964                           -1.3931263562 
#>                                   theta                                   theta 
#>                           -1.0679788694                            0.4390664046
fit_model("probit", method = "SNVA", param_type = "CP_trans")
#> $`lower bound`
#> [1] -13035.3136
#> 
#> $par
#>                             (Intercept)                                     trt 
#>                           -3.7445370729                            0.6055055795 
#> nsx(log(y), df = 3, intercept = FALSE)1 nsx(log(y), df = 3, intercept = FALSE)2 
#>                            2.6601098819                            5.0025483695 
#> nsx(log(y), df = 3, intercept = FALSE)3                                   theta 
#>                            2.9728743169                           -1.8488522015 
#>                                   theta                                   theta 
#>                           -1.5640410826                            0.2936106644
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
#>  PH           795.5  797.1  805.2  808.0  812.2  813.1     5
#>  PH     (2L)  507.8  513.5  523.7  523.4  528.5  545.4     5
#>  PH     (4L)  364.0  373.2  383.4  389.0  393.8  396.7     5
#>  PO          1150.9 1169.0 1185.8 1191.1 1193.3 1224.4     5
#>  PO     (2L)  716.8  720.5  742.7  743.1  757.6  775.5     5
#>  PO     (4L)  513.9  533.3  540.1  546.1  549.7  557.6     5
#>  probit      1556.9 1581.4 1598.9 1594.3 1620.0 1641.7     5
#>  probit (2L)  917.1  949.6  950.3  956.0  963.5  965.2     5
#>  probit (4L)  655.9  667.3  679.8  681.2  685.6  709.0     5
#> 
#> Method: GVA
#> -----------
#> Unit: milliseconds
#>         expr    min     lq   mean median     uq    max neval
#>  PH           248.3  251.2  252.3  253.9  254.0  254.0     5
#>  PH     (2L)  170.9  171.1  172.4  171.1  171.6  177.2     5
#>  PH     (4L)  139.9  140.6  142.3  141.0  143.3  146.4     5
#>  PO          1006.4 1007.8 1014.7 1014.5 1021.0 1023.8     5
#>  PO     (2L)  599.5  599.8  608.8  601.0  611.9  631.8     5
#>  PO     (4L)  419.2  422.0  431.0  437.4  438.3  438.3     5
#>  probit      1826.4 1830.4 1841.2 1833.4 1837.6 1878.2     5
#>  probit (2L) 1063.8 1067.2 1072.3 1072.1 1073.2 1085.1     5
#>  probit (4L)  720.9  724.5  727.3  727.4  728.4  735.1     5
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
#>  PH           440.8  441.2  446.6  442.7  450.2  458.4     5
#>  PH     (2L)  300.1  300.6  308.4  305.0  306.7  329.4     5
#>  PH     (4L)  234.2  236.6  245.3  244.7  251.2  259.9     5
#>  PO          3755.4 3761.6 3779.1 3767.3 3794.0 3817.2     5
#>  PO     (2L) 2718.2 2723.8 2797.3 2730.3 2896.5 2917.5     5
#>  PO     (4L) 1661.8 1675.4 1723.1 1699.8 1773.4 1804.8     5
#>  probit      4025.7 4025.9 4039.3 4040.7 4048.6 4055.4     5
#>  probit (2L) 2510.9 2556.5 2754.4 2767.2 2925.8 3011.6     5
#>  probit (4L) 1653.2 1677.2 1717.8 1685.1 1715.9 1857.5     5
#> 
#> Method: SNVA (CP_trans)
#> -----------------------
#> Unit: milliseconds
#>         expr    min     lq   mean median     uq    max neval
#>  PH           635.9  638.9  654.2  663.3  664.6  668.6     5
#>  PH     (2L)  332.4  339.6  343.3  340.4  344.1  359.8     5
#>  PH     (4L)  297.7  297.7  301.9  298.9  300.0  315.0     5
#>  PO          5751.8 5770.4 5790.6 5783.6 5803.8 5843.5     5
#>  PO     (2L) 2227.3 2594.0 2543.1 2613.8 2614.7 2665.7     5
#>  PO     (4L) 3225.5 3298.8 3321.4 3311.8 3385.0 3385.7     5
#>  probit      8828.4 8830.0 8917.9 8892.9 8928.4 9109.7     5
#>  probit (2L) 3593.4 3711.3 3758.5 3729.6 3878.8 3879.1     5
#>  probit (4L) 2547.2 2553.4 2580.6 2578.6 2599.3 2624.5     5
```
