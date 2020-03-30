
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
#>  PH           741.7  742.7  748.3  743.5  745.6  767.9     5
#>  PH     (2L)  463.8  474.1  481.6  485.4  487.9  496.7     5
#>  PH     (4L)  347.4  355.8  365.0  367.1  368.9  385.6     5
#>  PO          1094.6 1098.6 1107.6 1108.4 1118.0 1118.5     5
#>  PO     (2L)  678.7  682.9  693.0  691.2  703.5  708.8     5
#>  PO     (4L)  497.4  504.2  507.7  512.2  512.2  512.7     5
#>  probit      1484.2 1492.9 1494.7 1496.4 1499.5 1500.8     5
#>  probit (2L)  905.9  906.7  908.6  908.0  910.9  911.2     5
#>  probit (4L)  641.8  646.4  646.7  647.5  648.7  649.1     5
#> 
#> Method: GVA
#> -----------
#> Unit: milliseconds
#>         expr    min     lq   mean median     uq    max neval
#>  PH           242.2  242.4  244.6  244.9  245.5  248.1     5
#>  PH     (2L)  165.5  166.1  167.1  167.4  167.6  169.1     5
#>  PH     (4L)  136.4  136.4  137.0  137.0  137.2  137.9     5
#>  PO           969.6  976.0  977.8  976.3  980.8  986.3     5
#>  PO     (2L)  584.8  585.9  587.0  587.6  587.7  589.0     5
#>  PO     (4L)  408.9  410.0  410.6  410.4  410.8  412.7     5
#>  probit      2984.4 2996.7 3002.7 3002.9 3010.4 3019.2     5
#>  probit (2L) 1739.9 1740.8 1743.6 1742.6 1747.3 1747.6     5
#>  probit (4L) 1183.8 1184.3 1191.5 1193.3 1194.5 1201.8     5
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
#>  PH           429.5  431.2  433.3  434.3  435.0  436.6     5
#>  PH     (2L)  302.6  328.5  343.3  331.1  348.2  406.0     5
#>  PH     (4L)  236.7  238.7  243.7  240.6  244.0  258.6     5
#>  PO          4674.6 4678.1 4965.8 4819.2 4822.7 5834.5     5
#>  PO     (2L) 2729.1 2890.7 3031.2 2989.4 3091.2 3455.6     5
#>  PO     (4L) 1887.9 1890.1 1949.7 1919.6 1999.1 2051.7     5
#>  probit      5513.2 5574.1 5674.2 5691.1 5700.1 5892.5     5
#>  probit (2L) 3223.6 3254.4 3322.9 3310.7 3314.2 3511.5     5
#>  probit (4L) 2138.3 2165.8 2231.6 2189.6 2213.2 2451.2     5
#> 
#> Method: SNVA (CP_trans)
#> -----------------------
#> Unit: milliseconds
#>         expr    min     lq   mean median     uq    max neval
#>  PH           614.8  645.6  716.9  674.2  764.2  885.8     5
#>  PH     (2L)  333.3  338.0  346.8  345.4  352.8  364.3     5
#>  PH     (4L)  295.2  297.3  351.5  298.2  321.1  545.9     5
#>  PO          5626.6 5652.3 5775.5 5691.9 5718.9 6187.7     5
#>  PO     (2L) 2269.9 2294.5 2423.4 2444.6 2456.7 2651.4     5
#>  PO     (4L) 2456.7 2459.5 2553.9 2542.1 2593.8 2717.2     5
#>  probit      8531.0 8616.7 8775.9 8767.4 8771.1 9193.5     5
#>  probit (2L) 5628.4 5641.3 5887.6 5770.1 5804.3 6593.7     5
#>  probit (4L) 2951.2 3082.0 3218.2 3108.7 3369.6 3579.5     5
```
