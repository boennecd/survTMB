#' @importFrom stats cov2cor
#' @export
print.MGSM_ADFit <- function(x, ...){
  cat(sprintf(
    "\nMGSM estimated with method %s with link %s from call:\n",
    sQuote(x$method), sQuote(x$link)),
    paste0("  ", deparse(x$ADFun_cl), collapse = "\n"), "\n",
    paste0("  ", deparse(x$fit_cl)  , collapse = "\n"), "\n\n", sep = "")

  cat("Estimated fixed effects:\n")
  is_fix <- !grepl("^theta", names(x$params))
  fix_par <- x$params[is_fix]
  print(fix_par)

  cat("\nEstimated random effect covariance matrix (correlation matrix) is:\n")
  vcov_par <- x$params[!is_fix]
  vcov_m   <- as.matrix(.theta_to_cov(vcov_par))
  colnames(vcov_m) <- rownames(vcov_m) <- x$rng_names

  cor_m <- cov2cor(vcov_m)
  diag(cor_m) <- sqrt(diag(vcov_m))

  to_print <- cbind(vcov_m, NA_real_, cor_m)
  colnames(to_print) <- c(colnames(vcov_m), "     ", colnames(cor_m))
  print(to_print, na.print = "")
  cat("(standard deviations are in the diagonal of the correlation matrix)\n")

  cat(sprintf(
    "\nEstimated %s is %.2f\n\n",
    if(x$is_va) "lower bound" else "log-likelihood", -x$optim$value))

  invisible(list(fix_par = fix_par, vcov = vcov))
}

#' @export
print.MGSM_ADFun <- function(x, ...){
  cat(sprintf(
    "\nMGSM objective function with link %s from call:\n",
    sQuote(x$link)),
    paste0("  ", deparse(x$cl), collapse = "\n"), "\n\n", sep = "")

  yn <- function(z)
    if(z) "Yes" else "No"


  cat("The following is available:\n")
  to_print <- c(
    `Laplace approximation`  = yn(!is.null(x$laplace)),
    `GVA`                    = yn(!is.null(x$gva)),
    `SNVA`                   = yn(!is.null(x$snva)),
    `Dense Hessian with VA`  = yn(x$dense_hess),
    `Sparse Hessian with VA` = yn(x$sparse_hess))

  for(i in seq_along(to_print))
    cat(sprintf("%-23s %3s\n", names(to_print)[i], to_print[i]))
  cat("\n")
}

