#' @export
fit_mgsm <- function(object, method, optim = object$opt_func, ...){
  method <- method[1]
  stopifnot(
    inherits(object, "GSM_ADFun"),
    is.character(method), method %in% c(
      .laplace_char, .gva_char, .snva_char))

  optim_args <- object[[tolower(method)]]
  stopifnot(!is.null(optim_args))

  # fit model
  dots <- list(...)
  optim_args[names(dots)] <- dots
  fit <- do.call(optim, optim_args)

  # get parameters
  is_va <- method %in% c(.gva_char, .snva_char)
  if(is_va){
    params    <- optim_args$get_params(fit$par)
    va_params <- tail(fit$par, -length(params))

  } else if(method == .laplace_char){
    params <- fit$par
    va_params <- NULL

  } else
    stop(sprintf("Unkown method: %s", sQuote(method)))

  structure(list(
    params = params, va_params = va_params, link = object$link,
    ADFun_cl = object$cl, fit_cl = match.call(), method = method,
    optim = fit, is_va = is_va, fix_names = colnames(object$X),
    rng_names = colnames(object$Z)), class = "GSM_ADFit")
}
