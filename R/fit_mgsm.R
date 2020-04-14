#' Fit Mixed Generalized Survival Model using an Objective Function
#'
#' @description
#' Fits a mixed generalized survival model using an objective function.
#'
#' @param object an object with class \code{MGSM_ADFun}.
#' @param method character with the approximation method to use.
#' @param optim general-purpose optimization function with an interface like
#'              \code{\link{optim}}.
#' @param ... additional arguments passed to \code{optim}.
#'
#' @details
#' The advantage of using this function is that the output is separated,
#' named, and a \code{MGSM_ADFit} class is returned. Thus, other S3 function
#' in this package can be used.
#'
#' Other functions than \code{\link{optim}} can be used to estimate the
#' model as long as the function has a similar interface.
#'
#' See the README at \url{https://github.com/boennecd/survTMB} for
#' further information and examples.
#'
#' @return
#' An \code{MGSM_ADFit} object. The elements are
#' \item{params}{estimated model parameters.}
#' \item{va_params}{estimated variational parameters if a variational approximation is used.}
#' \item{link}{character with the link function which is used in the model.}
#' \item{ADFun_cl}{matched call when the objective function was constructed.}
#' \item{fit_cl}{matched call.}
#' \item{method}{character with the used approximation method.}
#' \item{optim}{returned object from \code{optim}.}
#' \item{is_va}{logical for whether a variational approximation is used.}
#' \item{fix_names}{character vector with fixed effect names.}
#' \item{rng_names}{character vector with random effect names.}
#'
#' @examples
#' library(survTMB)
#' if(require(coxme)){
#'   # construct function with a random intercept and a proportional hazard
#'   # link function
#'   func <- make_mgsm_ADFun(
#'     Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
#'     df = 3L, data = eortc, link = "PH",
#'     do_setup = c("Laplace", "GVA", "SNVA"), n_threads = 1L,
#'     dense_hess = TRUE, sparse_hess = TRUE)
#'   print(func)
#'
#'   lfit <- fit_mgsm(func, "Laplace")
#'   print(lfit)
#'
#'   gfit <- fit_mgsm(func, "GVA")
#'   print(gfit)
#'
#'   sfit <- fit_mgsm(func, "SNVA")
#'   print(sfit)
#'
#'   # compute Hessians
#'   gpar <- with(gfit, c(params, va_params))
#'   dhess <- func$gva$he   (gpar) # dense
#'   shess <- func$gva$he_sp(gpar) # sparse
#'
#'   # they match
#'   stopifnot(all.equal(dhess, as.matrix(shess), check.attributes = FALSE))
#' }
#'
#' @seealso
#' \code{\link{make_mgsm_ADFun}}
#'
#' @export
fit_mgsm <- function(object, method, optim = object$opt_func, ...){
  method <- method[1]
  stopifnot(
    inherits(object, "MGSM_ADFun"),
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
    rng_names = colnames(object$Z)), class = "MGSM_ADFit")
}
