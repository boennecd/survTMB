#' Create mgsm Object to Pass to psqn Method
#' @inheritParams make_mgsm_ADFun
#'
#' @param do_setup character vector indicating which approximation to setup.
#' See \code{\link{make_mgsm_ADFun}}.
#'
#' @importFrom stats setNames
#' @export
make_mgsm_psqn_obj <- function(
  formula, data, df, tformula = NULL, Z, cluster,
  do_setup = c("GVA", "SNVA"), n_nodes = 20L,
  link = c("PH", "PO", "probit"),
  theta = NULL, beta = NULL, opt_func = .opt_default, n_threads = 1L,
  skew_start = -.0001, kappa = .MGSM_default_kappa){
  #####
  # checks
  stopifnot(
    is.function(opt_func),
    all(do_setup %in% c(.gva_char, .snva_char)))

  #####
  # get arguments to pass on
  cluster <- substitute(cluster)
  args_pass <- mgsm_setup(
    formula = formula, data = data, df = df, tformula = tformula, Z = Z,
    cluster = cluster, n_nodes = n_nodes,
    link = link, theta = theta, beta = beta, opt_func = opt_func,
    n_threads = n_threads, skew_start = skew_start, kappa = kappa)

  #####
  # create C++ object and return
  # split data
  cluster_data <- tapply(
    1:length(args_pass$grp), args_pass$grp, function(indices)
      list(tobs  = args_pass$tobs [indices],
           event = args_pass$event[indices],
           X     = args_pass$X    [indices, , drop = FALSE],
           XD    = args_pass$XD   [indices, , drop = FALSE],
           Z     = args_pass$Z    [indices, , drop = FALSE]),
      simplify = FALSE)

  # find starting values for VA parameters
  data_ad_func <- list(
    tobs = args_pass$tobs, event = args_pass$event, X = args_pass$X,
    XD = args_pass$XD, Z = args_pass$Z, grp = args_pass$grp,
    link = args_pass$link, grp_size = args_pass$grp_size,
    n_threads = args_pass$n_threads)
  params <- list(eps = args_pass$eps, kappa = args_pass$kappa,
                 b = args_pass$b, theta = args_pass$theta)
  if(length(args_pass$skew_start) == 1L && args_pass$n_rng > 1L)
    args_pass$skew_start <- rep(args_pass$skew_start, args_pass$n_rng)

  va_start <- .get_MGSM_VA_start(
    n_rng = args_pass$n_rng, params = params, data_ad_func = data_ad_func,
    opt_func = opt_func, skew_start = args_pass$skew_start, is_cp = FALSE)
  names(va_start) <- mgsm_get_snva_names(
    args_pass$n_rng, args_pass$n_grp, "DP")

  # create the C++ object
  setup_atomic_cache(
    n_nodes = args_pass$n_nodes, type = .snva_char,
    link = args_pass$link)
  cpp_ptr <- psqn_get_mgsm_funcs(
    data = cluster_data, eps = args_pass$eps, kappa = args_pass$kappa,
    b = args_pass$b, theta = args_pass$theta, theta_va = va_start,
    n_nodes = args_pass$n_nodes, link = args_pass$link,
    max_threads = n_threads)

  # return the object with starting values
  par <- c(args_pass$b, args_pass$theta, va_start)

  n_threads <- args_pass$n_threads
  fn <- function(x, ...)
    eval_psqn_mgsm(x, ptr = cpp_ptr, n_threads = n_threads)
  gr <- function(x, ...)
    grad_psqn_mgsm(x, ptr = cpp_ptr, n_threads = n_threads)

  out <- structure(list(
    par = par, cpp_ptr = cpp_ptr, y = args_pass$y, event = args_pass$event,
    X = args_pass$X, XD = args_pass$XD, Z = args_pass$Z,
    grp = args_pass$grp, terms = args_pass$terms, cl = match.call(),
    link = args_pass$link, fn = fn, gr = gr, n_threads = n_threads),
    class = "mgsm_psqn")

  rm(list = setdiff(ls(), c("cpp_ptr", "out", "n_threads")))
  out
}

#' Optimize mgsm Object Using the psqn Method
#'
#' @param object object of class \code{"mgsm_psqn"}.
#' @param par initial values for the parameters to be optimized over. A
#' default value is used if this is \code{NULL}.
#' @param rel_eps relative convergence threshold.
#' @param max_it maximum number of iterations.
#' @param n_threads number of threads to use.
#' @param c1,c2,use_bfgs,cg_tol,strong_wolfe,trace arguments passed to
#' \code{\link{psqn}}.
#'
#' @export
optim_mgsm_psqn <- function(object, par = NULL,
                            rel_eps = sqrt(.Machine$double.eps),
                            max_it = 1000L, n_threads = object$n_threads,
                            c1 = 1e-4, c2 = .9, use_bfgs = TRUE,
                            cg_tol = .5, strong_wolfe = TRUE, trace = 0L){
  #####
  # chekcs
  stopifnot(
    is.numeric(par) || is.null(par),
    inherits(object, "mgsm_psqn"),
    is.numeric(rel_eps), length(rel_eps) == 1L, is.finite(rel_eps),
    is.integer(max_it), length(max_it) == 1L, is.finite(n_threads),
    is.integer(n_threads), length(n_threads) == 1L, n_threads > 0,
    is.integer(trace), length(trace) == 1L, is.finite(trace),
    is.numeric(c1), length(c1) == 1L, is.finite(c1),
    is.numeric(c2), length(c2) == 1L, is.finite(c2),
    is.numeric(cg_tol), length(cg_tol) == 1L, is.finite(cg_tol),
    is.logical(use_bfgs), length(use_bfgs) == 1L, is.finite(use_bfgs),
    is.logical(strong_wolfe), length(strong_wolfe) == 1L,
    is.finite(strong_wolfe))

  #####
  # optimize
  if(is.null(par))
    par <- object$par
  names(par) <- names(object$par)

  out <- psqn_optim_mgsm(val = par, ptr = object$cpp_ptr, rel_eps = rel_eps,
                         max_it = max_it, n_threads = n_threads, c1 = c1,
                         c2 = c2, use_bfgs = use_bfgs, trace = trace,
                         cg_tol = cg_tol, strong_wolfe = strong_wolfe)

  # sort out variational parameters from model parameters
  is_val <- grepl("^g\\d+:", names(par))
  out$params    <- out$par[!is_val]
  out$va_params <- out$par[ is_val]

  out
}
