.laplace_char <- "Laplace"
.gva_char     <- "GVA"
.snva_char    <- "SNVA"

.MGSM_defaul_eps <- .Machine$double.xmin^(1/5)
.MGSM_default_kappa <- 100
.skew_boundary <- 0.99527

#' Maps Between Centralized and Direct Parameters
#'
#' @description
#' Maps between centralized parameters (CP) and direct parameters (DP) for
#' the multivariate skew-normal distribution.
#'
#' @param mu mean of the random variable.
#' @param Sigma covariance matrix of the random variable.
#' @param gamma pearson skewness coefficients of the random variable.
#' @param xi location parameter of the random variable.
#' @param Psi scale parameter of the random variable.
#' @param alpha skewness parameter of the random variable.
#'
#' @details
#' We say that a random variable is skew-normal distributed if its density
#' is given by
#'
#' \deqn{2\phi^{(K)}(y - \xi; \Psi)\Phi(\alpha^\top diag(\psi)^{-1}(y - \xi))}
#'
#' where \eqn{\psi} is the square root of the diagonal of \eqn{\Psi}. The
#' two functions maps between the above DP to the CP where
#' \eqn{\mu = E(Y)}, \eqn{\Sigma = Var(Y)}, and
#' \eqn{\gamma_i =  E((Y_i - E(Y_i))^3) / Var(Y_i)^{3/2}}.
#'
#' @examples
#' mu <- 1:3
#' Sig <- diag(3)
#' Sig[lower.tri(Sig)] <- Sig[upper.tri(Sig)] <- .25
#' gamma <- ((1:3) - 2) / 10
#'
#' dps <- cp_to_dp(mu = mu, Sigma = Sig, gamma = gamma)
#' dps
#' cps <- dp_to_cp(xi = dps$xi, Psi = dps$Psi, alpha = dps$alpha)
#'
#' stopifnot(all.equal(mu   , cps$mu))
#' stopifnot(all.equal(Sig  , cps$Sigma))
#' stopifnot(all.equal(gamma, cps$gamma))
#'
#' @export
cp_to_dp <- function(mu, Sigma, gamma){
  K <- length(mu)
  stopifnot(NROW(Sigma) == K, NCOL(Sigma) == K,
            length(gamma) == K)

  if(!is.matrix(Sigma) && length(Sigma) == 1L)
    Sigma <- as.matrix(Sigma)

  sign_gamma <- ifelse(gamma > 0, 1, -1)
  gamma <- abs(gamma)
  cf <- (2 * gamma / (4 - pi))^(1/3)

  delta <- sign_gamma * cf / sqrt(1 + cf^2) * sqrt(pi / 2)
  sig_bar <- sqrt(1 - 2 / pi  * delta^2)
  psi <- sqrt(diag(Sigma)) / sig_bar

  k <- delta * psi
  xi <- mu - sqrt(2 / pi) * k
  Psi <- Sigma + 2 / pi * outer(k, k)
  rho <- solve(Psi, k)
  denom <- drop(1 - k %*% rho)
  if(denom > 0.)
    alpha <- psi * rho / sqrt(denom)
  else {
    warning("cp_to_dp: invalid gamma parameter")
    alpha_max <- 20. # TODO: some better value?
    alpha <- sign(rho) * alpha_max
  }

  list(xi = xi, Psi = Psi, alpha = alpha)
}

#' @rdname cp_to_dp
#' @export
dp_to_cp <- function(xi, Psi, alpha){
  if(!is.matrix(Psi) && length(Psi) == 1L)
    Psi <- as.matrix(Psi)
  rho <- alpha / sqrt(diag(Psi))
  k <- drop(Psi %*% rho)
  k <- k / sqrt(drop(1 + rho %*% k))

  mu <- xi + sqrt(2 / pi) * k

  Sigma <- Psi - 2 / pi * outer(k, k)

  gamma <-  sqrt(2 / pi) * k / sqrt(diag(Psi))
  gamma <- (4 - pi) / 2 * gamma^3 / (1 - gamma^2)^(3 / 2)

  list(mu = mu, Sigma = Sigma, gamma = gamma)
}

#' @importFrom lbfgs lbfgs
.opt_default <- function(
  par, fn, gr, ..., control = list(
    reltol = sqrt(.Machine$double.eps), maxit = 500L, trace = 0L)){
  if(length(par) > 1000L){
    # Current statues codes are (https://github.com/chokkan/liblbfgs/blob/7fc787678e4a7f02eaef1c21b36b9bc3bcc0d39b/include/lbfgs.h#L75-L146)
    #    LBFGS_SUCCESS = 0
    #    LBFGS_CONVERGENCE = 0
    #    https://github.com/chokkan/liblbfgs/blob/7fc787678e4a7f02eaef1c21b36b9bc3bcc0d39b/lib/lbfgs.c#L650-L651
    #    "Success: met stopping criteria (ftol)."
    #    LBFGS_STOP 1
    #    /** The initial variables already minimize the objective function. */
    #    LBFGS_ALREADY_MINIMIZED 2
    #
    #    /** Unknown error. */
    #    LBFGSERR_UNKNOWNERROR = -1024
    #    /** Logic error. */
    #    LBFGSERR_LOGICERROR -1023
    #    /** Insufficient memory. */
    #    LBFGSERR_OUTOFMEMORY -1022
    #    /** The minimization process has been canceled. */
    #    LBFGSERR_CANCELED -1021
    #    /** Invalid number of variables specified. */
    #    LBFGSERR_INVALID_N -1020
    #    /** Invalid number of variables (for SSE) specified. */
    #    LBFGSERR_INVALID_N_SSE -1019
    #    /** The array x must be aligned to 16 (for SSE). */
    #    LBFGSERR_INVALID_X_SSE -1018
    #    /** Invalid parameter lbfgs_parameter_t::epsilon specified. */
    #    LBFGSERR_INVALID_EPSILON -1017
    #    /** Invalid parameter lbfgs_parameter_t::past specified. */
    #    LBFGSERR_INVALID_TESTPERIOD -1016
    #    /** Invalid parameter lbfgs_parameter_t::delta specified. */
    #    LBFGSERR_INVALID_DELTA -1015
    #    /** Invalid parameter lbfgs_parameter_t::linesearch specified. */
    #    LBFGSERR_INVALID_LINESEARCH -1014
    #    /** Invalid parameter lbfgs_parameter_t::max_step specified. */
    #    LBFGSERR_INVALID_MINSTEP -1013
    #    /** Invalid parameter lbfgs_parameter_t::max_step specified. */
    #    LBFGSERR_INVALID_MAXSTEP -1012
    #    /** Invalid parameter lbfgs_parameter_t::ftol specified. */
    #    LBFGSERR_INVALID_FTOL -1011
    #    /** Invalid parameter lbfgs_parameter_t::wolfe specified. */
    #    LBFGSERR_INVALID_WOLFE -1010
    #    /** Invalid parameter lbfgs_parameter_t::gtol specified. */
    #    LBFGSERR_INVALID_GTOL -1009
    #    /** Invalid parameter lbfgs_parameter_t::xtol specified. */
    #    LBFGSERR_INVALID_XTOL -1008
    #    /** Invalid parameter lbfgs_parameter_t::max_linesearch specified. */
    #    LBFGSERR_INVALID_MAXLINESEARCH -1007
    #    /** Invalid parameter lbfgs_parameter_t::orthantwise_c specified. */
    #    LBFGSERR_INVALID_ORTHANTWISE -1006
    #    /** Invalid parameter lbfgs_parameter_t::orthantwise_start specified. */
    #    LBFGSERR_INVALID_ORTHANTWISE_START -1005
    #    /** Invalid parameter lbfgs_parameter_t::orthantwise_end specified. */
    #    LBFGSERR_INVALID_ORTHANTWISE_END -1004
    #    /** The line-search step went out of the interval of uncertainty. */
    #    LBFGSERR_OUTOFINTERVAL -1003
    #    /** A logic error occurred; alternatively, the interval of uncertainty
    #        became too small. */
    #    LBFGSERR_INCORRECT_TMINMAX -1002
    #    /** A rounding error occurred; alternatively, no line-search step
    #        satisfies the sufficient decrease and curvature conditions. */
    #    LBFGSERR_ROUNDING_ERROR -1001
    #    /** The line-search step became smaller than lbfgs_parameter_t::min_step. */
    #    LBFGSERR_MINIMUMSTEP -1000
    #    /** The line-search step became larger than lbfgs_parameter_t::max_step. */
    #    LBFGSERR_MAXIMUMSTEP -999
    #    /** The line-search routine reaches the maximum number of evaluations. */
    #    LBFGSERR_MAXIMUMLINESEARCH -998
    #    /** The algorithm routine reaches the maximum number of iterations. */
    #    LBFGSERR_MAXIMUMITERATION -997
    #    /** Relative width of the interval of uncertainty is at most
    #        lbfgs_parameter_t::xtol. */
    #    LBFGSERR_WIDTHTOOSMALL -996
    #    /** A logic error (negative line-search step) occurred. */
    #    LBFGSERR_INVALIDPARAMETERS -995
    #    /** The current search direction increases the objective function value. */
    #    LBFGSERR_INCREASEGRADIENT -994
    delta <- if(!is.null(control$reltol))
      control$reltol else sqrt(.Machine$double.eps)
    max_iterations <- if(!is.null(control$maxit))
      control$maxit else formals(lbfgs)$max_iterations
    invisible <- if(!is.null(control$trace))
      as.integer(control$trace == 0) else 1L
    out <- lbfgs(
      call_eval = fn, call_grad = gr, vars = par, invisible = invisible,
      m = 6L, epsilon = 1e-5, delta = delta, past = 50L,
      max_iterations = max_iterations, max_linesearch = 100L,
      linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_WOLFE")
    names(out$par) <- names(par)
    out$optimizer <- "lbfgs"
    out$ok <- out$convergence >= 0
    return(out)
  }

  cl <- match.call()
  cl[[1L]] <- quote(stats::optim)
  if(is.null(cl$method))
    cl$method <- "BFGS"

  out <- eval(cl, parent.frame())
  out$optimizer <- "optim"
  out$ok <- out$convergence %in% 0:1
  out
}



#' Maps Between a Covariance Matrix and Its Log-Cholesky Parametrization
#'
#' @description
#' Let \eqn{C} be a \eqn{K}-dimensional covariance matrix. Let
#' \eqn{C = LL^\top} where \eqn{L} is the Cholesky decomposition.
#' Then the two functions map between \eqn{C} and
#' a vector containing the non-zero elements of L where the diagonal
#' entries have been log transformed.
#'
#' @param cov the covariance matrix.
#' @param theta numeric vector with log-cholesky parametrization.
#'
#' @examples
#' set.seed(1)
#' Sigma <- drop(rWishart(1, 4, diag(4)))
#' log_chol <- cov_to_theta(Sigma)
#' log_chol
#' stopifnot(all.equal(Sigma, theta_to_cov(log_chol)))
#'
#' @export
cov_to_theta <- function(cov){
  n_rng <- NCOL(cov)
  ch <- t(chol(cov))
  diag(ch) <- log(diag(ch))
  keep <- lower.tri(ch, TRUE)
  lower_tri <- ch[keep]
  names(lower_tri) <-
    outer(1:n_rng, 1:n_rng, function(x, y) paste0("L", x, ".", y))[keep]
  lower_tri
}

#' @rdname cov_to_theta
#' @importFrom utils head tail
#' @export
theta_to_cov <- function(theta){
  dim <- .5 * (sqrt(8 * length(theta) + 1) - 1)
  if(dim < 2L)
    return(exp(2 * theta))

  L <- matrix(0, dim, dim)
  L[lower.tri(L, TRUE)] <- theta
  diag(L) <- exp(diag(L))
  tcrossprod(L)
}

.set_use_own_VA_method <- with(new.env(), {
  .use_own_VA_method <- TRUE
  function(use_own_VA_method){
    stopifnot(
      is.logical(use_own_VA_method),
      length(use_own_VA_method) == 1L,
      !is.na(use_own_VA_method))
    .use_own_VA_method <<- use_own_VA_method
    invisible(use_own_VA_method)
  }
})

.get_use_own_VA_method <- with(
  environment(.set_use_own_VA_method),
  function()
    .use_own_VA_method)

#' Free Memory Used by TMB
#' @description
#' Essentially a wrapper around \code{\link{FreeADFun}}.
#'
#' @param obj object of type MGSM_ADFun.
#'
#' @importFrom TMB FreeADFun
#' @seealso \code{\link{make_mgsm_ADFun}}
#' @export
free_laplace <- function(obj){
  stopifnot(inherits(obj, "MGSM_ADFun"))

  lap <- obj$laplace
  if(!is.null(lap))
    FreeADFun(lap)

  invisible()
}

#' Construct Objective Functions with Derivatives for a Mixed Generalized
#' Survival Model
#'
#' @description
#' Constructs an objective function for a particular generalized survival
#' model applied to a given data set.
#'
#' @param formula two-sided \code{\link{formula}} where the left-hand side is a
#'                \code{\link{Surv}} object and the right-hand side is the
#'                fixed effects.
#' @param data \code{\link{data.frame}} with variables used in the model.
#' @param df integer with the degrees of freedom used for the baseline
#'           spline.
#' @param tformula \code{\link{formula}} with baseline survival function.
#'                 The time variable must be the same
#'                 symbol as used in the left-hand-side of \code{formula}.
#'                 \code{NULL} implies that \code{df} is passed to
#'                 \code{\link{nsx}}.
#' @param Z one-sided \code{\link{formula}} where the right-hand side are
#'          the random effects.
#' @param cluster vector with integers or factors for group identifiers
#'                (one for each observation).
#' @param do_setup character vector indicating which approximation to setup.
#'                 See 'Details'.
#' @param n_nodes integer with the number of nodes to use in (adaptive)
#'                Gauss-Hermite quadrature.
#' @param link character specifying the link function.
#' @param param_type characters for the parameterization used with the SNVA.
#'                   See 'Details'.
#' @param theta starting values for covariance matrix.
#' @param beta starting values for fixed effect coefficients.
#' @param opt_func general optimization function to use. It
#'                 needs to have an interface like \code{\link{optim}}.
#' @param n_threads integer with number of threads to use.
#' @param skew_start starting value for the Pearson's moment coefficient of
#'                   skewness parameter when a SNVA is used. Currently,
#'                   a somewhat arbitrary value.
#' @param dense_hess logical for whether to make dense Hessian computation
#'                   available. Memory and computation time is saved if it is
#'                   \code{FALSE}.
#' @param sparse_hess logical for whether to make sparse Hessian computation
#'                    available. Memory and computation time is saved if it is
#'                    \code{FALSE}.
#' @param kappa numeric scalar with the penalty in the relaxed problem
#' ensuring the monotonicity of the survival curve.
#' @param dtformula \code{\link{formula}} with the derivative of the
#' baseline survival function.
#'
#' @details
#' Possible link functions for \code{link} are:
#' \code{"PH"} for proportional hazard,
#' \code{"PO"} for proportional odds, and \code{"probit"} for a probit link
#' function.
#'
#' The available estimation methods for \code{do_setup} are:
#' \code{"Laplace"} for a Laplace approximation,
#' \code{"GVA"} for a Gaussian variational approximation (GVA),
#' and \code{"SNVA"} for a skew-normal variational approximation (SNVA).
#'
#' The parameterizations for the SNVA are selected by \code{param_type}.
#' Possible arguments are: \code{"DP"} for direct parameterization,
#' \code{"CP"} for centralized parameterization, and
#' \code{"CP_trans"} for transformed centralized parameterization where the
#' skew parameters are transformed by a logistic function to be in the
#' appropriate range.
#'
#' See the README \url{https://github.com/boennecd/survTMB} for
#' further information and examples.
#'
#' @return
#' An object of class \code{MGSM_ADFun}. The elements are:
#' \item{laplace}{object to perform a Laplace approximation if \code{do_setup} contains \code{"Laplace"}.}
#' \item{gva}{object to perform a GVA if \code{do_setup} contains \code{"GVA"}.}
#' \item{snva}{object to perform a SNVA if \code{do_setup} contains \code{"SNVA"}.}
#' \item{y}{numeric vector with outcomes.}
#' \item{event}{numeric vector with event indicators.}
#' \item{X}{fixed effect design matrix.}
#' \item{XD}{derivative of fixed effect design matrix with respect to time.}
#' \item{Z}{Random effect design matrix.}
#' \item{grp}{integer vector with group identifier.}
#' \item{terms}{\code{\link{list}} with \code{\link{terms.object}}s.}
#' \item{link}{character with the link function.}
#' \item{cl}{matched call.}
#' \item{dense_hess}{dense_hess argument.}
#' \item{sparse_hess}{sparse_hess argument.}
#'
#' @examples
#' library(survTMB)
#' if(require(coxme)){
#'   # construct function with a random intercept and a proportional hazard
#'   # link function
#'   func <- make_mgsm_ADFun(
#'     Surv(y, uncens) ~ trt, cluster = as.factor(center), Z = ~ 1,
#'     df = 3L, data = eortc, link = "PH",
#'     do_setup = c("Laplace", "GVA", "SNVA"), n_threads = 1L)
#'   print(func)
#' }
#'
#' @seealso
#' \code{\link{fit_mgsm}}
#'
#' @importFrom TMB MakeADFun
#' @importFrom stats model.frame model.response terms model.matrix lm lm.fit predict qnorm sd
#' @importFrom rstpm2 nsx
#' @importFrom survival coxph frailty Surv
#' @importFrom Matrix sparseMatrix
#' @export
make_mgsm_ADFun <- function(
  formula, data, df, tformula = NULL, Z, cluster,
  do_setup = c("Laplace", "GVA", "SNVA"), n_nodes = 20L,
  param_type = c("DP", "CP_trans", "CP"), link = c("PH", "PO", "probit"),
  theta = NULL, beta = NULL, opt_func = .opt_default, n_threads = 1L,
  skew_start = -.0001, dense_hess = FALSE,
  sparse_hess = FALSE, kappa = .MGSM_default_kappa,
  dtformula = NULL){
  #####
  # checks
  param_type <- param_type[1]
  stopifnot(
    param_type %in% c("DP", "CP_trans", "CP"),
    is.logical(dense_hess), length(dense_hess) == 1L, !is.na(dense_hess),
    is.logical(sparse_hess), length(sparse_hess) == 1L, !is.na(sparse_hess),
    is.function(opt_func),
    all(do_setup %in% c(.laplace_char, .gva_char, .snva_char)))

  #####
  # get arguments to pass on
  cluster <- substitute(cluster)
  args_pass <- mgsm_setup(
    formula = formula, data = data, df = df, tformula = tformula, Z = Z,
    cluster = cluster, n_nodes = n_nodes,
    link = link, theta = theta, beta = beta, opt_func = opt_func,
    n_threads = n_threads, skew_start = skew_start, kappa = kappa,
    dtformula = dtformula)

  #####
  # setup ADFun object for the Laplace approximation
  data_ad_func <- list(
    tobs = args_pass$tobs, event = args_pass$event, X = args_pass$X,
    XD = args_pass$XD, Z = args_pass$Z, grp = args_pass$grp,
    link = args_pass$link, grp_size = args_pass$grp_size,
    n_threads = args_pass$n_threads)

  # assign parameter list
  params <- list(eps = args_pass$eps, kappa = args_pass$kappa,
                 b = args_pass$b, theta = args_pass$theta)

  laplace_out <- if(.laplace_char %in% do_setup)
    .get_laplace_func(data_ad_func, params, args_pass$n_rng,
                      args_pass$n_grp, args_pass$inits)
  else
    NULL

  #####
  # setup ADFun object for the GVA
  gva_out <- if(.gva_char %in% do_setup)
    .get_gva_func(args_pass$n_rng, args_pass$n_grp, params, data_ad_func,
                  args_pass$n_nodes, dense_hess,
                  sparse_hess, args_pass$inits, opt_func)
  else
    NULL

  #####
  # setup ADFun object for the SNVA
  snva_out <- if(.snva_char %in% do_setup)
    .get_snva_out(args_pass$n_rng, args_pass$n_grp, params,
                  args_pass$skew_start, param_type,
                  data_ad_func, args_pass$n_nodes, dense_hess,
                  sparse_hess, opt_func, gva_obj = gva_out,
                  inits = args_pass$inits)
  else
    NULL

  structure(
    list(laplace = laplace_out, gva = gva_out, snva = snva_out,
         y = args_pass$y, event = args_pass$event, X = args_pass$X,
         XD = args_pass$XD, Z = args_pass$Z, grp = args_pass$grp,
         terms = args_pass$terms, cl = match.call(),
         link = args_pass$link, opt_func = opt_func,
         dense_hess = dense_hess, sparse_hess = sparse_hess),
    class = "MGSM_ADFun")
}

mgsm_setup <- function(
  formula, data, df, tformula, Z, cluster, n_nodes,
  link, theta, beta, opt_func, n_threads, skew_start, kappa,
  dtformula){
  link <- link[1]
  stopifnot(
    is.integer(df), df > 0L, inherits(formula, "formula"),
    is.null(tformula) || inherits(tformula, "formula"),
    inherits(Z, "formula"), is.data.frame(data),
    !missing(cluster),
    is.integer(n_nodes), length(n_nodes) == 1L, n_nodes > 0L,
    link %in% c("PH", "PO", "probit"),
    is.null(theta) || (is.numeric(theta) && all(is.finite(theta))),
    is.null(beta) || (is.numeric(beta) && all(is.finite(beta))),
    is.function(opt_func),
    is.integer(n_threads) && n_threads > 0L && length(n_threads) == 1L,
    is.numeric(skew_start), length(skew_start) == 1L,
    is.numeric(kappa), length(kappa) == 1L, is.finite(kappa),
    is.null(dtformula) || (!is.null(tformula) &&
                             inherits(dtformula, "formula")))
  eval(bquote(stopifnot(
    .(-.skew_boundary) < skew_start && skew_start < .(.skew_boundary))))

  eps <- .MGSM_defaul_eps

  #####
  # get the cluster variable
  grp <- eval(cluster, data)
  stopifnot(is.factor(grp) || is.integer(grp) || is.character(grp),
            isTRUE(length(grp) == NROW(data)))
  if(!is.factor(grp))
    grp <- as.factor(grp)
  grp <- as.integer(grp)
  n_grp <- length(unique(grp))

  #####
  # change order of data set and assign group size variables
  ord <- order(grp)
  data <- data[ord, ]
  grp <- grp[ord]
  grp_size <- table(grp)

  #####
  # get design matrices and outcomes
  need_theta <- is.null(theta)
  need_beta  <- is.null(beta)

  gsm_output <- gsm(formula = formula, tformula = tformula, data = data,
                    df = df, link = link, n_threads = n_threads,
                    do_fit = need_beta, opt_func = opt_func, eps = eps,
                    kappa = kappa, dtformula = dtformula)

  mt_X <- gsm_output$mt_Z
  X <- gsm_output$Z
  n_x_fix <- NCOL(X)
  is_fix <- seq_len(n_x_fix)

  # add time-varying baseline
  mt_b <- gsm_output$mt_X
  X <- cbind(X, gsm_output$X)
  XD <- cbind(matrix(0., NROW(X), n_x_fix), gsm_output$XD)
  colnames(XD) <- colnames(X)

  # get outcome
  y <- gsm_output$y
  stopifnot(inherits(y, "Surv"), isTRUE(attr(y, "type") == "right"))
  event <- y[, 2]
  tobs  <- y[, 1]

  # get random effect covariates
  formula_Z <- Z
  mf_Z <- model.frame(formula_Z, data = data)
  mt_Z <- terms(mf_Z)
  Z <- model.matrix(mt_Z, mf_Z)
  n_rng <- NCOL(Z)
  stopifnot(NCOL(X) > 0, NCOL(Z) > 0, NROW(X) == NROW(Z),
            length(y) == NROW(X))

  #####
  # get starting values
  inits <- list()
  if(need_beta){
    gsm_est <- gsm_output$fit

    inits$coef <- numeric(NCOL(X))
    inits$coef[-is_fix] <- gsm_est$beta
    inits$coef[ is_fix] <- gsm_est$gamma
  } else {
    stopifnot(length(beta) == NCOL(X))
    inits$coef <- beta
  }
  names(inits$coef) <- colnames(X)
  rm(gsm_output)

  # the user may have provided values
  theta <- if(!need_theta){
    stopifnot(all(dim(theta) == n_rng))
    theta

  } else local({
    sds <- apply(Z, 2, sd)
    sds[sds == 0] <- 1
    diag((.1 / sds)^2, length(sds))
  })
  theta <- cov_to_theta(theta)
  theta <- setNames(theta, paste0("theta:", names(theta)))
  beta <- inits$coef

  list(tobs = tobs, event = event, X = X, XD = XD, Z = Z,
       grp = grp - 1L, link = link, grp_size = grp_size,
       n_threads = n_threads, eps = eps, kappa = kappa, b = beta,
       theta = theta, terms = list(
         X = mt_X, Z = mt_Z, baseline = mt_b), n_rng = n_rng, n_grp = n_grp,
       inits = inits, n_nodes = n_nodes, y = y, skew_start = skew_start)
}

.tmb_set_n_threads <- function(n, is_laplace = FALSE){
  out <- if(is_laplace || !.get_use_own_VA_method())
    set_n_threads(n) else set_n_threads(-1L)
  out
}

.get_MGSM_VA_start <- function(
  n_rng, params, data_ad_func, opt_func, skew_start = NULL, is_cp = NULL){
  # set the initial values
  grp_end <- cumsum(data_ad_func$grp_size)
  grp_start <- c(1L, head(grp_end, -1) + 1L)
  b <- params$b
  sig <- theta_to_cov(params$theta)
  if(!is.matrix(sig))
    sig <- as.matrix(sig)
  sig_inv <- solve(sig)
  chol_sig_inv <- chol(sig_inv)

  is_snva <- !is.null(skew_start) && !is.null(is_cp)
  theta_VA <- mapply(function(istart, iend){
    # get the data we need
    idx <- istart:iend
    n <- length(idx)
    X  <- t(data_ad_func$X [idx, ])
    XD <- t(data_ad_func$XD[idx, ])
    Z  <- t(data_ad_func$Z[idx, ])
    X_arg <- matrix(nrow = 0, ncol = n)
    y <- data_ad_func$event[idx]

    # get the offset
    if(length(b) > 0){
      offset_eta  <- drop(b %*% X)
      offset_etaD <- drop(b %*% XD)
    } else
      offset_eta <- offset_etaD <- numeric(n)

    # make Taylor approximation
    opt_obj <- get_gsm_pointer(
      X = X_arg, XD = X_arg, Z = Z, y = y, eps = params$eps,
      kappa = params$kappa, link = data_ad_func$link,
      n_threads = 1L, offset_eta = offset_eta, offset_etaD = offset_etaD)

    fn <- function(x, ...)
      -gsm_eval_ll(ptr = opt_obj, beta = numeric(), gamma = x) +
      sum((chol_sig_inv %*% x)^2) / 2
    gr <- function(x, ...)
      -gsm_eval_grad(ptr = opt_obj, beta = numeric(), gamma = x) +
      drop(sig_inv %*% x)
    he <- function(x, ...)
      -gsm_eval_hess(ptr = opt_obj, beta = numeric(), gamma = x) + sig_inv

    opt_ret <- opt_func(numeric(NCOL(sig)), fn, gr)
    mu <- opt_ret$par
    sig_use <- solve(he(mu))

    if(is_snva)
      # SNVA
      if(is_cp){
        # centralized parameters
        get_skew_trans <- function(x)
          log((.skew_boundary + x) / (.skew_boundary - x))

        dp_pars <- cp_to_dp(mu = mu, Sigma = sig_use, gamma = skew_start)
        cp_pars <- dp_to_cp(xi = dp_pars$xi, Psi = dp_pars$Psi,
                            alpha = dp_pars$alpha)
        c(cp_pars$mu, cov_to_theta(cp_pars$Sigma),
          get_skew_trans(cp_pars$gamma))

      } else {
        # direct parameters
        dp_pars <- cp_to_dp(mu = mu, Sigma = sig_use, gamma = skew_start)
        c(dp_pars$xi, cov_to_theta(dp_pars$Psi), dp_pars$alpha)

      }
    else
      # GVA
      c(mu, cov_to_theta(sig_use))
  }, istart = grp_start, iend = grp_end)
  c(theta_VA)
}

.get_laplace_func <- function(data_ad_func, params, n_rng, n_grp, inits) {
  setup_atomic_cache(
    n_nodes = 15L, type = .laplace_char, link = data_ad_func$link)

  # get Laplace AD function
  .tmb_set_n_threads(data_ad_func$n_threads, is_laplace = TRUE)
  data_ad_func <- c(
    list(app_type = .laplace_char), data_ad_func)
  # TODO: initialize random effects in a smarter way...
  params$u <- matrix(0., n_rng, n_grp)
  adfunc_laplace <- MakeADFun(
    data = data_ad_func, parameters = params, DLL = "survTMB",
    silent = TRUE, random = "u")

  # we make a wrapper object to account for the eps and kappa and allow the
  # user to change these
  n_threads <- data_ad_func$n_threads
  eps <- adfunc_laplace$par["eps"]
  kappa <- adfunc_laplace$par["kappa"]
  fn <- adfunc_laplace$fn
  gr <- adfunc_laplace$gr
  he <- adfunc_laplace$he
  get_x <- function(x)
    c(eps = eps, kappa = kappa, x)

  out <- adfunc_laplace[
    !names(adfunc_laplace) %in% c("par", "fn", "gr", "he")]

  par <- adfunc_laplace$par[-(1:2)]
  names(par)[seq_along(inits$coef)] <- names(inits$coef)

  out <- c(list(
      par = par,
      fn = function(x, ...){
        .tmb_set_n_threads(n_threads, is_laplace = TRUE)
        fn(get_x(x))
      },
      gr = function(x, ...){
        .tmb_set_n_threads(n_threads, is_laplace = TRUE)
        gr(get_x(x))[-(1:2)]
      },
      he = function(x, ...){
        .tmb_set_n_threads(n_threads, is_laplace = TRUE)
        he(get_x(x))[-(1:2), -(1:2), drop = FALSE]
      },
      # function to set penalty parameters
      update_pen = function(eps, kappa){
        p_env <- parent.env(environment())
        if(!missing(eps))
          assign("eps"  , eps  , p_env)
        if(!missing(kappa))
          assign("kappa", kappa, p_env)
        invisible(with(p_env, c(eps = eps, kappa = kappa)))
      },
      adfun = adfunc_laplace
    ), out)

  # clean up and return
  rm(list = setdiff(ls(), c(
    "out", "n_threads", "fn", "gr", "he", "get_x", "eps", "kappa")))
  out
}

# VA util funcs
.get_par_va <- function(params)
  with(params, c(c(eps = eps, kappa = kappa), b, theta, theta_VA))

.eval_hess_sparse <- function(ptr, par){
  out <- VA_funcs_eval_hess_sparse(ptr, par)
  Matrix::sparseMatrix(
    i = out$row_idx + 1L, j = out$col_idx + 1L, x = out$val,
    symmetric = TRUE)
}

.get_gva_func <- function(n_rng, n_grp, params, data_ad_func, n_nodes,
                          dense_hess, sparse_hess, inits, opt_func) {
  setup_atomic_cache(
    n_nodes = n_nodes, type = .gva_char, link = data_ad_func$link)

  # set names
  theta_VA <- .get_MGSM_VA_start(
    n_rng = n_rng, params = params, data_ad_func = data_ad_func,
    opt_func = opt_func)
  names(theta_VA) <- theta_VA_names <-
    mgsm_get_gva_names(n_rng, n_grp, params$theta)

  get_gva_out <- function(theta_VA){
    adfunc_VA <- local({
      data_ad_func <- c(
        list(app_type = .gva_char), data_ad_func,
        list(n_nodes = n_nodes, dense_hess = dense_hess,
             sparse_hess = sparse_hess))
      params$theta_VA <- theta_VA

      if(.get_use_own_VA_method()){
        data_ad_func$param_type <- "DP"
        within(list(), {
          ptr <- get_VA_funcs(data = data_ad_func, parameters = params)
          fn <- function(par)
            VA_funcs_eval_lb(ptr, par)
          gr <- function(par)
            drop(VA_funcs_eval_grad(ptr, par))
          he <- function(par)
            VA_funcs_eval_hess(ptr, par)
          he_sp <- function(x, ...)
            .eval_hess_sparse(ptr, x)
          par <- .get_par_va(params)
        })

      } else {
        .tmb_set_n_threads(data_ad_func$n_threads, is_laplace = FALSE)
        within(MakeADFun(
          data = data_ad_func, parameters = params, DLL = "survTMB",
          silent = TRUE),
          he_sp <- function(x, ...)
            stop())
      }
    })

    # we make a wrapper object to account for the eps and kappa and allow the
    # user to change these
    gva_out <- with(new.env(), {
      n_threads <- data_ad_func$n_threads
      eps <- adfunc_VA$par["eps"]
      kappa <- adfunc_VA$par["kappa"]
      fn    <- adfunc_VA$fn
      gr    <- adfunc_VA$gr
      he    <- adfunc_VA$he
      he_sp <- adfunc_VA$he_sp
      get_x <- function(x)
        c(eps = eps, kappa = kappa, x)

      out <- adfunc_VA[
        !names(adfunc_VA) %in% c("par", "fn", "gr", "he", "he_sp")]

      par <- adfunc_VA$par[-(1:2)]
      names(par)[seq_along(inits$coef)] <- names(inits$coef)
      idx_va <- (length(par) - length(theta_VA_names) + 1):length(par)
      names(par)[idx_va] <-
        theta_VA_names

      c(list(
        par = par,
        fn = function(x, ...){
          .tmb_set_n_threads(n_threads, is_laplace = FALSE)
          fn(get_x(x))
        },
        gr = function(x, ...){
          .tmb_set_n_threads(n_threads, is_laplace = FALSE)
          gr(get_x(x))[-(1:2)]
        },
        he = function(x, ...){
          .tmb_set_n_threads(n_threads, is_laplace = FALSE)
          he(get_x(x))[-(1:2), -(1:2), drop = FALSE]
        },
        he_sp = function(x, ...){
          he_sp(get_x(x))[-(1:2), -(1:2), drop = FALSE]
        },
        # function to set penalty parameters
        update_pen = function(eps, kappa){
          p_env <- parent.env(environment())
          if(!missing(eps))
            assign("eps"  , eps  , p_env)
          if(!missing(kappa))
            assign("kappa", kappa, p_env)
          invisible(with(p_env, c(eps = eps, kappa = kappa)))
        },
        # function to get parameters
        get_params = function(x)
          x[-idx_va],
        control = list(maxit = 1000L),
        adfun = adfunc_VA
      ), out)
    })
  }

  get_gva_out(theta_VA)
}

mgsm_get_gva_names <- function(n_rng, n_grp, theta){
  theta_VA_names <- c(paste0("mu", 1:n_rng), names(theta))
  c(outer(
    theta_VA_names, paste0("g", 1:n_grp), function(x, y)
      paste0(y, ":", x)))
}

.get_snva_out <- function(
  n_rng, n_grp, params, skew_start, param_type, data_ad_func,
  n_nodes, dense_hess, sparse_hess, opt_func, gva_obj, inits) {
  # setup cache
  setup_atomic_cache(
    n_nodes = n_nodes, type = .snva_char, link = data_ad_func$link)

  #####
  # setup starting values for VA parameters
  # get GVA object if needed
  if(is.null(gva_obj))
    gva_obj <- .get_gva_func(
      n_rng, n_grp, params, data_ad_func, n_nodes, dense_hess = FALSE,
      sparse_hess = FALSE, inits, opt_func)

  # set the initial values
  n_mu     <- n_rng
  n_rho    <- n_rng
  n_Lambda <- (n_rng * (n_rng + 1L)) / 2L
  n_p_grp  <- n_mu + n_Lambda + n_rho
  theta_VA <- rep(NA_real_, n_p_grp * n_grp)

  # find gva solution
  gva_opt <- with(gva_obj, opt_func(par, fn, gr,
                                    control = list(maxit = 1000L)))
  beta  <- params$b     <-
    gva_opt$par[1:length(params$b)]
  theta <- params$theta <-
    gva_opt$par[1:length(params$theta) + length(beta)]
  gva_va_vals <- gva_opt$par[
    -seq_len(length(beta) + length(theta))]

  if(length(skew_start) == 1L && n_rng > 1L)
    skew_start <- rep(skew_start, n_rng)

  n_p_grp_gva <-  n_p_grp - n_rho
  theta_VA <- if(param_type == "DP"){
    n_lower_tri <- (n_rng * (n_rng - 1L)) / 2L

    vapply(1:n_grp, function(i){
      # setup mean and VA variance
      gva_par <- gva_va_vals[(i - 1L) * n_p_grp_gva + 1:n_p_grp_gva]
      Sig <- theta_to_cov(gva_par[-seq_len(n_mu)])

      dp_pars <- cp_to_dp(mu = gva_par[1:n_mu], Sigma = Sig,
                          gamma = skew_start)

      c(dp_pars$xi, cov_to_theta(dp_pars$Psi), dp_pars$alpha)
    }, numeric(n_p_grp), USE.NAMES = FALSE)

  } else {
    stopifnot(param_type == "CP_trans")
    get_skew_trans <- function(x)
      log((.skew_boundary + x) / (.skew_boundary - x))

    vapply(1:n_grp, function(i){
      # setup mean and VA variance
      gva_par <- gva_va_vals[(i - 1L) * n_p_grp_gva + 1:n_p_grp_gva]

      Sig <- theta_to_cov(gva_par[-seq_len(n_mu)])
      dp_pars <- cp_to_dp(mu = gva_par[1:n_mu], Sigma = Sig,
                          gamma = skew_start)

      cp_pars <- dp_to_cp(xi = dp_pars$xi, Psi = dp_pars$Psi,
                          alpha = dp_pars$alpha)

      out <- c(cp_pars$mu, cov_to_theta(cp_pars$Sigma),
               get_skew_trans(cp_pars$gamma))
    }, numeric(n_p_grp), USE.NAMES = FALSE)

  }

  # set names
  theta_VA_names <- mgsm_get_snva_names(n_rng, n_grp, param_type)
  theta_VA <- structure(c(theta_VA), names = theta_VA_names)

  adfunc_VA <- local({
    data_ad_func <- c(
      list(app_type = .snva_char), data_ad_func,
      list(n_nodes = n_nodes, param_type = param_type,
           dense_hess = dense_hess, sparse_hess = sparse_hess))
    params$theta_VA <- theta_VA

    if(.get_use_own_VA_method())
      within(list(), {
        ptr <- get_VA_funcs(data = data_ad_func, parameters = params)
        fn <- function(par)
          VA_funcs_eval_lb(ptr, par)
        gr <- function(par)
          drop(VA_funcs_eval_grad(ptr, par))
        he <- function(par)
          VA_funcs_eval_hess(ptr, par)
        he_sp <- function(x, ...)
          .eval_hess_sparse(ptr, x)
        par <- .get_par_va(params)
      })
    else {
      .tmb_set_n_threads(data_ad_func$n_threads, is_laplace = FALSE)
      within(MakeADFun(
        data = data_ad_func, parameters = params, DLL = "survTMB",
        silent = TRUE),
        he_sp <- function(x, ...)
          stop())
    }
  })

  # we make a wrapper object to account for the eps and kappa and allow the
  # user to change these
  beta <- params$b
  theta <- params$theta
  get_snva_out <- function(theta_VA){
    with(new.env(), {
      n_threads <- data_ad_func$n_threads
      eps <- adfunc_VA$par["eps"]
      kappa <- adfunc_VA$par["kappa"]
      fn <- adfunc_VA$fn
      gr <- adfunc_VA$gr
      he <- adfunc_VA$he
      he_sp <- adfunc_VA$he_sp
      get_x <- function(x)
        c(eps = eps, kappa = kappa, x)

      out <- adfunc_VA[
        !names(adfunc_VA) %in% c("par", "fn", "gr", "he", "he_sp")]

      par <- adfunc_VA$par[-(1:2)]
      names(par)[seq_along(beta)] <- names(beta)
      idx_va <- (length(par) - length(theta_VA_names) + 1):length(par)
      names(par)[idx_va] <- theta_VA_names

      c(list(
        par = par,
        fn = function(x, ...){
          .tmb_set_n_threads(n_threads, is_laplace = FALSE)
          fn(get_x(x))
        },
        gr = function(x, ...){
          .tmb_set_n_threads(n_threads, is_laplace = FALSE)
          gr(get_x(x))[-(1:2)]
        },
        he = function(x, ...){
          .tmb_set_n_threads(n_threads, is_laplace = FALSE)
          he(get_x(x))[-(1:2), -(1:2), drop = FALSE]
        },
        he_sp = function(x, ...){
          he_sp(get_x(x))[-(1:2), -(1:2), drop = FALSE]
        },
        # function to set penalty parameters
        update_pen = function(eps, kappa){
          p_env <- parent.env(environment())
          if(!missing(eps))
            assign("eps"  , eps  , p_env)
          if(!missing(kappa))
            assign("kappa", kappa, p_env)
          invisible(with(p_env, c(eps = eps, kappa = kappa)))
        },
        # function to get parameters
        get_params = function(x)
          x[-idx_va],
        control = list(maxit = 1000L),
        adfun = adfunc_VA
      ), out)
    })
  }

  # find initial va parameters
  func <- get_snva_out(theta_VA)
  coefs_start <- c(params$b, params$theta)
  do_drop <- seq_along(coefs_start)
  get_x <- function(x)
    c(coefs_start, x)

  fn <- function(x)
    func$fn(get_x(x))
  gr <- function(x)
    func$gr(get_x(x))[-do_drop]

  opt_out <- opt_func(theta_VA, fn = fn, gr = gr,
                      control = list(maxit = 1000L))
  get_snva_out(opt_out$par)
}

mgsm_get_snva_names <- function(n_rng, n_grp, param_type){
  keep <- which(lower.tri(diag(n_rng)))
  L <- outer(
    1:n_rng, 1:n_rng, function(x, y) paste0("L", x, ".", y))[keep]

  skew_name <- switch(param_type,
                      DP = "alpha",
                      `CP_trans` = "skew_trans",
                      stop("unknown 'param_type'"))
  proto <- c(paste0("mu", 1:n_rng), paste0("log_sd", 1:n_rng), L,
             paste0(skew_name, 1:n_rng))

  c(outer(proto, paste0("g", 1:n_grp), function(x, y)
    paste0(y, ":", x)))
}
