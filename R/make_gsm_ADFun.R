.laplace_char <- "Laplace"
.gva_char     <- "GVA"
.snva_char    <- "SNVA"

.MGSM_defaul_eps <- .Machine$double.eps^(1/2)
.MGSM_default_kappa <- 1e8

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
    reltol = sqrt(.Machine$double.eps), maxit = 100L, trace = 0L)){
  if(length(par) > 1000L){
    # Current statues codes are (https://github.com/chokkan/liblbfgs/blob/7fc787678e4a7f02eaef1c21b36b9bc3bcc0d39b/include/lbfgs.h#L75-L146)
    #    LBFGS_SUCCESS = 0
    #    LBFGS_CONVERGENCE = 0
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
      m = 6L, epsilon = 0, delta = delta, past = 1L,
      max_iterations = max_iterations)
    names(out$par) <- names(par)
    return(out)
  }

  cl <- match.call()
  cl[[1L]] <- quote(stats::optim)
  if(is.null(cl$method))
    cl$method <- "BFGS"

  out <- eval(cl, parent.frame())
  out
}

#' @importFrom survival survfit
Shat <- function(obj){
  ## predicted survival for individuals (adjusted for covariates)
  newobj <- survfit(obj, se.fit = FALSE)
  surv <- newobj$surv
  rr <- try(predict(obj, type = "risk"), silent = TRUE)
  if (inherits(rr, "try-error"))
    rr <- 1
  surv2 <- surv[match(obj$y[, ncol(obj$y) - 1], newobj$time)]

  surv2^rr
}

#' Maps Between a Covariance Matrix and Its Log-Cholesky Parametrization
#'
#' @description
#' Let \eqn{C} be a \eqn{K}-dimensional covariance matrix. Let
#' \eqn{C = LL^\top} where \eqn{L} is the Cholesky decomposition
#' and \eqn{L = diag(\psi) H} where  \eqn{H} is a lower triangular matrix
#' with ones in the diagonal. Then the two functions map between \eqn{C} and
#' a vector where the first \eqn{K} elements are the log of \eqn{psi} and
#' the remaining \eqn{K (K - 1) / 2} elements are the non-zero and
#' non-diagonal entries of \eqn{H}.
#'
#' @param cov the covariance matrix.
#' @param theta numeric vector where the
#'              first \eqn{K} elements are the log of \eqn{psi} and the
#'              remaining \eqn{K (K - 1) / 2} elements are the non-zero and
#'              non-diagonal entries of \eqn{H}.
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
  log_sd <- structure(log(diag(ch)), names = paste0("log_sd", 1:n_rng))
  keep <- lower.tri(ch)
  lower_tri <-(diag(diag(ch)^(-1), n_rng) %*% ch)[keep]
  if(n_rng > 1L)
    names(lower_tri) <-
      outer(1:n_rng, 1:n_rng, function(x, y) paste0("L", x, ".", y))[keep]
  c(log_sd, lower_tri)
}

#' @rdname cov_to_theta
#' @importFrom utils head tail
#' @export
theta_to_cov <- function(theta){
  dim <- .5 * (sqrt(8 * length(theta) + 1) - 1)
  if(dim < 2L)
    return(exp(2 * theta))

  L <- diag(dim)
  L[lower.tri(L)] <- tail(theta, -dim)
  L <- diag(exp(head(theta, dim))) %*% L
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
#' \item{gva}{object to perform a GVA if \code{do_setup} contains \code{"GVA"} or \code{"SNVA"}.}
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
#' @importFrom stats model.frame model.response terms model.matrix lm lm.fit predict qnorm
#' @importFrom rstpm2 nsx
#' @importFrom survival coxph frailty Surv
#' @importFrom Matrix sparseMatrix
#' @export
make_mgsm_ADFun <- function(
  formula, data, df, Z, cluster, do_setup = c("Laplace", "GVA", "SNVA"),
  n_nodes = 20L, param_type = c("DP", "CP_trans", "CP"),
  link = c("PH", "PO", "probit"), theta = NULL, beta = NULL,
  opt_func = .opt_default, n_threads = 1L,
  skew_start = -.0001, dense_hess = FALSE,
  sparse_hess = FALSE){
  link <- link[1]
  param_type <- param_type[1]
  skew_boundary <- 0.99527
  stopifnot(
    is.integer(df), df > 0L, inherits(formula, "formula"),
    inherits(Z, "formula"), is.data.frame(data),
    !missing(cluster),
    all(do_setup %in% c(.laplace_char, .gva_char, .snva_char)),
    is.integer(n_nodes), length(n_nodes) == 1L, n_nodes > 0L,
    link %in% c("PH", "PO", "probit"),
    param_type %in% c("DP", "CP_trans", "CP"),
    is.integer(n_threads) && n_threads > 0L && length(n_threads) == 1L,
    is.numeric(skew_start), length(skew_start) == 1L,
    is.logical(dense_hess), length(dense_hess) == 1L, !is.na(dense_hess),
    is.logical(sparse_hess), length(sparse_hess) == 1L, !is.na(sparse_hess))
  eval(bquote(stopifnot(
    .(-skew_boundary) < skew_start && skew_start < .(skew_boundary))))

  #####
  # get the cluster variable
  cluster <- substitute(cluster)
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
  # get design matrices and outcome
  mf_X <- model.frame(formula, data = data)
  mt_X <- terms(mf_X)
  X <- model.matrix(mt_X, mf_X)
  n_x_fix <- NCOL(X)

  y <- model.response(mf_X)
  stopifnot(inherits(y, "Surv"), isTRUE(attr(y, "type") == "right"))
  event <- y[, 2]
  tobs  <- y[, 1]

  formula_Z <- Z
  mf_Z <- model.frame(formula_Z, data = data)
  mt_Z <- terms(mf_Z)
  Z <- model.matrix(mt_Z, mf_Z)
  n_rng <- NCOL(Z)
  stopifnot(NCOL(X) > 0, NCOL(Z) > 0, NROW(X) == NROW(Z),
            length(y) == NROW(X))

  #####
  # add time-varying baseline
  time_var <- formula[[2L]][[2L]]
  formula_b <- eval(
    bquote(~ nsx(log(.(time_var)), df = .(df), intercept = FALSE) - 1))
  mt_b <- terms(model.frame(formula_b, data = data[event > 0, ]))
  X <- cbind(X, model.matrix(mt_b, data))

  #####
  # approximate derivatives w/ finite difference
  is_fix <- 1:n_x_fix
  XD <- X
  XD[,  is_fix] <- 0.
  XD[, -is_fix] <- local({
    dt <- .Machine$double.eps^(1/3)
    dat_m1 <- dat_p1 <- eval(bquote(data.frame(.(time_var))), data)
    x <- dat_m1[[1L]]
    dat_p1[[1L]] <- x * exp( dt)
    dat_m1[[1L]] <- x * exp(-dt)
    (model.matrix(mt_b, dat_p1) - model.matrix(mt_b, dat_m1)) / 2 / dt / x
  })

  #####
  # get starting values w/ coxph and then a lm fit
  need_theta <- is.null(theta)
  need_beta  <- is.null(beta)

  stopifnot(isTRUE(attr(mt_Z, "intercept") == 1L),
            colnames(Z)[1] == "(Intercept)")

  inits <- list()

  if(need_beta){
    inits$link_Hat <- local({
        cox_fit <- coxph(formula, data, model = TRUE)
        S_hat <- Shat(cox_fit)
        if(link == "PH")
          pmax(log(.Machine$double.eps) / 4, log(-log(S_hat)))
        else if(link == "PO")
          log((1 - S_hat) / S_hat)
        else if(link == "probit")
          -qnorm(S_hat)
        else
          stop(sprintf("%s not implemented", sQuote(link)))
      })

    inits$coef <- local({
      keep <- event > 0
      lm.fit(x = X[keep, , drop = FALSE],
             y = inits$link_Hat[keep])$coefficients
    })
  }

  if(need_theta)
    inits$theta <- local({
      cox_frm <- eval(bquote(
        update(formula, . ~ . + frailty(.(cluster),
                                        distribution = "gaussian"))))

      cox_fit <- coxph(cox_frm, data)

      cox_fit$history[[1L]]$theta
    })

  #####
  # setup ADFun object for the Laplace approximation
  data_ad_func <- list(
    tobs = tobs, event = event, X = X, XD = XD, Z = Z, grp = grp - 1L,
    link = link, grp_size = grp_size, n_threads = n_threads)

  # the user may have provided values
  theta <- if(!need_theta){
    stopifnot(all(dim(theta) == n_rng))
    theta

  } else local({
    out <- matrix(0., n_rng, n_rng)
    out[1L] <- inits$theta
    if(n_rng > 1L)
      diag(out)[-1L] <- rep(1, n_rng - 1L)
    out
  })

  theta <- cov_to_theta(theta)

  beta <- if(!is.null(beta)){
    stopifnot(length(beta) == NCOL(X))
    beta
  } else
    inits$coef
  names(beta) <- colnames(X)

  # assign parameter list
  params = list(
    eps = .MGSM_defaul_eps, kappa = .MGSM_default_kappa, b = beta,
    theta = theta)

  laplace_out <- if(.laplace_char %in% do_setup)
    .get_laplace_func(data_ad_func, params, n_rng, n_grp, inits)
  else
    NULL

  #####
  # setup ADFun object for the GVA
  gva_out <- if(.gva_char %in% do_setup || .snva_char %in% do_setup)
    .get_gva_func(n_rng, n_grp, params, data_ad_func, n_nodes, dense_hess,
                  sparse_hess, inits, opt_func)
  else
    NULL

  #####
  # setup ADFun object for the SNVA
  snva_out <- if(.snva_char %in% do_setup)
    .get_snva_out(n_rng, n_grp, gva_out, params, skew_start, param_type,
                  skew_boundary, data_ad_func, n_nodes, dense_hess,
                  sparse_hess, opt_func)
  else
    NULL

  structure(
    list(laplace = laplace_out, gva = gva_out, snva = snva_out, y = y,
         event = event, X = X, XD = XD, Z = Z, grp = grp, terms = list(
           X = mt_X, Z = mt_Z, baseline = mt_b), cl = match.call(),
         link = link, opt_func = opt_func, dense_hess = dense_hess,
         sparse_hess = sparse_hess),
    class = "MGSM_ADFun")
}

.get_laplace_func <- function(data_ad_func, params, n_rng, n_grp, inits) {
  # get Laplace AD function
  adfunc_laplace <- local({
    data_ad_func <- c(
      list(app_type = .laplace_char), data_ad_func)
    # TODO: initialize random effects in a smarter way...
    params$u <- matrix(0., n_rng, n_grp)
    MakeADFun(
      data = data_ad_func, parameters = params, DLL = "survTMB",
      silent = TRUE, random = "u")
  })

  # we make a wrapper object to account for the eps and kappa and allow the
  # user to change these
  with(new.env(), {
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

    c(list(
      par = par,
      fn = function(x, ...){ fn(get_x(x))                               },
      gr = function(x, ...){ gr(get_x(x))[-(1:2)]                       },
      he = function(x, ...){ he(get_x(x))[-(1:2), -(1:2), drop = FALSE] },
      # function to set penalty parameters
      update_pen = function(eps, kappa){
        p_env <- parent.env(environment())
        if(!missing(eps))
          assign("eps"  , eps  , p_env)
        if(!missing(kappa))
          assign("kappa", kappa, p_env)
        invisible(with(p_env, c(eps = eps, kappa = kappa)))
      }
    ), out)
  })
}

# VA util funcs
.get_par_va <- function(params)
  with(params, {
    names(b)        <- rep("b", length(b))
    names(theta)    <- rep("theta", length(theta))
    names(theta_VA) <- rep("theta_VA", length(theta_VA))
    c(eps = eps, kappa = kappa, b, theta, theta_VA)
  })
.eval_hess_sparse <- function(ptr, par){
  out <- VA_funcs_eval_hess_sparse(ptr, par)
  Matrix::sparseMatrix(
    i = out$row_idx + 1L, j = out$col_idx + 1L, x = out$val,
    symmetric = TRUE)
}

.get_gva_func <- function(n_rng, n_grp, params, data_ad_func, n_nodes,
                          dense_hess, sparse_hess, inits, opt_func) {
  # setup cache
  setup_atomic_cache(
    n_nodes = n_nodes, type = .gva_char, link = data_ad_func$link)

  # set the initial values
  n_mu     <- n_rng
  n_Lambda <- (n_rng * (n_rng + 1L)) / 2L
  n_p_grp  <- n_mu + n_Lambda
  theta_VA <- rep(NA_real_, n_p_grp * n_grp)

  # set means to zero
  idx_mean <- sapply(1:n_grp - 1L, function(i) i * n_p_grp + 1:n_mu)
  theta_VA[idx_mean] <- 0.
  # set parameters for Lambda
  idx_Lambda <- sapply(1:n_grp - 1L, function(i)
    i * n_p_grp + n_mu + 1:n_Lambda)
  theta_VA[idx_Lambda] <- params$theta
  # set names
  theta_VA_names <- c(paste0("mu", 1:n_rng), names(params$theta))
  theta_VA_names <- c(outer(
    theta_VA_names, paste0("g", 1:n_grp), function(x, y)
    paste0(y, ":", x)))
  names(theta_VA) <- theta_VA_names

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

      } else
        within(MakeADFun(
          data = data_ad_func, parameters = params, DLL = "survTMB",
          silent = TRUE),
          he_sp <- function(x, ...)
            stop())
    })

    # we make a wrapper object to account for the eps and kappa and allow the
    # user to change these
    gva_out <- with(new.env(), {
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
        fn = function(x, ...){ fn(get_x(x))                               },
        gr = function(x, ...){ gr(get_x(x))[-(1:2)]                       },
        he = function(x, ...){ he(get_x(x))[-(1:2), -(1:2), drop = FALSE] },
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
        control = list(maxit = 1000L)
      ), out)
    })
  }

  # find initial VA params
  func <- get_gva_out(theta_VA)
  coefs_start <- c(params$b, params$theta)

  # find parameters for each group
  par <- c(coefs_start, theta_VA)
  do_drop <- seq_along(coefs_start)
  get_x <- function(x){
    par[-do_drop] <- x
    par
  }

  fn <- function(x)
    func$fn(get_x(x))
  gr <- function(x)
    func$gr(get_x(x))[-do_drop]

  opt_out <- opt_func(theta_VA, fn = fn, gr = gr,
                      control = list(maxit = 1000L,
                                     reltol = .Machine$double.eps^(1/5)))
  get_gva_out(opt_out$par)
}

.get_snva_out <- function(n_rng, n_grp, gva_out, params, skew_start,
                          param_type, skew_boundary, data_ad_func, n_nodes,
                          dense_hess, sparse_hess, opt_func) {
  # setup cache
  setup_atomic_cache(
    n_nodes = n_nodes, type = .snva_char, link = data_ad_func$link)

  # set the initial values
  n_mu     <- n_rng
  n_rho    <- n_rng
  n_Lambda <- (n_rng * (n_rng + 1L)) / 2L
  n_p_grp  <- n_mu + n_Lambda + n_rho
  theta_VA <- rep(NA_real_, n_p_grp * n_grp)

  #####
  # setup starting values for VA parameters
  # find GVA solution
  beta <- params$b
  theta <- params$theta
  gva_opt <- with(gva_out, opt_func(par, fn, gr,
                                    control = list(maxit = 1000L)))
  gva_va_vals <- gva_opt$par[-seq_len(length(beta) + length(theta))]

  params$b     <- gva_opt$par[1:length(beta)]
  params$theta <- gva_opt$par[1:length(theta) + length(beta)]

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
      log((skew_boundary + x) / (skew_boundary - x))

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
  theta_VA_names <- local({
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
  })
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
    else
      within(MakeADFun(
        data = data_ad_func, parameters = params, DLL = "survTMB",
        silent = TRUE),
        he_sp <- function(x, ...)
          stop())
  })

  # we make a wrapper object to account for the eps and kappa and allow the
  # user to change these
  get_snva_out <- function(theta_VA){
    with(new.env(), {
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
        fn = function(x, ...){ fn(get_x(x))                               },
        gr = function(x, ...){ gr(get_x(x))[-(1:2)]                       },
        he = function(x, ...){ he(get_x(x))[-(1:2), -(1:2), drop = FALSE] },
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
        control = list(maxit = 1000L)
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
