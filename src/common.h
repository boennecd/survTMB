/* Common arguments used in approximations.
 *
 * Args:
 *   tobs: observed times.
 *   event: zero/one variables for whether the event is observed.
 *   X: fixed effect design matrix.
 *   XD: derivative of X wrt time.
 *   Z: random effect design matrix.
 *   grp: integer vector with group identifier.
 *   eps: small tolerance for positivity constraint on hazards.
 *   kappa: one-sided L2 penalty for positivity constraint.
 *   b: fixed effect coefficients.
 *   theta: covariance matrix parameters. See get_vcov_from_trian.
 */
#ifndef COMMON_ARGS
#define COMMON_ARGS(TYPE)                                      \
  parallel_accumulator<TYPE> &result,                          \
  vector<TYPE> const &tobs, vector<TYPE> const &event,         \
  matrix<TYPE> const &X, matrix<TYPE> const &XD,               \
 matrix<TYPE> const &Z, vector<int> const &grp,                \
 TYPE const &eps, TYPE const &kappa, vector<TYPE> const &b,    \
  vector<TYPE> const &theta, std::string const &link,          \
  vector<int> &grp_size
#endif
