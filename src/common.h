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
#define COMMON_ARGS(TYPE, ACCUMLATOR)                          \
  ACCUMLATOR<TYPE> &result,                                    \
  vector<TYPE> const &tobs, vector<TYPE> const &event,         \
  matrix<TYPE> const &X, matrix<TYPE> const &XD,               \
 matrix<TYPE> const &Z, vector<int> const &grp,                \
 TYPE const &eps, TYPE const &kappa, vector<TYPE> const &b,    \
  vector<TYPE> const &theta, std::string const &link,          \
  vector<int> const &grp_size
#endif

#ifndef COMMON_CALL
#define COMMON_CALL                                              \
  result, tobs, event, X, XD, Z, grp, eps, kappa, b, theta,      \
  link, grp_size
#endif

#ifndef SETUP_DATA
#define SETUP_DATA                                             \
  /* common data objects and parameters */                     \
  DATA_STRING(app_type);                                       \
  DATA_VECTOR(tobs);                                           \
                                                               \
  DATA_VECTOR(event);                                          \
  DATA_MATRIX(X);                                              \
  DATA_MATRIX(XD);                                             \
  DATA_MATRIX(Z);                                              \
  DATA_IVECTOR(grp);                                           \
  DATA_STRING(link);                                           \
  DATA_IVECTOR(grp_size);                                      \
  DATA_INTEGER(n_threads);                                     \
                                                               \
  /* They are marked as parameter such that the user can       \
   * change them later */                                      \
  PARAMETER(eps);                                              \
  PARAMETER(kappa);                                            \
                                                               \
  PARAMETER_VECTOR(b);                                         \
  PARAMETER_VECTOR(theta);
#endif

#ifndef SETUP_DATA_CHECK
#define SETUP_DATA_CHECK                                       \
  /* checks */                                                 \
  {                                                            \
    unsigned const n = tobs.size();                            \
    auto check_rows = [n](matrix<Type> const &x,               \
                          char const *msg){                    \
      if(x.rows() != n)                                        \
        error(msg);                                            \
    };                                                         \
    check_rows(X , "invalid 'X'");                             \
    check_rows(XD, "invalid 'XD'");                            \
    check_rows(Z , "invalid 'Z'");                             \
    if(n != event.size())                                      \
      error("invalid 'event'");                                \
    if(n != grp.size())                                        \
      error("invalid 'grp'");                                  \
                                                               \
    if(b.size()  != X.cols())                                  \
      error("invalid 'b'");                                    \
    if(XD.cols() != X.cols())                                  \
      error("invalid 'XD' (# columns)");                       \
                                                               \
    if(Z.cols() != survTMB::get_rng_dim(theta))                \
      error("invalid 'Z' (# columns: %d %d)", Z.cols(),        \
            survTMB::get_rng_dim(theta));                      \
                                                               \
    std::size_t grp_size_sum = 0;                              \
    for(int i = 0; i < grp_size.size(); ++i)                   \
      grp_size_sum += grp_size[i];                             \
    if(n != grp_size_sum)                                      \
      error("invalid 'grp_size'");                             \
                                                               \
    int g_max = 0;                                             \
    for(int i = 0; i < grp.size(); ++i)                        \
      g_max = std::max(grp[i], g_max);                         \
    if(g_max + 1L != grp_size.size())                          \
      error(                                                   \
        "invalid 'grp_size' (does not match with number of groups)"); \
                                                               \
    for(int i = 1L; i < grp.size(); ++i){                      \
      if(grp[i - 1L] > grp[i])                                 \
        error("'grp' is not sorted");                          \
      if(grp[i - 1L] - grp[i] < -1L)                           \
        error("too big gap in 'grp'");                         \
    }                                                          \
  }
#endif

