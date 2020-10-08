#define INCLUDE_RCPP
#include "get-x.h"
#include "utils.h"
#include "joint-utils.h"
#include "snva-utils.h"
#include "convert-eigen-arma.h"
#include <limits>
#include "bases-wrapper.h"

namespace {
constexpr size_t const INT_NS = 0L,
                     INT_POLY = 1L;

/**
 abstract base class to use for the integration in the hazard and to
 evaluate the basis functions.
 */
template<class Type>
class splines_n_cum_haz_base {
public:
  /** returns true if there is a term for each basis. */
  virtual bool has_b() const = 0;
  virtual bool has_g() const = 0;
  virtual bool has_m() const = 0;
  virtual bool has_g_surv() const = 0;
  virtual bool has_m_surv() const = 0;

  /** evalutes the bases and return by reference. */
  virtual void eval_b(arma::vec&, double const) const = 0;
  virtual void eval_g(arma::vec&, double const) const = 0;
  virtual void eval_m(arma::vec&, double const) const = 0;
  virtual void eval_g_surv(arma::vec&, double const) const = 0;
  virtual void eval_m_surv(arma::vec&, double const) const = 0;

  /** returns the input vector to pass along to eval_int. */
  virtual vector<Type> get_x(
      Type const, Type const, vector<Type> const&, vector<Type> const&,
      matrix<Type> const&, vector<Type> const&, vector<Type> const&,
      matrix<Type> const&) const = 0;

  /** evalutes the integral in the survival function.  */
  virtual Type eval_int(vector<Type> const&) = 0;

  virtual ~splines_n_cum_haz_base() = default;
};

/** class for a general basis. */
template<class Type, class Basis>
class splines_n_cum_haz_T final : public splines_n_cum_haz_base<Type> {
  using cum_integral =
    fastgl::joint::snva_integral
  <typename Type::value_type, Basis, Basis, Basis>;

  bool const use_log;
  std::unique_ptr<Basis> b, g, m, g_surv, m_surv;
  cum_integral cum_haz;

  bool has_b() const {
    return static_cast<bool>(b);
  }
  bool has_g() const {
    return static_cast<bool>(g);
  }
  bool has_m() const {
    return static_cast<bool>(m);
  }
  bool has_g_surv() const {
    return static_cast<bool>(g_surv);
  }
  bool has_m_surv() const {
    return static_cast<bool>(m_surv);
  }

  arma::vec to_arma_vec(vector<Type> const &x){
    arma::vec out(x.size());
    for(int i = 0; i < x.size(); ++i)
      out[i] = asDouble(x[i]);
    return out;
  }

public:
  splines_n_cum_haz_T(
    size_t const n_nodes, vector<Type> const &scoefs,
    vector<Type> const &gcoefs, vector<Type> const &mcoefs,
    size_t const n_y, bool const use_log, vector<Type> const& mcoefs_surv,
    vector<Type> const& gcoefs_surv):
  use_log(use_log),
  b(get_basis<Basis>(to_arma_vec(scoefs))),
  g(get_basis<Basis>(to_arma_vec(gcoefs))),
  m(get_basis<Basis>(to_arma_vec(mcoefs))),
  g_surv(get_basis<Basis>(to_arma_vec(gcoefs_surv))),
  m_surv(get_basis<Basis>(to_arma_vec(mcoefs_surv))),
  cum_haz("cum haz integral", n_nodes, b.get(), g_surv.get(), m_surv.get(),
          n_y, use_log)
  {
    if((!has_g() and has_g_surv()) or (
        has_g_surv() and gcoefs.size() != gcoefs_surv.size()))
      throw std::invalid_argument(
          "splines_n_cum_haz_T(): invalid gcoefs_surv");
    else if(
        (!has_m() and has_m_surv()) or (
          has_m_surv() and mcoefs.size() != mcoefs_surv.size()))
      throw std::invalid_argument(
          "splines_n_cum_haz_T(): invalid mcoefs_surv");
  }

  void eval_b(arma::vec &out, double x) const {
    if(use_log)
      x = log(x);

    b->operator()(out, x);
  }
  void eval_g(arma::vec &out, double const x) const {
    g->operator()(out, x);
  }
  void eval_m(arma::vec &out, double const x) const {
    m->operator()(out, x);
  }
  void eval_g_surv(arma::vec &out, double const x) const {
    g_surv->operator()(out, x);
  }
  void eval_m_surv(arma::vec &out, double const x) const {
    m_surv->operator()(out, x);
  }

  vector<Type> get_x(
      Type const lb, Type const ub, vector<Type> const &omega,
      vector<Type> const &alpha, matrix<Type> const &B,
      vector<Type> const &U, vector<Type> const &k,
      matrix<Type> const &lambda) const {
    return cum_haz.get_x(lb, ub, omega, alpha, B, U, k, lambda);
  }

  Type eval_int(vector<Type> const &params) {
    vector<Type> y(1L);
    cum_haz(params, y);
    return y[0];
  }
};

template<class Type>
size_t get_basis_dim(vector<Type> const &coefs, size_t const basis_type){
  if(basis_type == INT_NS)
    return coefs.size();
  else if(basis_type == INT_POLY)
    return coefs.size() / 2L;

  throw std::runtime_error("'basis_type' not impemented");
  return 0L;
}

/** class for a natural cubic spline. */
using splines::ns;
template<class Type>
using splines_n_cum_haz_ns = splines_n_cum_haz_T<Type, ns>;

/** class for orthogonal polynomials. */
using pol = poly::orth_poly;
template<class Type>
using splines_n_cum_haz_pol = splines_n_cum_haz_T<Type, pol>;

template<class Type>
class VA_worker {
  Rcpp::List data, parameters;

  const DATA_MATRIX(markers);
  const DATA_IVECTOR(n_markers);
  const DATA_VECTOR(m_time);
  const DATA_VECTOR(mcoefs);
  const DATA_VECTOR(gcoefs);
  const DATA_VECTOR(mcoefs_surv);
  const DATA_VECTOR(gcoefs_surv);
  const DATA_MATRIX(X);

  const DATA_VECTOR(tstart);
  const DATA_VECTOR(tstop);
  const DATA_VECTOR(outcomes);
  const DATA_VECTOR(scoefs);
  const DATA_MATRIX(Z);

  const DATA_INTEGER(basis_type);
  const DATA_INTEGER(n_threads);
  const DATA_LOGICAL(sparse_hess);
  const DATA_LOGICAL(use_log);
  const DATA_INTEGER(n_nodes);

  const PARAMETER_MATRIX(gamma);
  const PARAMETER_MATRIX(B);
  const PARAMETER_VECTOR(Psi);
  const PARAMETER_VECTOR(Sigma);
  const PARAMETER_VECTOR(delta);
  const PARAMETER_VECTOR(omega);
  const PARAMETER_VECTOR(alpha);
  const PARAMETER_VECTOR(va_par);

public:
#ifdef _OPENMP
  std::size_t const n_blocks = n_threads;
#else
  std::size_t const n_blocks = 1L;
#endif

public:
  std::size_t const n_y = markers.rows();

  using cum_base_T = splines_n_cum_haz_base<Type>;
  std::vector<std::unique_ptr<cum_base_T> >
    get_splines_n_cum_ints()
    const {
    using output_T = std::unique_ptr<cum_base_T>;

    std::vector<output_T> out;
    out.reserve(n_blocks);
    for(size_t i = 0; i < n_blocks; ++i)
      if     (basis_type == INT_NS)
        out.emplace_back(new splines_n_cum_haz_ns<Type>(
            n_nodes, scoefs, gcoefs, mcoefs, n_y, use_log, mcoefs_surv,
            gcoefs_surv));
      else if(basis_type == INT_POLY)
        out.emplace_back(new splines_n_cum_haz_pol<Type>(
            n_nodes, scoefs, gcoefs, mcoefs, n_y, use_log, mcoefs_surv,
            gcoefs_surv));
      else
        throw std::invalid_argument("'basis_type' not implemented");

    return out;
  }

private:
  cum_base_T &get_splines_n_cum_haz(
      std::vector<std::unique_ptr<cum_base_T> > &splines_n_cum_ints)
  const {
#ifdef _OPENMP
    return *splines_n_cum_ints[omp_get_thread_num()];
#else
    return *splines_n_cum_ints[0L];
#endif
  }

public:
  std::size_t const dim_m = get_basis_dim(mcoefs, basis_type),
                    dim_g = get_basis_dim(gcoefs, basis_type),
                    n_fix = gamma.rows(),
                    dim_b = get_basis_dim(scoefs, basis_type),
                        K = n_y * dim_m,
                      n_Z = Z.rows(),
             n_obs_marker = markers.cols(),
                 n_groups = outcomes.size(),
                   n_pars =
                     gamma.rows() * gamma.cols() + B.rows() * B.cols() +
                     Psi.size() + Sigma.size() + delta.size() +
                     omega.size() + alpha.size() + va_par.size();

  VA_worker(Rcpp::List data, Rcpp::List parameters):
    data(data), parameters(parameters) {
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
    data = Rcpp::List();
    parameters = Rcpp::List();

    /* checks */
    if((size_t)n_markers.size() != n_groups)
      throw std::invalid_argument("VA_worker: invalid n_markers");
    else if((size_t)tstart.size() != n_groups)
      throw std::invalid_argument("VA_worker: invalid tstart");
    else if((size_t)tstop.size() != n_groups)
      throw std::invalid_argument("VA_worker: invalid tstop");
    else if((size_t)outcomes.size() != n_groups)
      throw std::invalid_argument("VA_worker: invalid outcomes");
    else if((size_t)Z.cols() != n_groups)
      throw std::invalid_argument("VA_worker: invalid Z");

    {
      size_t sum_n_markers = 0L;
      for(unsigned i = 0; i < n_markers.size(); ++i)
        sum_n_markers += n_markers[i];
      if(sum_n_markers != (size_t)markers.cols())
        throw std::invalid_argument("VA_worker: invalid number of markers");
    }

    if((size_t)m_time.size() != n_obs_marker)
      throw std::invalid_argument("VA_worker: invalid m_time");

    else if((size_t)B.cols() != n_y or (size_t)B.rows() != dim_g)
      throw std::invalid_argument("VA_worker: invalid B");
    else if(gamma.rows() != X.rows() or (size_t)gamma.cols() != n_y)
      throw std::invalid_argument("VA_worker: invalid gamma");
    else if((size_t)X.cols() != n_obs_marker)
      throw std::invalid_argument("VA_worker: invalid X");
    else if(survTMB::get_rng_dim(Psi) != K)
      throw std::invalid_argument("VA_worker: invalid Psi");
    else if(survTMB::get_rng_dim(Sigma) != n_y)
      throw std::invalid_argument("VA_worker: invalid Sigma");
    else if(delta.size() != Z.rows())
      throw std::invalid_argument("VA_worker: invalid delta");
    else if((size_t)omega.size() != dim_b)
      throw std::invalid_argument("VA_worker: invalid omega");
    else if((size_t)alpha.size() != n_y)
      throw std::invalid_argument("VA_worker: invalid alpha");
  }

  template<class Tout>
  vector<Tout> get_concatenated_args() const {
    vector<Tout> out(n_pars);

    Type *o = &out[0];
    auto add_mat = [&](matrix<Type> const &X){
      for(size_t j = 0; j < (size_t)X.cols(); ++j)
        for(size_t i = 0; i < (size_t)X.rows(); ++i)
          *o++ = Tout(X(i, j));
    };
    auto add_vec = [&](vector<Type> const &x){
      for(size_t i = 0; i < (size_t)x.size(); ++i)
        *o++ = Tout(x[i]);
    };

    add_mat(gamma);
    add_mat(B);

    add_vec(Psi);
    add_vec(Sigma);
    add_vec(delta);
    add_vec(omega);
    add_vec(alpha);
    add_vec(va_par);

    return out;
  }

  Type operator()
    (vector<Type> &args,
     std::vector<std::unique_ptr<cum_base_T> > &splines_n_cum_ints)
    const {
    if((size_t)args.size() != n_pars)
      throw std::invalid_argument("VA_worker::operator(): invalid args");

    /* assign the parameters from args */
    Type *a = &args[0];

    auto set_mat = [&](size_t const nr, size_t const nc){
      matrix<Type> out(nr, nc);
      for(size_t j = 0; j < nc; ++j)
        for(size_t i = 0; i < nr; ++i)
          out(i, j) = *a++;
      return out;
    };
    auto set_vec = [&](size_t const n){
      vector<Type> out(n);
      for(size_t i = 0; i < n; ++i)
        out[i] = *a++;
      return out;
    };
    auto set_vcov = [&](size_t const n){
      return survTMB::get_vcov_from_trian(set_vec((n * (n + 1L)) / 2L));
    };

    matrix<Type> agamma = set_mat(gamma.rows(), gamma.cols()),
                     aB = set_mat(B.rows(), B.cols()),
                   aPsi = set_vcov(K),
                 aSigma = set_vcov(n_y);
    vector<Type> adelta = set_vec(delta.size()),
                 aomega = set_vec(omega.size()),
                 aalpha = set_vec(alpha.size());

    auto const ava_par = ([&](){
      return GaussHermite::SNVA::SNVA_MD_theta_DP_to_DP(
        a, va_par.size(), K);
    })();

    /* get the cumulative hazard integral object */
    auto &my_int = get_splines_n_cum_haz(splines_n_cum_ints);

    /* vector needed to be passed to the integral object */
    vector<Type> cum_int_arg = ([&](){
      vector<Type> U(K), k(K);
      matrix<Type> Lambda(K, K);

      return my_int.get_x(Type(1.) /* lb */, Type(2.) /* ub */, aomega,
                          aalpha, aB, U, k, Lambda);
    })();

    /* evaluates the cumulative hazard integral */
    bool const is_in_parallel = CppAD::thread_alloc::in_parallel(),
                   has_b      = my_int.has_b(),
                   has_g      = my_int.has_g(),
                   has_m      = my_int.has_m(),
                   has_g_surv = my_int.has_g_surv(),
                   has_m_surv = my_int.has_m_surv();
    auto comp_cum_haz = [&](
      Type const lb, Type const ub, vector<Type> const &U,
      vector<Type> const &k, matrix<Type> const &Lambda){
      cum_int_arg[0L] = lb;
      cum_int_arg[1L] = ub;

      Type *x = &cum_int_arg[
        2L + has_b * dim_b + n_y + has_g_surv * dim_g * n_y];
      if(has_m_surv){
        for(size_t i = 0; i < K; ++i)
          *x++ = U[i];
        for(size_t i = 0; i < K; ++i)
          *x++ = k[i];
        for(size_t j = 0; j < K; ++j)
          for(size_t i = 0; i < K; ++i)
            *x++ = Lambda(i, j);
      }

      return my_int.eval_int(cum_int_arg);
    };

    /* compute other intermediaries */
    Type log_det_sigma;
    matrix<Type> sigma_inv = atomic::matinvpd(aSigma, log_det_sigma);
    Type log_det_psi;
    matrix<Type> psi_inv = atomic::matinvpd(aPsi, log_det_psi);

    Type const one(1.),
              zero(0.),
              half(.5),
            two_pi(2. / M_PI),
       sqrt_two_pi(sqrt(two_pi)),
             small(std::numeric_limits<double>::epsilon());

    /* evaluate the lower bound */
    survTMB::accumulator_mock<Type> result;
    arma::vec b_wrk(has_b ? dim_b : 0L),
              g_wrk(has_g ? dim_g : 0L),
              m_wrk(has_m ? dim_m : 0L);
    size_t marker_idx(0L);
    for(size_t g = 0; g < n_groups; ++g){
      /* is this our cluster? */
      if(is_in_parallel and !is_my_region(*result.obj)){
        result.obj->parallel_region();
        marker_idx += n_markers[g];
        continue;
      }

      /* get VA parameters */
      vector<Type> const  &va_mu = ava_par.va_mus[g],
                         &va_rho = ava_par.va_rhos[g];
      matrix<Type> const &Lambda = ava_par.va_lambdas[g];

      /* the term we add in the end */
      Type term(0.);

      /* compute quantities needed later */
      vector<Type> k = ([&](){
        vector<Type> out = Lambda * va_rho;
        Type denom = sqrt(one + vec_dot(va_rho, out));
        for(int j = 0; j < out.size(); ++j)
          out[j] /= denom;
        return out;
      })();
      vector<Type> const U_mean = va_mu + sqrt_two_pi * k;

      /* evaluate fixed time-invariant effect */
      vector<Type> fix_invariant(n_y);
      for(size_t j = 0; j < n_y; ++j){
        fix_invariant[j] = zero;
        for(size_t i = 0; i < n_fix; i++)
          fix_invariant[j] += X(i, marker_idx) * agamma(i, j);
      }

      /* terms from the survival outcome */
      {
        Type surv_term(0.);
        Type const lb = tstart[g],
                   ub = tstop [g];
        Type z_dot_d(0.);
        for(int j = 0; j < Z.rows(); ++j)
          z_dot_d += Z(j, g) * adelta[j];

        Type const alpha_fix_invariant = vec_dot(aalpha, fix_invariant);

        if(asDouble(outcomes[g]) > 0.){
          /* add the jump term */
          surv_term += z_dot_d + alpha_fix_invariant;

          double const double_ub = asDouble(tstop [g]);
          if(has_b){
            my_int.eval_b(b_wrk, double_ub);

            for(size_t j = 0; j < b_wrk.n_elem; ++j)
              surv_term += Type(b_wrk[j]) * aomega[j];
          }

          if(has_g_surv){
            my_int.eval_g_surv(g_wrk, double_ub);
            vector<Type> const type_g_wrk = vec_eigen_arma<Type>(g_wrk);
            for(size_t j = 0; j < n_y; ++j)
              for(size_t i = 0; i < dim_g; ++i)
                surv_term += aalpha[j] * type_g_wrk[i] * aB(i, j);
          }

          if(has_m_surv){
            my_int.eval_m_surv(m_wrk, double_ub);
            vector<Type> const type_m_wrk = vec_eigen_arma<Type>(m_wrk);

            Type const *u_ptr = &U_mean[0L];
            for(size_t j = 0; j < n_y; ++j)
              for(size_t i = 0; i < dim_m; ++i, ++u_ptr)
                surv_term += aalpha[j] * type_m_wrk[i] * *u_ptr;
          }
        }

        /* add term from survival probability */
        Type const f1 = exp(z_dot_d + alpha_fix_invariant),
                   f2 = comp_cum_haz(lb, ub, va_mu, k, Lambda);
        surv_term -= f1 * f2;
        term += surv_term;
      }

      /* terms from marker */
      {
        Type marker_term(0.);
        size_t const marker_end = marker_idx + n_markers[g];
        for(; marker_idx < marker_end; ++marker_idx){
          vector<Type> residual(n_y);
          for(size_t i = 0; i < n_y; ++i)
            residual[i] = markers(i, marker_idx) - fix_invariant[i];

          double const obs_time_i = asDouble(m_time[marker_idx]);
          if(has_g){
            my_int.eval_g(g_wrk, obs_time_i);
            vector<Type> const type_g_wrk = vec_eigen_arma<Type>(g_wrk);

            for(size_t j = 0; j < n_y; ++j)
              for(size_t i = 0; i < dim_g; ++i)
                residual[j] -= aB(i, j) * type_g_wrk[i];
          }

          if(has_m){
            my_int.eval_m(m_wrk, obs_time_i);
            vector<Type> const type_m_wrk = vec_eigen_arma<Type>(m_wrk);

            {
              Type const *mu_i = &U_mean[0L];
              for(size_t j = 0; j < n_y; ++j)
                for(size_t i = 0; i < dim_m; ++i, ++mu_i)
                  residual[j] -= *mu_i * type_m_wrk[i];
            }

            /* TODO: Do something smarter... */
            matrix<Type> M(n_y, K);
            M.setZero();
            for(size_t i = 0; i < n_y; ++i){
              Type const *m_i = &type_m_wrk[0L];
              size_t const j_end =  (i + 1L) * dim_m;
              for(size_t j = i * dim_m; j < j_end; ++j, ++m_i)
                M(i, j) = *m_i;
            }

            matrix<Type> tmp_mat = M *  Lambda * M.transpose();
            vector<Type> tmp_vec = M * k;
            marker_term -=
              half * (mat_mult_trace(tmp_mat, sigma_inv)
                        - two_pi * quad_form_sym(tmp_vec, sigma_inv));
          }

          marker_term -= half * quad_form_sym(residual, sigma_inv);
        }

        term += marker_term;

        /* terms from the random effect prior */
        /* TODO: can be done smarter */
        Type prior_term = - half * (
          quad_form_sym(U_mean, psi_inv)
          + mat_mult_trace(Lambda, psi_inv)
          - two_pi * quad_form_sym(k, psi_inv));

        term += prior_term;
      }

      /* add entropy terms*/
      {
        Type misc_term = half * atomic::logdet(Lambda);
        auto const llt_mat = Lambda.llt();
        vector<Type> const va_rho_scaled =
          (llt_mat.matrixU() * va_rho.matrix()).array() + small;
        misc_term -= GaussHermite::SNVA::entropy_term(
          vec_dot(va_rho_scaled, va_rho_scaled), n_nodes);
        term += misc_term;
      }

      result += term;
    }

    if(!is_my_region(*result.obj))
      /* only have to add one more term so just return */
      return result;

    Type norm_constant =
      - Type(double(n_obs_marker * n_y) / 2. * log(2 * M_PI))
      - Type(double(n_groups) / 2.) * log_det_psi
      - Type(double(n_obs_marker) / 2.) * log_det_sigma
      + Type(double(K * n_groups) / 2.)
      - Type(double(n_groups) * M_LN2);

    result += norm_constant;

    return result;
  }
};

class VA_func {
  using ADd   = CppAD::AD<double>;
  using ADdd  = CppAD::AD<ADd>;
  using ADddd = CppAD::AD<ADdd>;
  template<class Type>
  using ADFun = CppAD::ADFun<Type>;

  size_t n_pars;

public:
  size_t get_n_pars() const {
    return n_pars;
  }

  std::vector<std::unique_ptr<splines_n_cum_haz_base<ADd> > >
    splines_n_cum_ints_ADd;
  std::vector<std::unique_ptr<ADFun<double> > > funcs;

  VA_func(Rcpp::List data, Rcpp::List parameters){
    {
      /* to compute function and gradient */
      VA_worker<ADd> w(data, parameters);
      splines_n_cum_ints_ADd = w.get_splines_n_cum_ints();
      funcs.resize(w.n_blocks);
      n_pars = w.n_pars;
      vector<ADd> args = w.get_concatenated_args<ADd>();

#ifdef _OPENMP
#pragma omp parallel for if(w.n_blocks > 1L) firstprivate(args)
#endif
      for(unsigned i = 0; i < w.n_blocks; ++i){
        funcs[i].reset(new ADFun<double>());

        CppAD::Independent(args);
        vector<ADd> y(1);
        y[0] = w(args, splines_n_cum_ints_ADd);

        funcs[i]->Dependent(args, y);
        funcs[i]->optimize();
      }
    }
  }
};
} // namesapce

// [[Rcpp::export(rng = false)]]
SEXP get_joint_funcs
  (Rcpp::List data, Rcpp::List parameters){
  shut_up();

  unsigned const n_threads(data["n_threads"]);
  setup_parallel_ad setup_ADd(n_threads);

  return Rcpp::XPtr<VA_func>(new VA_func(data, parameters));
}

// [[Rcpp::export(rng = false)]]
double joint_funcs_eval_lb(SEXP p, SEXP par){
  shut_up();

  Rcpp::XPtr<VA_func > ptr(p);
  std::vector<std::unique_ptr<CppAD::ADFun<double> > > &funcs = ptr->funcs;
  vector<double> parv = get_vec<double>(par);
  if((size_t)parv.size() != ptr->get_n_pars())
    throw std::invalid_argument("joint_funcs_eval_lb: invalid par");

  unsigned const n_blocks = ptr->funcs.size();
  double out(0);
#ifdef _OPENMP
/* #pragma omp parallel for if(n_blocks > 1L) firstprivate(parv) reduction(+:out) */
#pragma omp parallel for ordered if(n_blocks > 1L) firstprivate(parv)  schedule(static, 1)
#endif
  for(unsigned i = 0; i < n_blocks; ++i){
    double const term = funcs[i]->Forward(0, parv)[0];
#ifdef _OPENMP
    /* TODO: replace with a reduction */
#pragma omp ordered
#endif
    out += term;
  }

  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector joint_funcs_eval_grad(SEXP p, SEXP par){
  shut_up();

  Rcpp::XPtr<VA_func > ptr(p);
  std::vector<std::unique_ptr<CppAD::ADFun<double> > > &funcs = ptr->funcs;
  vector<double> parv = get_vec<double>(par);
  if((size_t)parv.size() != ptr->get_n_pars())
    throw std::invalid_argument("joint_funcs_eval_grad: invalid par");

  unsigned const n_blocks = ptr->funcs.size();
  vector<double> grad(parv.size());
  grad.setZero();

#ifdef _OPENMP
#pragma omp parallel for ordered if(n_blocks > 1L) firstprivate(parv) schedule(static, 1)
#endif
  for(unsigned i = 0; i < n_blocks; ++i){
    funcs[i]->Forward(0, parv);
    vector<double> w(1);
    w[0] = 1;

    vector<double> grad_i = funcs[i]->Reverse(1, w);
#ifdef _OPENMP
    /* TODO: replace with a reduction */
#pragma omp ordered
#endif
    grad += grad_i;
  }

  std::size_t const n = grad.size();
  Rcpp::NumericVector out(n);
  for(unsigned i = 0; i < n; ++i)
    out[i] = grad[i];

  return out;
}
