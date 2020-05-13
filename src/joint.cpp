#define INCLUDE_RCPP
#include "get-x.h"
#include "utils.h"
#include "joint-utils.h"
#include "snva-utils.h"
#include "splines.h"
#include "convert-eigen-arma.h"
#include <limits>

namespace {

template<class Type>
class VA_worker {
  Rcpp::List data, parameters;

  const DATA_MATRIX(markers);
  const DATA_IVECTOR(n_markers);
  const DATA_VECTOR(m_time);
  const DATA_VECTOR(mknots);

  const DATA_VECTOR(tstart);
  const DATA_VECTOR(tstop);
  const DATA_VECTOR(outcomes);
  const DATA_VECTOR(sknots);
  const DATA_MATRIX(Z);

  const DATA_INTEGER(n_threads);
  const DATA_LOGICAL(sparse_hess);
  const DATA_INTEGER(n_nodes);

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

private:
  std::size_t const n_y = markers.rows(),
                  dim_m = mknots.size(),
                  dim_b = sknots.size(),
                      K = n_y * dim_m,
                    n_Z = Z.rows(),
               n_groups = outcomes.size(),
                 n_pars = B.rows() * B.cols() + Psi.size() + Sigma.size() +
                   delta.size() + omega.size() + alpha.size() +
                   va_par.size();

public:
  using cum_integral =
    fastgl::joint::snva_integral
    <typename Type::value_type, splines::ns, splines::ns>;

  /* small class to hold spline bases and object to compute the cumulative
   * hazard */
  class splines_n_cum_haz {
    splines::ns get_basis(vector<Type> const &knots) const {
      arma::vec bk(2L);
      bk[0L] = asDouble(knots[0L]);
      bk[1L] = asDouble(knots[knots.size() - 1L]);

      arma::vec ik(knots.size() - 2L);
      for(int i = 1L; i < knots.size() - 1L; ++i)
        ik[i - 1L] = asDouble(knots[i]);

      return splines::ns(bk, ik, true);
    }

  public:
    splines::ns b, m;
    cum_integral cum_haz;

    splines_n_cum_haz(
      size_t const n_nodes, vector<Type> const &sknots,
      vector<Type> const &mknots, size_t const n_y):
      b(get_basis(sknots)),
      m(get_basis(mknots)),
      cum_haz("cum haz integral", n_nodes, b, m, n_y) { }
  };

  std::vector<splines_n_cum_haz> get_splines_n_cum_ints() const {
    std::vector<splines_n_cum_haz> out;
    out.reserve(n_blocks);
    for(size_t i = 0; i < n_blocks; ++i)
      out.emplace_back(n_nodes, sknots, mknots, n_y);

    return out;
  }

private:
  splines_n_cum_haz &get_splines_n_cum_haz(
      std::vector<splines_n_cum_haz> &splines_n_cum_ints) const {
#ifdef _OPENMP
    return splines_n_cum_ints[omp_get_thread_num()];
#else
    return splines_n_cum_ints[0L];
#endif
  }

public:
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

    if(m_time.size() != markers.cols())
      throw std::invalid_argument("VA_worker: invalid m_time");

    else if(mknots.size() < 2L)
      throw std::invalid_argument("VA_worker: invalid mknots");
    else if(sknots.size() < 2L)
      throw std::invalid_argument("VA_worker: invalid sknots");

    else if((size_t)B.cols() != n_y or (size_t)B.rows() != dim_m)
      throw std::invalid_argument("VA_worker: invalid B");
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
    for(size_t j = 0; j < (size_t)B.cols(); ++j)
      for(size_t i = 0; i < (size_t)B.rows(); ++i)
        *o++ = Tout(B(i, j));

    auto add_vec = [&](vector<Type> const &x){
      for(size_t i = 0; i < (size_t)x.size(); ++i)
        *o++ = Tout(x[i]);
    };
    add_vec(Psi);
    add_vec(Sigma);
    add_vec(delta);
    add_vec(omega);
    add_vec(alpha);
    add_vec(va_par);

    return out;
  }

  Type operator()
    (vector<Type> &args, std::vector<splines_n_cum_haz> &splines_n_cum_ints)
    const {
    if((size_t)args.size() != n_pars)
      throw std::invalid_argument("VA_worker::operator(): invalid args");

    /* assign the parameters from args */
    Type *a = &args[0];

    matrix<Type> aB(B.rows(), B.cols());
    for(size_t j = 0; j < (size_t)B.cols(); ++j)
      for(size_t i = 0; i < (size_t)B.rows(); ++i)
        aB(i, j) = *a++;

    auto set_vec = [&](size_t const n){
      vector<Type> out(n);
      for(size_t i = 0; i < n; ++i)
        out[i] = *a++;
      return out;
    };
    auto set_vcov = [&](size_t const n){
      return survTMB::get_vcov_from_trian(set_vec((n * (n + 1L)) / 2L));
    };
    matrix<Type> aPsi = set_vcov(K),
               aSigma = set_vcov(n_y);
    vector<Type> adelta = set_vec(delta.size()),
                 aomega = set_vec(omega.size()),
                 aalpha = set_vec(alpha.size());

    /* TODO: delete */
    /*auto print_mat = [&](std::string const &name, matrix<Type> const &x){
      Rcpp::Rcout << name << ":\n";
      for(int i = 0; i < x.rows(); ++i){
        for(int j = 0; j < x.cols(); ++j)
          Rcpp::Rcout << asDouble(x(i, j)) << '\t';
        Rcpp::Rcout << '\n';
      }
    };
    auto print_vec = [&](std::string const &name, vector<Type> const &x){
      Rcpp::Rcout << name << ": ";
      for(int i = 0; i < x.size(); ++i)
        Rcpp::Rcout << asDouble(x[i]) << '\t';
      Rcpp::Rcout << '\n';
    };*/

    auto const ava_par = ([&](){
      vector<Type> theta_va = set_vec(va_par.size());
      return  GaussHermite::SNVA::SNVA_MD_theta_DP_to_DP(theta_va, K);
    })();

    /* get the cumulative hazard integral object */
    auto &my_int = get_splines_n_cum_haz(splines_n_cum_ints);

    /* vector needed to be passed to the integral object */
    vector<Type> cum_int_arg = ([&](){
      vector<Type> U(K), k(K);
      matrix<Type> Lambda(K, K);

      return my_int.cum_haz.get_x(
        Type(1.) /* lb */, Type(2.) /* ub */, aomega, aalpha, U, k, Lambda);
    })();

    /* evaluates the cumulative hazard integral */
    auto comp_cum_haz = [&](
      Type const lb, Type const ub, vector<Type> const &U,
      vector<Type> const &k, matrix<Type> const &Lambda){
      cum_int_arg[0L] = lb;
      cum_int_arg[1L] = ub;

      Type *x = &cum_int_arg[2L + omega.size() + alpha.size()];
      for(size_t i = 0; i < K; ++i)
        *x++ = U[i];
      for(size_t i = 0; i < K; ++i)
        *x++ = k[i];
      for(size_t j = 0; j < K; ++j)
        for(size_t i = 0; i < K; ++i)
          *x++ = Lambda(i, j);

      vector<Type> y(1L);
      my_int.cum_haz(cum_int_arg, y);
      return y[0L];
    };

    /* compute other intermediaries */
    Type log_det_sigma;
    matrix<Type> sigma_inv = atomic::matinvpd(aSigma, log_det_sigma);
    Type log_det_psi;
    matrix<Type> psi_inv = atomic::matinvpd(aPsi, log_det_psi);

    Type const one(1.),
       sqrt_two_pi(sqrt(2. / M_PI)),
              half(.5),
            two_pi(2. / M_PI),
             small(std::numeric_limits<double>::epsilon());

    GaussHermite::HermiteData<Type> const &xw =
      GaussHermite::GaussHermiteDataCached<Type>(n_nodes);

    /* evaluate the lower bound */
    survTMB::accumulator_mock<Type> result;
    bool const is_in_parallel = CppAD::thread_alloc::in_parallel();
    arma::vec b_wrk(my_int.b.get_n_basis()),
              m_wrk(my_int.m.get_n_basis());
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

      /* terms from the survival outcome */
      {
        Type surv_term(0.);
        Type const lb = tstart[g],
                   ub = tstop [g];
        Type z_dot_d(0.);
        for(int j = 0; j < Z.rows(); ++j)
          z_dot_d += Z(j, g) * adelta[j];

        if(asDouble(outcomes[g]) > 0.){
          /* add the jump term */
          my_int.b(b_wrk, log(asDouble(ub)));
          for(size_t j = 0; j < b_wrk.n_elem; ++j)
            surv_term += Type(b_wrk[j]) * aomega[j];
          surv_term += z_dot_d;

          my_int.m(m_wrk, asDouble(ub));
          vector<Type> const type_m_wrk = vec_eigen_arma<Type>(m_wrk);
          for(int j = 0; j < B.cols(); ++j){
            int const inc = j * B.rows();
            for(int i = 0; i < B.rows(); ++i)
              surv_term += aalpha[j] * type_m_wrk[i] * (
                aB(i, j) + U_mean[i + inc]);
          }
        }

        /* add term from survival probability */
        vector<Type> U = va_mu;
        for(int j = 0; j < B.cols(); ++j){
          int const inc = j * B.rows();
          for(int i = 0; i < B.rows(); ++i)
            U[i + inc] += aB(i, j);
        }

        Type const f1 = exp(z_dot_d),
                   f2 = comp_cum_haz(lb, ub, U, k, Lambda);
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
            residual[i] = markers(i, marker_idx);

          my_int.m(m_wrk, asDouble(m_time[marker_idx]));
          vector<Type> type_m_wrk = vec_eigen_arma<Type>(m_wrk);
          Type const *mu_i = &U_mean[0L];
          for(size_t j = 0; j < n_y; ++j)
            for(size_t i = 0; i < dim_m; ++i)
              residual[j] -= (*mu_i++ + aB(i, j)) * type_m_wrk[i];

          marker_term -= half * quad_form_sym(residual, sigma_inv);

          /* TODO: Do something smarter... */
          matrix<Type> M(n_y, K);
          M.setZero();
          for(size_t i = 0; i < n_y; ++i){
            Type *m_i = &type_m_wrk[0L];
            size_t const j_end =  (i + 1L) * dim_m;
            for(size_t j = i * dim_m; j < j_end; ++j)
              M(i, j) = *m_i++;
          }

          matrix<Type> tmp_mat = M * Lambda * M.transpose();
          vector<Type> tmp_vec = M * k;
          marker_term -=
            half * (mat_mult_trace(tmp_mat, sigma_inv)
            - two_pi * quad_form_sym(tmp_vec, sigma_inv));

          term += marker_term;
        }

        /* terms from the random effect prior*/
        /* TODO: can be done smarter */
        Type prior_term = -half * (
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
          vec_dot(va_rho_scaled, va_rho_scaled), xw);
        term += misc_term;
      }

      result += term;
    }

    if(!is_my_region(*result.obj))
      /* only have to add one more term so just return */
      return result;

    Type norm_constant =
      - Type(double(markers.cols() * n_y) / 2. * log(2 * M_PI))
      - Type(double(n_groups) / 2.) * log_det_psi
      - Type(double(markers.cols()) / 2.) * log_det_sigma
      + Type(double(K * n_groups) / 2.)
      - Type(double(n_groups) * M_LN2);

    result += norm_constant;

    return result;
  }
};

struct VA_func {
  using ADd   = CppAD::AD<double>;
  using ADdd  = CppAD::AD<ADd>;
  using ADddd = CppAD::AD<ADdd>;
  template<class Type>
  using ADFun = CppAD::ADFun<Type>;

  template<class Type>
  using splines_n_cum_haz = typename VA_worker<Type>::splines_n_cum_haz;
  std::vector<splines_n_cum_haz<ADd> >splines_n_cum_ints_ADd;

  std::vector<std::unique_ptr<ADFun<double> > > funcs;

  VA_func(Rcpp::List data, Rcpp::List parameters){
    {
      /* to compute function and gradient */
      VA_worker<ADd> w(data, parameters);
      splines_n_cum_ints_ADd = w.get_splines_n_cum_ints();
      funcs.resize(w.n_blocks);
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
  setup_parallel_ad<          double  > setup_ADd(n_threads);
  return Rcpp::XPtr<VA_func>(new VA_func(data, parameters));
}

// [[Rcpp::export(rng = false)]]
double joint_funcs_eval_lb(SEXP p, SEXP par){
  shut_up();

  Rcpp::XPtr<VA_func > ptr(p);
  std::vector<std::unique_ptr<CppAD::ADFun<double> > > &funcs = ptr->funcs;
  vector<double> parv = get_vec<double>(par);

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
