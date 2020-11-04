#define INCLUDE_RCPP
#include "pedigree-Laplace.h"
#include "get-x.h"
#include "utils.h"
#include "pnorm-log.h"
#include "laplace-util.h"

namespace {

/* object to hold the data for a given cluster */
template<class Type>
struct cluster_data {
  Rcpp::List data;

  const DATA_VECTOR(event);
  const DATA_MATRIX(X);
  const DATA_MATRIX(XD);
  const DATA_MATRIX(Z);

  std::vector<matrix<Type> > const cor_mats = ([&](){
    Rcpp::List list_w_cor_mats = data["cor_mats"];

    std::vector<matrix<Type> > out;
    out.reserve(list_w_cor_mats.size());
    for(auto x : list_w_cor_mats)
      out.emplace_back(get_mat<Type>(x));

    return out;
  })();

  R_len_t const n_members = event.size();

  /** not thread-safe because of the R interaction! */
  cluster_data(Rcpp::List data): data(data) {
    if(event.size() != n_members)
      throw std::runtime_error("cluster_data<Type>: invalid event");
    else if(X.cols()  != n_members)
      throw std::runtime_error("cluster_data<Type>: invalid X");
    else if(XD.cols() != n_members)
      throw std::runtime_error("cluster_data<Type>: invalid XD");
    else if(Z.cols() != n_members)
      throw std::runtime_error("cluster_data<Type>: invalid Z");
    for(auto &V : cor_mats)
      if(V.rows() != n_members or V.cols() != n_members)
        throw std::runtime_error("cluster_data<Type>: invalid cor_mats");

    data = R_NilValue;
  }

  void check(vector<Type> const &omega, vector<Type> const &beta,
             vector<Type> const &log_sds, vector<Type> const &rng_mode){
    if(rng_mode.size() != n_members)
      throw std::runtime_error("cluster_data<Type>: invalid va_par");
    if(X.rows() != omega.size())
      throw std::runtime_error("cluster_data<Type>: invalid c_data (X)");
    else if(Z.rows() != beta.size())
      throw std::runtime_error("cluster_data<Type>: invalid c_data (Z)");
    else if(cor_mats.size() != static_cast<size_t>(log_sds.size()))
      throw std::runtime_error(
          "cluster_data<Type>: invalid c_data (cor_mats)");
    else if(XD.rows() != omega.size())
      throw std::runtime_error("cluster_data<Type>: invalid c_data (XD)");
  }
};

} // namespace

namespace survTMB {

template<class Type>
void pedigree_laplace
  (parallel_accumulator<Type> &out, SEXP dat, vector<Type> const &omega,
   vector<Type> const &beta, vector<Type> const &log_sds,
   vector<Type> const &rng_modes){
  // setup
  Rcpp::List data(dat),
            c_dat(static_cast<SEXP>(data["c_dat"]));

  DATA_SCALAR(eps);
  DATA_SCALAR(kappa);
  DATA_STRING(link);

  Type const two(2.),
         eps_log = log(eps);

  vector<Type> const var = ([&](){
    vector<Type> out(log_sds.size());
    for(int i = 0; i < out.size(); ++i)
      out[i] = exp(two * log_sds[i]);
    return out;
  })();

  // compute the output cluster by cluster
  Type const * mi = &rng_modes[0L];
  for(auto di : c_dat){
    cluster_data<Type> cdi(Rcpp::List(static_cast<SEXP>(di)));
    R_len_t const n_members = cdi.n_members;
    if(!is_my_region(*out.obj)){
      // not this threads region
      mi += n_members;
      out.obj->parallel_region();
      continue;
    }

    Type new_terms(0.);

    // perform the computation. Start with the conditional density terms
    vector<Type> c_mode(n_members);
    cdi.check(omega, beta, log_sds, c_mode);
    for(R_len_t i = 0; i < n_members; ++i, ++mi){
      c_mode[i] = *mi;
      vector<Type> const x  = cdi.X .col(i),
                         xd = cdi.XD.col(i),
                         z  = cdi.Z .col(i);
      Type const eta = vec_dot(omega, x) + vec_dot(beta, z) + c_mode[i],
                etaD = vec_dot(omega, xd),
               event = cdi.event[i];

      if     (link == "PH")
        new_terms -= survTMB::laplace_PH_terms<Type>(
          eta, etaD, event, eps, eps_log, kappa);
      else if(link == "PO")
        new_terms -= survTMB::laplace_PO_terms<Type>(
          eta, etaD, event, eps, eps_log, kappa);
      else if(link == "probit")
        new_terms -= survTMB::laplace_probit_terms<Type>(
          eta, etaD, event, eps, eps_log, kappa);
      else
        throw std::runtime_error("link not implemented");
    }

    // handle terms from the unconditional density of the random
    // effects
    matrix<Type> sigma(n_members, n_members);
    sigma.setZero();
    for(R_len_t i = 0; i < var.size(); ++i)
      for(R_len_t j = 0; j < n_members; ++j)
        for(R_len_t k = 0; k < n_members; ++k)
          sigma(k, j) += var[i] * cdi.cor_mats[i](k, j);

    density::MVNORM_t<Type> norm_dist(sigma, true);
    new_terms += norm_dist(c_mode); // this is the negative log density
    out += new_terms;
  }
}

using ADd   = CppAD::AD<double>;
using ADdd  = CppAD::AD<CppAD::AD<double> >;
using ADddd = CppAD::AD<CppAD::AD<CppAD::AD<double> > >;

template void pedigree_laplace<double>
  (parallel_accumulator<double> &out, SEXP dat,
   vector<double> const &omega, vector<double> const &beta,
   vector<double> const &log_sds, vector<double> const &rng_modes);
template void pedigree_laplace<ADd   >
  (parallel_accumulator<ADd   > &out, SEXP dat,
   vector<ADd   > const &omega, vector<ADd   > const &beta,
   vector<ADd   > const &log_sds, vector<ADd   > const &rng_modes);
template void pedigree_laplace<ADdd  >
  (parallel_accumulator<ADdd  > &out, SEXP dat,
   vector<ADdd  > const &omega, vector<ADdd  > const &beta,
   vector<ADdd  > const &log_sds, vector<ADdd  > const &rng_modes);
template void pedigree_laplace<ADddd >
  (parallel_accumulator<ADddd > &out, SEXP dat,
   vector<ADddd > const &omega, vector<ADddd > const &beta,
   vector<ADddd > const &log_sds, vector<ADddd > const &rng_modes);
} // namespace survTMB
