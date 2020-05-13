#ifndef JOINT_UTILS_H
#define JOINT_UTILS_H

#define INCLUDE_RCPP
#include "tmb_includes.h"

#include "fastgl.h"
#include "pnorm-log.h"

namespace fastgl {
namespace joint {

/*
  This function performs an approximation of

 \begin{align*}
 y &= \int_l^u g(o; \vec g) d o \\
 g(o; \vec g) &= 2\exp(
 \vec\omega^\top\vec b(o)
 + M_{\vec\alpha}(o)\vec U
 +\frac 12 M_{\vec\alpha}(o)\Lambda_i
 M_{\vec\alpha}(o)^\top)\Phi(
 M_{\vec\alpha}(o)
 \vec k) \\
 \vec g &= (\vec \omega^\top, \vec\alpha^\top, \vec U^\top,
 \vec k^\top,  (\text{vec} \Lambda)^\top
 )^\top\in\mathbb{R}^{e + r + K(2 + K)} \\
 \vec\omega &\in\mathbb{R}^e \\
 \vec\alpha &\in\mathbb{R}^r \\
 \vec U,\vec k &\in\mathbb{R}^K \\
 \Lambda &\in \mathbb{R}^{K\times K} \\
 M_{\vec\alpha}(t) &= \vec\alpha^\top (I \otimes \vec m(t)^\top) \\
 \vec m(t) &:\, (0,\infty) \rightarrow \mathbb{R}^s \\
 \vec b(t) &:\, (0,\infty) \rightarrow \mathbb{R}^e
 \end{align*}
 */
template<class Type, class B, class M>
class snva_integral : public CppAD::atomic_base<Type> {
  size_t const n_nodes;
  std::vector<QuadPair<double> > const& xw =
    GLPairsCached<double>(n_nodes);
  std::vector<QuadPair<Type> > const& xw_type =
    GLPairsCached<Type>(n_nodes);

  B const b;
  M const m;

  size_t const dim_omega = b.get_n_basis(),
               dim_alpha,
               dim_m     = m.get_n_basis(),
               dim_U     = dim_alpha * dim_m;
  /* objects needed for function evaluation */
  mutable arma::vec fbi = arma::vec(b.get_n_basis()),
                    fmi = arma::vec(m.get_n_basis());
  mutable vector<double> fma = vector<double>(dim_U),
                      fomega = vector<double>(dim_omega),
                      falpha = vector<double>(dim_alpha),
                          fU = vector<double>(dim_U),
                          fk = vector<double>(dim_U);
  mutable matrix<double> fLambda = matrix<double>(dim_U, dim_U);

  /* objects needed for reverse evaluation */
  mutable vector<Type> rbi = vector<Type>(b.get_n_basis()),
                       rmi = vector<Type>(m.get_n_basis()),
                       rma = vector<Type>(dim_U),
                    romega = vector<Type>(dim_omega),
                    ralpha = vector<Type>(dim_alpha),
                 alpha_rhs = vector<Type>(dim_U),
                        rU = vector<Type>(dim_U),
                        rk = vector<Type>(dim_U);
  mutable matrix<Type> rLambda = matrix<Type>(dim_U, dim_U);

public:
  snva_integral(char const *name, size_t const n_nodes,
                B const &b, M const &m, size_t const dim_alpha):
  CppAD::atomic_base<Type>(name), n_nodes(n_nodes), b(b), m(m),
  dim_alpha(dim_alpha) {
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
  }

  virtual bool forward(std::size_t p, std::size_t q,
                       const CppAD::vector<bool> &vx,
                       CppAD::vector<bool> &vy,
                       const CppAD::vector<Type> &tx,
                       CppAD::vector<Type> &ty){
    if(q > 0)
      return false;

    /* get parameters and integral bounds */
    double lb, ub;
    {
      size_t i(0L);
      lb = asDouble(tx[i++]),
      ub = asDouble(tx[i++]);

      for(size_t j = 0; j < dim_omega; ++j)
        fomega[j] = asDouble(tx[i++]);
      for(size_t j = 0; j < dim_alpha; ++j)
        falpha[j] = asDouble(tx[i++]);
      for(size_t j = 0; j < dim_U; ++j)
        fU[j] = asDouble(tx[i++]);
      for(size_t j = 0; j < dim_U; ++j)
        fk[j] = asDouble(tx[i++]);
      for(size_t j = 0; j < dim_U; ++j)
        for(size_t k = 0; k < dim_U; ++k)
          fLambda(k, j) = asDouble(tx[i++]);
    }

    /* perform an approximation of the integral */
    double const d1 = (ub - lb) / 2.,
                 d2 = (ub + lb) / 2.;

    double out(0.);
    for(auto const &xwi : xw){
      /* evalutes splines and related objects */
      {
        double const node = d1 * xwi.x + d2;
        m(fmi, node);
        b(fbi, log(node));

        {
          size_t i(0L);
          for(size_t j = 0; j < dim_alpha; ++j)
            for(size_t k = 0; k < m.get_n_basis(); ++k)
              fma[i++] = falpha[j] * fmi[k];
        }
      }

      /* evalute integrand */
      double dot_o_b(0.);
      for(size_t i = 0; i < dim_omega; ++i)
        dot_o_b += fomega[i] * fbi[i];

      double const v1 = pnorm_log((fma * fk).sum()),
                   v2 = dot_o_b +
                     (fma * fU).sum() + .5 * (fma * (fLambda * fma)).sum();

      out += xwi.weight * exp(v1 + v2);
    }

    out *= d1 * 2;
    ty[0L] = Type(out);

    /* set variable flags */
    if (vx.size() > 0) {
      bool anyvx = false;
      for (std::size_t i = 0; i < vx.size(); i++)
        anyvx |= vx[i];
      for (std::size_t i = 0; i < vy.size(); i++)
        vy[i] = anyvx;
    }

    return true;
  }

  virtual bool reverse(std::size_t q, const CppAD::vector<Type> &tx,
                       const CppAD::vector<Type> &ty,
                       CppAD::vector<Type> &px,
                       const CppAD::vector<Type> &py){
    if(q > 0)
      return false;

    /* get parameters and integral bounds */
    double lb, ub;
    {
      size_t i(0L);
      lb = asDouble(tx[i++]),
      ub = asDouble(tx[i++]);

      for(size_t j = 0; j < dim_omega; ++j)
        romega[j] = tx[i++];
      for(size_t j = 0; j < dim_alpha; ++j)
        ralpha[j] = tx[i++];
      for(size_t j = 0; j < dim_U; ++j)
        rU[j] = tx[i++];
      for(size_t j = 0; j < dim_U; ++j)
        rk[j] = tx[i++];
      for(size_t j = 0; j < dim_U; ++j)
        for(size_t k = 0; k < dim_U; ++k)
          rLambda(k, j) = tx[i++];
    }

    /* perform an approximation of the integral */
    Type const d1((ub - lb) / 2.);
    double const d2 = (ub + lb) / 2.;

    for(size_t i = 0; i < px.size(); ++i)
      px[i] = Type(0.);

    Type const ZERO(0.), ONE(1.), HALF(.5);

    for(auto xwi : xw_type){
      /* evalutes splines and related objects */
      {
        double const node = asDouble(d1 * xwi.x) + d2;
        m(fmi, node);
        b(fbi, log(node));

        for(size_t i = 0; i < fbi.size(); ++i)
          rbi[i] = Type(fbi[i]);
        for(size_t i = 0; i < fmi.size(); ++i)
          rmi[i] = Type(fmi[i]);

        {
          size_t i(0L);
          for(size_t j = 0; j < dim_alpha; ++j)
            for(size_t k = 0; k < dim_m; ++k)
              rma[i++] = ralpha[j] * rmi[k];
        }
      }

      /* evaluate intermediary constants */
      vector<Type> lambda_ma = rLambda * rma;
      Type const ma_k = (rma * rk).sum(),
                   v1 = pnorm_log(ma_k),
                   v2 = (romega * rbi).sum() + (rma * rU).sum() +
                     HALF * (rma * lambda_ma).sum(),
                   v3 = dnorm(ma_k, ZERO, ONE, 1L),
                    g = xwi.weight * exp(v1 + v2),
                    h = xwi.weight * exp(v3 + v2);

      /* add terms to gradient */
      {
        size_t i(2L);
        /* omega */
        for(size_t j = 0; j < dim_omega; ++j)
          px[i++] += g * rbi[j];
        /* alpha */
        for(size_t j = 0; j < dim_U; ++j){
          alpha_rhs[j]  = rU[j];
          alpha_rhs[j] += lambda_ma[j];
          alpha_rhs[j] *= g;
          alpha_rhs[j] += h * rk[j];
        }
        {
          size_t s(0L);
          for(size_t j = 0; j < dim_alpha; ++j){
            Type term(0.);
            for(size_t k = 0; k < dim_m; ++k)
              term += rmi[k] * alpha_rhs[s++];
            px[i++] += term;
          }
        }
        /* U */
        for(size_t j = 0; j < dim_U; ++j)
          px[i++] += g * rma[j];
        /* K */
        for(size_t j = 0; j < dim_U; ++j)
          px[i++] += h * rma[j];
        /* Lambda */
        Type const mult = HALF * g;
        for(size_t j = 0; j < dim_U; ++j){
          size_t const dj = j * dim_U;
          size_t dk(0L);
          for(size_t k = 0; k < j; ++k, dk += dim_U){
            Type const inc = mult * rma[j] * rma[k];
            px[i + k + dj] += inc;
            px[i + j + dk] += inc;
          }
          px[i + j + dj] +=  mult * rma[j] * rma[j];

        }
      }
    }

    Type const mult = Type(ub - lb) * py[0L];
    for(size_t i = 2L; i < px.size(); ++i)
      px[i] *= mult;

    return true;
  }

  virtual bool rev_sparse_jac(size_t q, const CppAD::vector<bool>& rt,
                              CppAD::vector<bool>& st) {
    bool anyrt = false;
    for (std::size_t i = 0; i < rt.size(); i++)
      anyrt |= rt[i];
    for (std::size_t i = 0; i < st.size(); i++)
      st[i] = anyrt;
    return true;
  }

  template<class T>
  CppAD::vector<AD<T> > get_x
    (AD<T> const lb, AD<T> const ub,
     vector<AD<T> > const &omega, vector<AD<T> > const &alpha,
     vector<AD<T> > const &U, vector<AD<T> > const &k,
     matrix<AD<T> > const &Lambda) const {
    size_t const n_ele = 2L + dim_omega + dim_alpha + dim_U * (2L + dim_U);
#ifdef DO_CHECKS
    if(omega.size() != (int)dim_omega)
      throw std::runtime_error("get_x: invalid omega");
    if(alpha.size() != (int)dim_alpha)
      throw std::runtime_error("get_x: invalid alpha");
    if(U.size() != (int)dim_U)
      throw std::runtime_error("get_x: invalid U");
    if(k.size() != (int)dim_U)
      throw std::runtime_error("get_x: invalid k");
    if(Lambda.rows() != (int)dim_U or Lambda.cols() != (int)dim_U)
      throw std::runtime_error("get_x: invalid Lambda");
#endif
    CppAD::vector<AD<Type> > tx(n_ele);

    size_t i(0L);
    tx[i++] = lb;
    tx[i++] = ub;
    for(int j = 0; j < omega.size(); ++j)
      tx[i++] = omega[j];
    for(int j = 0; j < alpha.size(); ++j)
      tx[i++] = alpha[j];
    for(int j = 0; j < U.size(); ++j)
      tx[i++] = U[j];
    for(int j = 0; j < k.size(); ++j)
      tx[i++] = k[j];
    for(int j = 0; j < Lambda.cols(); ++j)
      for(int k = 0; k < Lambda.rows(); ++k)
        tx[i++] = Lambda(k, j);

    return tx;
  }
};

} // namespace joint
} // namespace fastgl

#endif
