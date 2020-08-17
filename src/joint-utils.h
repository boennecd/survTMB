#ifndef JOINT_UTILS_H
#define JOINT_UTILS_H

#define INCLUDE_RCPP
#include "tmb_includes.h"

#include "fastgl.h"
#include "pnorm-log.h"
#include "memory.h"

namespace fastgl {
namespace joint {

/*
  This function performs an approximation of

 \begin{align*}
 y &= \int_l^u g(o; \vec g) d o \\
 g(o; \vec g) &= 2\exp(
 \vec\omega^\top\vec b(o)
 + G_{\vec\alpha}(o)\vec b
 + M_{\vec\alpha}(o)\vec U
 +\frac 12 M_{\vec\alpha}(o)\Lambda_i
 M_{\vec\alpha}(o)^\top)\Phi(
 M_{\vec\alpha}(o)
 \vec k) \\
 \vec g &= (\vec \omega^\top, \vec\alpha^\top, \vec b^\top,
 \vec U^\top, \vec k^\top,  (\text{vec} \Lambda)^\top
 )^\top\in\mathbb{R}^{e + r + c + K(2 + K)} \                     \
 \vec\omega &\in\mathbb{R}^e \                                    \
 \vec\alpha &\in\mathbb{R}^r \                                    \
 \vec B &\in\mathbb{R}^c \                                        \
 \vec U,\vec k &\in\mathbb{R}^K \                                 \
 \Lambda &\in \mathbb{R}^{K\times K} \                            \
 G_{\vec\alpha}(t) &= \vec\alpha^\top (I \otimes \vec g(t)^\top) \\
 M_{\vec\alpha}(t) &= \vec\alpha^\top (I \otimes \vec m(t)^\top) \\
 \vec b(t) &:\, (0,\infty) \rightarrow \mathbb{R}^e \             \
 \vec g(t) &:\, (0,\infty) \rightarrow \mathbb{R}^l \             \
 \vec m(t) &:\, (0,\infty) \rightarrow \mathbb{R}^s \             \
 \end{align*}
 */
template<class Type, class B, class G, class M>
class snva_integral : public CppAD::atomic_base<Type> {
  size_t const n_nodes;
  std::vector<QuadPair<double> > const& xw =
    GLPairsCached<double>(n_nodes);
  std::vector<QuadPair<Type> > const& xw_type =
    GLPairsCached<Type>(n_nodes);

  std::unique_ptr<B> const b;
  std::unique_ptr<G> const g;
  std::unique_ptr<M> const m;

  bool const has_b = static_cast<bool>(b),
             has_g = static_cast<bool>(g),
             has_m = static_cast<bool>(m);

  size_t const dim_omega = has_b ? b->get_n_basis() : 0L,
               dim_alpha,
               dim_g     = has_g ? g->get_n_basis() : 0L,
               dim_B     = dim_g * dim_alpha,
               dim_m     = has_m ? m->get_n_basis() : 0L,
               dim_U     = dim_alpha * dim_m;
  bool const use_log;

  /* objects needed for function evaluation */
  mutable arma::vec fbi = arma::vec(dim_omega),
                    fgi = arma::vec(dim_g),
                    fmi = arma::vec(dim_m);
  mutable vector<double> fma = vector<double>(dim_U),
                         fga = vector<double>(dim_B),
                      fomega = vector<double>(dim_omega),
                      falpha = vector<double>(dim_alpha),
                          fB = vector<double>(dim_B),
                          fU = vector<double>(dim_U),
                          fk = vector<double>(dim_U);
  mutable matrix<double> fLambda = matrix<double>(dim_U, dim_U);

  /* objects needed for reverse evaluation */
  mutable vector<Type> rbi = vector<Type>(dim_omega),
                       rgi = vector<Type>(dim_g),
                       rmi = vector<Type>(dim_m),
                       rma = vector<Type>(dim_U),
                       rga = vector<Type>(dim_B),
                    romega = vector<Type>(dim_omega),
                    ralpha = vector<Type>(dim_alpha),
                 alpha_rhs = vector<Type>(dim_U),
                        rB = vector<Type>(dim_B),
                        rU = vector<Type>(dim_U),
                        rk = vector<Type>(dim_U);
  mutable matrix<Type> rLambda = matrix<Type>(dim_U, dim_U);

public:
  snva_integral(char const *name, size_t const n_nodes,
                B const *b_in, G const *g_in, M const *m_in,
                size_t const dim_alpha, bool const use_log):
  CppAD::atomic_base<Type>(name), n_nodes(n_nodes),
  b(b_in ? new B(*b_in) : nullptr),
  g(g_in ? new G(*g_in) : nullptr),
  m(m_in ? new M(*m_in) : nullptr),
  dim_alpha(dim_alpha),
  use_log(use_log) {
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

      auto set_vec = [&](vector<double> &x){
        for(int j = 0; j < x.size(); ++j)
          x[j] = asDouble(tx[i++]);
      };

      if(has_b)
        set_vec(fomega);
      set_vec(falpha);
      if(has_g)
        set_vec(fB);
      if(has_m){
        set_vec(fU);
        set_vec(fk);

        for(size_t j = 0; j < dim_U; ++j)
          for(size_t k = 0; k < dim_U; ++k)
            fLambda(k, j) = asDouble(tx[i++]);
      }
    }

    /* perform an approximation of the integral */
    double const d1 = (ub - lb) / 2.,
                 d2 = (ub + lb) / 2.;

    double out(0.);
    for(auto const &xwi : xw){
      /* evalutes splines and related objects */
      double const node = d1 * xwi.x + d2;
      if(has_m){
        m->operator()(fmi, node);

        size_t i(0L);
        for(size_t j = 0; j < dim_alpha; ++j)
          for(size_t k = 0; k < dim_m; ++k)
            fma[i++] = falpha[j] * fmi[k];
      }

      if(has_g){
        g->operator()(fgi, node);

        size_t i(0L);
        for(size_t j = 0; j < dim_alpha; ++j)
          for(size_t k = 0; k < dim_g; ++k)
            fga[i++] = falpha[j] * fgi[k];
      }

      if(has_b)
        if(use_log)
          b->operator()(fbi, log(node));
        else
          b->operator()(fbi,     node );

      /* evalute integrand */
      double const v1 = pnorm_log(vec_dot(fma, fk)),
                   v2 = vec_dot(fomega, fbi) +
                     vec_dot(fga, fB) +
                     vec_dot(fma, fU) +
                     .5 * quad_form_sym(fma, fLambda);

      out += xwi.weight * exp(v1 + v2);
    }

    out *= has_m ? d1 * 2 : d1 * 4;

    ty[0L] = Type(out);

    /* set variable flags */
    if (vx.size() > 0) {
      bool anyvx = false;
      for (std::size_t i = 0; i < vx.size() and !anyvx; i++)
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

      auto set_vec = [&](vector<Type> &x){
        for(int j = 0; j < x.size(); ++j)
          x[j] = tx[i++];
      };

      if(has_b)
        set_vec(romega);
      set_vec(ralpha);
      if(has_g)
        set_vec(rB);
      if(has_m){
        set_vec(rU);
        set_vec(rk);

        for(size_t j = 0; j < dim_U; ++j)
          for(size_t k = 0; k < dim_U; ++k)
            rLambda(k, j) = tx[i++];
      }
    }

    /* perform an approximation of the integral */
    Type const d1((ub - lb) / 2.);
    double const d2 = (ub + lb) / 2.;

    for(size_t i = 0; i < px.size(); ++i)
      px[i] = Type(0.);

    Type const ZERO(0.), ONE(1.), HALF(.5);

    for(auto xwi : xw_type){
      /* evalutes splines and related objects */
      double const node = asDouble(d1 * xwi.x) + d2;
      if(has_m){
        m->operator()(fmi, node);
        for(size_t i = 0; i < fmi.size(); ++i)
          rmi[i] = Type(fmi[i]);

        size_t i(0L);
        for(size_t j = 0; j < dim_alpha; ++j)
          for(size_t k = 0; k < dim_m; ++k)
            rma[i++] = ralpha[j] * rmi[k];
      }

      if(has_g){
        g->operator()(fgi, node);
        for(size_t i = 0; i < fgi.size(); ++i)
          rgi[i] = Type(fgi[i]);

        size_t i(0L);
        for(size_t j = 0; j < dim_alpha; ++j)
          for(size_t k = 0; k < dim_g; ++k)
            rga[i++] = ralpha[j] * rgi[k];
      }

      if(has_b){
        if(use_log)
          b->operator()(fbi, log(node));
        else
          b->operator()(fbi,     node );
        for(size_t i = 0; i < fbi.size(); ++i)
          rbi[i] = Type(fbi[i]);
      }

      /* evaluate intermediary constants */
      vector<Type> lambda_ma = rLambda * rma;
      Type const ma_k = vec_dot(rma, rk),
                   v1 = pnorm_log(ma_k),
                   v2 = vec_dot(romega, rbi) +
                     vec_dot(rga, rB) +
                     vec_dot(rma, rU) +
                     HALF * vec_dot(rma, lambda_ma),
                   v3 = dnorm(ma_k, ZERO, ONE, 1L),
            integrand = xwi.weight * exp(v1 + v2),
                    h = xwi.weight * exp(v3 + v2);

      /* add terms to gradient */
      {
        size_t i(2L);
        /* omega */
        for(size_t j = 0; j < dim_omega; ++j)
          px[i++] += integrand * rbi[j];
        /* alpha */
        for(size_t j = 0; j < dim_U; ++j){
          alpha_rhs[j]  = rU[j];
          alpha_rhs[j] += lambda_ma[j];
          alpha_rhs[j] *= integrand;
          alpha_rhs[j] += h * rk[j];
        }
        {
          size_t sm(0L), sg(0L);
          for(size_t j = 0; j < dim_alpha; ++j){
            Type term_m(0.);
            for(size_t k = 0; k < dim_m; ++k)
              term_m += rmi[k] * alpha_rhs[sm++];

            Type term_g(0.);
            for(size_t k = 0; k < dim_g; ++k)
              term_g += rgi[k] * rB[sg++];
            term_g *= integrand;

            px[i++] += term_m + term_g;
          }
        }
        /* B */
        for(size_t j = 0; j < dim_B; ++j)
          px[i++] += integrand * rga[j];
        /* U */
        for(size_t j = 0; j < dim_U; ++j)
          px[i++] += integrand * rma[j];
        /* K */
        for(size_t j = 0; j < dim_U; ++j)
          px[i++] += h * rma[j];
        /* Lambda */
        Type const mult = HALF * integrand;
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

    Type const mult =
      has_m ? Type(ub - lb) * py[0L] : 2 * Type(ub - lb) * py[0L];
    for(size_t i = 2L; i < px.size(); ++i)
      px[i] *= mult;

    return true;
  }

  virtual bool rev_sparse_jac(size_t q, const CppAD::vector<bool>& rt,
                              CppAD::vector<bool>& st) {
    bool anyrt = false;
    for (std::size_t i = 0; i < rt.size() and !anyrt; i++)
      anyrt |= rt[i];
    for (std::size_t i = 0; i < st.size(); i++)
      st[i] = anyrt;
    return true;
  }

  template<class T>
  CppAD::vector<AD<T> > get_x
    (AD<T> const lb, AD<T> const ub,
     vector<AD<T> > const &omega, vector<AD<T> > const &alpha,
     matrix<AD<T> > const &b_arg, vector<AD<T> > const &U,
     vector<AD<T> > const &k, matrix<AD<T> > const &Lambda) const {
    size_t const n_ele =
      2L + dim_omega + dim_alpha + dim_B + dim_U * (2L + dim_U);
#ifdef DO_CHECKS
    if(has_b and omega.size() != (int)dim_omega)
      throw std::runtime_error("get_x: invalid omega");
    if(alpha.size() != (int)dim_alpha)
      throw std::runtime_error("get_x: invalid alpha");
    if(b_arg.cols() != (int)dim_alpha or
         (has_g and b_arg.rows() != (int)dim_g))
      throw std::runtime_error("get_x: invalid b_arg");
    if(has_m and U.size() != (int)dim_U)
      throw std::runtime_error("get_x: invalid U");
    if(has_m and k.size() != (int)dim_U)
      throw std::runtime_error("get_x: invalid k");
    if(has_m and
         (Lambda.rows() != (int)dim_U or Lambda.cols() != (int)dim_U))
      throw std::runtime_error("get_x: invalid Lambda");
#endif
    CppAD::vector<AD<Type> > tx(n_ele);

    size_t i(0L);
    tx[i++] = lb;
    tx[i++] = ub;
    auto add_vec = [&](vector<AD<T> > const &x){
      for(int j = 0; j < x.size(); ++j)
        tx[i++] = x[j];
    };
    auto add_mat = [&](matrix<AD<T> > const &X){
      for(int j = 0; j < X.cols(); ++j)
        for(int k = 0; k < X.rows(); ++k)
          tx[i++] = X(k, j);
    };

    if(has_b)
      add_vec(omega);
    add_vec(alpha);
    if(has_g)
      add_mat(b_arg);
    if(has_m){
      add_vec(U);
      add_vec(k);
      add_mat(Lambda);
    }

    return tx;
  }
};

} // namespace joint
} // namespace fastgl

#endif
