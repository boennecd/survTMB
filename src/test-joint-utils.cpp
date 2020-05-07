#include "testthat-wrap.h"
#include "joint-utils.h"
#include "splines.h"
#include <vector>

using namespace fastgl::joint;

context("joint-utils unit tests") {
  test_that("eval_snva_integral gives the correct result") {
    using ADd = AD<double>;
    constexpr size_t const dim_o(3L),
                           dim_a(2L),
                           dim_m(3L),
                               K = dim_a * dim_m,
                         n_nodes(30L);

    vector<ADd> omega(dim_o);
    omega << -0.74, 0.76, 0.20;

    vector<ADd> alpha(dim_a);
    alpha << 0.7, 0.6;

    vector<ADd> U(K),
                k(K);
    U << 0.18, -0.66, -0.02, -1.49, -0.99, 0.65;
    k << -0.155666233880193, -0.202626997173659, 0.480913001875679,
         -0.185234121879782, 0.138273358586317, -0.056526844705098;

    matrix<ADd> Lambda(K, K);
    Lambda << 1.08, 0.12, -0.36, -0.48, 0.36, -0.12, 0.12, 0.36, 0, 0.12,
              0, -0.12, -0.36, 0, 0.84, 0.12, 0.12, 0.12, -0.48, 0.12, 0.12,
              0.84, -0.12, 0.24, 0.36, 0, 0.12, -0.12, 0.84, -0.12, -0.12,
              -0.12, 0.12, 0.24, -0.12, 0.6;

    ADd const lb(1.2), ub(9.5);

    arma::vec m_bk(2L);
    m_bk[0L] = 0.;
    m_bk[1L] = 10;
    arma::vec m_ik(1L);
    m_ik[0L] = 5;
    splines::ns m(m_bk, m_ik, true);

    arma::vec b_bk(2L);
    b_bk[0L] = 0.;
    b_bk[1L] = 2.30258509299405;
    arma::vec b_ik(2L);
    b_ik[0L] = 0.767528364331349;
    b_ik[1L] = 1.5350567286627;
    splines::ns b(b_bk, b_ik, false);

    snva_integral<double, splines::ns, splines::ns> func(
          "snva_integral", n_nodes, b, m, dim_a);

    auto x = func.get_x(lb, ub, omega, alpha, U, k, Lambda);
    {
      CppAD::vector<ADd> y(1L);
      CppAD::Independent(x);
      func(x, y);
      CppAD::ADFun<double> afunc(x, y);
      afunc.optimize();

      expect_equal(5.83072289287952, asDouble(y[0L]));

      CppAD::vector<double> xx(x.size());
      for(size_t i = 0; i < x.size(); ++i)
        xx[i] = asDouble(x[i]);

      auto yy = afunc.Forward(0, xx);
      expect_equal(5.83072289287952, yy[0L]);

      constexpr size_t n_grad_ele = 53L;
      constexpr double const grad[n_grad_ele] = {
        1.06189406698313, 2.29990047787843, 1.12303952006957, -0.539084307993634,
        -2.63809475605053, 0.867304891575122, 1.6389584379032, 0.552593787296199,
        0.743404192778676, 1.40482151820274, 0.473651817682456, 0.734949096311874,
        1.34505017304729, 0.356446809637724, 0.629956368267321, 1.15290014832625,
        0.305525836832335, 0.108269380727632, 0.114400655605564, 0.0251879067710015,
        0.0928023263379704, 0.098057704804769, 0.0215896343751442, 0.114400655605564,
        0.233102582120057, 0.0653672205918621, 0.098057704804769, 0.199802213245763,
        0.0560290462215961, 0.0251879067710015, 0.0653672205918621, 0.185308170946502,
        0.0215896343751442, 0.0560290462215961, 0.158835575097001, 0.0928023263379704,
        0.098057704804769, 0.0215896343751442, 0.0795448511468318, 0.0840494612612306,
        0.0185054008929807, 0.098057704804769, 0.199802213245763, 0.0560290462215961,
        0.0840494612612306, 0.17125903992494, 0.0480248967613681, 0.0215896343751442,
        0.0560290462215961, 0.158835575097001, 0.0185054008929807, 0.0480248967613681,
        0.136144778654573 };

      vector<double> w(1L);
      w[0L] = 1.;
      auto dx = afunc.Reverse(1, w);
      expect_true(dx.size() == n_grad_ele + 2L);

      {
        size_t i = 0L;
        /* omega  */
        for(size_t j = 0; j < dim_o; ++j, ++i)
          expect_equal(grad[i], dx[i + 2L]);
        /* alpha */
        for(size_t j = 0; j < dim_a; ++j, ++i)
          expect_equal(grad[i], dx[i + 2L]);
        /* U */
        for(size_t j = 0; j < K; ++j, ++i)
          expect_equal(grad[i], dx[i + 2L]);
        /* k */
        for(size_t j = 0; j < K; ++j, ++i)
          expect_equal(grad[i], dx[i + 2L]);
        /* Lambda */
        for(size_t j = 0; j < K * K; ++j, ++i)
          expect_equal(grad[i], dx[i + 2L]);
      }

      /* change upper and lower bounds and compute agian */
      xx[0L] = 2.4;
      xx[1L] = 6.7;

      yy = afunc.Forward(0, xx);
      expect_equal(1.97650616568221, yy[0L]);
      constexpr double const grad2[n_grad_ele] = {
        0.810859563385157, 0.775789550032837, -0.148800927381389, -0.290605976451132,
        -1.75197653501335, 0.52211309198903, 0.541534259100077, -0.151034243578911,
        0.447525507419168, 0.46417222208578, -0.129457923067638, 0.462987273833512,
        0.481704805733632, -0.137187065918621, 0.396846234714439, 0.41288983348597,
        -0.117588913644532, 0.0744002796404626, 0.0698790585360406, -0.0146263142437222,
        0.0637716682632537, 0.0598963358880348, -0.0125368407803333,
        0.0698790585360406, 0.0747181299338131, -0.0226422144324493,
        0.0598963358880348, 0.0640441113718398, -0.0194076123706708,
        -0.0146263142437222, -0.0226422144324493, 0.0144534802620029,
        -0.0125368407803333, -0.0194076123706708, 0.012388697367431,
        0.0637716682632537, 0.0598963358880348, -0.0125368407803333,
        0.0546614299399317, 0.0513397164754584, -0.010745863526, 0.0598963358880348,
        0.0640441113718398, -0.0194076123706708, 0.0513397164754584,
        0.0548949526044341, -0.0166350963177178, -0.0125368407803333,
        -0.0194076123706708, 0.012388697367431, -0.010745863526, -0.0166350963177178,
        0.010618883457798 };
      dx = afunc.Reverse(1, w);
      expect_true(dx.size() == n_grad_ele + 2L);
      {
        size_t i = 0L;
        /* omega  */
        for(size_t j = 0; j < dim_o; ++j, ++i)
          expect_equal(grad2[i], dx[i + 2L]);
        /* alpha */
        for(size_t j = 0; j < dim_a; ++j, ++i)
          expect_equal(grad2[i], dx[i + 2L]);
        /* U */
        for(size_t j = 0; j < K; ++j, ++i)
          expect_equal(grad2[i], dx[i + 2L]);
        /* k */
        for(size_t j = 0; j < K; ++j, ++i)
          expect_equal(grad2[i], dx[i + 2L]);
        /* Lambda */
        for(size_t j = 0; j < K * K; ++j, ++i)
          expect_equal(grad2[i], dx[i + 2L]);
      }
    }
  }
}
