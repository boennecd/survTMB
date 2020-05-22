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
                           dim_g(4L),
                           dim_B = dim_g * dim_a,
                         n_nodes(30L);

    vector<ADd> omega(dim_o);
    omega << -0.74, 0.76, 0.20;

    vector<ADd> alpha(dim_a);
    alpha << 0.7, 0.6;

    vector<ADd> U(K),
                k(K);
    U << -0.79, -0.67, 0.05, -0.71, -0.13, -0.33;
    k << -0.155666233880193, -0.202626997173659, 0.480913001875679,
         -0.185234121879782, 0.138273358586317, -0.056526844705098;

    matrix<ADd> B(dim_g, dim_a),
           Lambda(K, K);
    B << 0.97, -0.78, 0.01, -0.86, -0.07, 0.98, 0.02, -0.03;
    Lambda <<  1.08, 0.12, -0.36, -0.48, 0.36, -0.12, 0.12, 0.36, 0, 0.12,
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

    arma::vec g_bk(2L);
    g_bk[0L] = 0.;
    g_bk[1L] = 10;
    arma::vec g_ik(2L);
    g_ik[0L] = 3. + 1. / 3.;
    g_ik[1L] = 6. + 2. / 3.;
    splines::ns g(g_bk, g_ik, true);

    arma::vec b_bk(2L);
    b_bk[0L] = 0.;
    b_bk[1L] = 2.30258509299405;
    arma::vec b_ik(2L);
    b_ik[0L] = 0.767528364331349;
    b_ik[1L] = 1.5350567286627;
    splines::ns b(b_bk, b_ik, false);

    snva_integral<double, splines::ns, splines::ns, splines::ns> func(
          "snva_integral", n_nodes, &b, &g, &m, dim_a);

    auto x = func.get_x(lb, ub, omega, alpha, B, U, k, Lambda);
    {
      CppAD::vector<ADd> y(1L);
      CppAD::Independent(x);
      func(x, y);
      CppAD::ADFun<double> afunc(x, y);
      afunc.optimize();

      constexpr double const intgral_val = 7.26194858630117;
      expect_equal(intgral_val, asDouble(y[0L]));

      CppAD::vector<double> xx(x.size());
      for(size_t i = 0; i < x.size(); ++i)
        xx[i] = asDouble(x[i]);

      auto yy = afunc.Forward(0, xx);
      expect_equal(intgral_val, yy[0L]);

      constexpr size_t n_grad_ele = 61L;
      constexpr double const grad[n_grad_ele] = {
         1.14724243199794, 2.95735564711031, 0.61440889004403, -0.683696240047505,
        -0.545434183517209, 1.32971805754415, 0.916240709055219, 1.53837667227965,
        0.137091513534909, 1.13975833503785, 0.785349179190188, 1.3186085762397,
        0.117507011601351, 0.939852768581435, 2.11722385520291, 0.209187744618701,
        0.805588087355515, 1.81476330445964, 0.179303781101744, 0.804763180771597,
        1.77128412626126, 0.0538008668494866, 0.68979701208994, 1.51824353679537,
        0.0461150287281313, 0.118257598027757, 0.125133891434628, 0.0134464198397455,
        0.101363655452363, 0.107257621229681, 0.0115255027197819, 0.125133891434628,
        0.313145617644883, 0.0095611312927022, 0.107257621229681, 0.2684105294099,
        0.00819525539374475, 0.0134464198397455, 0.0095611312927022,
        0.216884681710166, 0.0115255027197819, 0.00819525539374475, 0.185901155751571,
        0.101363655452363, 0.107257621229681, 0.0115255027197819, 0.0868831332448826,
        0.0919351039111553, 0.00987900233124163, 0.107257621229681, 0.2684105294099,
        0.00819525539374475, 0.0919351039111553, 0.230066168065629, 0.00702450462320979,
        0.0115255027197819, 0.00819525539374475, 0.185901155751571, 0.00987900233124163,
        0.00702450462320979, 0.159343847787061 };

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
        /* B */
        for(size_t j = 0; j < dim_B; ++j, ++i)
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
      expect_equal(2.57545288795293, yy[0L]);
      constexpr double const grad2[n_grad_ele] = {
        0.992292690709001, 1.06060442669076, -0.328799704440544, -0.209138559359073,
        -1.23315658694422, 0.854059238136823, 0.448065505222882, 0.392597786372208,
        -0.215301458456684, 0.732050775545848, 0.384056147333899, 0.336512388319035,
        -0.184544107248586, 0.635340452881667, 0.720719656265072, -0.255759215333374,
        0.544577531041429, 0.617759705370062, -0.219222184571464, 0.56552663022551,
        0.643315977120891, -0.231236549558922, 0.484737111621866, 0.55141369467505,
        -0.198202756764791, 0.0864009273176511, 0.0865190883641926, -0.0242352561257939,
        0.0740579377008438, 0.0741592185978794, -0.0207730766792519,
        0.0865190883641926, 0.101582765509422, -0.0383072038003378, 0.0741592185978794,
        0.0870709418652192, -0.0328347461145753, -0.0242352561257939,
        -0.0383072038003378, 0.0229951847286293, -0.0207730766792519,
        -0.0328347461145753, 0.0197101583388251, 0.0740579377008438,
        0.0741592185978794, -0.0207730766792519, 0.063478232315009, 0.063565044512468,
        -0.0178054942965017, 0.0741592185978794, 0.0870709418652192,
        -0.0328347461145753, 0.063565044512468, 0.0746322358844736, -0.0281440680982074,
        -0.0207730766792519, -0.0328347461145753, 0.0197101583388251,
        -0.0178054942965017, -0.0281440680982074, 0.0168944214332787 };
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
        /* B */
        for(size_t j = 0; j < dim_B; ++j, ++i)
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

  test_that("eval_snva_integral gives the correct result without b") {
    using ADd = AD<double>;
    constexpr size_t const dim_o(0L),
                           dim_a(2L),
                           dim_m(3L),
                               K = dim_a * dim_m,
                           dim_g(4L),
                           dim_B = dim_g * dim_a,
                         n_nodes(30L);

    vector<ADd> omega(dim_o);

    vector<ADd> alpha(dim_a);
    alpha << 0.7, 0.6;

    vector<ADd> U(K),
                k(K);
    U << -0.79, -0.67, 0.05, -0.71, -0.13, -0.33;
    k << -0.155666233880193, -0.202626997173659, 0.480913001875679,
         -0.185234121879782, 0.138273358586317, -0.056526844705098;

    matrix<ADd> B(dim_g, dim_a),
    Lambda(K, K);
    B << 0.97, -0.78, 0.01, -0.86, -0.07, 0.98, 0.02, -0.03;
    Lambda <<  1.08, 0.12, -0.36, -0.48, 0.36, -0.12, 0.12, 0.36, 0, 0.12,
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

    arma::vec g_bk(2L);
    g_bk[0L] = 0.;
    g_bk[1L] = 10;
    arma::vec g_ik(2L);
    g_ik[0L] = 3. + 1. / 3.;
    g_ik[1L] = 6. + 2. / 3.;
    splines::ns g(g_bk, g_ik, true);;

    snva_integral<double, splines::ns, splines::ns, splines::ns> func(
        "snva_integral", n_nodes, nullptr, &g, &m, dim_a);

    auto x = func.get_x(lb, ub, omega, alpha, B, U, k, Lambda);
    {
      CppAD::vector<ADd> y(1L);
      CppAD::Independent(x);
      func(x, y);
      CppAD::ADFun<double> afunc(x, y);
      afunc.optimize();

      constexpr double const intgral_val = 6.00523840539583;
      expect_equal(intgral_val, asDouble(y[0L]));

      CppAD::vector<double> xx(x.size());
      for(size_t i = 0; i < x.size(); ++i)
        xx[i] = asDouble(x[i]);

      auto yy = afunc.Forward(0, xx);
      expect_equal(intgral_val, yy[0L]);

      constexpr size_t n_grad_ele = 58L;
      constexpr double const grad[n_grad_ele] = {
        -0.615681785182354, -0.90689670151562, 1.20528124625157, 0.860834305998444,
        1.21360048030537, 0.0144390525796759, 1.03309821107277, 0.737857976570095,
        1.04022898311889, 0.0123763307825794, 0.908936699277026, 1.73357348825448,
        0.0709680923166048, 0.779088599380308, 1.48592013278956, 0.0608297934142327,
        0.784474922392819, 1.46443035153139, -0.0273961962840262, 0.672407076336702,
        1.25522601559833, -0.0234824539577367, 0.119235936396173, 0.120753983224687,
        0.00541547061720333, 0.10220223119672, 0.103503414192589, 0.00464183195760285,
        0.120753983224687, 0.253865283546395, -0.00490692209026202, 0.103503414192589,
        0.217598814468338, -0.00420593322022458, 0.00541547061720333,
        -0.00490692209026202, 0.154975676808193, 0.00464183195760286,
        -0.00420593322022458, 0.132836294407023, 0.10220223119672, 0.103503414192589,
        0.00464183195760286, 0.087601912454331, 0.0887172121650761, 0.00397871310651673,
        0.103503414192589, 0.217598814468338, -0.00420593322022458, 0.0887172121650761,
        0.18651326954429, -0.00360508561733536, 0.00464183195760285,
        -0.00420593322022458, 0.132836294407023, 0.00397871310651673,
        -0.00360508561733536, 0.113859680920305 };

      vector<double> w(1L);
      w[0L] = 1.;
      auto dx = afunc.Reverse(1, w);
      expect_true(dx.size() == n_grad_ele + 2L);

      {
        size_t i = 0L;
        /* omega  */
        /*for(size_t j = 0; j < dim_o; ++j, ++i)
          expect_equal(grad[i], dx[i + 2L]);*/
        /* alpha */
        for(size_t j = 0; j < dim_a; ++j, ++i)
          expect_equal(grad[i], dx[i + 2L]);
        /* B */
        for(size_t j = 0; j < dim_B; ++j, ++i)
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
    }
  }

  test_that("eval_snva_integral gives the correct result without g") {
    using ADd = AD<double>;
    constexpr size_t const dim_o(3L),
                           dim_a(2L),
                           dim_m(3L),
                               K = dim_a * dim_m,
                           dim_g(0L),
                         n_nodes(30L);

    vector<ADd> omega(dim_o);
    omega << -0.74, 0.76, 0.20;

    vector<ADd> alpha(dim_a);
    alpha << 0.7, 0.6;

    vector<ADd> U(K),
    k(K);
    U << -0.79, -0.67, 0.05, -0.71, -0.13, -0.33;
    k << -0.155666233880193, -0.202626997173659, 0.480913001875679,
         -0.185234121879782, 0.138273358586317, -0.056526844705098;

    matrix<ADd> B(dim_g, dim_a),
           Lambda(K, K);
    Lambda <<  1.08, 0.12, -0.36, -0.48, 0.36, -0.12, 0.12, 0.36, 0, 0.12,
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

    snva_integral<double, splines::ns, splines::ns, splines::ns> func(
        "snva_integral", n_nodes, &b, nullptr, &m, dim_a);

    auto x = func.get_x(lb, ub, omega, alpha, B, U, k, Lambda);
    {
      CppAD::vector<ADd> y(1L);
      CppAD::Independent(x);
      func(x, y);
      CppAD::ADFun<double> afunc(x, y);
      afunc.optimize();

      constexpr double const intgral_val = 6.48324263010285;
      expect_equal(intgral_val, asDouble(y[0L]));

      CppAD::vector<double> xx(x.size());
      for(size_t i = 0; i < x.size(); ++i)
        xx[i] = asDouble(x[i]);

      auto yy = afunc.Forward(0, xx);
      expect_equal(intgral_val, yy[0L]);

      constexpr size_t n_grad_ele = 53L;
      constexpr double const grad[n_grad_ele] = {
        1.16645901569117, 2.59408381254936, 0.78928104422283, -2.12101474056299,
        -0.28954490600768, 0.938203095875693, 1.85634472217966, 0.322929109101604,
        0.804174082179166, 1.59115261901114, 0.276796379229946, 0.800696769637818,
        1.54647968567511, 0.164804151310031, 0.686311516832415, 1.32555401629295,
        0.141260701122883, 0.11908952496286, 0.124184958582097, 0.0186722013297797,
        0.102076735682451, 0.106444250213226, 0.016004743996954, 0.124184958582097,
        0.269491981467769, 0.029567159181209, 0.106444250213226, 0.230993126972374,
        0.0253432792981792, 0.0186722013297797, 0.029567159181209, 0.188065042561797,
        0.016004743996954, 0.0253432792981792, 0.161198607910112, 0.102076735682451,
        0.106444250213226, 0.016004743996954, 0.0874943448706723, 0.0912379287541935,
        0.0137183519973892, 0.106444250213226, 0.230993126972374, 0.0253432792981792,
        0.0912379287541935, 0.197994108833463, 0.0217228108270107, 0.016004743996954,
        0.0253432792981792, 0.161198607910112, 0.0137183519973892, 0.0217228108270107,
        0.138170235351524 };

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
        /* B */
        /*for(size_t j = 0; j < dim_B; ++j, ++i)
          expect_equal(grad[i], dx[i + 2L]);*/
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
    }
  }

  test_that("eval_snva_integral gives the correct result without m") {
    using ADd = AD<double>;
    constexpr size_t const dim_o(3L),
                           dim_a(2L),
                           dim_m(0L),
                               K = dim_a * dim_m,
                           dim_g(4L),
                           dim_B = dim_g * dim_a,
                         n_nodes(30L);

    vector<ADd> omega(dim_o);
    omega << -0.74, 0.76, 0.20;

    vector<ADd> alpha(dim_a);
    alpha << 0.7, 0.6;

    vector<ADd> U(K),
                k(K);

    matrix<ADd> B(dim_g, dim_a),
    Lambda(K, K);
    B << 0.97, -0.78, 0.01, -0.86, -0.07, 0.98, 0.02, -0.03;

    ADd const lb(1.2), ub(9.5);

    arma::vec g_bk(2L);
    g_bk[0L] = 0.;
    g_bk[1L] = 10;
    arma::vec g_ik(2L);
    g_ik[0L] = 3. + 1. / 3.;
    g_ik[1L] = 6. + 2. / 3.;
    splines::ns g(g_bk, g_ik, true);

    arma::vec b_bk(2L);
    b_bk[0L] = 0.;
    b_bk[1L] = 2.30258509299405;
    arma::vec b_ik(2L);
    b_ik[0L] = 0.767528364331349;
    b_ik[1L] = 1.5350567286627;
    splines::ns b(b_bk, b_ik, false);

    snva_integral<double, splines::ns, splines::ns, splines::ns> func(
        "snva_integral", n_nodes, &b, &g, nullptr, dim_a);

    auto x = func.get_x(lb, ub, omega, alpha, B, U, k, Lambda);
    {
      CppAD::vector<ADd> y(1L);
      CppAD::Independent(x);
      func(x, y);
      CppAD::ADFun<double> afunc(x, y);
      afunc.optimize();

      constexpr double const intgral_val = 20.9588272700042;
      expect_equal(intgral_val, asDouble(y[0L]));

      CppAD::vector<double> xx(x.size());
      for(size_t i = 0; i < x.size(); ++i)
        xx[i] = asDouble(x[i]);

      auto yy = afunc.Forward(0, xx);
      expect_equal(intgral_val, yy[0L]);

      constexpr size_t n_grad_ele = 13L;
      constexpr double const grad[n_grad_ele] = {
        4.15687817968065, 8.47365360644718, 1.39583552811081, 5.40135344816136,
        -2.52633036440621, 4.16872715454956, 3.0847235343567, 4.22364619736634,
        0.104503613597825, 3.57319470389962, 2.64404874373432, 3.62026816917115,
        0.089574525940993 };

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
        /* B */
        for(size_t j = 0; j < dim_B; ++j, ++i)
          expect_equal(grad[i], dx[i + 2L]);
        /* U */
        /*for(size_t j = 0; j < K; ++j, ++i)
          expect_equal(grad[i], dx[i + 2L]);*/
        /* k */
        /*for(size_t j = 0; j < K; ++j, ++i)
          expect_equal(grad[i], dx[i + 2L]);*/
        /* Lambda */
        /*for(size_t j = 0; j < K * K; ++j, ++i)
          expect_equal(grad[i], dx[i + 2L]);*/
      }
    }
  }
}
