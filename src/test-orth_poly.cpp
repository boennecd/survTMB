#include "testthat-wrap.h"
#include "orth_poly.h"

context("testing orth_poly") {
  test_that("orth_poly gives the same as poly in R") {
    /*
     set.seed(1)
     dput(x <- round(rnorm(6), 2))
     obj <- poly(x, degree = 3)
     dput(attr(obj, "coefs"))
     dput(unclass(obj))
    */
    arma::vec const x = { -0.63, 0.18, -0.84, 1.6, 0.33, -0.82 },
                alpha = { -0.03, 0.673718350183412, 0.388455148829439 },
                norm2 = { 1, 6, 4.4708, 2.48205707947988, 0.141037779245163 };
    arma::mat basis =
      { { -0.283764870210735, 0.0993177045737573, -0.383082574784492,
          0.770894564072497, 0.170258922126441, -0.373623745777468, 0.0235472844805747,
          -0.538774146147203, 0.305295082237885, 0.485387375757513, -0.551505552689155,
          0.276049956360386, 0.78636696277879, 0.160707894213561, -0.375908167424547,
          0.0573748191776599, -0.396941602521703, -0.231599906223761 } };
    basis.reshape(6L, 3L);

    arma::mat Xout;
    poly::orth_poly const obj =
      poly::orth_poly::get_poly_basis(x, 3L, Xout);

    expect_true(obj.alpha.n_elem == alpha.n_elem);
    for(size_t i = 0; i < alpha.n_elem; ++i)
      expect_equal(obj.alpha[i], alpha[i]);

    expect_true(obj.norm2.n_elem == norm2.n_elem);
    for(size_t i = 0; i < norm2.n_elem; ++i)
      expect_equal(obj.norm2[i], norm2[i]);

    expect_true(basis.n_cols + 1L == Xout.n_cols);
    expect_true(basis.n_rows == Xout.n_rows);
    for(size_t j = 1; j < Xout.n_cols; ++j)
      for(size_t i = 0; i < Xout.n_rows; ++i)
        expect_equal(Xout.at(i, j), basis.at(i, j - 1));

    for(size_t i = 0; i < Xout.n_rows; ++i){
      arma::vec const b = obj(x[i]);
      expect_true(b.n_elem == Xout.n_cols);
      for(size_t j = 0; j < Xout.n_cols; ++j)
        expect_equal(Xout.at(i, j), b[j]);
    }
  }
}
