#include "testthat-wrap.h"
#include "orth_poly.h"

context("testing orth_poly") {
  test_that("orth_poly gives the same as poly in R") {
    /*
     set.seed(1)
     dput(x <- round(rnorm(6), 2))
     obj <- poly(x, degree = 3)
     attr(obj, "coefs")$norm2 <- attr(obj, "coefs")$norm2 / NROW(obj)
     dput(attr(obj, "coefs"))
     dput(cbind(1, predict(obj, x)))
    */
    arma::vec const x = { -0.63, 0.18, -0.84, 1.6, 0.33, -0.82 },
                alpha = { -0.03, 0.673718350183412, 0.388455148829439 },
                norm2 = { 0.166666666666667, 1, 0.745133333333333, 0.413676179913314,
                          0.0235062965408605 };
    arma::mat basis =
      { { 1, 1, 1, 1, 1, 1, -0.695079138943395, 0.243277698630188,
          -0.938356837573584, 1.88829832746289, 0.417047483366037, -0.915187532942137,
          0.0576788318055649, -1.31972174466434, 0.747817172463846, 1.18895139819447,
          -1.35090719440005, 0.676181536600509, 1.92619780939021, 0.393652338460405,
          -0.920783200334851, 0.14053903106972, -0.97230438386083, -0.567301594724647 } };
    basis.reshape(6L, 4L);

    arma::mat Xout;
    poly::orth_poly const obj =
      poly::orth_poly::get_poly_basis(x, 3L, Xout);

    expect_true(obj.alpha.n_elem == alpha.n_elem);
    for(size_t i = 0; i < alpha.n_elem; ++i)
      expect_equal(obj.alpha[i], alpha[i]);

    expect_true(obj.norm2.n_elem == norm2.n_elem);
    for(size_t i = 0; i < norm2.n_elem; ++i)
      expect_equal(obj.norm2[i], norm2[i]);

    expect_true(basis.n_cols == Xout.n_cols);
    expect_true(basis.n_rows == Xout.n_rows);
    for(size_t j = 0; j < Xout.n_cols; ++j)
      for(size_t i = 0; i < Xout.n_rows; ++i)
        expect_equal(Xout.at(i, j), basis.at(i, j));

    for(size_t i = 0; i < Xout.n_rows; ++i){
      arma::vec const b = obj(x[i]);
      expect_true(b.n_elem == Xout.n_cols);
      for(size_t j = 0; j < Xout.n_cols; ++j)
        expect_equal(Xout.at(i, j), b[j]);
    }
  }
}
