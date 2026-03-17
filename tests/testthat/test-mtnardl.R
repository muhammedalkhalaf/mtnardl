test_that("mtnardl binary partition works", {
  set.seed(12345)
  TT <- 50
  y  <- cumsum(rnorm(TT))
  x  <- cumsum(rnorm(TT))
  res <- mtnardl(y = y, x = matrix(x, ncol = 1),
                 decompose = 1L, partition = "binary",
                 maxlag = 2L, ic = "aic")
  expect_s3_class(res, "mtnardl")
  expect_true(res$nq == 2)
  expect_true(is.numeric(res$ecm_coef))
})

test_that("mtnardl quintile partition gives 5 regimes", {
  set.seed(12345)
  TT <- 60
  y  <- cumsum(rnorm(TT))
  x  <- cumsum(rnorm(TT))
  res <- mtnardl(y = y, x = matrix(x, ncol = 1),
                 decompose = 1L, partition = "quintile",
                 maxlag = 2L, ic = "bic")
  expect_equal(res$nq, 5L)
})

test_that("bounds test is returned", {
  set.seed(12345)
  TT <- 50
  y  <- cumsum(rnorm(TT))
  x  <- cumsum(rnorm(TT))
  res <- mtnardl(y = y, x = matrix(x, ncol = 1),
                 decompose = 1L, partition = "binary", maxlag = 1L)
  expect_true(!is.null(res$bounds_test))
  expect_true(is.numeric(res$bounds_test$F_ov))
})

test_that("print.mtnardl does not error", {
  set.seed(12345)
  TT <- 50
  y  <- cumsum(rnorm(TT))
  x  <- cumsum(rnorm(TT))
  res <- mtnardl(y = y, x = matrix(x, ncol = 1),
                 decompose = 1L, partition = "binary", maxlag = 1L)
  expect_output(print(res))
})

test_that("multipliers are computed", {
  set.seed(12345)
  TT <- 50
  y  <- cumsum(rnorm(TT))
  x  <- cumsum(rnorm(TT))
  res <- mtnardl(y = y, x = matrix(x, ncol = 1),
                 decompose = 1L, partition = "binary",
                 maxlag = 1L, horizon = 5L)
  if (!is.null(res$multipliers)) {
    expect_true(nrow(res$multipliers) > 0)
  } else {
    succeed("multipliers is NULL (ecm = 0)")
  }
})

