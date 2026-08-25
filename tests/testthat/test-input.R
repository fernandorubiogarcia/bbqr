test_that("invalid responses are rejected", {
  s <- sim_binary()
  d <- s$data

  # continuous response
  expect_error(bbqr(x1 ~ x2 + x3, data = d, ndraw = 200), "binary")

  # 0/1/2 coding
  d2 <- d; d2$y[1:5] <- 2
  expect_error(bbqr(s$formula, d2, ndraw = 200), "binary")

  # constant response
  d3 <- d; d3$y <- 1
  expect_error(bbqr(s$formula, d3, ndraw = 200), "only one value")

  # three-level factor
  d4 <- d; d4$y <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  expect_error(bbqr(s$formula, d4, ndraw = 200), "two levels")

})

test_that("rows with missing values are dropped, as model.frame does", {
  s <- sim_binary()
  d <- s$data
  n0 <- nrow(d)
  d$x1[1:3] <- NA
  fit <- bbqr(s$formula, d, ndraw = 300, burn = 50)
  expect_equal(fit$n, n0 - 3L)
})

test_that("MCMC control arguments are validated", {
  s <- sim_binary()
  expect_error(bbqr(s$formula, s$data, quantile = 0, ndraw = 200), "between")
  expect_error(bbqr(s$formula, s$data, quantile = 1, ndraw = 200), "between")
  expect_error(bbqr(s$formula, s$data, ndraw = NULL), "ndraw")
  expect_error(bbqr(s$formula, s$data, ndraw = 100, keep = 3), "divisible")
  expect_error(bbqr(s$formula, s$data, ndraw = 100, keep = 200), "larger")
  expect_error(bbqr(s$formula, s$data, ndraw = 100, burn = 100), "burn")
  expect_error(bbqr(s$formula, s$data, ndraw = 100, burn = -1), "burn")
})

test_that("beta_init is matched by name and length-checked", {
  s <- sim_binary()
  # named, deliberately out of order
  init <- c(x3 = 0.1, x1 = 2, x2 = -2)
  fit <- bbqr(s$formula, s$data, ndraw = 400, burn = 100, beta_init = init)
  expect_s3_class(fit, "bbqr")

  expect_error(bbqr(s$formula, s$data, ndraw = 400, beta_init = c(1, 2)),
               "length")
})

test_that("prior() validates and reports method-specific defaults", {
  # The package default is the fully proper "v5" hierarchy for every penalty
  for (pen in c("alasso", "lasso", "none"))
    expect_identical(prior(pen)$model, "v5", info = pen)
  # under which the sigma prior is the mean-one reference for both penalties
  # that sample sigma ...
  expect_equal(c(prior("alasso")$a, prior("alasso")$b), c(2, 2))
  expect_equal(c(prior("lasso")$a, prior("lasso")$b), c(2, 2))
  # ... and the lasso's delta prior is proper.
  expect_equal(prior("lasso")$alpha_delta, 2)
  expect_equal(prior("lasso")$beta_delta, 2)

  # The documented values of the source papers are reached with model = "v3"
  expect_equal(prior("alasso", model = "v3")$a, 1)
  expect_equal(prior("lasso", model = "v3")$a, 0.1)
  # (1, 0) is the improper flat prior p(delta) ~ 1 of the 2013 paper
  expect_equal(prior("lasso", model = "v3")$alpha_delta, 1)
  expect_equal(prior("lasso", model = "v3")$beta_delta, 0)
  # The adaptive lasso needs a proper prior on delta under every model: with k
  # penalties the marginal posterior of delta flattens, so a flat prior is not
  # integrable.
  expect_equal(prior("alasso")$alpha_delta, 2)
  expect_equal(prior("alasso")$beta_delta, 2)
  expect_equal(prior("alasso", model = "v3")$beta_delta, 2)
  expect_equal(prior("none")$beta_var, 100)

  expect_error(prior("alasso", a = -1), "positive")
  # beta_delta may be exactly zero for the single-penalty lasso only, and
  # only under the published hierarchy
  expect_equal(prior("lasso", model = "v3", beta_delta = 0)$beta_delta, 0)
  expect_error(prior("lasso", beta_delta = 0), "strictly positive")
  expect_error(prior("lasso", beta_delta = -1), "non-negative")
  expect_equal(prior("alasso", alpha_delta = 2, beta_delta = 2)$beta_delta, 2)
  expect_error(prior("alasso", alpha_delta = 1, beta_delta = 0), "improper")
  expect_error(prior("alasso", a = 0), "positive")
  expect_error(prior("alasso", target_acc = 1.2), "between")
  expect_error(prior("nope"))

  expect_s3_class(prior("lasso"), "bbqr.prior")
})

test_that("a custom prior changes the fit", {
  s <- sim_binary()
  set.seed(8)
  a <- bbqr(s$formula, s$data, penalty = "none", ndraw = 600, burn = 100,
            prior = prior("none", beta_var = 0.01))
  set.seed(8)
  b <- bbqr(s$formula, s$data, penalty = "none", ndraw = 600, burn = 100,
            prior = prior("none", beta_var = 1000))
  # a much tighter prior must shrink the slopes
  expect_lt(sum(abs(coef(a)[-1])), sum(abs(coef(b)[-1])))
})

test_that("a formula without an intercept still fits", {
  s <- sim_binary()
  fit <- bbqr(y ~ 0 + x1 + x2 + x3, data = s$data, ndraw = 400, burn = 100)
  expect_equal(fit$nvar, 3L)
})
