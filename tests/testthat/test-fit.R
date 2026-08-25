test_that("all three penalties run and return a well-formed object", {
  s <- sim_binary()
  for (pen in c("none", "lasso", "alasso")) {
    fit <- bbqr(s$formula, data = s$data, penalty = pen,
                ndraw = 600, burn = 100)
    expect_s3_class(fit, "bbqr")
    expect_equal(fit$penalty, pen)
    expect_equal(dim(fit$beta), c(500L, 3L))
    expect_length(fit$beta0, 500L)
    expect_equal(colnames(fit$beta), c("x1", "x2", "x3"))
    expect_true(all(is.finite(fit$beta)))
    expect_true(all(is.finite(fit$beta0)))
    expect_length(coef(fit), 4L)
    expect_equal(names(coef(fit))[1], "(Intercept)")
  }
})

test_that("thinning and burn-in produce the right number of draws", {
  s <- sim_binary()
  fit <- bbqr(s$formula, data = s$data, ndraw = 1000, burn = 200, keep = 5)
  # 1000/5 = 200 rows returned by Fortran, less floor(200/5) = 40 burn rows
  expect_equal(nrow(fit$beta), 160L)
  expect_equal(fit$ndraw_kept, 160L)
})

test_that("each sampler recovers the direction of the true coefficients", {
  s <- sim_binary(n = 400, seed = 11)
  for (pen in c("none", "lasso", "alasso")) {
    fit <- bbqr(s$formula, data = s$data, penalty = pen,
                ndraw = 3000, burn = 500)
    expect_gt(cos_sim(coef(fit)[-1], s$beta), 0.9)
  }
})

test_that("the sign of each non-zero coefficient is recovered", {
  s <- sim_binary(n = 400, seed = 12)
  fit <- bbqr(s$formula, data = s$data, ndraw = 3000, burn = 500)
  cf <- coef(fit)[-1]
  expect_gt(cf[["x1"]], 0)
  expect_lt(cf[["x2"]], 0)
})

test_that("results are reproducible under set.seed", {
  s <- sim_binary()
  set.seed(123); a <- bbqr(s$formula, data = s$data, ndraw = 500, burn = 100)
  set.seed(123); b <- bbqr(s$formula, data = s$data, ndraw = 500, burn = 100)
  expect_identical(a$beta, b$beta)
  expect_identical(a$beta0, b$beta0)
})

test_that("the low-level fitters agree with the bbqr() wrapper", {
  s <- sim_binary()
  set.seed(5); a <- bbqr(s$formula, data = s$data, penalty = "alasso",
                         ndraw = 500, burn = 100)
  set.seed(5); b <- bbqr_alasso(s$formula, data = s$data,
                                ndraw = 500, burn = 100)
  expect_equal(a$beta, b$beta)

  set.seed(5); c1 <- bbqr(s$formula, data = s$data, penalty = "none",
                          ndraw = 500, burn = 100)
  set.seed(5); c2 <- bbqr_none(s$formula, data = s$data,
                               ndraw = 500, burn = 100)
  expect_equal(c1$beta, c2$beta)
})

test_that("a two-level factor response is accepted", {
  s <- sim_binary()
  d <- s$data
  d$y <- factor(ifelse(d$y == 1, "yes", "no"), levels = c("no", "yes"))
  fit <- bbqr(s$formula, data = d, ndraw = 400, burn = 100)
  expect_s3_class(fit, "bbqr")
  expect_equal(fit$ylevels, c("no", "yes"))
})
