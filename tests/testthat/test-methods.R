test_that("the back-transform preserves the linear predictor exactly", {
  set.seed(31)
  n <- 50; p <- 4
  Xraw <- matrix(stats::rnorm(n * p, mean = c(0, 10, -5, 100),
                              sd = c(1, 0.01, 50, 3)), n, p, byrow = TRUE)
  x_mu <- colMeans(Xraw)
  x_sd <- apply(Xraw, 2L, stats::sd)
  Xstd <- sweep(sweep(Xraw, 2L, x_mu, "-"), 2L, x_sd, "/")

  B  <- matrix(stats::rnorm(7 * p), 7, p)   # draws on the standardized scale
  b0 <- stats::rnorm(7)

  bt <- bbqr:::.bbqr_backtransform(B, b0, list(x_mu = x_mu, x_sd = x_sd))

  # eta_std[d, i] and eta_raw[d, i] must agree for every draw and observation
  eta_std <- outer(b0, rep(1, n)) + B %*% t(Xstd)
  eta_raw <- outer(bt$beta0, rep(1, n)) + bt$beta %*% t(Xraw)
  expect_equal(eta_std, eta_raw)
})

test_that("a NULL scaler leaves the draws untouched", {
  B <- matrix(1:6, 3, 2); b0 <- c(1, 2, 3)
  bt <- bbqr:::.bbqr_backtransform(B, b0, NULL)
  expect_identical(bt$beta, B)
  expect_identical(bt$beta0, b0)
})

test_that("standardize rescues covariates on very different scales", {
  s <- sim_binary(n = 400, seed = 21)
  # x1 blown up 100-fold and x2 shrunk and shifted: the fixed priors are then
  # badly matched to the raw scale, which is exactly when standardizing pays.
  d <- s$data
  d$x1 <- d$x1 * 100
  d$x2 <- d$x2 / 50 + 7
  truth <- s$beta / c(100, 1 / 50, 1)   # truth in the rescaled units

  set.seed(3); raw <- coef(bbqr(s$formula, d, ndraw = 2500, burn = 500))
  set.seed(3); std <- coef(bbqr(s$formula, d, ndraw = 2500, burn = 500,
                                standardize = TRUE))

  expect_gt(cos_sim(std[-1], truth), 0.99)
  expect_gt(cos_sim(std[-1], truth), cos_sim(raw[-1], truth))
})

test_that("standardize is near-neutral when covariates share a scale", {
  s <- sim_binary(n = 400, seed = 25)
  set.seed(4); raw <- coef(bbqr(s$formula, s$data, ndraw = 2500, burn = 500))
  set.seed(4); std <- coef(bbqr(s$formula, s$data, ndraw = 2500, burn = 500,
                                standardize = TRUE))
  expect_gt(cos_sim(raw[-1], std[-1]), 0.98)
})

test_that("summary() reports intervals and a selection flag", {
  s <- sim_binary(n = 400, seed = 22)
  sm <- summary(bbqr(s$formula, s$data, ndraw = 2000, burn = 400))
  expect_s3_class(sm, "bbqr.summary")
  tab <- sm$coefficients
  expect_equal(nrow(tab), 4L)
  expect_true(all(tab$lower <= tab$mean))
  expect_true(all(tab$mean <= tab$upper))
  expect_type(tab$selected, "logical")
  # x1 and x2 are genuinely non-zero and should be picked up at n = 400
  expect_true(tab["x1", "selected"])
  expect_true(tab["x2", "selected"])
  expect_equal(sm$level, 0.95)

  sm90 <- summary(bbqr(s$formula, s$data, ndraw = 1000, burn = 200),
                  level = 0.90)
  expect_equal(sm90$level, 0.90)
  expect_error(summary(bbqr(s$formula, s$data, ndraw = 300, burn = 50),
                       level = 1.5))
})

test_that("predict() returns valid probabilities, indices and labels", {
  s <- sim_binary(n = 300, seed = 23)
  fit <- bbqr(s$formula, s$data, ndraw = 1500, burn = 300)

  p <- predict(fit, newdata = s$data, type = "response")
  expect_length(p, 300L)
  expect_true(all(p >= 0 & p <= 1))

  eta <- predict(fit, newdata = s$data, type = "link")
  expect_true(all(is.finite(eta)))

  cl <- predict(fit, newdata = s$data, type = "class")
  expect_true(all(cl %in% c(0L, 1L)))
  expect_identical(cl, as.integer(eta > 0))
  # should beat guessing the majority class
  expect_gt(mean(cl == s$data$y), max(table(s$data$y)) / 300)

  expect_error(predict(fit), "newdata")
  expect_error(predict(fit, newdata = s$data[, c("y", "x1")]))
})

test_that("predict(type = 'response') matches the ALD survival function", {
  s <- sim_binary(n = 200, seed = 24)
  fit <- bbqr(s$formula, s$data, quantile = 0.3, ndraw = 800, burn = 200)
  eta <- predict(fit, newdata = s$data, type = "link")
  p   <- predict(fit, newdata = s$data, type = "response")
  tau <- 0.3
  sig <- mean(fit$sigma)
  ref <- ifelse(eta > 0,
                1 - tau * exp(-(1 - tau) * eta / sig),
                (1 - tau) * exp(tau * eta / sig))
  expect_equal(p, ref)
  # the two branches agree at eta = 0
  expect_equal(1 - tau * exp(0), (1 - tau) * exp(0))
})

test_that("multiple quantiles give a bbqr.list with one fit per level", {
  s <- sim_binary()
  fit <- bbqr(s$formula, s$data, quantile = c(0.25, 0.5, 0.75),
              ndraw = 500, burn = 100)
  expect_s3_class(fit, "bbqr.list")
  expect_length(fit, 3L)
  expect_equal(vapply(fit, function(z) z$quantile, numeric(1)),
               c(0.25, 0.5, 0.75), ignore_attr = TRUE)
  expect_equal(dim(coef(fit)), c(4L, 3L))
  expect_length(summary(fit), 3L)
})

test_that("print methods run without error", {
  s <- sim_binary()
  fit  <- bbqr(s$formula, s$data, ndraw = 400, burn = 100)
  mfit <- bbqr(s$formula, s$data, quantile = c(0.4, 0.6),
               ndraw = 400, burn = 100)
  expect_output(print(fit), "Bayesian binary quantile regression")
  expect_output(print(summary(fit)), "credible intervals")
  expect_output(print(mfit), "2 quantiles")
  expect_output(print(prior("alasso")), "alasso")
})

test_that("plot() runs for each type", {
  s <- sim_binary()
  fit <- bbqr(s$formula, s$data, ndraw = 400, burn = 100)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_invisible(plot(fit, which = "x1"))
  expect_silent(plot(fit, which = c("x1", "x2"), type = "trace"))
  expect_silent(plot(fit, which = 1, type = "density"))
  expect_error(plot(fit, which = "nope"))
})
