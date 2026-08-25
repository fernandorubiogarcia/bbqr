test_that("adaptive-lasso defaults implement the V5 appendix hierarchy", {
  pr <- prior("alasso")
  expect_identical(pr$model, "v5")
  expect_equal(pr$alpha_delta, 2)
  expect_equal(pr$beta_delta, 2)
  expect_equal(c(pr$alpha_omega, pr$beta_omega), c(2, 2))
  expect_true(is.na(pr$b0_var))

  s <- sim_binary()
  fit <- bbqr_alasso(s$formula, s$data, ndraw = 600, burn = 100)
  expect_identical(fit$model, "v5")
  expect_true(fit$derived)
  expect_equal(fit$anchor, "sigma1")
  expect_true(all(fit$sigma == 1))
  expect_true(all(is.finite(fit$delta) & fit$delta > 0))
  expect_true(all(is.finite(fit$omega) & fit$omega > 0))
  expect_true(all(is.finite(fit$lambdasq) & fit$lambdasq > 0))
  ## The intercept prior is the tau-calibrated reference: (0, 8) at the median.
  expect_equal(unname(fit$b0_prior[c("mean", "var")]), c(0, 8))

  ## The appendix also covers a sampled sigma, one argument away.
  free <- bbqr_alasso(s$formula, s$data, ndraw = 600, burn = 100,
                      anchor = "free")
  expect_gt(stats::sd(free$sigma), 0)
  expect_true(free$derived)
})

test_that("adaptive-lasso rejects an improper prior on delta", {
  ## beta_delta = 0 is the flat prior, under which the joint posterior is not
  ## integrable in delta. It is admissible for the single-penalty lasso only,
  ## and there only under the published hierarchy: "v5" means every component
  ## is proper.
  expect_error(prior("alasso", beta_delta = 0), "improper")
  expect_silent(prior("lasso", model = "v3", beta_delta = 0))
  expect_error(prior("lasso", beta_delta = 0), "strictly positive")

  ## A saved prior object edited by hand must be caught at fit time too.
  bad <- prior("alasso")
  bad$beta_delta <- 0
  s <- sim_binary()
  expect_error(bbqr_alasso(s$formula, s$data, ndraw = 200, prior = bad),
               "improper")
})

test_that("the delta hyperprior reaches the sampler", {
  ## A rate large enough to dominate the likelihood must pull delta down
  ## relative to a diffuse one. If alpha_delta/beta_delta were dropped on the
  ## way to the kernel, these two chains would agree.
  s <- sim_binary()
  tight <- bbqr_alasso(s$formula, s$data, ndraw = 3000, burn = 600,
                       prior = prior("alasso", alpha_delta = 2, beta_delta = 50))
  loose <- bbqr_alasso(s$formula, s$data, ndraw = 3000, burn = 600,
                       prior = prior("alasso", alpha_delta = 2, beta_delta = 2))
  expect_lt(mean(tight$delta), mean(loose$delta))
  expect_lt(mean(tight$delta), 0.5)
})

test_that("delta is chain-length stable under a proper prior", {
  ## Under the flat prior delta drifts without bound, so the posterior mean
  ## grows with the number of draws. A proper prior must not do that.
  s <- sim_binary()
  short <- bbqr_alasso(s$formula, s$data, ndraw = 2000,  burn = 400)
  long  <- bbqr_alasso(s$formula, s$data, ndraw = 16000, burn = 3200)
  expect_lt(abs(mean(long$delta) - mean(short$delta)),
            2 * max(mean(short$delta), mean(long$delta)))
  expect_lt(mean(long$delta), 20)
})
