test_that("anchor = 'sigma1' holds the scale at sigma_fixed", {
  s <- sim_binary()
  fit <- bbqr(s$formula, data = s$data, anchor = "sigma1",
              ndraw = 600, burn = 100)
  expect_true(all(fit$sigma == 1))

  fit2 <- bbqr(s$formula, data = s$data, anchor = "sigma1", sigma_fixed = 4,
               ndraw = 600, burn = 100)
  expect_true(all(fit2$sigma == 4))
})

test_that("anchor = 'free' actually samples the scale", {
  s <- sim_binary()
  fit <- bbqr(s$formula, data = s$data, anchor = "free",
              ndraw = 600, burn = 100)
  expect_gt(stats::sd(fit$sigma), 0)
})

## The three projection anchors rescale the coefficient vector after it is
## drawn. Under the default model = "v5" that is a departure from the derived
## posterior, so the fit warns and records derived = FALSE; the geometry the
## tests below check holds regardless, and the warning is asserted rather than
## silenced so that a fit which stopped warning would be noticed.
test_that("anchor = 'beta1' pins the first slope at exactly 1", {
  s <- sim_binary()
  expect_warning(fit <- bbqr(s$formula, data = s$data, anchor = "beta1",
                             ndraw = 600, burn = 100),
                 "derived posterior")
  expect_true(all(fit$beta[, "x1"] == 1))
  expect_false(fit$derived)
})

test_that("anchor = 'norm1' puts (beta0, beta) on the unit sphere", {
  s <- sim_binary()
  expect_warning(fit <- bbqr(s$formula, data = s$data, anchor = "norm1",
                             ndraw = 600, burn = 100),
                 "derived posterior")
  nrm <- sqrt(fit$beta0^2 + rowSums(fit$beta^2))
  expect_true(all(abs(nrm - 1) < 1e-8))
})

test_that("anchor = 'normslopes' normalises the slopes, not the full vector", {
  s <- sim_binary()
  expect_warning(fit <- bbqr(s$formula, data = s$data, anchor = "normslopes",
                             ndraw = 600, burn = 100),
                 "derived posterior")
  expect_true(all(abs(sqrt(rowSums(fit$beta^2)) - 1) < 1e-8))
  # the intercept is carried along by the rescaling, not pinned, so the joint
  # norm must generally differ from 1
  joint <- sqrt(fit$beta0^2 + rowSums(fit$beta^2))
  expect_false(all(abs(joint - 1) < 1e-8))
})

test_that("the two norm anchors are genuinely different constraints", {
  s <- sim_binary()
  ## Under the published hierarchy the projections are not flagged, which is
  ## the historical contract and keeps this comparison warning-free.
  set.seed(2); a <- bbqr(s$formula, s$data, anchor = "norm1",
                         prior = v3_prior(), ndraw = 800, burn = 200)
  set.seed(2); b <- bbqr(s$formula, s$data, anchor = "normslopes",
                         prior = v3_prior(), ndraw = 800, burn = 200)
  expect_false(isTRUE(all.equal(coef(a), coef(b))))
  expect_true(a$derived && b$derived)
})

test_that("the default anchor is sigma = 1 for every penalty", {
  s <- sim_binary()
  ## The likelihood identifies beta only up to a positive scale, so something
  ## has to close the ray. Holding sigma at 1 does it outright: the coefficient
  ## scale is then stated rather than left to the prior, the intervals do not
  ## inherit the width of an unidentified direction, and the tau-calibrated
  ## intercept prior of model = "v5" is exact there. "free" -- sigma sampled
  ## from its full conditional, as the adaptive-lasso and lasso derivations
  ## also allow -- stays one argument away.
  for (pen in c("alasso", "lasso", "none")) {
    fit <- bbqr(s$formula, s$data, penalty = pen, ndraw = 300, burn = 50)
    expect_equal(fit$anchor, "sigma1", info = pen)
    expect_true(all(fit$sigma == 1), info = pen)
    expect_true(fit$derived, info = pen)
  }
  ## The low-level fitters agree with the wrapper on the default.
  expect_equal(bbqr_alasso(s$formula, s$data, ndraw = 300, burn = 50)$anchor,
               "sigma1")
  expect_equal(bbqr_lasso(s$formula, s$data, ndraw = 300, burn = 50)$anchor,
               "sigma1")
})

test_that("an unknown anchor name is rejected", {
  s <- sim_binary()
  expect_error(bbqr(s$formula, s$data, anchor = "nonsense", ndraw = 300))
  expect_error(bbqr_lasso(s$formula, s$data, anchor = "nonsense", ndraw = 300))
})

test_that("penalty = 'none' has no shrinkage layer, the others do", {
  s <- sim_binary()
  none <- bbqr(s$formula, s$data, penalty = "none", ndraw = 400, burn = 100)
  las  <- bbqr(s$formula, s$data, penalty = "lasso", ndraw = 400, burn = 100)
  ala  <- bbqr(s$formula, s$data, penalty = "alasso", ndraw = 400, burn = 100)

  expect_null(none$lambdasq)
  expect_true(is.na(none$accept))

  # one shared penalty for the lasso, one per slope for the adaptive lasso
  expect_true(is.numeric(las$lambdasq) && is.null(dim(las$lambdasq)))
  expect_equal(dim(ala$lambdasq), c(300L, 3L))
})
