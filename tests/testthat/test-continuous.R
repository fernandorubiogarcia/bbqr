## Continuous-response kernels (appendix V5C).
##
## The load-bearing checks are (1) that the sampler recovers the same quantile
## surface as a frequentist fit, and (2) that the y-standardization round-trip
## is exact -- that transformation is new code with no binary counterpart, and
## a sign or scale error in it would be invisible in every other test.

sim_cont <- function(n = 400, seed = 7, y_shift = 40, y_scale = 5) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y <- y_shift + 6 * d$x1 - 4.5 * d$x2 + rnorm(n, 0, y_scale)
  d
}

test_that("the unpenalised kernel recovers the quantile surface", {
  skip_on_cran()
  skip_if_not_installed("quantreg")
  d <- sim_cont()
  for (tau in c(0.25, 0.5, 0.75)) {
    set.seed(1)
    f <- cbqr(y ~ x1 + x2 + x3, d, penalty = "none", quantile = tau,
              ndraw = 6000, keep = 2)
    est <- c(mean(f$beta0), colMeans(f$beta))
    ## The reference is the frequentist estimator on the same data, not the
    ## data-generating truth: at n = 400 with residual sd 5 the two differ by
    ## more than the sampler does.
    ref <- unname(stats::coef(quantreg::rq(y ~ x1 + x2 + x3, tau = tau, data = d)))
    expect_lt(max(abs(est - ref)), 0.6)
  }
})

test_that("y standardization round-trips exactly", {
  skip_on_cran()
  d <- sim_cont()
  set.seed(11)
  fT <- cbqr(y ~ x1 + x2 + x3, d, penalty = "none", ndraw = 6000, keep = 2,
             standardize = TRUE)
  set.seed(11)
  fF <- cbqr(y ~ x1 + x2 + x3, d, penalty = "none", ndraw = 6000, keep = 2,
             standardize = FALSE)
  expect_equal(c(mean(fT$beta0), colMeans(fT$beta)),
               c(mean(fF$beta0), colMeans(fF$beta)),
               tolerance = 0.05, ignore_attr = TRUE)
  ## sigma is an INVERSE scale, so it must divide by y_sd on the way back. If
  ## it multiplied instead, this comparison would be off by y_sd^2.
  expect_equal(mean(fT$sigma), mean(fF$sigma), tolerance = 0.05)
  expect_true(fT$standardized$y)
  expect_false(fF$standardized$y)
  expect_equal(fT$standardized$y_sd, stats::sd(d$y))
})

test_that("all three penalties run and agree on a well-identified design", {
  skip_on_cran()
  d <- sim_cont()
  fits <- lapply(c("none", "lasso", "alasso"), function(pen) {
    set.seed(3)
    cbqr(y ~ x1 + x2 + x3, d, penalty = pen, ndraw = 6000, keep = 2)
  })
  est <- vapply(fits, function(f) colMeans(f$beta), numeric(3))
  expect_lt(max(abs(est[, 1] - est[, 2])), 0.3)
  expect_lt(max(abs(est[, 1] - est[, 3])), 0.3)
  expect_true(all(vapply(fits, function(f) all(is.finite(f$sigma)), logical(1))))
})

test_that("the layer defaults for q are the published ones, and differ", {
  ## Not a style choice: the published lasso penalty is sigma-free (q = 0) and
  ## the published adaptive penalty is a half power (q = 1). See Remark 5C.
  expect_equal(bbqr:::.cbqr_default_q("lasso"), 0)
  expect_equal(bbqr:::.cbqr_default_q("alasso"), 1)
  skip_on_cran()
  d <- sim_cont()
  set.seed(5); fl <- cbqr(y ~ x1 + x2, d, penalty = "lasso",  ndraw = 3000, keep = 2)
  set.seed(5); fa <- cbqr(y ~ x1 + x2, d, penalty = "alasso", ndraw = 3000, keep = 2)
  expect_equal(fl$q, 0)
  expect_equal(fa$q, 1)
})

test_that("q = 2 takes the slice branch and still mixes", {
  skip_on_cran()
  d <- sim_cont()
  set.seed(9)
  f <- cbqr(y ~ x1 + x2 + x3, d, penalty = "alasso", q = 2,
            ndraw = 6000, keep = 2)
  expect_true(all(is.finite(f$sigma)))
  ## A stuck slice sampler would return long runs of the identical value; the
  ## exact-Gamma branch at q = 1 gives ~100% unique draws, so anything below
  ## half is a failure rather than a tuning problem.
  expect_gt(length(unique(f$sigma)) / length(f$sigma), 0.5)
})

test_that("the continuous clamp-site contract is 16 names and is not the binary one", {
  expect_length(bbqr:::.CBQR_CLIP_SITES, 16L)
  expect_length(bbqr:::.BBQR_CLIP_SITES, 18L)
  ## The binary vector carries two sites that guard the deleted y* draw.
  expect_setequal(setdiff(bbqr:::.BBQR_CLIP_SITES, bbqr:::.CBQR_CLIP_SITES),
                  c("ystar", "v4_sd"))
  ## Feeding the binary names to the continuous resolver must fail loudly, not
  ## silently misalign every index.
  expect_error(
    bbqr:::.cbqr_clip_vector(stats::setNames(rep(TRUE, 16), bbqr:::.BBQR_CLIP_SITES[1:16])),
    "DIFFERENT contract")
})

test_that("continuous entry points reject binary and malformed input", {
  d <- sim_cont(n = 60)
  expect_error(cbqr(I(y > 40) ~ x1, d, ndraw = 100), "must be numeric and continuous")
  expect_error(cbqr(y ~ x1, d, penalty = "none", q = 2, ndraw = 100), "no effect")
  expect_error(cbqr(y ~ x1, d, penalty = "lasso", clip = TRUE, ndraw = 100),
               "only available for penalty")
  expect_error(cbqr(y ~ x1, d, ndraw = 100, prior = list(nonsense = 1)),
               "Unknown prior component")
  expect_error(cbqr(y ~ x1, d, ndraw = 100, prior = list(a = -1)), "positive")
  d2 <- d; d2$y <- rep(c(0, 1), length.out = nrow(d))
  expect_error(cbqr(y ~ x1, d2, ndraw = 100), "fewer than three distinct")
})

test_that("multi-tau returns an ordered list of fits", {
  skip_on_cran()
  d <- sim_cont()
  m <- cbqr(y ~ x1 + x2, d, penalty = "none", quantile = c(0.25, 0.5, 0.75),
            ndraw = 3000, keep = 2)
  expect_s3_class(m, "cbqr.list")
  expect_length(m, 3L)
  b0 <- vapply(m, function(f) mean(f$beta0), numeric(1))
  expect_true(all(diff(b0) > 0))   # intercept increases with tau
})

test_that("no anchor arguments leak into the continuous interface", {
  ## sigma is identified once y is observed, so an anchor is a
  ## misspecification rather than a convention; the arguments must be absent,
  ## not silently accepted and ignored.
  a <- names(formals(cbqr))
  expect_false(any(c("anchor", "sigma_fixed") %in% a))
  expect_true(all(c("q", "standardize") %in% a))
})

## ---------------------------------------------------------------------------
## Methods. The predict test is the important one: before predict.cbqr existed,
## the inherited binary method returned P(y = 1 | x) from the ALD distribution
## function, which on a continuous fit saturates at 1 and returns a vector of
## ones with no error. A regression there would be silent.
## ---------------------------------------------------------------------------

test_that("predict returns the fitted quantile, not a probability", {
  skip_on_cran()
  skip_if_not_installed("quantreg")
  d <- sim_cont()
  set.seed(21)
  f <- cbqr(y ~ x1 + x2 + x3, d, penalty = "none", ndraw = 4000, keep = 2)
  p <- predict(f, newdata = utils::head(d, 20))
  ## The binary method returns P(y = 1 | x), so it is confined to [0, 1] and
  ## saturates at 1 here. Both facts are checked rather than a range bound on
  ## the linear predictor, which spans roughly 40 +/- 22 on this design and
  ## makes a tight interval a test of the simulation rather than of predict().
  expect_false(all(p >= 0 & p <= 1))
  expect_gt(diff(range(p)), 5)          # not saturated at a constant
  expect_gt(stats::cor(p, utils::head(d$y, 20)), 0.8)
  ref <- stats::predict(quantreg::rq(y ~ x1 + x2 + x3, tau = 0.5, data = d),
                        newdata = utils::head(d, 20))
  expect_lt(max(abs(p - ref)), 1.0)
})

test_that("predict(interval = TRUE) brackets the point estimate", {
  skip_on_cran()
  d <- sim_cont()
  set.seed(22)
  f <- cbqr(y ~ x1 + x2, d, penalty = "none", ndraw = 3000, keep = 2)
  nd <- utils::head(d, 10)
  m <- predict(f, nd, interval = TRUE)
  expect_equal(colnames(m), c("fit", "lower", "upper"))
  expect_true(all(m[, "lower"] < m[, "fit"]))
  expect_true(all(m[, "fit"] < m[, "upper"]))
  expect_equal(as.vector(m[, "fit"]), predict(f, nd), tolerance = 0.02)
})

test_that("print and summary work and report on the original y scale", {
  skip_on_cran()
  d <- sim_cont()
  set.seed(23)
  f <- cbqr(y ~ x1 + x2, d, penalty = "alasso", ndraw = 3000, keep = 2)
  out <- utils::capture.output(print(f))
  expect_true(any(grepl("continuous response", out)))
  expect_true(any(grepl("q = 1 (layer default)", out, fixed = TRUE)))
  ## No identification line: there is no anchor in the continuous model.
  expect_false(any(grepl("Identification", out)))

  s <- summary(f)
  expect_s3_class(s, "cbqr.summary")
  expect_true("sigma" %in% rownames(s$sigma))
  sout <- utils::capture.output(print(s))
  expect_true(any(grepl("ORIGINAL y scale", sout)))
  ## The binary summary's "identified only up to a positive scale" caveat is
  ## false here and must not appear.
  expect_false(any(grepl("up to a positive scale", sout)))
})

test_that("a non-default q is flagged in the printed output", {
  skip_on_cran()
  d <- sim_cont()
  set.seed(24)
  f <- cbqr(y ~ x1 + x2, d, penalty = "alasso", q = 2, ndraw = 2000, keep = 2)
  expect_true(any(grepl("NON-DEFAULT", utils::capture.output(print(f)))))
})

test_that("plot accepts sigma as a parameter and rejects unknown names", {
  skip_on_cran()
  d <- sim_cont()
  set.seed(25)
  f <- cbqr(y ~ x1 + x2, d, penalty = "none", ndraw = 2000, keep = 2)
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(plot(f, which = "sigma", type = "density"))
  expect_error(plot(f, which = "not_a_parameter"), "Unknown parameter")
})

test_that("the continuous adaptive lasso ships the loosened omega hyperprior", {
  ## Gamma(0.1, 0.1), not the binary hierarchy's Gamma(2, 2). Experiment C1's
  ## 6 x 2 factorial measured the tighter default as BINDING: the posterior
  ## wants ~23% of omega's mass below 1e-3 and Gamma(2,2) permits none of it,
  ## costing 5.8% RMSE and 2.3 points of coverage at n = 1000. See the comment
  ## on .cbqr_prior for the full decomposition.
  pr <- bbqr:::.cbqr_prior("alasso")
  expect_equal(c(pr$alpha_omega, pr$beta_omega), c(0.1, 0.1))

  ## Propriety is what must not regress: Proposition 1C needs both strictly
  ## positive and imposes no lower bound, so the loosened prior is exactly as
  ## proper as the reference. A zero here would be an improper hierarchy.
  expect_true(pr$alpha_omega > 0 && pr$beta_omega > 0)

  ## The change is scoped. The lasso layer's own scale hyperprior was NOT
  ## covered by that factorial and must stay at the reference value, and the
  ## binary constructor must be untouched -- the binary sensitivity study
  ## reached a different conclusion about what the binding costs.
  expect_equal(c(bbqr:::.cbqr_prior("lasso")$alpha_tau,
                 bbqr:::.cbqr_prior("lasso")$beta_tau), c(2, 2))
  pb <- prior("alasso", model = "v5")
  expect_equal(c(pb$alpha_omega, pb$beta_omega), c(2, 2))

  ## And it must still fit.
  d <- sim_cont(n = 200)
  f <- cbqr_alasso(y ~ x1 + x2 + x3, d, ndraw = 400, burn = 100)
  expect_true(all(is.finite(f$omega) & f$omega > 0))
  expect_true(all(is.finite(f$beta)))
})
