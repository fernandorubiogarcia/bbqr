## model = "v3" / "v4" / "v5": the three adaptive-lasso hierarchies of
## appendices V3, V4 and V5.
##
## The three differ only in the prior on omega and the prior on the intercept,
## and both differences are carried as ARGUMENT VALUES rather than as branches
## on a model identifier:
##
##   v3   alpha_omega = beta_omega = 0        p(omega) prop 1/omega, flat beta0
##   v4   alpha_omega = beta_omega = 2        proper omega,          flat beta0
##   v5   v4 plus b0_prec = 1/V00 > 0         proper omega and beta0
##
## so "does the switch work" is really "do these numbers reach the kernel".

test_that("the version presets are the ones the appendices specify", {
  p3 <- prior("alasso", model = "v3")
  p4 <- prior("alasso", model = "v4")
  p5 <- prior("alasso", model = "v5")

  ## v3 is the historical specification and is retained only as a baseline.
  expect_identical(c(p3$alpha_omega, p3$beta_omega), c(0, 0))
  expect_identical(c(p3$a, p3$b), c(1, 1))
  expect_true(is.infinite(p3$b0_var))

  ## v4 needs a > 1: a flat intercept against a free sigma otherwise leaves
  ## int_0^eps sigma^(a-2) dsigma divergent, and a = 1 sits on that boundary.
  expect_identical(c(p4$alpha_omega, p4$beta_omega), c(2, 2))
  expect_identical(c(p4$a, p4$b), c(2, 2))
  expect_true(is.infinite(p4$b0_var))

  ## v5's intercept prior is tau-dependent, so it cannot be resolved until the
  ## fit knows the quantile; NA is the sentinel for "use the reference".
  expect_identical(c(p5$alpha_omega, p5$beta_omega), c(2, 2))
  expect_true(is.na(p5$b0_var) && is.na(p5$b0_mean))

  ## delta is Gamma(2, 2) in all three: a flat prior on delta is a separate
  ## propriety failure and not a competing setting.
  for (p in list(p3, p4, p5)) {
    expect_identical(p$alpha_delta, 2)
    expect_identical(p$beta_delta, 2)
  }
})

test_that("the package default is v5 at sigma = 1, and says so", {
  ## A user who calls bbqr() with nothing but a formula must get a posterior
  ## that exists. The published hierarchy stays one argument away and is
  ## labelled as improper wherever it is printed.
  s <- sim_binary()
  fit <- bbqr(s$formula, data = s$data, ndraw = 300, burn = 50)
  expect_identical(fit$model, "v5")
  expect_identical(fit$anchor, "sigma1")
  expect_true(fit$derived)
  expect_output(print(fit), 'model = "v5"')
  expect_output(print(summary(fit)), 'model = "v5"')
  expect_output(print(prior("alasso")), "proper omega and intercept")
  expect_output(print(prior("lasso")), "proper tau_h, delta and intercept")
  expect_output(print(prior("alasso", model = "v3")), "IMPROPER")
  expect_output(print(prior("lasso", model = "v3")), "IMPROPER")

  ## A projection anchor under the default is flagged on the printed object.
  expect_warning(pj <- bbqr(s$formula, data = s$data, anchor = "norm1",
                            ndraw = 300, burn = 50), "derived posterior")
  expect_output(print(pj), "not a derived posterior")

  ## An object from a build before the version split has no `model` and the
  ## print method says nothing about a hierarchy rather than guessing.
  old <- fit; old$model <- NULL
  expect_false(any(grepl("Hierarchy", capture.output(print(old)))))
})

test_that("the tau-calibrated reference intercept prior is the ALD moment match", {
  for (tau in c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9)) {
    r <- bbqr:::.bbqr_b0_reference(tau)
    ## mean = -theta, the negated mixture constant; var = Var(ALD(0,1,tau)).
    expect_equal(r$mean, -(1 - 2 * tau) / (tau * (1 - tau)))
    expect_equal(r$var, (1 - 2 * tau + 2 * tau^2) / (tau^2 * (1 - tau)^2))
  }
  ## At the median the prior is centred, and its variance is that of the
  ## Laplace(scale 2) error.
  expect_equal(bbqr:::.bbqr_b0_reference(0.5), list(mean = 0, var = 8))
  ## Away from the median it is NOT centred. Centring at zero would put the
  ## prior predictive mean of Pr(y = 1) at 0.73 rather than 0.5 at tau = 0.1.
  expect_true(bbqr:::.bbqr_b0_reference(0.1)$mean < -8)
  expect_error(bbqr:::.bbqr_b0_reference(0), "strictly between")
  expect_error(bbqr:::.bbqr_b0_reference(1), "strictly between")
})

test_that("a prior cannot name one hierarchy and carry another", {
  ## Each of these would produce a fit labelled with a model it does not
  ## implement, which is worse than an error.
  expect_error(prior("alasso", model = "v4", alpha_omega = 0, beta_omega = 0),
               "strictly positive")
  expect_error(prior("alasso", model = "v3", alpha_omega = 2, beta_omega = 2),
               "alpha_omega = beta_omega = 0")
  expect_error(prior("alasso", model = "v5", b0_var = Inf), "cannot be Inf")
  expect_error(prior("alasso", model = "v4", b0_var = 8), "flat intercept")
  expect_error(prior("lasso", model = "v4"), "omega exists only")
  expect_error(prior("none", model = "v4"), "omega exists only")
  expect_error(prior("alasso", model = "v9"), "should be one of")
})

test_that("a prior serialized before the version split is read as v3", {
  ## Such an object described the V3 hierarchy, so the V3 values are filled in
  ## rather than erroring: the numbers reproduce what it actually meant.
  old <- prior("alasso", model = "v3")
  old$model <- old$alpha_omega <- old$beta_omega <- NULL
  old$b0_mean <- old$b0_var <- NULL
  p <- bbqr:::.bbqr_resolve_prior(old, "alasso")
  expect_identical(p$model, "v3")
  expect_identical(c(p$alpha_omega, p$beta_omega), c(0, 0))
  expect_true(is.infinite(p$b0_var))

  ## But a hand-edited object that names a proper model while carrying V3's
  ## log-uniform prior on omega is refused.
  bad <- prior("alasso", model = "v4")
  bad$alpha_omega <- 0
  expect_error(bbqr:::.bbqr_resolve_prior(bad, "alasso"), "does not exist")
})

test_that("v4 and v5 refuse every device their appendices exclude", {
  s <- sim_binary()
  f <- function(...) bbqr_alasso(s$formula, data = s$data, ndraw = 200,
                                 burn = 40, ...)
  for (m in c("v4", "v5")) {
    p <- prior("alasso", model = m)
    ## No clamp, no smooth barrier, no cancelling GIG root: the appendices say
    ## in terms that none of these is part of the target.
    expect_error(f(prior = p, clip = "v7clip"), "admits no numerical clamps")
    expect_error(f(prior = p, clip = "v4"), "admits no numerical clamps")
    expect_error(f(prior = p, kappa_lam = 1e-6), "kappa_lam = 0")
    expect_error(f(prior = p, sampler = "v4_invgaus"), "gig_half")
    ## The post-draw projections are not in any derivation either, but they
    ## are WARNED rather than refused: the anchor comparison is a deliberate
    ## part of the study, so the run proceeds and the departure is recorded as
    ## derived = FALSE. The projection is an exact reparameterisation of the
    ## likelihood but not of the prior, and under a shrinkage prior the penalty
    ## depends on |beta_j|, so rescaling mid-sweep changes it every sweep.
    for (a in c("beta1", "norm1", "normslopes")) {
      expect_warning(fit <- f(prior = p, anchor = a), "derived posterior")
      expect_false(fit$derived)
    }
    ## The two conventions the derivations do cover are accepted and derived.
    for (a in c("free", "sigma1")) {
      fit <- suppressWarnings(f(prior = p, anchor = a))
      expect_s3_class(fit, "bbqr")
      expect_true(fit$derived)
    }
  }
  ## v3 keeps them all reachable, because reproducing the historical build
  ## requires them.
  expect_s3_class(f(prior = prior("alasso", model = "v3"), clip = "v4",
                    sampler = "v4_invgaus"), "bbqr")
  expect_s3_class(f(prior = prior("alasso", model = "v3"), anchor = "norm1"),
                  "bbqr")
})

test_that("the a > 1 warning fires exactly where the theory says it should", {
  s <- sim_binary()
  p_bad <- prior("alasso", model = "v4", a = 1, b = 1)
  ## Free sigma: the intercept integral carries a factor 1/sigma, so a = 1 is
  ## the divergent boundary.
  expect_warning(bbqr_alasso(s$formula, data = s$data, ndraw = 200, burn = 40,
                             prior = p_bad, anchor = "free"),
                 "improper")
  ## Anchored sigma: that factor is constant and the condition is vacuous.
  expect_warning(bbqr_alasso(s$formula, data = s$data, ndraw = 200,
                             burn = 40, prior = p_bad, anchor = "sigma1"),
                 regexp = NA)
  ## The v4 reference satisfies it.
  expect_warning(bbqr_alasso(s$formula, data = s$data, ndraw = 200,
                             burn = 40, prior = prior("alasso", model = "v4"),
                             anchor = "free"),
                 regexp = NA)
})

test_that("the three hierarchies produce three different chains", {
  s <- sim_binary()
  run <- function(m) {
    set.seed(11)
    bbqr_alasso(s$formula, data = s$data, ndraw = 1500, burn = 300,
                prior = prior("alasso", model = m))
  }
  f3 <- run("v3"); f4 <- run("v4"); f5 <- run("v5")
  for (f in list(f3, f4, f5)) {
    expect_true(all(is.finite(f$beta)))
    expect_true(all(is.finite(f$omega)))
    expect_true(all(f$omega > 0))
  }
  expect_false(isTRUE(all.equal(f3$beta, f4$beta)))
  expect_false(isTRUE(all.equal(f4$beta, f5$beta)))
  expect_identical(c(f3$model, f4$model, f5$model), c("v3", "v4", "v5"))

  ## The resolved intercept prior is recorded, because under the reference
  ## calibration it depends on tau and so cannot be reconstructed from the
  ## prior object alone.
  expect_identical(unname(f4$b0_prior[["prec"]]), 0)
  expect_equal(unname(f5$b0_prior[["var"]]), 8)
  expect_equal(unname(f5$b0_prior[["prec"]]), 1 / 8)
  expect_identical(unname(f5$omega_prior), c(2, 2))

  ## No clamp is enabled under v4 or v5, so none can have bound.
  expect_false(any(f4$clip))
  expect_false(any(f5$clip))
})

test_that("a tight intercept prior pins beta0 at its mean", {
  ## The sharpest available check that b0_mean and b0_prec reach the kernel and
  ## enter the update the derivation says they do.
  s <- sim_binary()
  set.seed(7)
  f <- bbqr_alasso(s$formula, data = s$data, ndraw = 800, burn = 200,
                   prior = prior("alasso", model = "v5",
                                 b0_mean = 3, b0_var = 1e-8))
  expect_equal(mean(f$beta0), 3, tolerance = 1e-3)
  expect_lt(sd(f$beta0), 1e-3)
})

test_that("omega avoids the origin under v4 and v5 but not under v3", {
  ## V3's posterior does not exist, and the divergence sits on the wedge where
  ## omega, lambda_j^2 and s_j vanish together. The observable signature is
  ## omega mass piling up near zero, which the Gamma(2, 2) prior of v4 and v5
  ## removes because its density vanishes linearly at the origin.
  ##
  ## Note this is a diagnostic, not a proof: a Gibbs sampler on an improper
  ## posterior is null-recurrent and runs without complaint, so a short chain
  ## can miss the drift entirely. The chain here is long enough on this design
  ## to see it, and the test is one-sided for that reason.
  skip_on_cran()
  s <- sim_binary(n = 200, beta = c(1.5, -1, 0, 0, 0, 0))
  run <- function(m) {
    set.seed(4321)
    bbqr_alasso(s$formula, data = s$data, ndraw = 6000, burn = 1500,
                prior = prior("alasso", model = m))$omega
  }
  expect_gt(min(run("v3")), 0)          # still a valid draw, just tiny
  expect_gt(min(run("v4")), 1e-4)
  expect_gt(min(run("v5")), 1e-4)
})


## ---------------------------------------------------------------------------
## The lasso and the unpenalised kernel under the same v3 / v5 split.
## ---------------------------------------------------------------------------

test_that("the lasso carries the same defect as the adaptive lasso", {
  ## Its published tau_h conditional is Gamma(delta, lambda^2) -- no additive
  ## constants -- which is the signature of p(tau_h) prop 1/tau_h, the same
  ## log-uniform hyperprior that leaves the adaptive lasso improper in omega.
  ## With g = lambda^2*tau_h ~ Gamma(delta,1), u_j = lambda^2*s_j/2 ~ Exp(1)
  ## and q_j = beta_j/sqrt(s_j) ~ N(0,1) the prior measure factorises as
  ## (d tau_h / tau_h) times proper factors, and along tau_h -> 0 one has
  ## lambda^2 -> Inf, s_j -> 0 and beta_j = O(sqrt(tau_h)), so the likelihood
  ## is bounded below while the measure is not integrable.
  p3 <- prior("lasso", model = "v3")
  p5 <- prior("lasso", model = "v5")
  expect_identical(c(p3$alpha_tau, p3$beta_tau), c(0, 0))
  expect_identical(c(p5$alpha_tau, p5$beta_tau), c(2, 2))
  ## v3 keeps the published flat delta prior and a = b = 0.1; v5 makes every
  ## component proper. Note a = 0.1 is well inside the region where a flat
  ## intercept against a free sigma diverges, not merely on its boundary.
  expect_identical(c(p3$alpha_delta, p3$beta_delta), c(1, 0))
  expect_identical(c(p5$alpha_delta, p5$beta_delta), c(2, 2))
  expect_identical(c(p3$a, p3$b), c(0.1, 0.1))
  expect_true(is.infinite(p3$b0_var))
  expect_true(is.na(p5$b0_var))
  ## A prior may not name v5 while carrying the published hyperprior.
  expect_error(prior("lasso", model = "v5", alpha_tau = 0, beta_tau = 0),
               "strictly positive")
  expect_error(prior("lasso", model = "v3", alpha_tau = 2, beta_tau = 2),
               "alpha_tau = beta_tau = 0")
})

test_that("the unpenalised kernel needs only the intercept", {
  ## beta_j ~ N(0, beta_var) is already proper and sigma is fixed at 1, so
  ## there is no unidentified scale ray and the flat intercept is admissible
  ## given class overlap. A proper intercept prior removes that last condition.
  p3 <- prior("none", model = "v3")
  p5 <- prior("none", model = "v5")
  expect_true(is.infinite(p3$b0_var))
  expect_true(is.na(p5$b0_var))
  expect_identical(p3$beta_var, p5$beta_var)
})

test_that("v3 and v5 give different chains in all three kernels", {
  s <- sim_binary(n = 200, beta = c(1.5, -1, 0, 0))
  run <- function(pen, m) {
    set.seed(23)
    suppressWarnings(bbqr(s$formula, data = s$data, penalty = pen,
                          anchor = "sigma1", ndraw = 1200, burn = 250,
                          prior = prior(pen, model = m)))
  }
  for (pen in c("alasso", "lasso", "none")) {
    a <- run(pen, "v3"); b <- run(pen, "v5")
    expect_true(all(is.finite(a$beta)), info = pen)
    expect_true(all(is.finite(b$beta)), info = pen)
    expect_true(all(is.finite(b$beta0)), info = pen)
    expect_false(isTRUE(all.equal(a$beta, b$beta)), info = pen)
    expect_true(b$derived, info = pen)
    expect_equal(unname(b$b0_prior[["var"]]), 8)
  }
})

test_that("a tight intercept prior pins beta0 in every kernel", {
  ## The sharpest check that b0_mean and b0_prec reach each Fortran kernel in
  ## the right argument positions -- a misalignment there would pass a logical
  ## where a double is expected and corrupt silently rather than error.
  s <- sim_binary()
  for (pen in c("alasso", "lasso", "none")) {
    set.seed(31)
    f <- bbqr(s$formula, data = s$data, penalty = pen, anchor = "sigma1",
              ndraw = 700, burn = 150,
              prior = prior(pen, model = "v5", b0_mean = 2.5, b0_var = 1e-8))
    expect_equal(mean(f$beta0), 2.5, tolerance = 1e-3, info = pen)
    expect_lt(sd(f$beta0), 1e-3)
  }
})

test_that("the lasso's tau_h prior is live and keeps it off the origin", {
  skip_on_cran()
  s <- sim_binary(n = 200, beta = c(1.5, -1, 0, 0, 0))
  run <- function(m) {
    set.seed(5)
    bbqr(s$formula, data = s$data, penalty = "lasso", anchor = "sigma1",
         ndraw = 4000, burn = 800, prior = prior("lasso", model = m))$tauhyper
  }
  a <- run("v3"); b <- run("v5")
  ## Under the published prior tau_h roams over many orders of magnitude; the
  ## Gamma(2,2) prior vanishes linearly at the origin and stops it.
  expect_gt(min(b), min(a))
  expect_gt(min(b), 1e-3)
})
