## Geweke's joint-distribution test of the sampler itself.
##
##   Geweke, J. (2004). Getting It Right: Joint Distribution Tests of Posterior
##   Simulators. Journal of the American Statistical Association 99, 799-804.
##
## Every other test in this suite checks a consequence of the sampler being
## right: that constraints hold, that coefficients are recovered, that draws are
## finite. None of them checks the thing itself -- that the transition kernel
## leaves the posterior invariant. A kernel with a wrong full conditional can
## satisfy all of them and still converge to the wrong distribution.
##
## The test draws the same joint p(beta0, beta, y) = p(beta0, beta) p(y | .) two
## ways:
##
##   marginal-conditional  theta ~ prior, then y ~ p(y | theta).  Exact, iid.
##   successive-conditional  a two-block Gibbs chain,
##                             y     ~ p(y | theta)      exact, closed form
##                             theta ~ K(theta -> . ; y)  the package's sweep
##     whose invariant law is p(theta, y) if and only if K leaves p(theta | y)
##     invariant.
##
## So the two agree iff the sampler is correct, and the comparison is a z test
## on moments of the two samples. p(y | theta) is written out here from the ALD
## rather than taken from predict(), so the test does not lean on the code it is
## testing.
##
## Scope: penalty = "none" at anchor = "sigma1", where theta = (beta0, beta) is
## the whole non-latent state and beta_init/beta0_init can set it. The BINARY
## lasso and adaptive-lasso kernels carry lambda^2, omega and delta across
## sweeps with no way to initialise them, so the same construction needs a
## full-state entry point; benchmarks/geweke_joint.R records what that would
## take. Their CONTINUOUS counterparts do have that entry point, and are tested
## in test-geweke-continuous.R.
##
## Sweeps per iteration: the persistent state is (beta0, beta, z), and every
## call resets z = 1 before its loop, so ONE sweep per call is not the
## stationary kernel -- it is "reset z, then sweep". z is redrawn from its full
## conditional inside sweep 1, so the reset is a first-sweep perturbation that
## decays: max |z| runs 5.26, 2.95, 1.65, 1.71 at L = 1, 3, 10, 40. A genuine
## kernel error would not decay with L. L = 20 is used here.

ald_cdf <- function(u, tau)
  ifelse(u <= 0, tau * exp((1 - tau) * u), 1 - (1 - tau) * exp(-tau * u))

## P(y = 1 | eta) = P(eps > -eta) for eps ~ ALD(tau, sigma = 1)
p_y1 <- function(eta, tau) 1 - ald_cdf(-eta, tau)

## Numerical standard error of a mean from an autocorrelated chain.
nse_batch <- function(x, nb = 40) {
  bs <- floor(length(x) / nb)
  stats::var(colMeans(matrix(x[seq_len(nb * bs)], nrow = bs))) / nb
}

geweke_z <- function(tau, M, L = 20, n = 10, beta_var = 4, b0_var = 4,
                     seed = 20260823) {
  set.seed(seed)
  X <- matrix(round(stats::rnorm(n * 2), 3), n, 2)
  vn <- c("V1", "V2")
  FORM <- y ~ V1 + V2
  PR <- prior("none", beta_var = beta_var, b0_mean = 0, b0_var = b0_var)

  draw_prior <- function()
    list(b0 = stats::rnorm(1, 0, sqrt(b0_var)), b = stats::rnorm(2, 0, sqrt(beta_var)))
  draw_y <- function(th)
    stats::rbinom(n, 1, p_y1(as.numeric(th$b0 + X %*% th$b), tau))
  gstat <- function(th, y) c(th$b0, th$b[1], th$b[2],
                             th$b0^2, th$b[1]^2, th$b[2]^2, sum(y))

  MC <- t(vapply(seq_len(M), function(i) {
    th <- draw_prior(); gstat(th, draw_y(th))
  }, numeric(7)))

  th <- draw_prior()
  SC <- matrix(NA_real_, M, 7)
  for (t in seq_len(M)) {
    y <- draw_y(th)
    SC[t, ] <- gstat(th, y)
    if (length(unique(y)) < 2L) next          # needs both classes; theta unchanged
    dd <- data.frame(y = y); dd[vn] <- X
    fit <- bbqr_none(FORM, dd, quantile = tau, ndraw = L, burn = 0, keep = 1,
                     prior = PR, anchor = "sigma1",
                     beta_init = stats::setNames(th$b, vn), beta0_init = th$b0)
    th <- list(b0 = fit$beta0[L], b = as.numeric(fit$beta[L, ]))
  }
  (colMeans(MC) - colMeans(SC)) /
    sqrt(apply(MC, 2, stats::var) / M + apply(SC, 2, nse_batch))
}

test_that("the unpenalised kernel passes Geweke's joint-distribution test", {
  skip_on_cran()
  ## |z| > 4 on any of seven moments. Under a correct sampler each z is
  ## standard normal, so the seven together clear this about 9997 times in
  ## 10000; a kernel that targets the wrong law does not clear it at all.
  z <- geweke_z(tau = 0.5, M = 5000)
  expect_lt(max(abs(z)), 4)
})

test_that("the joint-distribution test also passes away from the median", {
  skip_on_cran()
  ## tau != 1/2 is where the ALD is asymmetric and where the intercept prior is
  ## calibrated, so a sign or 1-tau slip that cancels at the median shows here.
  expect_lt(max(abs(geweke_z(tau = 0.30, M = 5000))), 4)
  expect_lt(max(abs(geweke_z(tau = 0.75, M = 5000))), 4)
})

test_that("the joint-distribution test detects a sampler that is wrong", {
  skip_on_cran()
  ## A test that never fails proves nothing. Give the sampler a prior different
  ## from the one the truth simulator drew from -- the same footprint a wrong
  ## full conditional leaves -- and the statistic must move.
  set.seed(11)
  n <- 10
  X <- matrix(round(stats::rnorm(n * 2), 3), n, 2)
  vn <- c("V1", "V2"); M <- 5000; L <- 20; tau <- 0.5
  mc_var <- 5; fit_var <- 4                     # 25% apart
  PR <- prior("none", beta_var = fit_var, b0_mean = 0, b0_var = 4)

  dp <- function() list(b0 = stats::rnorm(1, 0, 2), b = stats::rnorm(2, 0, sqrt(mc_var)))
  dy <- function(th) stats::rbinom(n, 1, p_y1(as.numeric(th$b0 + X %*% th$b), tau))
  g  <- function(th, y) c(th$b0, th$b[1], th$b[2], th$b0^2, th$b[1]^2, th$b[2]^2, sum(y))

  MC <- t(vapply(seq_len(M), function(i) { th <- dp(); g(th, dy(th)) }, numeric(7)))
  th <- dp(); SC <- matrix(NA_real_, M, 7)
  for (t in seq_len(M)) {
    y <- dy(th); SC[t, ] <- g(th, y)
    if (length(unique(y)) < 2L) next
    dd <- data.frame(y = y); dd[vn] <- X
    fit <- bbqr_none(y ~ V1 + V2, dd, quantile = tau, ndraw = L, burn = 0,
                     keep = 1, prior = PR, anchor = "sigma1",
                     beta_init = stats::setNames(th$b, vn), beta0_init = th$b0)
    th <- list(b0 = fit$beta0[L], b = as.numeric(fit$beta[L, ]))
  }
  z <- (colMeans(MC) - colMeans(SC)) /
       sqrt(apply(MC, 2, stats::var) / M + apply(SC, 2, nse_batch))
  expect_gt(max(abs(z)), 4)
})
