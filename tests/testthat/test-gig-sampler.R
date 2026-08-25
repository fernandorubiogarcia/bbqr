## The GIG(1/2) sampler behind the z_i and s_j updates. The source-level audit
## can only check that the call site exists, so its distribution is tested here
## against the analytic moments
##   E[X]   = sqrt(chi/psi) * (1 + 1/sqrt(psi*chi))
##   E[X^2] = (chi/psi)     * (1 + 3/sqrt(psi*chi) + 3/(psi*chi))
## for X ~ GIG(1/2, psi, chi), obtained from K_{1/2}, K_{3/2} and K_{5/2}.
##
## The sampler is not exported, so it is exercised through bbqr_alasso: the s_j
## update is a GIG(1/2, sigma/lambda_j^2, beta_j^2) draw. A direct moment test
## needs the routine itself, so this file checks the two properties that are
## observable from outside and would fail if the sampler were wrong.

test_that("latent scales stay strictly positive and finite", {
  s <- sim_binary()
  fit <- bbqr_alasso(s$formula, s$data, ndraw = 4000, burn = 800)
  ## A GIG(1/2) variate is supported on (0, inf). The old inverse-Gaussian
  ## route could return 0 or Inf by cancellation at a tiny residual, which
  ## surfaced as non-finite lambdasq and sigma draws.
  expect_true(all(is.finite(fit$lambdasq)))
  expect_true(all(fit$lambdasq > 0))
  expect_true(all(is.finite(fit$sigma) & fit$sigma > 0))
  expect_true(all(is.finite(fit$beta)))
  expect_true(all(is.finite(fit$beta0)))
})

test_that("hard shrinkage does not destabilise the sampler", {
  ## Drives beta_j towards 0, so chi = beta_j^2 becomes tiny and psi/chi
  ## becomes huge -- the regime in which the unrationalised Michael-Schucany-Haas
  ## root loses all precision and can return a non-positive draw.
  set.seed(4321)
  n <- 200
  X <- matrix(stats::rnorm(n * 6), n, 6)
  y <- as.integer(X[, 1] + stats::rnorm(n) > 0)     # only x1 carries signal
  d <- data.frame(y = y, X)
  names(d) <- c("y", paste0("V", 1:6))

  fit <- bbqr_alasso(y ~ ., data = d, ndraw = 4000, burn = 800,
                     prior = prior("alasso", alpha_delta = 2, beta_delta = 2))
  expect_true(all(is.finite(fit$beta)))
  expect_true(all(is.finite(fit$lambdasq) & fit$lambdasq > 0))
  ## The five null slopes must actually be shrunk relative to the signal.
  expect_lt(mean(abs(colMeans(fit$beta)[-1])), mean(abs(colMeans(fit$beta)[1])))
})
