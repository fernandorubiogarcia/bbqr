## Geweke's joint-distribution test for the continuous kernels -- including the
## penalised ones, which the binary test could not reach.
##
##   Geweke, J. (2004). Getting It Right: Joint Distribution Tests of Posterior
##   Simulators. JASA 99(467), 799-804.
##
## test-geweke-joint.R covers bbqr_none only, because theta = (beta0, beta) is
## the whole non-latent state there and nothing else can be initialised. cbqr
## takes sigma_init, delta_init and omega_init/tau_init in its prior, so the
## penalised hierarchies -- the ones the propriety corrections are actually
## about -- can be tested the same way.
##
## The hierarchies below were read off the full conditionals in the Fortran and
## checked against them term by term. For the adaptive lasso, kernel step (7)
## draws lambda_j^2 ~ InvGamma(delta + 1, sigma^q s_j / 2 + omega), which is an
## InvGamma(delta, omega) prior updated by the Exp(sigma^q / (2 lambda_j^2))
## term in step (5); step (8) draws omega ~ Gamma(alpha_om + k delta,
## beta_om + sum 1/lambda_j^2), which is Gamma(alpha_om, beta_om) updated by k
## InvGamma(delta, omega) factors. The lasso differs in two ways that matter:
## lambda^2 is a single global scale with a Gamma (not InvGamma) prior, and
## s_j's rate multiplies by sigma^q where the adaptive lasso divides.
##
## Sweeps per iteration. z, s and lambda^2 reset to 1 on every call, so a short
## run is not the stationary kernel. The perturbation decays in L:
##   cbqr_none    L = 1, 3, 10, 20, 40 -> 6.78, 2.84, 1.86, 1.25, 1.93
##   cbqr_alasso  L = 1, 5, 15, 30, 60 -> 12.69, 7.20, 3.00, 2.41, 2.25
##   cbqr_lasso   L = 1, 15, 30, 60    -> 15.35, 2.67, 0.95, 1.61
## A wrong kernel would not decay. Quadrupling M at fixed L left the adaptive
## lasso at 1.39 -> 1.95 with the largest z moving to a different statistic,
## where a bias would have doubled it: the residual is noise.
##
## Statistics. sigma, delta, omega and tau_h have Gamma marginals under the
## joint, so their raw moments are safe. beta is a scale mixture whose mixing
## variance has infinite expectation whenever delta <= 1, which the delta prior
## gives positive mass, so beta and y enter only through bounded functionals --
## plus log|beta_j|, which is what actually sees an error in the beta scale.

ald_r <- function(m, tau, sigma) {
  p <- stats::runif(m)
  ifelse(p <= tau, log(p / tau) / ((1 - tau) * sigma),
                   -log((1 - p) / (1 - tau)) / (tau * sigma))
}
gw_nse <- function(x, nb = 40) {
  bs <- floor(length(x) / nb)
  stats::var(colMeans(matrix(x[seq_len(nb * bs)], nrow = bs))) / nb
}

## One driver for all three continuous penalties. `mc_*` overrides let the
## truth simulator disagree with the sampler on purpose, which is how the test
## is shown to have power rather than merely to pass.
gw_cont <- function(penalty, tau = 0.5, M = 3000, L = 60, n = 12,
                    a = 4, b = 4, ad = 4, bd = 2, ah = 2, bh = 2, bv = 4,
                    b0v = 4, q = NULL, mc_a = NULL, mc_ad = NULL,
                    mc_bv = NULL, mc_tau = NULL, seed = 20260823) {
  if (is.null(mc_a))   mc_a   <- a
  if (is.null(mc_ad))  mc_ad  <- ad
  if (is.null(mc_bv))  mc_bv  <- bv
  if (is.null(mc_tau)) mc_tau <- tau
  if (is.null(q)) q <- switch(penalty, alasso = 1, lasso = 0, none = 0)
  set.seed(seed)
  X <- matrix(round(stats::rnorm(n * 2), 3), n, 2)
  vn <- c("V1", "V2")

  ## prior draw and the prior list handed to the sampler, per penalty
  dp <- function() {
    sig <- stats::rgamma(1, mc_a, rate = b)
    out <- list(b0 = stats::rnorm(1, 0, sqrt(b0v)), sig = sig)
    if (penalty == "none") {
      out$b <- stats::rnorm(2, 0, sqrt(mc_bv))
    } else {
      out$del <- stats::rgamma(1, mc_ad, rate = bd)
      out$hyp <- stats::rgamma(1, ah, rate = bh)
      s <- if (penalty == "alasso") {
        lam <- 1 / stats::rgamma(2, shape = out$del, rate = out$hyp)
        stats::rexp(2, rate = sig^q / (2 * lam))
      } else {
        lam <- stats::rgamma(1, shape = out$del, rate = out$hyp)
        stats::rexp(2, rate = sig^q * lam / 2)
      }
      out$b <- stats::rnorm(2, 0, sqrt(s))
    }
    out
  }
  mkpr <- function(th) {
    pr <- list(a = a, b = b, b0_mean = 0, b0_var = b0v, sigma_init = th$sig)
    if (penalty == "none") return(c(pr, list(beta_var = bv)))
    pr <- c(pr, list(alpha_delta = ad, beta_delta = bd, delta_init = th$del))
    if (penalty == "alasso")
      c(pr, list(alpha_omega = ah, beta_omega = bh, omega_init = th$hyp))
    else
      c(pr, list(alpha_tau = ah, beta_tau = bh, tau_init = th$hyp))
  }
  dy <- function(th) as.numeric(th$b0 + X %*% th$b) + ald_r(n, mc_tau, th$sig)
  gs <- function(th, y) {
    base <- c(th$b0, th$b0^2, tanh(th$b[1]), tanh(th$b[2]),
              log(abs(th$b[1])), log(abs(th$b[2])),
              th$sig, th$sig^2, log(th$sig), mean(y > 0), mean(tanh(y)))
    if (penalty == "none") base
    else c(base, th$del, log(th$del), th$hyp, log(th$hyp))
  }
  fitter <- switch(penalty, none = cbqr_none, lasso = cbqr_lasso,
                   alasso = cbqr_alasso)
  nstat <- length(gs(dp(), rep(0, n)))

  MC <- t(vapply(seq_len(M), function(i) { th <- dp(); gs(th, dy(th)) },
                 numeric(nstat)))
  th <- dp(); SC <- matrix(NA_real_, M, nstat)
  for (t in seq_len(M)) {
    y <- dy(th); SC[t, ] <- gs(th, y)
    dd <- data.frame(y = y); dd[vn] <- X
    args <- list(y ~ V1 + V2, dd, quantile = tau, ndraw = L, burn = 0, keep = 1,
                 prior = mkpr(th), standardize = FALSE,
                 beta_init = stats::setNames(th$b, vn), beta0_init = th$b0)
    if (penalty != "none") args$q <- q
    fit <- try(do.call(fitter, args), silent = TRUE)
    if (inherits(fit, "try-error")) next
    nxt <- list(b0 = fit$beta0[L], b = as.numeric(fit$beta[L, ]),
                sig = fit$sigma[L])
    if (penalty == "alasso") { nxt$del <- fit$delta[L]; nxt$hyp <- fit$omega[L] }
    if (penalty == "lasso")  { nxt$del <- fit$delta[L]; nxt$hyp <- fit$tau_hyper[L] }
    if (all(is.finite(unlist(nxt)))) th <- nxt
  }
  (colMeans(MC) - colMeans(SC)) /
    sqrt(apply(MC, 2, stats::var) / M + apply(SC, 2, gw_nse))
}

## |z| > 4.5 on any of eleven to fifteen moments. Under a correct sampler that
## happens about once in ten thousand runs; the observed controls sit at 1.8 to
## 2.9, and the injected errors below at 15 to 27.
test_that("the continuous unpenalised kernel passes Geweke's test", {
  skip_on_cran()
  expect_lt(max(abs(gw_cont("none", tau = 0.5, M = 4000, L = 20))), 4.5)
})

test_that("the continuous adaptive-lasso kernel passes Geweke's test", {
  skip_on_cran()
  ## The V5 correction lives in this hierarchy, so this is the one that matters.
  expect_lt(max(abs(gw_cont("alasso", tau = 0.5, M = 3000, L = 60))), 4.5)
})

test_that("the continuous lasso kernel passes Geweke's test", {
  skip_on_cran()
  expect_lt(max(abs(gw_cont("lasso", tau = 0.5, M = 3000, L = 60))), 4.5)
})

test_that("the continuous test also passes away from the median", {
  skip_on_cran()
  expect_lt(max(abs(gw_cont("alasso", tau = 0.30, M = 3000, L = 60))), 4.5)
  expect_lt(max(abs(gw_cont("alasso", tau = 0.75, M = 3000, L = 60))), 4.5)
})

test_that("the continuous test detects a sampler that is wrong", {
  skip_on_cran()
  ## Each of these is a hierarchy the sampler is not the posterior for.
  expect_gt(max(abs(gw_cont("alasso", M = 3000, L = 60, mc_ad = 5))), 4.5)
  expect_gt(max(abs(gw_cont("alasso", M = 3000, L = 60, mc_a  = 5))), 4.5)
  expect_gt(max(abs(gw_cont("lasso",  M = 3000, L = 60, mc_tau = 0.55))), 4.5)
  expect_gt(max(abs(gw_cont("none",   M = 4000, L = 20, mc_bv = 8))), 4.5)
})
