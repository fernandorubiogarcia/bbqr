## The z/s draw has two implementations of the same GIG(1/2) conditional.
##
## They are the SAME expression: with t = mu*nu/lambda, the textbook form used
## by the v4 routine is
##     q = mu*(2 + t - sqrt(4t + t^2)) / 2
## and rgig_half's is
##     q = 2*mu / (2 + t + sqrt(4t + t^2)),
## which are equal because (2+t-s)(2+t+s) = (2+t)^2 - (4t+t^2) = 4. rgig_half is
## the rationalized form of the same algebra.
##
## So the switch is not two different laws, it is two different floating-point
## paths to one law, and the arms it defines only separate where the textbook
## form loses precision -- which is at LARGE t, where it subtracts two nearly
## equal quantities and eventually returns q = 0, hence a latent of Inf. That is
## the failure v4's guards were written to catch, and it is why "historical
## routine + v4 guards" is the only historical configuration worth running.
##
## These tests therefore pin three things: the switch is wired and recorded, the
## two paths agree where the expression is well conditioned, and they part
## company where it is not.

test_that("the sampler argument is validated and recorded", {
  s <- sim_binary(n = 80)
  expect_error(bbqr_alasso(s$formula, data = s$data, sampler = "rgig"),
               "should be one of")
  expect_error(bbqr(s$formula, data = s$data, penalty = "lasso",
                    sampler = "v4_invgaus"),
               "only available for penalty")

  for (nm in c("gig_half", "v4_invgaus")) {
    fit <- bbqr_alasso(s$formula, data = s$data, ndraw = 200, burn = 40,
                       prior = v3_prior(), sampler = nm, clip = "v4")
    expect_identical(fit$sampler, nm)
  }
  ## The default must stay the exact routine. A silent switch to the historical
  ## one would move every number the package produces.
  expect_identical(
    bbqr_alasso(s$formula, data = s$data, ndraw = 200, burn = 40)$sampler,
    "gig_half")
})

test_that("the two routines implement the same law", {
  s <- sim_binary(n = 120)
  run <- function(nm) {
    set.seed(4242L)
    bbqr_alasso(s$formula, data = s$data, ndraw = 500, burn = 100,
                prior = v3_prior(), sampler = nm, clip = "v4")
  }
  a <- run("gig_half"); b <- run("v4_invgaus")
  ## Same seed, same guards, and on data this well behaved t never reaches the
  ## regime where the textbook form loses digits -- so the chains should agree
  ## to within rounding. A LARGE discrepancy here would mean the historical path
  ## is wired up wrong, most likely by inverting the draw the wrong way round,
  ## which is the mistake this catches.
  expect_equal(colMeans(a$beta), colMeans(b$beta), tolerance = 1e-6)
  expect_true(all(is.finite(b$beta)))
  expect_gt(cos_sim(colMeans(b$beta), s$beta), 0.5)
})

test_that("the textbook form is the one that loses the digits", {
  ## A pure-arithmetic check of the two expressions, kept here because it is the
  ## whole reason the switch exists. No RNG and no fit: just the algebra.
  q_v4  <- function(mu, t) mu*(2 + t - sqrt(4*t + t*t))/2
  q_gig <- function(mu, t) 2*mu/(2 + t + sqrt(4*t + t*t))

  ## Well conditioned: the two agree to machine precision, which is the identity
  ## asserted in the comment above.
  for (t in 10^c(-6, -2, 1, 4))
    expect_equal(q_v4(1, t), q_gig(1, t), tolerance = 1e-8)

  ## Ill conditioned: the textbook form cancels to nothing and the latent, which
  ## is its reciprocal, goes to Inf. The rationalized form stays finite and
  ## accurate. This is the divergence the ablation is measuring.
  expect_gt(abs(q_v4(1, 1e8) - q_gig(1, 1e8)) / q_gig(1, 1e8), 0.1)
  expect_identical(q_v4(1, 1e12), 0)
  expect_true(is.infinite(1 / q_v4(1, 1e12)))
  expect_true(is.finite(1 / q_gig(1, 1e12)))
})

test_that("selecting the routine leaves the clip contract alone", {
  s <- sim_binary(n = 120)
  fit <- bbqr_alasso(s$formula, data = s$data, ndraw = 300, burn = 60,
                     prior = v3_prior(), sampler = "v4_invgaus", clip = "v4")
  expect_identical(names(fit$clip_hits), bbqr:::.BBQR_CLIP_SITES)
  ## The v4 argument floors bound psi and chi below by eps_small, which caps t at
  ## nu/eps_small and so keeps the textbook form out of its worst regime. They
  ## are what makes the historical routine survivable at all, and they must stay
  ## reachable and counted on this path.
  expect_true(all(fit$clip[c("v4_z_psi", "v4_z_chi",
                             "v4_s_lam", "v4_s_psi", "v4_s_chi")]))
  ## And the bounds on the sampled latents are the backstop for when it does not
  ## survive: an Inf from a cancelled q is caught here.
  expect_true(all(fit$clip[c("v4_zbound", "v4_sbound")]))
})

test_that("the three arms of the decomposition are distinct configurations", {
  s <- sim_binary(n = 120)
  run <- function(nm, cl) {
    set.seed(4242L)
    bbqr_alasso(s$formula, data = s$data, ndraw = 400, burn = 80,
                prior = v3_prior(), anchor = "free", sampler = nm, clip = cl)
  }
  v6      <- run("gig_half",   "v6")
  gig_v4  <- run("gig_half",   "v4")
  invg_v4 <- run("v4_invgaus", "v4")

  ## Arm 1 -> 2 moves the guards; arm 2 -> 3 moves the routine. The provenance
  ## fields must distinguish all three, since that is what the analysis joins on.
  expect_identical(v6$sampler,      "gig_half")
  expect_identical(gig_v4$sampler,  "gig_half")
  expect_identical(invg_v4$sampler, "v4_invgaus")
  expect_false(any(v6$clip))
  expect_true(all(gig_v4$clip[bbqr:::.BBQR_CLIP_PRESETS$v4]))
  expect_identical(gig_v4$clip, invg_v4$clip)
})
