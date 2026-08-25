## The per-site clamp switches. These replace the single `use_boundaries` flag,
## so the first thing to pin down is that the two corners still mean what they
## used to: all sites on is the old TRUE, all sites off is the old FALSE.

## Same seed, so any difference between two fits is the clamp configuration
## rather than a different position in one RNG stream.
paired <- function(expr, seed = 4242L) {
  set.seed(seed)
  force(expr)
}

test_that("the site vector is the Fortran contract's length", {
  sites <- bbqr:::.BBQR_CLIP_SITES
  expect_length(sites, 18L)
  expect_false(anyDuplicated(sites) > 0)

  s <- sim_binary()
  fit <- bbqr_alasso(s$formula, data = s$data, ndraw = 200, burn = 40)
  ## The kernel writes one counter per site; a mismatch here means the R site
  ## table and the C_* parameters in QRb_AL_mcmc.f95 have drifted apart.
  expect_length(fit$clip_hits, length(sites))
  expect_identical(names(fit$clip_hits), sites)
  expect_identical(names(fit$clip), sites)
})

test_that("use_boundaries addresses the eight gated sites and no others", {
  sites <- bbqr:::.BBQR_CLIP_SITES
  gated <- bbqr:::.BBQR_CLIP_GATED

  ## This is the back-compatibility contract. The v4 guards were applied
  ## unconditionally in v4, never gated by use_boundaries, so widening the site
  ## table must not give use_boundaries = TRUE a reach it never had -- otherwise
  ## re-running the v7 driver against this build would silently change its
  ## `_clip` arm from eight clamps to fifteen.
  on  <- bbqr:::.bbqr_resolve_clip(NULL, use_boundaries = TRUE)
  off <- bbqr:::.bbqr_resolve_clip(NULL, use_boundaries = FALSE)
  expect_identical(names(which(on)), gated)
  expect_false(any(off))

  ## A length-1 logical clip means the same thing, so it has the same reach.
  expect_identical(bbqr:::.bbqr_resolve_clip(TRUE, FALSE), on)
  expect_identical(bbqr:::.bbqr_resolve_clip(FALSE, TRUE), off)

  ## An explicit length-15 logical is taken literally, which is the only way to
  ## ask for a gated clamp and a v4 guard at the same time.
  all_on <- bbqr:::.bbqr_resolve_clip(rep(TRUE, length(sites)), FALSE)
  expect_true(all(all_on))
})

test_that("clip = NULL reproduces use_boundaries, both ways", {
  s <- sim_binary()
  for (cfg in list(list(FALSE, "v6"), list(TRUE, "v7clip"))) {
    a <- paired(bbqr_alasso(s$formula, data = s$data, ndraw = 400, burn = 80,
                            prior = v3_prior(),
                            anchor = "sigma1", use_boundaries = cfg[[1]]))
    b <- paired(bbqr_alasso(s$formula, data = s$data, ndraw = 400, burn = 80,
                            prior = v3_prior(),
                            anchor = "sigma1", clip = cfg[[2]]))
    expect_identical(as.numeric(a$beta), as.numeric(b$beta))
    expect_identical(a$beta0, b$beta0)
  }
})

test_that("the v7clip preset is the old all-clamps-on set", {
  s <- sim_binary()
  a <- paired(bbqr_alasso(s$formula, data = s$data, ndraw = 400, burn = 80,
                          prior = v3_prior(),
                          anchor = "sigma1", use_boundaries = TRUE))
  b <- paired(bbqr_alasso(s$formula, data = s$data, ndraw = 400, burn = 80,
                          prior = v3_prior(),
                          anchor = "sigma1", clip = "v7clip"))
  expect_identical(as.numeric(a$beta), as.numeric(b$beta))
})

test_that("v4's clamp set is disjoint from the gated set", {
  p <- bbqr:::.BBQR_CLIP_PRESETS
  expect_length(intersect(p$v7clip, p$v4), 0L)
  expect_length(p$v6, 0L)
  ## v4 ran use_boundaries = FALSE, so its guards and the gated clamps were
  ## never on together. "v6" is not "v4 minus clamps": it lacks v4's guards too.
  expect_true(all(c(p$v7clip, p$v4) %in% bbqr:::.BBQR_CLIP_SITES))
})

test_that("sites unreachable under anchor = 'sigma1' are inert", {
  s <- sim_binary()
  sites <- bbqr:::.BBQR_CLIP_SITES
  base <- paired(bbqr_alasso(s$formula, data = s$data, ndraw = 400, burn = 80,
                             prior = v3_prior(),
                             anchor = "sigma1", clip = "v7clip"))
  for (inert in c("sigma", "v4_rate_sig")) {
    cl <- sites %in% c(bbqr:::.BBQR_CLIP_PRESETS$v7clip, inert)
    alt <- paired(bbqr_alasso(s$formula, data = s$data, ndraw = 400, burn = 80,
                              prior = v3_prior(),
                              anchor = "sigma1", clip = cl))
    expect_identical(as.numeric(base$beta), as.numeric(alt$beta),
                     info = paste(inert, "should be inert under sigma1"))
    expect_equal(unname(alt$clip_hits[[inert]]), 0)
  }
})

test_that("v4_sd covers the initialisation sweep only", {
  ## v4 floored the truncation sd at its L80, in the initialisation sweep, and
  ## left the main-loop assignment at its L97 bare -- by then v4's hard z bound
  ## already kept phisq*z/sigma positive. Flooring both would be a stronger
  ## intervention than v4 applied, so the scope is pinned here.
  ##
  ## The floor is unreachable at sigma = 1: phisq = 2/(tau(1-tau)) is 8 at the
  ## median and z starts at 1, so the quantity is 8, far above eps_small = 1e-6.
  ## Holding sigma at 1e8 makes it 8e-8 and the floor bites -- exactly once per
  ## observation, in that one sweep. If the main-loop floor were still present
  ## the count would run well past n.
  s <- sim_binary(n = 150)
  fit <- bbqr_alasso(s$formula, data = s$data, ndraw = 500, burn = 100,
                     prior = v3_prior(),
                     anchor = "sigma1", sigma_fixed = 1e8, clip = "v4_sd")
  expect_equal(unname(fit$clip_hits[["v4_sd"]]), 150)
})

test_that("counters are recorded even when the switch is off", {
  s <- sim_binary(n = 300, beta = c(1.5, -1, 0, 0.5))
  sites <- bbqr:::.BBQR_CLIP_SITES
  off <- paired(bbqr_alasso(s$formula, data = s$data, ndraw = 600, burn = 120,
                            prior = v3_prior(),
                            anchor = "sigma1", clip = "v6"))
  ## Nothing was clamped, yet the run still reports which clamps would have
  ## fired -- that is what lets one baseline arm screen every site.
  expect_true(all(off$clip == FALSE))
  expect_true(any(off$clip_hits > 0))

  ## A site that would have fired must change the chain once switched on, and a
  ## site that would not must leave it alone. This is the whole premise of
  ## reading the counters as a screen.
  for (site in sites[!sites %in% c("sigma", "v4_rate_sig")]) {
    one <- paired(bbqr_alasso(s$formula, data = s$data, ndraw = 600, burn = 120,
                              prior = v3_prior(),
                              anchor = "sigma1", clip = site))
    changed <- !identical(as.numeric(off$beta), as.numeric(one$beta))
    expect_equal(changed, off$clip_hits[[site]] > 0,
                 info = sprintf("%s: would-fire=%g but changed=%s",
                                site, off$clip_hits[[site]], changed))
  }
})

test_that("a named clip overrides use_boundaries per site", {
  r <- bbqr:::.bbqr_resolve_clip(c(delta = FALSE), use_boundaries = TRUE)
  expect_false(r[["delta"]])
  ## The baseline is the gated eight, so everything else gated is on and the v4
  ## guards stay off -- this is "v7's clamp set minus delta", the leave-one-out
  ## arm, not "everything minus delta".
  expect_identical(names(which(r)), setdiff(bbqr:::.BBQR_CLIP_GATED, "delta"))

  ## A v4 guard can still be named explicitly, on top of that baseline.
  r3 <- bbqr:::.bbqr_resolve_clip(c(v4_z_psi = TRUE), use_boundaries = TRUE)
  expect_identical(names(which(r3)), c(bbqr:::.BBQR_CLIP_GATED, "v4_z_psi"))

  ## A character clip is absolute, not an override: it names exactly what is on.
  r2 <- bbqr:::.bbqr_resolve_clip("v4", use_boundaries = TRUE)
  expect_identical(unname(which(r2)),
                   unname(which(bbqr:::.BBQR_CLIP_SITES %in%
                                  bbqr:::.BBQR_CLIP_PRESETS$v4)))
})

test_that("bad clip specifications are rejected", {
  s <- sim_binary()
  expect_error(bbqr:::.bbqr_resolve_clip(c(nope = TRUE), FALSE), "unknown clip site")
  expect_error(bbqr:::.bbqr_resolve_clip(c(TRUE, FALSE), FALSE), "length 1 or 18")
  expect_error(bbqr:::.bbqr_resolve_clip(c(delta = NA), FALSE), "missing values")
  expect_error(bbqr:::.bbqr_clip_preset("no_such_preset"), "unknown clip site")
  ## The lasso and unpenalised kernels still take one all-or-nothing switch, so
  ## accepting `clip` there and dropping it would mislabel an ablation arm.
  expect_error(bbqr(s$formula, data = s$data, penalty = "lasso", ndraw = 100,
                    burn = 20, clip = c(delta = FALSE)),
               "only available for penalty")
})
