## The anchor x penalty grid. Every combination is offered except
## penalty = "none" with anchor = "free", which is not identified.
##
## Under the default model = "v5" the three projection anchors warn that the
## chain does not target a derived posterior, so the grid is run through
## `projected()`, which captures that warning where it is expected and lets
## any other warning through.
projected <- function(expr, anchor) {
  if (anchor %in% c("beta1", "norm1", "normslopes"))
    expect_warning(out <- expr, "derived posterior")
  else
    out <- expr
  out
}

test_that("every offered anchor/penalty combination runs and is finite", {
  s <- sim_binary()
  for (p in c("none", "lasso", "alasso")) {
    for (a in .bbqr_allowed_anchors(p)) {
      fit <- projected(bbqr(s$formula, s$data, penalty = p, anchor = a,
                            ndraw = 600, burn = 100), a)
      expect_s3_class(fit, "bbqr")
      expect_equal(fit$anchor, a)
      expect_true(all(is.finite(fit$beta)),
                  label = sprintf("finite beta for %s/%s", p, a))
      expect_true(all(is.finite(fit$beta0)))
      expect_true(all(is.finite(fit$sigma)))
      ## derived = TRUE exactly for the two conventions the derivations cover.
      expect_equal(fit$derived, a %in% c("sigma1", "free"),
                   info = sprintf("%s/%s", p, a))
    }
  }
})

test_that("penalty = 'none' with anchor = 'free' is refused, not silently run", {
  s <- sim_binary()
  expect_false("free" %in% .bbqr_allowed_anchors("none"))
  expect_error(bbqr(s$formula, s$data, penalty = "none", anchor = "free",
                    ndraw = 300, burn = 50))
  expect_error(bbqr_none(s$formula, s$data, anchor = "free", ndraw = 300))
})

test_that("each anchor imposes its constraint under every penalty", {
  s <- sim_binary()
  for (p in c("none", "lasso", "alasso")) {
    ok <- .bbqr_allowed_anchors(p)

    if ("sigma1" %in% ok) {
      f1 <- bbqr(s$formula, s$data, penalty = p, anchor = "sigma1",
                 ndraw = 600, burn = 100)
      expect_true(all(f1$sigma == 1), label = sprintf("sigma1 under %s", p))
    }
    if ("beta1" %in% ok) {
      f2 <- projected(bbqr(s$formula, s$data, penalty = p, anchor = "beta1",
                           ndraw = 600, burn = 100), "beta1")
      expect_true(all(f2$beta[, 1L] == 1), label = sprintf("beta1 under %s", p))
    }
    if ("norm1" %in% ok) {
      f3 <- projected(bbqr(s$formula, s$data, penalty = p, anchor = "norm1",
                           ndraw = 600, burn = 100), "norm1")
      expect_true(all(abs(sqrt(f3$beta0^2 + rowSums(f3$beta^2)) - 1) < 1e-8),
                  label = sprintf("norm1 under %s", p))
    }
    if ("normslopes" %in% ok) {
      f4 <- projected(bbqr(s$formula, s$data, penalty = p,
                           anchor = "normslopes", ndraw = 600, burn = 100),
                      "normslopes")
      expect_true(all(abs(sqrt(rowSums(f4$beta^2)) - 1) < 1e-8),
                  label = sprintf("normslopes under %s", p))
    }
  }
})

test_that("the norm anchors let sigma absorb the rescaling", {
  s <- sim_binary()
  # Under sigma1 the scale is pinned; under a norm anchor it must move, because
  # the projection pushes the scale factor into it.
  pin <- bbqr(s$formula, s$data, penalty = "none", anchor = "sigma1",
              ndraw = 600, burn = 100)
  mov <- projected(bbqr(s$formula, s$data, penalty = "none", anchor = "norm1",
                        ndraw = 600, burn = 100), "norm1")
  expect_true(all(pin$sigma == 1))
  expect_gt(stats::sd(mov$sigma), 0)
})

test_that("the projection preserves the linear predictor it rescales", {
  # (beta0, beta, sigma) -> (c beta0, c beta, sigma/c) leaves sign(beta0 + x'b)
  # alone, so a norm anchor must not change which side of the threshold a point
  # falls on relative to its own fit.
  s <- sim_binary(n = 300, seed = 41)
  fit <- projected(bbqr(s$formula, s$data, penalty = "none", anchor = "norm1",
                        ndraw = 1500, burn = 300), "norm1")
  eta <- predict(fit, newdata = s$data, type = "link")
  expect_true(all(is.finite(eta)))
  # classification from a unit-norm fit still tracks the response
  expect_gt(mean((eta > 0) == (s$data$y == 1)), 0.6)
})
