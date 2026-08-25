#' Lasso binary quantile regression
#'
#' Low-level sampler for binary quantile regression with a standard lasso
#' penalty: one shrinkage parameter shared by all slopes. This implements the
#' hierarchical sampler of Benoit, Al-Hamzawi and Yu (2013). Most users should
#' call [bbqr()] with `penalty = "lasso"` instead; this function is exported
#' for direct access to the sampler.
#'
#' The default `anchor = "sigma1"` holds the ALD inverse-scale at
#' `sigma_fixed`; `anchor = "free"` samples it from its full conditional,
#' which is what the published algorithm does. Which hierarchy is sampled --
#' the published one, whose posterior does not exist, or the corrected
#' `"v5"` default -- is set by `model` on the [prior()] object.
#'
#' @inheritParams bbqr
#' @return An object of class `"bbqr"`, as documented in [bbqr()].
#' @references
#' Benoit, D. F., Al-Hamzawi, R. and Yu, K. (2013). Bayesian lasso binary
#' quantile regression. *Computational Statistics*, 28(6), 2861--2873.
#' \doi{10.1007/s00180-013-0407-8}
#' @seealso [bbqr()], and [cbqr_lasso()] for the same penalty layer
#'   fitted to an observed continuous response
#' @examples
#' set.seed(1)
#' n <- 150
#' X <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
#' y <- as.numeric(X %*% c(1, -1, 0) + rnorm(n) > 0)
#' d <- data.frame(y = y, X)
#' fit <- bbqr_lasso(y ~ x1 + x2 + x3, data = d, ndraw = 400, burn = 100)
#' coef(fit)
#' @export
bbqr_lasso <- function(formula, data = NULL, quantile = 0.5,
                       ndraw = 10000, burn = NULL, keep = 1,
                       prior = NULL,
                       anchor = c("sigma1", "free", "beta1", "norm1",
                                  "normslopes"),
                       sigma_fixed = 1, standardize = FALSE,
                       use_boundaries = FALSE, beta_init = NULL,
                       beta0_init = 0) {
  anchor <- match.arg(anchor)
  md <- .bbqr_model_data(formula, data)
  .bbqr_check_mcmc(quantile, ndraw, keep)
  burn <- .bbqr_check_burn(burn, ndraw)
  pr <- .bbqr_resolve_prior(prior, "lasso")
  sw <- .bbqr_anchor_switches(anchor, "lasso")

  std <- .bbqr_maybe_standardize(md, standardize)
  binit <- .bbqr_init_beta(beta_init, md$names, md$nvar)
  out_rows <- ndraw / keep

  ## Resolve the tau-calibrated reference intercept prior, and carry it as a
  ## PRECISION so that the flat prior is the exact point b0_prec = 0.
  b0_mean <- pr$b0_mean
  b0_var  <- pr$b0_var
  if (is.na(b0_var) || is.na(b0_mean)) {
    ref <- .bbqr_b0_reference(quantile)
    if (is.na(b0_mean)) b0_mean <- ref$mean
    if (is.na(b0_var))  b0_var  <- ref$var
  }
  b0_prec <- if (is.infinite(b0_var)) 0 else 1 / b0_var

  ## Derivation fidelity. The projection anchors rescale the coefficient vector
  ## AFTER it is drawn. That map leaves the likelihood exactly invariant --
  ## beta -> beta/c, beta0 -> beta0/c, sigma -> sigma*c preserves sigma*beta --
  ## but it does NOT leave the prior invariant, and under a shrinkage prior the
  ## penalty depends on |beta_j|, so rescaling mid-sweep changes the effective
  ## penalty every iteration. The chain then targets no written-down model.
  ## This is a warning rather than an error because the anchor comparison is a
  ## deliberate part of the study; the departure is recorded on the fit as
  ## `derived = FALSE` so no row can be mistaken for a derived posterior.
  derived <- TRUE
  if (!identical(pr$model, "v3") && !anchor %in% c("free", "sigma1")) {
    derived <- FALSE
    warning("model = ", sQuote(pr$model), " with anchor = ", sQuote(anchor),
            ": the post-draw projection is an exact reparameterisation of the ",
            "likelihood but not of the prior, so this chain does not target a ",
            "derived posterior. Reported as derived = FALSE.", call. = FALSE)
  }

  fn <- .Fortran("qrb_l_mcmc",
    as.integer(md$n), as.integer(md$nvar), as.integer(ndraw), as.integer(keep),
    as.integer(md$y), as.double(quantile), as.double(std$X),
    as.double(pr$a), as.double(pr$b), as.double(pr$delta_init),
    as.double(pr$tau_init), as.double(pr$mh_sd0), as.double(pr$target_acc),
    as.integer(burn), as.logical(use_boundaries),
    as.logical(sw$fix_sigma), as.double(sigma_fixed),
    as.logical(sw$constrain_beta_norm), as.logical(sw$fix_beta1),
    as.logical(sw$norm_slopes_only),
    as.double(pr$alpha_delta), as.double(pr$beta_delta),
    as.double(pr$alpha_tau), as.double(pr$beta_tau),
    as.double(b0_mean), as.double(b0_prec),
    betadraw     = double(out_rows * md$nvar),
    beta0draw    = double(out_rows),
    sigmadraw    = double(out_rows),
    lambda2draw  = double(out_rows),
    tauhyperdraw = double(out_rows),
    deltadraw    = double(out_rows),
    accdraw      = double(out_rows),
    beta_init    = binit,
    beta0_init   = as.double(beta0_init),
    PACKAGE = "bbqr")

  out <- .bbqr_assemble(fn, md, std$scaler, quantile, "lasso", anchor,
                        out_rows, ndraw, burn, keep, pr,
                        extra = list(
                          lambdasq  = fn$lambda2draw,
                          tauhyper  = fn$tauhyperdraw,
                          delta     = fn$deltadraw))
  ## Provenance: which hierarchy was sampled, whether it is the derived
  ## target, and the intercept prior actually used (tau-dependent under the
  ## reference calibration, so it is recorded rather than reconstructed).
  out$model    <- pr$model
  out$derived  <- derived
  out$b0_prior <- c(mean = b0_mean, var = b0_var, prec = b0_prec)
  out$tau_prior <- c(shape = pr$alpha_tau, rate = pr$beta_tau)
  out
}
