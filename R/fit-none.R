#' Unpenalised binary quantile regression
#'
#' Low-level sampler for binary quantile regression with a vague normal prior
#' on the slopes and no shrinkage layer, implementing Benoit and Van den Poel
#' (2012). Most users should call [bbqr()] with `penalty = "none"` instead.
#'
#' The default `anchor = "sigma1"` reproduces the published algorithm, which
#' holds the ALD scale at 1. `"free"` is deliberately absent from the choices
#' here, unlike the penalised fitters: with neither a penalty nor an anchor the
#' likelihood is flat along the scale ray, so the model is not identified and
#' only the vague prior stops the chain drifting.
#'
#' @inheritParams bbqr
#' @return An object of class `"bbqr"`, as documented in [bbqr()].
#' @references
#' Benoit, D. F. and Van den Poel, D. (2012). Binary quantile regression: a
#' Bayesian approach based on the asymmetric Laplace distribution.
#' *Journal of Applied Econometrics*, 27(7), 1174--1188.
#' \doi{10.1002/jae.1216}
#' @param burn Number of initial iterations to discard. Defaults to
#'   `ndraw / 10`. This sampler has no penalty hyperparameter and no
#'   Metropolis step, so nothing adapts over it.
#' @seealso [bbqr()], and [cbqr_none()] for the same penalty layer
#'   fitted to an observed continuous response
#' @examples
#' set.seed(1)
#' n <- 150
#' X <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
#' y <- as.numeric(X %*% c(1, -1, 0) + rnorm(n) > 0)
#' d <- data.frame(y = y, X)
#' fit <- bbqr_none(y ~ x1 + x2 + x3, data = d, ndraw = 400, burn = 100)
#' coef(fit)
#' @export
bbqr_none <- function(formula, data = NULL, quantile = 0.5,
                      ndraw = 10000, burn = NULL, keep = 1,
                      prior = NULL,
                      anchor = c("sigma1", "beta1", "norm1", "normslopes"),
                      standardize = FALSE,
                      use_boundaries = FALSE, beta_init = NULL,
                      beta0_init = 0) {
  anchor <- match.arg(anchor)
  md <- .bbqr_model_data(formula, data)
  .bbqr_check_mcmc(quantile, ndraw, keep)
  burn <- .bbqr_check_burn(burn, ndraw)
  pr <- .bbqr_resolve_prior(prior, "none")
  sw <- .bbqr_anchor_switches(anchor, "none")

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

  fn <- .Fortran("qrb_bqr_mcmc",
    as.integer(md$n), as.integer(md$nvar), as.integer(ndraw), as.integer(keep),
    as.integer(md$y), as.double(quantile), as.double(std$X),
    as.double(pr$beta_var), as.logical(use_boundaries),
    as.logical(sw$constrain_beta_norm), as.logical(sw$fix_beta1),
    as.logical(sw$norm_slopes_only),
    as.double(b0_mean), as.double(b0_prec),
    betadraw   = double(out_rows * md$nvar),
    beta0draw  = double(out_rows),
    sigmadraw  = double(out_rows),
    beta_init  = binit,
    beta0_init = as.double(beta0_init),
    PACKAGE = "bbqr")

  ## This sampler has no Metropolis step, so there is no acceptance rate.
  fn$accdraw <- rep(NA_real_, out_rows)

  out <- .bbqr_assemble(fn, md, std$scaler, quantile, "none", anchor,
                        out_rows, ndraw, burn, keep, pr, extra = list())
  ## Provenance: which hierarchy was sampled, whether it is the derived
  ## target, and the intercept prior actually used (tau-dependent under the
  ## reference calibration, so it is recorded rather than reconstructed).
  out$model    <- pr$model
  out$derived  <- derived
  out$b0_prior <- c(mean = b0_mean, var = b0_var, prec = b0_prec)
  out
}
