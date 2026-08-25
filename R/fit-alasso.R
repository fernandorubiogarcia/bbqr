#' Adaptive-lasso binary quantile regression
#'
#' Low-level sampler for binary quantile regression with an adaptive-lasso
#' penalty, in which each slope receives its own shrinkage parameter. Most
#' users should call [bbqr()] with `penalty = "alasso"` instead; this function
#' is exported for direct access to the sampler.
#' With the defaults -- `model = "v5"` on the prior, `anchor = "sigma1"`, no
#' clamps -- every transition implements the corresponding full conditional
#' in the V5 appendix: the ALD inverse-scale is held at `sigma_fixed`, the
#' `omega` draw is `Gamma(alpha_omega + k delta, beta_omega + sum(1 /
#' lambdasq))`, the intercept carries its \eqn{\tau}-calibrated normal prior,
#' and `delta` is updated by a random-walk Metropolis step on the log scale
#' under a proper `Gamma(alpha_delta, beta_delta)` prior. That last prior is
#' required rather than optional: see [prior()]. `anchor = "free"` samples
#' the inverse-scale from its full conditional instead, which the appendix
#' also covers.
#'
#' @inheritParams bbqr
#' @param clip Which numerical clamps to apply, for isolating which individual
#'   clamp a result depends on. There are eighteen sites: the eight the single
#'   `use_boundaries` flag used to gate together (`"ystar"`, `"z"`, `"s"`,
#'   `"sigma"`, `"lambdasq"`, `"omega"`, `"delta"`, `"mh_sd"`) and ten that
#'   restore guards the Experiment v4 build applied unconditionally and later
#'   builds removed
#'   (`"v4_sd"`, `"v4_z_psi"`, `"v4_z_chi"`, `"v4_zbound"`, `"v4_s_lam"`,
#'   `"v4_s_psi"`, `"v4_s_chi"`, `"v4_sbound"`, `"v4_rate_sig"`,
#'   `"v4_rate_om"`).
#'
#'   `use_boundaries` addresses the first eight sites and only those: the v4
#'   guards were applied unconditionally rather than gated by it, so
#'   `use_boundaries = TRUE` does not switch them on. Reaching them requires an
#'   explicit `clip`.
#'
#'   Accepts three forms. A **named logical** vector overrides individual sites
#'   on top of `use_boundaries`, so `clip = c(delta = FALSE)` with
#'   `use_boundaries = TRUE` means every gated clamp except `delta`. A
#'   **character** vector names exactly the sites that are on and ignores
#'   `use_boundaries` entirely, which is what makes a named arm reproducible;
#'   the preset names `"none"` (nothing; `"v6"` is a synonym naming the build
#'   that first shipped it), `"v7clip"` (the eight gated clamps) and `"v4"`
#'   (the ten v4 guards, gated clamps off) are accepted here. An **unnamed
#'   logical** of length 1 or 18 is taken in site order; a length-1 value
#'   behaves exactly like `use_boundaries`, so it too reaches only the eight
#'   gated sites. `NULL`, the default, is `use_boundaries` throughout.
#'
#'   `clip = "v4"` reproduces v4's clamp configuration, not a v4 rerun: v4 also
#'   drew `z` and `s` with a hand-rolled inverse-Gaussian rather than the
#'   current GIG(1/2) routine, and ran an effectively flat `delta` prior. Those
#'   are a sampler and a prior, not clamps.
#' @param sampler Which routine draws the latents `z` and `s`. Both full
#'   conditionals are GIG(1/2) and two routines implement that draw, which are
#'   not numerically equivalent:
#'
#'   \describe{
#'     \item{`"gig_half"`}{The default. `rgig_half`, whose `q` is rationalized
#'       to avoid cancellation, with a limiting branch at `chi == 0`.}
#'     \item{`"v4_invgaus"`}{The routine the v4 build shipped: draw the
#'       reciprocal of the latent as inverse Gaussian from the textbook `q`
#'       expression and invert. That expression subtracts two nearly equal
#'       quantities when `chi` is small, so it can return a latent at zero or
#'       `Inf` -- the failure v4's guards were written to contain.}
#'   }
#'
#'   The switch exists so that the v4-to-v6 difference can be attributed. `clip`
#'   settles the guards and `sampler` settles the routine; holding one fixed
#'   while moving the other is what separates the two mechanisms. Both routines
#'   consume one normal and one uniform per draw, so two fits from the same seed
#'   that differ only in `sampler` stay aligned in the RNG stream and their
#'   paired difference is the routine alone. Reproducing v4's z and s draws
#'   takes `sampler = "v4_invgaus"` together with `clip = "v4"`, since v4 applied
#'   those argument floors unconditionally; `"v4_invgaus"` with the floors off is
#'   accepted but is not a configuration v4 ever ran.
#' @param kappa_lam Coefficient of an \eqn{\exp(-\kappa_\lambda / \lambda_j^2)}
#'   tilt on each \eqn{\lambda_j^2} conditional. The tilt is conjugate, so it
#'   only shifts the inverse-Gamma scale to
#'   \eqn{\sigma s_j / 2 + \omega + \kappa_\lambda}; it was tried as a barrier
#'   keeping \eqn{\lambda_j^2} away from zero. It is not part of any derived
#'   target and does not repair the \eqn{\omega} defect -- the normalising
#'   constant travels with the scale, which removes the very factor that made
#'   the origin integrable -- so it must be `0` (the default) under
#'   `model = "v4"` and `"v5"`. Under `model = "v3"` a positive value is
#'   accepted, for reproducing the runs that used it.
#' @section Derivation version:
#' Which of the three adaptive-lasso hierarchies is sampled is set by
#' `model` on the [prior()] object, not by an argument here. `"v4"` and `"v5"`
#' are derivation-faithful: their appendices state that no clamp, no
#' sampler-argument floor, and no `exp(-kappa_lam / lambda_j^2)` tilt is part
#' of the target, so a fit that names them while enabling one of those devices
#' is refused rather than recorded. Concretely, under `model = "v4"` or
#' `"v5"` this function requires `clip = "none"` (the default), `kappa_lam = 0`
#' and `sampler = "gig_half"`. `model = "v3"` leaves all three switches
#' reachable, because reproducing what the earlier builds actually ran needs
#' them.
#'
#' One further check is a warning rather than an error. Under `model = "v4"`
#' with a free `sigma`, a flat intercept is admissible only if `a > 1`, and the
#' historical default `a = 1` sits exactly on the divergent boundary. That
#' configuration is what Experiments v4 and v6 ran, so it stays available and
#' warns instead of failing.
#'
#' @return An object of class `"bbqr"`, as documented in [bbqr()], with two
#'   extra components: `clip`, the resolved logical(18) of switches, and
#'   `clip_hits`, the number of times each site's clamp actually bound.
#'   `clip_hits` is counted whether or not the switch was on, so a fit with
#'   every clamp off still reports how often each clamp would have fired.
#'   Denominators differ by site. `ystar`, `z`, `v4_zarg` and `v4_zbound` are
#'   out of `ndraw * n`; `s`, `lambdasq`, `v4_sarg` and `v4_sbound` out of
#'   `ndraw * p`; `mh_sd` out of `burn`; `delta` out of the accepted Metropolis
#'   proposals; `omega` and `v4_rate_om` out of `ndraw`. Multiply by the number
#'   of guards in a group where it holds more than one: `v4_zarg` has two,
#'   `v4_sarg` three, each rate group two. `v4_sd` is out of `n` alone, since
#'   v4 floored that standard deviation in the initialisation sweep only and
#'   left the main-loop assignment bare. `sigma` and `v4_rate_sig` stay at zero
#'   under `anchor = "sigma1"`, where nothing ever reaches them.
#' @references
#' Rubio Garcia, F. (2023). Bayesian adaptive lasso binary quantile
#' regression with hybrid resampling for classification of imbalanced data.
#' M.S. thesis, Wichita State University, Dept. of Mathematics, Statistics,
#' and Physics.
#' \url{https://soar.wichita.edu/entities/publication/a2f86232-4704-4ec2-b685-751e7b04ec42}
#' @seealso [bbqr()], and [cbqr_alasso()] for the same penalty layer
#'   fitted to an observed continuous response
#' @examples
#' set.seed(1)
#' n <- 150
#' X <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
#' y <- as.numeric(X %*% c(1, -1, 0) + rnorm(n) > 0)
#' d <- data.frame(y = y, X)
#' fit <- bbqr_alasso(y ~ x1 + x2 + x3, data = d, ndraw = 400, burn = 100)
#' coef(fit)
#' @export
bbqr_alasso <- function(formula, data = NULL, quantile = 0.5,
                        ndraw = 10000, burn = NULL, keep = 1,
                        prior = NULL,
                        anchor = c("sigma1", "free", "beta1", "norm1",
                                   "normslopes"),
                        sigma_fixed = 1, standardize = FALSE,
                        use_boundaries = FALSE, clip = NULL,
                        sampler = c("gig_half", "v4_invgaus"),
                        kappa_lam = 0,
                        beta_init = NULL, beta0_init = 0) {
  anchor <- match.arg(anchor)
  md <- .bbqr_model_data(formula, data)
  .bbqr_check_mcmc(quantile, ndraw, keep)
  burn <- .bbqr_check_burn(burn, ndraw)
  pr <- .bbqr_resolve_prior(prior, "alasso")
  sw <- .bbqr_anchor_switches(anchor, "alasso")
  clipv <- .bbqr_resolve_clip(clip, use_boundaries)
  sampler <- match.arg(sampler)
  zs <- .bbqr_resolve_sampler(sampler)

  ## Resolve the tau-calibrated reference intercept prior, and convert the
  ## variance to a precision so that the flat intercept is the exact point
  ## b0_prec = 0 rather than a large finite variance.
  b0_mean <- pr$b0_mean
  b0_var  <- pr$b0_var
  if (is.na(b0_var) || is.na(b0_mean)) {
    ref <- .bbqr_b0_reference(quantile)
    if (is.na(b0_mean)) b0_mean <- ref$mean
    if (is.na(b0_var))  b0_var  <- ref$var
  }
  b0_prec <- if (is.infinite(b0_var)) 0 else 1 / b0_var

  ## Derivation fidelity. V4 and V5 are defined by their appendices, which
  ## state that no clamp, argument floor, or exp(-kappa/lambda^2) tilt is part
  ## of the target. Running them with such a device on would sample a
  ## different posterior under the name of a derived one, so it is refused
  ## rather than silently recorded. V3 keeps every switch addressable, because
  ## reproducing what Experiments v4 and v6 actually ran requires them.
  if (!identical(pr$model, "v3")) {
    on_sites <- names(clipv)[clipv]
    if (length(on_sites))
      stop("model = ", sQuote(pr$model), " is a derivation-faithful ",
           "hierarchy and admits no numerical clamps, but these are on: ",
           paste(on_sites, collapse = ", "),
           ". Pass clip = \"none\" (the default), or use model = \"v3\" to ",
           "reproduce a guarded historical configuration.", call. = FALSE)
    if (!isTRUE(all.equal(as.double(kappa_lam), 0)))
      stop("model = ", sQuote(pr$model), " has kappa_lam = 0 by construction; ",
           "the exp(-kappa_lam / lambda_j^2) tilt is not part of the derived ",
           "target. Got kappa_lam = ", kappa_lam, ".", call. = FALSE)
    if (!identical(sampler, "gig_half"))
      stop("model = ", sQuote(pr$model), " requires sampler = \"gig_half\"; ",
           "\"v4_invgaus\" evaluates the same GIG(1/2) root in a form that ",
           "loses precision as beta_j -> 0 and is retained only to reproduce ",
           "the historical build under model = \"v3\".", call. = FALSE)
    ## The derivations cover exactly two identification conventions: sigma
    ## sampled from its full conditional, and sigma held at a fixed positive
    ## value. The remaining anchors rescale the coefficient vector AFTER it is
    ## drawn -- a projection onto the unit sphere, or a division pinning
    ## beta_1 = 1 -- and no such step appears in any of them. Under a proper
    ## intercept prior the conflict is direct: the projection moves beta0 to a
    ## value the prior was not consulted about, so the chain no longer targets
    ## the posterior the prior defines.
  }

  ## The projection anchors rescale the coefficient vector AFTER it is drawn.
  ## The map leaves the likelihood exactly invariant -- beta -> beta/c,
  ## beta0 -> beta0/c, sigma -> sigma*c preserves sigma*beta -- but not the
  ## prior, and under a shrinkage prior the penalty depends on |beta_j|, so
  ## rescaling mid-sweep changes the effective penalty every iteration and the
  ## chain targets no written-down model. Warned rather than refused, because
  ## the anchor comparison is a deliberate part of the study; the departure is
  ## recorded as `derived = FALSE` so no row reads as a derived posterior.
  derived <- TRUE
  if (!identical(pr$model, "v3") && !anchor %in% c("free", "sigma1")) {
    derived <- FALSE
    warning("model = ", sQuote(pr$model), " with anchor = ", sQuote(anchor),
            ": the post-draw projection is an exact reparameterisation of the ",
            "likelihood but not of the prior, so this chain does not target a ",
            "derived posterior. Reported as derived = FALSE.", call. = FALSE)
  }

  ## Propriety of the flat intercept under a free sigma needs a > 1; the
  ## historical default a = 1 sits exactly on the divergent boundary. This is
  ## a warning rather than an error because reproducing that configuration is
  ## a legitimate thing to want -- it is what Experiments v4 and v6 ran.
  if (identical(pr$model, "v4") && !sw$fix_sigma && pr$a <= 1)
    warning("model = \"v4\" with anchor = ", sQuote(anchor), " and a = ",
            pr$a, ": a flat intercept against a free sigma requires a > 1, ",
            "so this posterior is improper. The v4 reference is a = b = 2.",
            call. = FALSE)

  std <- .bbqr_maybe_standardize(md, standardize)
  binit <- .bbqr_init_beta(beta_init, md$names, md$nvar)
  out_rows <- ndraw / keep

  fn <- .Fortran("qrb_al_mcmc",
    as.integer(md$n), as.integer(md$nvar), as.integer(ndraw), as.integer(keep),
    as.integer(md$y), as.double(quantile), as.double(std$X),
    as.double(pr$a), as.double(pr$b), as.double(pr$delta_init),
    as.double(pr$omega_init), as.double(pr$mh_sd0), as.double(pr$target_acc),
    as.integer(burn), as.logical(clipv), as.integer(zs),
    as.logical(sw$fix_sigma), as.double(sigma_fixed),
    as.logical(sw$constrain_beta_norm), as.logical(sw$fix_beta1),
    as.logical(sw$norm_slopes_only),
    as.double(pr$alpha_delta), as.double(pr$beta_delta),
    as.double(kappa_lam),
    as.double(pr$alpha_omega), as.double(pr$beta_omega),
    as.double(b0_mean), as.double(b0_prec),
    betadraw     = double(out_rows * md$nvar),
    beta0draw    = double(out_rows),
    sigmadraw    = double(out_rows),
    lambdasqdraw = double(out_rows * md$nvar),
    omegadraw    = double(out_rows),
    deltadraw    = double(out_rows),
    accdraw      = double(out_rows),
    cliphits     = double(length(clipv)),
    beta_init    = binit,
    beta0_init   = as.double(beta0_init),
    PACKAGE = "bbqr")

  out <- .bbqr_assemble(fn, md, std$scaler, quantile, "alasso", anchor,
                        out_rows, ndraw, burn, keep, pr,
                        extra = list(
                          lambdasq = matrix(fn$lambdasqdraw, out_rows, md$nvar,
                                            dimnames = list(NULL, md$names)),
                          omega    = fn$omegadraw,
                          delta    = fn$deltadraw))
  ## Not per-draw quantities, so they bypass `extra` (which trims by
  ## retained-draw index) and are attached directly.
  out$clip      <- clipv
  out$clip_hits <- stats::setNames(fn$cliphits, names(clipv))
  ## Provenance, not a tuning record: the two routines are different samplers
  ## for the same conditional, so a result is not interpretable without it.
  out$sampler   <- sampler
  ## Likewise the derivation version and the resolved intercept prior.  The
  ## intercept prior is tau-dependent under the reference calibration, so the
  ## values that were actually used are recorded rather than reconstructed.
  out$model     <- pr$model
  out$derived   <- derived
  out$b0_prior  <- c(mean = b0_mean, var = b0_var, prec = b0_prec)
  out$omega_prior <- c(shape = pr$alpha_omega, rate = pr$beta_omega)
  out
}
