## Internal helpers shared by the three samplers. Not exported.

## Build the design matrix and response from a formula, splitting the intercept
## off from the slopes. The Fortran kernels treat the intercept separately (it
## is never penalised), so it must not sit in X.
.bbqr_model_data <- function(formula, data) {
  if (missing(formula) || is.null(formula))
    stop("'formula' has to be specified.", call. = FALSE)
  mf    <- stats::model.frame(formula = formula, data = data)
  Xfull <- stats::model.matrix(attr(mf, "terms"), data = mf)
  y     <- stats::model.response(mf)

  int_idx <- which(colnames(Xfull) == "(Intercept)")
  X <- if (length(int_idx)) Xfull[, -int_idx, drop = FALSE] else Xfull

  ylev <- NULL
  if (is.factor(y)) {
    if (nlevels(y) != 2L)
      stop("A factor response must have exactly two levels.", call. = FALSE)
    ylev <- levels(y)
    y <- as.integer(y) - 1L
  } else if (is.logical(y)) {
    y <- as.integer(y)
  }
  if (!all(y %in% c(0, 1)))
    stop("The response must be binary: 0/1, logical, or a two-level factor.",
         call. = FALSE)
  if (length(unique(y)) < 2L)
    stop("The response takes only one value; there is nothing to fit.",
         call. = FALSE)
  if (ncol(X) < 1L)
    stop("At least one predictor is required.", call. = FALSE)
  if (anyNA(X) || anyNA(y))
    stop("Missing values are not supported; remove or impute them first.",
         call. = FALSE)

  list(y = as.integer(y), X = X, names = colnames(X), ylevels = ylev,
       n = length(y), nvar = ncol(X), terms = attr(mf, "terms"))
}

## Shared argument checking for the MCMC controls.
.bbqr_check_mcmc <- function(quantile, ndraw, keep) {
  if (!is.numeric(quantile) || length(quantile) != 1L || !is.finite(quantile) ||
      quantile <= 0 || quantile >= 1)
    stop("'quantile' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  if (is.null(ndraw))
    stop("Number of MCMC draws ('ndraw') has to be specified.", call. = FALSE)
  if (!is.numeric(ndraw) || length(ndraw) != 1L || ndraw < 1)
    stop("'ndraw' must be a single positive integer.", call. = FALSE)
  if (!is.numeric(keep) || length(keep) != 1L || keep < 1)
    stop("'keep' must be a single positive integer.", call. = FALSE)
  if (keep > ndraw)
    stop("'keep' cannot be larger than 'ndraw'.", call. = FALSE)
  if (ndraw %% keep != 0)
    stop(sprintf("'ndraw' (%d) must be exactly divisible by 'keep' (%d).",
                 as.integer(ndraw), as.integer(keep)), call. = FALSE)
  invisible(TRUE)
}

## Default burn-in is a tenth of the chain. Must leave at least one draw.
.bbqr_check_burn <- function(burn, ndraw) {
  if (is.null(burn)) burn <- floor(ndraw / 10)
  if (!is.numeric(burn) || length(burn) != 1L || burn < 0 || burn >= ndraw)
    stop("'burn' must be a single number with 0 <= burn < ndraw.", call. = FALSE)
  as.integer(burn)
}

## Accept a user prior, or build the default for this penalty. Kept in a helper
## so the fitters can name their argument 'prior' without shadowing prior().
.bbqr_resolve_prior <- function(p, penalty) {
  if (is.null(p)) return(prior(penalty))
  if (!inherits(p, "bbqr.prior"))
    stop("'prior' must be an object created by prior().", call. = FALSE)
  if (p$penalty != penalty)
    stop("'prior' was built for penalty = ", sQuote(p$penalty),
         ", but this sampler is ", sQuote(penalty), ".", call. = FALSE)
  ## Catch serialized or hand-edited prior objects that would make the
  ## adaptive-lasso posterior improper. prior() already enforces this; the check
  ## is repeated here because a bbqr.prior can be saved, modified and reloaded.
  if (penalty == "alasso" && !(is.numeric(p$beta_delta) &&
                               length(p$beta_delta) == 1L &&
                               is.finite(p$beta_delta) && p$beta_delta > 0))
    stop("This prior has beta_delta <= 0, which leaves the adaptive-lasso ",
         "posterior improper; rebuild it with prior(\"alasso\").", call. = FALSE)

  ## A bbqr.prior serialized before the V3/V4/V5 distinction existed has none
  ## of the fields below. Such an object described the V3 hierarchy, so fill
  ## the V3 values in rather than erroring: the numbers reproduce what that
  ## object actually meant.
  ## Fields common to all three penalties.
  if (is.null(p$model))   p$model   <- "v3"
  if (is.null(p$b0_mean)) p$b0_mean <- 0
  if (is.null(p$b0_var))  p$b0_var  <- Inf
  if (is.null(p$alpha_tau)) p$alpha_tau <- 0
  if (is.null(p$beta_tau))  p$beta_tau  <- 0
  if (!p$model %in% c("v3", "v4", "v5"))
    stop("This prior has model = ", sQuote(as.character(p$model)),
         "; expected one of \"v3\", \"v4\", \"v5\".", call. = FALSE)
  if (p$model == "v5" && identical(p$b0_var, Inf))
    stop("This prior names model = \"v5\" but has a flat intercept; ",
         "rebuild it with prior(", encodeString(penalty, quote = "\""),
         ", model = \"v5\").", call. = FALSE)
  if (penalty == "lasso" && p$model == "v5" &&
      !(p$alpha_tau > 0 && p$beta_tau > 0))
    stop("This prior names model = \"v5\" for the lasso but carries ",
         "alpha_tau = ", p$alpha_tau, ", beta_tau = ", p$beta_tau,
         ", i.e. the published log-uniform prior on tau_h, under which the ",
         "joint posterior does not exist; rebuild it with ",
         "prior(\"lasso\", model = \"v5\").", call. = FALSE)

  if (penalty == "alasso") {
    if (is.null(p$alpha_omega)) p$alpha_omega <- 0
    if (is.null(p$beta_omega))  p$beta_omega  <- 0
    ## Same reasoning as the beta_delta check: a hand-edited prior could name a
    ## proper model while carrying V3's log-uniform prior on omega, which has
    ## no posterior.
    if (p$model != "v3" && !(p$alpha_omega > 0 && p$beta_omega > 0))
      stop("This prior names model = ", sQuote(p$model), " but carries ",
           "alpha_omega = ", p$alpha_omega, ", beta_omega = ", p$beta_omega,
           ", i.e. V3's log-uniform prior on omega, under which the joint ",
           "posterior does not exist; rebuild it with prior(\"alasso\", ",
           "model = ", encodeString(p$model, quote = "\""), ").",
           call. = FALSE)
  }
  p
}

## Centre and scale the columns of X, remembering the transformation so the
## draws can be mapped back to the original units.
.bbqr_maybe_standardize <- function(md, standardize) {
  if (!isTRUE(standardize)) return(list(X = md$X, scaler = NULL))
  nvar <- md$nvar
  x_mu <- stats::setNames(rep(0, nvar), md$names)
  x_sd <- stats::setNames(rep(1, nvar), md$names)
  X_work <- md$X
  for (j in seq_len(nvar)) {
    x_mu[j] <- mean(md$X[, j])
    sdj     <- stats::sd(md$X[, j])
    if (!is.finite(sdj) || isTRUE(all.equal(sdj, 0))) sdj <- 1
    x_sd[j] <- sdj
    X_work[, j] <- (md$X[, j] - x_mu[j]) / sdj
  }
  list(X = X_work, scaler = list(x_mu = x_mu, x_sd = x_sd))
}

## Map draws obtained on the standardized scale back to the original covariate
## units. On the standardized scale the linear predictor is
##   b0* + sum_j B*_j (x_j - mu_j) / sd_j,
## so the slopes divide by sd and the intercept absorbs the centring offset:
##   beta_j = B*_j / sd_j,   beta0 = b0* - sum_j beta_j mu_j.
## This is an exact identity: the linear predictor is unchanged.
.bbqr_backtransform <- function(B, b0, scaler) {
  if (is.null(scaler)) return(list(beta = B, beta0 = b0))
  B_orig <- sweep(B, 2L, scaler$x_sd, "/")
  list(beta  = B_orig,
       beta0 = b0 - as.vector(B_orig %*% scaler$x_mu))
}

## Align user-supplied starting values to the column order of X.
.bbqr_init_beta <- function(beta_init, names, nvar) {
  if (is.null(beta_init)) return(double(nvar))
  if (!is.null(names(beta_init))) {
    tmp <- stats::setNames(rep(0, nvar), names)
    common <- intersect(names, names(beta_init))
    tmp[common] <- beta_init[common]
    return(as.double(tmp))
  }
  if (length(beta_init) != nvar)
    stop(sprintf("'beta_init' has length %d but there are %d predictors.",
                 length(beta_init), nvar), call. = FALSE)
  as.double(beta_init)
}

## The eighteen numerical clamp sites in the adaptive-lasso kernel, in the order
## they are reached inside one sweep. This vector IS the contract with the
## Fortran layer: it must match the C_* index parameters at the top of
## QRb_AL_mcmc.f95 exactly, name for name and position for position.
##
## Sites 1-8 are the clamps that the single `use_boundaries` flag used to gate
## together. Sites 9-18, the `v4_*` group, restore guards that the v4 build
## applied unconditionally and that v5 deleted outright. v4 ran with
## use_boundaries = FALSE, so sites 1-8 were all OFF in v4 and sites 9-18 were
## all ON: that combination, and only that one, is v4's clamp configuration.
##
## Two sites are unreachable under anchor = "sigma1": `sigma`, and
## `v4_rate_sig`, whose floors sit inside the sigma draw. Their switches have
## no effect there and their counters stay at zero, so crossing them in an
## ablation design over that anchor only produces duplicate arms.
##
## Only the adaptive-lasso kernel takes per-site switches. The lasso and
## unpenalised kernels still take a single use_boundaries scalar, so `clip` is
## rejected for those penalties rather than silently ignored.
.BBQR_CLIP_SITES <- c("ystar", "z", "s", "sigma", "lambdasq", "omega", "delta",
                      "mh_sd", "v4_sd", "v4_z_psi", "v4_z_chi", "v4_zbound",
                      "v4_s_lam", "v4_s_psi", "v4_s_chi", "v4_sbound",
                      "v4_rate_sig", "v4_rate_om")

## Named clamp configurations worth referring to by name rather than by an
## eight- or twelve-vector of logicals.
##
##   v6      the current default: no clamping at all
##   v7clip  what the v7 `_clip` arm ran: the eight gated clamps, no v4 guards
##   v4      v4's clamp set: the seven v4 guards, gated clamps all off
##
## Note that `v4_sd` covers the initialisation sweep only. v4 floored the
## truncation standard deviation there and left the main-loop assignment bare,
## because its hard z bound already kept that quantity positive on every later
## sweep. Flooring both would be a stronger intervention than v4 applied.
##
## `v4` reproduces v4's CLAMPS. It is not a v4 rerun: v4 also drew z and s with
## a hand-rolled rinvgaus() rather than rgig_half(), and ran an effectively flat
## delta prior. Those are a sampler and a prior, not clamps. Note that v4's own
## source comments attribute these guards to rinvgaus cancellation returning
## z <= 0 or Inf -- a failure mode rgig_half was written to avoid -- so the
## guards may well never fire on the current sampler. The counters say whether
## they do.
## The eight sites the single `use_boundaries` flag gated. Kept as its own
## constant because it is the back-compatibility contract, not just a preset:
## use_boundaries only ever addresses these, whatever else the site table grows.
.BBQR_CLIP_GATED <- c("ystar", "z", "s", "sigma", "lambdasq", "omega", "delta",
                      "mh_sd")

.BBQR_CLIP_PRESETS <- list(
  ## "none" and "v6" are the same empty set. "v6" names the build that first
  ## shipped it; "none" says what it is, and is the name the derivation-faithful
  ## models v4 and v5 require.
  none   = character(0),
  v6     = character(0),
  v7clip = .BBQR_CLIP_GATED,
  v4     = c("v4_sd", "v4_z_psi", "v4_z_chi", "v4_zbound",
             "v4_s_lam", "v4_s_psi", "v4_s_chi", "v4_sbound",
             "v4_rate_sig", "v4_rate_om")
)

## ---- which routine draws the latents z and s -------------------------------
## Both full conditionals are GIG(1/2), and two routines implement that draw.
## They are not numerically equivalent, so which one ran belongs in a result's
## provenance rather than being an implementation detail.
##
##   "gig_half"     rgig_half: q rationalized to avoid cancellation, plus a
##                  chi == 0 limiting branch. The current default.
##   "v4_invgaus"   the routine the v4 build shipped: draw the reciprocal of the
##                  latent as inverse Gaussian from the textbook q formula and
##                  invert. The formula cancels for tiny chi, which is the
##                  failure v4's guards were written to contain.
##
## Kept as a named integer vector because the value crossing into Fortran is
## positional: the names are the R-side contract, the integers are ZS_GIG and
## ZS_V4INVGAUS in QRb_AL_mcmc.f95, and the two must stay in step.
##
## Only the adaptive-lasso kernel takes the switch. The lasso and unpenalised
## kernels still hard-code their own draw, so `sampler` is rejected for those
## penalties rather than silently ignored.
.BBQR_ZS_SAMPLERS <- c(gig_half = 0L, v4_invgaus = 1L)

.bbqr_resolve_sampler <- function(sampler) {
  if (is.null(sampler)) return(.BBQR_ZS_SAMPLERS[["gig_half"]])
  if (!is.character(sampler) || length(sampler) != 1L || is.na(sampler))
    stop("'sampler' must be a single string, one of ",
         paste(sQuote(names(.BBQR_ZS_SAMPLERS)), collapse = " or "), ".",
         call. = FALSE)
  if (!sampler %in% names(.BBQR_ZS_SAMPLERS))
    stop("unknown 'sampler' ", sQuote(sampler), "; use one of ",
         paste(sQuote(names(.BBQR_ZS_SAMPLERS)), collapse = ", "), ".",
         call. = FALSE)
  .BBQR_ZS_SAMPLERS[[sampler]]
}

## Expand a preset name, or a character vector of site names, into the logical
## vector the kernel expects. Anything not named is FALSE.
.bbqr_clip_preset <- function(x) {
  if (length(x) == 1L && x %in% names(.BBQR_CLIP_PRESETS))
    x <- .BBQR_CLIP_PRESETS[[x]]
  bad <- setdiff(x, .BBQR_CLIP_SITES)
  if (length(bad))
    stop("unknown clip site(s): ", paste(sQuote(bad), collapse = ", "), ".",
         call. = FALSE)
  stats::setNames(.BBQR_CLIP_SITES %in% x, .BBQR_CLIP_SITES)
}

## Resolve the user's `clip` argument into the logical(8) the kernel expects.
##
## `use_boundaries` supplies the baseline, so the historical two-state behaviour
## is unchanged: clip = NULL reproduces it exactly. A named `clip` overrides
## individual sites on top of that baseline, which is what makes an ablation arm
## readable -- clip = c(delta = FALSE) against use_boundaries = TRUE is "all
## clamps except delta".
##
## Note that two of the eight sites are unreachable under some anchors. "sigma"
## is only ever touched when sigma is sampled (anchor != "sigma1") or when a
## norm anchor rescales it, so under anchor = "sigma1" its switch has no effect
## and its counter stays at zero. That is reported, not corrected: an ablation
## design that crosses an inert switch produces duplicate arms.
.bbqr_resolve_clip <- function(clip, use_boundaries) {
  if (!is.logical(use_boundaries) || length(use_boundaries) != 1L ||
      is.na(use_boundaries))
    stop("'use_boundaries' must be a single TRUE or FALSE.", call. = FALSE)
  sites <- .BBQR_CLIP_SITES
  ## `use_boundaries` is the baseline for the eight sites it historically gated,
  ## and ONLY those. The v4 guards were never gated by it -- they were applied
  ## unconditionally -- so use_boundaries = TRUE must not switch them on. If it
  ## did, re-running the v7 driver against this build would silently produce a
  ## different arm from the one v7 reported.
  out <- stats::setNames(
    isTRUE(use_boundaries) & sites %in% .BBQR_CLIP_GATED, sites)
  if (is.null(clip)) return(out)

  ## A character `clip` is absolute, not an override: it names exactly the sites
  ## that are on, so clip = "v4" means v4's set and nothing else regardless of
  ## use_boundaries. That is what makes a named arm reproducible.
  if (is.character(clip)) return(.bbqr_clip_preset(clip))

  if (!is.logical(clip) || anyNA(clip))
    stop("'clip' must be a logical vector with no missing values.",
         call. = FALSE)
  if (!is.null(names(clip)) && any(nzchar(names(clip)))) {
    bad <- setdiff(names(clip), sites)
    if (length(bad))
      stop("unknown clip site(s): ", paste(sQuote(bad), collapse = ", "),
           ". Available sites are ", paste(sQuote(sites), collapse = ", "), ".",
           call. = FALSE)
    if (anyDuplicated(names(clip)))
      stop("'clip' names each site at most once; duplicated: ",
           paste(sQuote(unique(names(clip)[duplicated(names(clip))])),
                 collapse = ", "), ".", call. = FALSE)
    out[names(clip)] <- clip
    return(out)
  }
  ## A scalar means the same thing as use_boundaries, so it obeys the same
  ## restriction: it addresses the eight gated sites, not v4's guards.
  if (length(clip) == 1L)
    return(stats::setNames(clip & sites %in% .BBQR_CLIP_GATED, sites))
  if (length(clip) != length(sites))
    stop("an unnamed 'clip' must have length 1 or ", length(sites),
         " (got ", length(clip), "). For a partial specification, name the ",
         "sites: clip = c(delta = FALSE).", call. = FALSE)
  stats::setNames(clip, sites)
}

## Which anchors each kernel can impose.
##
## Every combination is available except penalty = "none" with anchor = "free".
## That one is deliberately absent rather than unimplemented: with no penalty
## and no anchor the likelihood is flat along the scale ray, so the model is not
## identified at all and only the vague normal prior stops the chain drifting.
## Leaving it out also keeps the unpenalised kernel a literal implementation of
## Benoit and Van den Poel (2012), which fixes sigma = 1.
.bbqr_allowed_anchors <- function(penalty) {
  switch(penalty,
         alasso = c("sigma1", "free", "beta1", "norm1", "normslopes"),
         lasso  = c("sigma1", "free", "beta1", "norm1", "normslopes"),
         none   = c("sigma1", "beta1", "norm1", "normslopes"))
}

## Map the user-facing 'anchor' onto the Fortran identification switches.
## Both norm anchors go through the same projection block and differ only in
## how the scaling constant is chosen: over the whole vector (norm1) or over
## the slopes alone (normslopes). In both cases beta0 and sigma are rescaled
## alongside, which is what makes the projection a genuine reparameterization
## of the likelihood rather than a change of model.
.bbqr_anchor_switches <- function(anchor, penalty) {
  allowed <- .bbqr_allowed_anchors(penalty)
  if (!anchor %in% allowed)
    stop(sprintf("anchor = \"%s\" is not available for penalty = \"%s\"; use %s.",
                 anchor, penalty,
                 paste(sprintf("\"%s\"", allowed), collapse = " or ")),
         call. = FALSE)
  list(fix_sigma           = anchor == "sigma1",
       constrain_beta_norm = anchor %in% c("norm1", "normslopes"),
       fix_beta1           = anchor == "beta1",
       norm_slopes_only    = anchor == "normslopes")
}

## Reshape the Fortran output, drop burn-in, undo standardization, and wrap the
## result in a "bbqr" object.
.bbqr_assemble <- function(fn, md, scaler, quantile, penalty, anchor,
                           out_rows, ndraw, burn, keep, pr, extra) {
  B  <- matrix(fn$betadraw, nrow = out_rows, ncol = md$nvar)
  b0 <- fn$beta0draw

  ## Burn-in is specified in iterations; the returned chain is already thinned.
  burn_rows <- floor(burn / keep)
  idx <- seq.int(burn_rows + 1L, out_rows)
  if (length(idx) < 1L)
    stop("'burn' leaves no draws after thinning; lower 'burn' or 'keep'.",
         call. = FALSE)

  B  <- B[idx, , drop = FALSE]
  b0 <- b0[idx]

  bt <- .bbqr_backtransform(B, b0, scaler)
  B  <- bt$beta
  b0 <- bt$beta0
  colnames(B) <- md$names

  trim <- function(v) if (is.null(v)) NULL else
    if (is.matrix(v)) v[idx, , drop = FALSE] else v[idx]

  out <- list(
    penalty      = penalty,
    anchor       = anchor,
    quantile     = quantile,
    names        = md$names,
    beta         = B,
    beta0        = b0,
    sigma        = trim(fn$sigmadraw),
    accept       = if (all(is.na(fn$accdraw))) NA_real_ else
                     utils::tail(fn$accdraw, 1L),
    ndraw        = as.integer(ndraw),
    burn         = as.integer(burn),
    keep         = as.integer(keep),
    ndraw_kept   = length(idx),
    n            = md$n,
    nvar         = md$nvar,
    ylevels      = md$ylevels,
    terms        = md$terms,
    standardized = !is.null(scaler),
    prior        = pr)

  for (nm in names(extra)) out[[nm]] <- trim(extra[[nm]])

  class(out) <- "bbqr"
  out
}

## Human-readable label for each penalty, used by print/summary.
.bbqr_penalty_label <- function(penalty) {
  switch(penalty,
         alasso = "adaptive lasso",
         lasso  = "lasso",
         none   = "no penalty")
}

## Human-readable label for each identification anchor.
.bbqr_anchor_label <- function(anchor) {
  switch(anchor,
         sigma1     = "sigma = 1",
         free       = "none (sigma sampled)",
         beta1      = "first slope = 1",
         norm1      = "||(beta0, beta)|| = 1",
         normslopes = "||beta|| = 1 (slopes only)")
}
