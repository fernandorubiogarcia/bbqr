## Continuous-response quantile regression: the three penalty layers.
##
## Implements derivations/Mathematical_Derivations_V5C.tex against the kernels
## QRc_AL_mcmc.f95, QRc_L_mcmc.f95 and QRc_BQR_mcmc.f95.
##
## This file is deliberately SELF-CONTAINED and touches nothing on the binary
## path. `bbqr()` and its fitters are unchanged; the continuous entry points are
## `cbqr()` and the three `cbqr_*()` wrappers. The two paths share only the
## generic helpers in internals.R that do not look at the response
## (.bbqr_check_mcmc, .bbqr_check_burn, .bbqr_init_beta, .bbqr_backtransform).
##
## Three things here have no binary counterpart and are not oversights:
##   * `q` -- the exponent on sigma in the local-scale prior. See Proposition 5C
##     and the per-layer defaults in .cbqr_default_q().
##   * `standardize` defaults to TRUE and scales y as well as X. Proposition 5C
##     makes response standardization part of the specification, not a
##     preprocessing convenience, because no proper hierarchy is exactly
##     response-scale equivariant.
##   * there is no `anchor` and no `sigma_fixed`. sigma is identified once y is
##     observed, so the anchors are misspecifications rather than conventions.

## ---------------------------------------------------------------------------
## The sixteen clamp sites in the continuous adaptive-lasso kernel.
##
## This is a SEPARATE contract from .BBQR_CLIP_SITES and must not be conflated
## with it: the continuous kernel has no y* draw, so `ystar` and `v4_sd` are
## gone and every later index shifts down by one or two. Passing the binary
## eighteen-name vector to the continuous kernel would misalign every site
## silently. It must match the C_* parameters at the top of QRc_AL_mcmc.f95
## name for name and position for position.
##
## Unlike the binary case, `sigma` and `v4_rate_sig` are LIVE here: the binary
## kernel reaches them only under a free sigma, which the package default never
## uses, whereas the continuous kernel samples sigma every sweep.
## ---------------------------------------------------------------------------
.CBQR_CLIP_SITES <- c("z", "s", "sigma", "lambdasq", "omega", "delta",
                      "mh_sd", "v4_z_psi", "v4_z_chi", "v4_zbound",
                      "v4_s_lam", "v4_s_psi", "v4_s_chi", "v4_sbound",
                      "v4_rate_sig", "v4_rate_om")

## Default exponent q, BY LAYER. These are not a style choice: each reproduces
## its own published specification, and the two published specifications differ.
##
##   alasso  q = 1  Alhamzawi, Yu & Benoit (2012) eq. (2.6): the penalty is
##                  (sigma^(1/2)/lambda_j)|beta_j|, a half power.
##   lasso   q = 0  Benoit, Al-Hamzawi & Yu (2013) eq. (13): the penalty is
##                  lambda|beta_j|, free of sigma.
##   none    q      irrelevant; there is no penalty layer to carry sigma.
##
## In the binary model sigma is anchored, so the two exponents give the same
## sampler and the difference is invisible. It is not invisible here. See
## Remark 5C. q = 2 is the response-scale equivariant choice for both layers.
.cbqr_default_q <- function(penalty)
  switch(penalty, alasso = 1, lasso = 0, none = 0)

## Model frame for a CONTINUOUS response. Same intercept-splitting contract as
## .bbqr_model_data -- the kernels never penalise the intercept, so it must not
## sit in X -- but y stays numeric and is checked for the opposite thing: a
## response that takes only two values is almost certainly a binary problem
## handed to the wrong entry point.
.cbqr_model_data <- function(formula, data) {
  if (missing(formula) || is.null(formula))
    stop("'formula' has to be specified.", call. = FALSE)
  mf    <- stats::model.frame(formula = formula, data = data)
  Xfull <- stats::model.matrix(attr(mf, "terms"), data = mf)
  y     <- stats::model.response(mf)

  int_idx <- which(colnames(Xfull) == "(Intercept)")
  X <- if (length(int_idx)) Xfull[, -int_idx, drop = FALSE] else Xfull

  if (is.factor(y) || is.logical(y))
    stop("The response must be numeric and continuous. For a binary response ",
         "use bbqr(), which fits the threshold model.", call. = FALSE)
  if (!is.numeric(y))
    stop("The response must be numeric.", call. = FALSE)
  if (ncol(X) < 1L)
    stop("At least one predictor is required.", call. = FALSE)
  if (anyNA(X) || anyNA(y))
    stop("Missing values are not supported; remove or impute them first.",
         call. = FALSE)
  if (!all(is.finite(y)))
    stop("The response contains non-finite values.", call. = FALSE)
  if (length(unique(y)) < 3L)
    stop("The response takes fewer than three distinct values. If it is ",
         "binary, use bbqr(); the continuous kernels assume an observed ",
         "response and will not identify a threshold model.", call. = FALSE)

  list(y = as.double(y), X = X, names = colnames(X), ylevels = NULL,
       n = length(y), nvar = ncol(X), terms = attr(mf, "terms"))
}

## Centre and scale y. Proposition 5C: the model cannot be made invariant to
## y -> c*y by any proper prior, so the scale is fixed by the specification
## instead. Returns the transformation so draws can be mapped back exactly.
.cbqr_standardize_y <- function(y, standardize) {
  if (!isTRUE(standardize)) return(list(y = y, y_mu = 0, y_sd = 1))
  y_mu <- mean(y)
  y_sd <- stats::sd(y)
  if (!is.finite(y_sd) || isTRUE(all.equal(y_sd, 0))) y_sd <- 1
  list(y = (y - y_mu) / y_sd, y_mu = y_mu, y_sd = y_sd)
}

## Map draws back to the original response units. On the standardized scale
##   (y - m)/s = b0* + sum_j B*_j x_j   =>   y = (m + s b0*) + sum_j (s B*_j) x_j,
## so slopes and intercept multiply by s and the intercept gains m. sigma is an
## INVERSE scale, so it divides by s. This is an exact identity.
.cbqr_backtransform_y <- function(B, b0, sigma, y_mu, y_sd) {
  list(beta  = B * y_sd,
       beta0 = y_mu + y_sd * b0,
       sigma = sigma / y_sd)
}

## Prior defaults for the continuous hierarchy.
##
## Differences from the binary defaults, each traceable to a proposition:
##   * b0_mean = 0 and b0_var = 100 rather than the tau-calibrated
##     (m00 = -theta, V00(tau)). That calibration is derived from
##     Pr(y = 1) = G_tau(beta0) and is meaningless for an intercept on the
##     scale of y (V5C section on intercept reference values).
##   * No a > 1 guard. Proposition 4C: the flat intercept needs only
##     a + n - 1 > 0, so the a = b = 2 reference and even a flat intercept are
##     admissible at every design and every n. The binary guard must not be
##     ported and is deliberately absent.
##   * alpha_omega = beta_omega = 0.1 for the adaptive lasso, against the
##     binary hierarchy's 2. Gamma(a, a) has mean 1 for every a, so this holds
##     the prior MEAN fixed and only widens it -- but the shape also decides
##     the density at the origin, and that is the operative difference:
##     shape > 1 sends the density to ZERO as omega -> 0, shape < 1 sends it to
##     infinity. omega is the adaptive scale, so small omega means strong
##     shrinkage, and Gamma(2,2) forbids exactly the region the data want.
##     Measured in Experiment C1's 6 x 2 factorial (3,200 paired scenarios per
##     cell, n = 1000): the posterior puts ~23% of omega's mass below 1e-3
##     under an improper reference, and Gamma(2,2) permits 0.0000 of it. The
##     prior is binding, in the direction that costs accuracy. Loosening to
##     Gamma(0.1, 0.1) moves RMSE 0.0895 -> 0.0843 (-5.8%), coverage
##     0.8008 -> 0.8234 and Winkler 0.7074 -> 0.6664, wins at every tau, and
##     never trades one quantile against another. Gamma(0.01, 0.01) is
##     marginally better still but sits at 22.7% mass below 1e-3, hard against
##     the degenerate boundary where the adaptive weights blow up; 0.1 takes
##     ~87% of the available gain with numerical margin left.
##     Propriety is unaffected: Proposition 1C requires only that the
##     hyperparameters be strictly positive -- there is no lower bound, and its
##     normalising bound does not involve alpha_omega or beta_omega at all.
##     Gamma(0.1, 0.1) is exactly as proper as Gamma(2, 2).
##     NOT ported to the lasso layer (alpha_tau/beta_tau below) or to the
##     binary hierarchy in prior(): the factorial covered the continuous
##     adaptive lasso only, and the binary sensitivity study reached a
##     different conclusion about what the binding costs.
.cbqr_prior <- function(penalty, user = NULL) {
  base <- list(
    a = 2, b = 2,
    b0_mean = 0, b0_var = 100,
    mh_sd0 = 0.5, target_acc = 0.3,
    sigma_init = 1,
    delta_init = 1, alpha_delta = 2, beta_delta = 2)
  pr <- switch(penalty,
    alasso = c(base, list(omega_init = 1, alpha_omega = 0.1, beta_omega = 0.1)),
    lasso  = c(base, list(tau_init = 1, alpha_tau = 2, beta_tau = 2)),
    none   = c(base, list(beta_var = 100)))
  if (!is.null(user)) {
    unknown <- setdiff(names(user), names(pr))
    if (length(unknown))
      stop("Unknown prior component(s) for penalty = ", sQuote(penalty), ": ",
           paste(sQuote(unknown), collapse = ", "), call. = FALSE)
    pr[names(user)] <- user
  }
  for (nm in c("a", "b", "b0_var")) {
    if (!is.numeric(pr[[nm]]) || length(pr[[nm]]) != 1L || pr[[nm]] <= 0)
      stop(sprintf("Prior component '%s' must be a single positive number.", nm),
           call. = FALSE)
  }
  pr
}

.cbqr_check_q <- function(q) {
  if (!is.numeric(q) || length(q) != 1L || !is.finite(q) || q < 0)
    stop("'q' must be a single finite number >= 0.", call. = FALSE)
  if (q > 0 && q < 1)
    warning("q in (0, 1) makes the sigma full conditional non-log-concave; ",
            "the slice sampler still runs but its reliability is not argued. ",
            "See the q /= 1 subsection of appendix V5C.", call. = FALSE)
  as.double(q)
}

## Shared assembly for the continuous fitters.
.cbqr_assemble <- function(fn, md, scaler, ystd, quantile, penalty, q,
                           out_rows, ndraw, burn, keep, pr, extra = list()) {
  burn_rows <- floor(burn / keep)
  idx <- seq.int(burn_rows + 1L, out_rows)
  if (length(idx) < 1L)
    stop("'burn' leaves no draws after thinning; lower 'burn' or 'keep'.",
         call. = FALSE)

  B     <- matrix(fn$betadraw, nrow = out_rows, ncol = md$nvar)[idx, , drop = FALSE]
  b0    <- fn$beta0draw[idx]
  sigma <- fn$sigmadraw[idx]

  ## Undo the X standardization first (it is a statement about the linear
  ## predictor), then the y standardization (a statement about the scale).
  bt <- .bbqr_backtransform(B, b0, scaler)
  yb <- .cbqr_backtransform_y(bt$beta, bt$beta0, sigma, ystd$y_mu, ystd$y_sd)
  B <- yb$beta; b0 <- yb$beta0; sigma <- yb$sigma
  colnames(B) <- md$names

  trim <- function(v) if (is.null(v)) NULL else
    if (is.matrix(v)) v[idx, , drop = FALSE] else v[idx]

  out <- c(list(
    family    = "continuous",
    penalty   = penalty,
    q         = q,
    quantile  = quantile,
    beta      = B,
    beta0     = b0,
    sigma     = sigma,
    ndraw     = ndraw,
    burn      = burn,
    keep      = keep,
    ndraw_kept = length(idx),
    ## MH acceptance rate of the delta step, taken at the last retained draw.
    ## NA for penalty = "none", which has no Metropolis step.
    accept    = if (is.null(extra$acc)) NA_real_ else utils::tail(extra$acc, 1L),
    n         = md$n,
    nvar      = md$nvar,
    names     = md$names,
    terms     = md$terms,
    prior     = pr,
    standardized = list(x = !is.null(scaler), y = ystd$y_sd != 1 || ystd$y_mu != 0,
                        y_mu = ystd$y_mu, y_sd = ystd$y_sd)),
    lapply(extra, trim))
  class(out) <- c("cbqr", "bbqr")
  out
}

## ---------------------------------------------------------------------------
## Adaptive lasso
## ---------------------------------------------------------------------------
#' Continuous quantile regression with the adaptive lasso
#'
#' Fits the continuous-response adaptive-lasso hierarchy for an observed response. The
#' response is observed, so there is no latent *response* draw and no
#' identification
#' anchor: sigma is identified by the data and is sampled every sweep.
#'
#' @param formula A two-sided formula. The intercept is split off and never
#'   penalised, so it must not appear among the penalised columns.
#' @param data Optional data frame.
#' @param quantile Quantile level(s) in (0, 1).
#' @param ndraw,burn,keep MCMC controls. `ndraw` must be divisible by `keep`.
#' @param prior Optional named list overriding the continuous prior
#'   defaults. Accepted names are `a`, `b`, `b0_mean`, `b0_var`, `mh_sd0`,
#'   `target_acc`, `sigma_init`, `delta_init`, `alpha_delta`, `beta_delta`,
#'   (note `mh_sd0` defaults to 0.5 here, not the binary 0.05: `delta` mixes
#'   over a wider range once `sigma` is live)
#'   plus `omega_init`, `alpha_omega` and `beta_omega`. This is **not** a [prior()] object:
#'   that constructor builds the binary hierarchy and is rejected here.
#' @param standardize Centre and scale **both** `X` and `y`. Defaults to `TRUE`,
#'   which is not the binary default: no proper hierarchy is
#'   exactly response-scale equivariant, so fixing the scale is part of the
#'   specification rather than a preprocessing convenience. Draws are mapped
#'   back to the original units exactly.
#' @param beta_init,beta0_init Starting values.
#' @param q Exponent on sigma in the local-scale prior rate
#'   `sigma^q / (2 lambda_j^2)`. Defaults to 1, reproducing Alhamzawi, Yu &
#'   Benoit (2012). `q = 2` is response-scale equivariant in the penalty block
#'   -- no proper hierarchy can be made exactly invariant to rescaling `y` --
#'   and draws sigma by slice sampling rather than conjugately. See
#'   `vignette("bbqr")`.
#' @param clip Optional clamp configuration: a single logical, a named logical
#'   vector over the sixteen continuous clamp sites, or a character vector of
#'   site names. **Not** the binary eighteen-site contract.
#' @param sampler `"gig_half"` (default) or `"v4_invgaus"`.
#' @return An object of class `cbqr`.
#' @seealso [cbqr()] for the common interface, [bbqr_alasso()] for the
#'   binary counterpart. [prior()] configures that one, not this one.
#' @references
#' Alhamzawi, R., Yu, K. and Benoit, D. F. (2012). Bayesian adaptive lasso
#' quantile regression. *Statistical Modelling*, 12(3), 279--297.
#' \doi{10.1177/1471082X1101200304}
#'
#' Benoit, D. F., Al-Hamzawi, R. and Yu, K. (2013). Bayesian lasso binary
#' quantile regression. *Computational Statistics*, 28(6), 2861--2873.
#' \doi{10.1007/s00180-013-0439-0}
#'
#' Kozumi, H. and Kobayashi, G. (2011). Gibbs sampling methods for Bayesian
#' quantile regression. *Journal of Statistical Computation and Simulation*,
#' 81(11), 1565--1578. \doi{10.1080/00949655.2010.496117}
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(120), x2 = rnorm(120))
#' d$y <- 1 + 2 * d$x1 - d$x2 + rnorm(120)
#' fit <- cbqr_alasso(y ~ x1 + x2, d, ndraw = 400, burn = 100)
#' round(coef(fit), 2)
#' @export
cbqr_alasso <- function(formula, data = NULL, quantile = 0.5,
                        ndraw = 10000, burn = NULL, keep = 1,
                        prior = NULL, q = NULL,
                        standardize = TRUE, clip = NULL, sampler = NULL,
                        beta_init = NULL, beta0_init = 0) {
  .bbqr_check_mcmc(quantile, ndraw, keep)
  burn <- .bbqr_check_burn(burn, ndraw)
  md   <- .cbqr_model_data(formula, data)
  pr   <- .cbqr_prior("alasso", prior)
  q    <- .cbqr_check_q(if (is.null(q)) .cbqr_default_q("alasso") else q)

  clipv <- if (is.null(clip)) stats::setNames(rep(FALSE, 16L), .CBQR_CLIP_SITES)
           else .cbqr_clip_vector(clip)
  zs <- if (is.null(sampler) || identical(sampler, "gig_half")) 0L else 1L

  std   <- .bbqr_maybe_standardize(md, standardize)
  ystd  <- .cbqr_standardize_y(md$y, standardize)
  binit <- .bbqr_init_beta(beta_init, md$names, md$nvar)
  out_rows <- ndraw / keep

  fn <- .Fortran("qrc_al_mcmc",
    as.integer(md$n), as.integer(md$nvar), as.integer(ndraw), as.integer(keep),
    as.double(ystd$y), as.double(quantile), as.double(std$X),
    as.double(pr$a), as.double(pr$b), as.double(q), as.double(pr$sigma_init),
    as.double(pr$delta_init), as.double(pr$omega_init),
    as.double(pr$mh_sd0), as.double(pr$target_acc),
    as.integer(burn), as.logical(clipv), as.integer(zs),
    as.double(pr$alpha_delta), as.double(pr$beta_delta),
    as.double(pr$alpha_omega), as.double(pr$beta_omega),
    as.double(pr$b0_mean), as.double(1 / pr$b0_var),
    betadraw     = double(out_rows * md$nvar),
    beta0draw    = double(out_rows),
    sigmadraw    = double(out_rows),
    lambdasqdraw = double(out_rows * md$nvar),
    omegadraw    = double(out_rows),
    deltadraw    = double(out_rows),
    accdraw      = double(out_rows),
    cliphits     = double(16L),
    beta_init    = binit,
    beta0_init   = as.double(beta0_init),
    PACKAGE = "bbqr")

  out <- .cbqr_assemble(fn, md, std$scaler, ystd, quantile, "alasso", q,
                        out_rows, ndraw, burn, keep, pr,
                        extra = list(
                          lambdasq = matrix(fn$lambdasqdraw, out_rows, md$nvar,
                                            dimnames = list(NULL, md$names)),
                          omega    = fn$omegadraw,
                          delta    = fn$deltadraw,
                          acc      = fn$accdraw))
  out$clip      <- clipv
  out$clip_hits <- stats::setNames(fn$cliphits, .CBQR_CLIP_SITES)
  out$sampler   <- if (zs == 0L) "gig_half" else "v4_invgaus"
  out$b0_prior  <- c(mean = pr$b0_mean, var = pr$b0_var, prec = 1 / pr$b0_var)
  out
}

## ---------------------------------------------------------------------------
## Lasso
## ---------------------------------------------------------------------------
#' Continuous quantile regression with the lasso
#'
#' Fits the continuous-response lasso hierarchy for an observed response: a single shared
#' `lambda^2` carried as a rate with a Gamma prior and a `tau_h` hyperlayer.
#' The prior rate here is `sigma^q * lambda^2 / 2`, with `lambda^2` a single
#' global scale multiplying; the adaptive-lasso rate is
#' `sigma^q / (2 * lambda_j^2)`, one scale per coefficient and dividing. The
#' two forms are different by construction, not a typo.
#'
#' This is Benoit, Al-Hamzawi & Yu (2013)'s hierarchy, **not** the adaptive one
#' with `k = 1`, and the two are not reachable from one another.
#'
#' @param formula A two-sided formula. The intercept is split off and never
#'   penalised, so it must not appear among the penalised columns.
#' @param data Optional data frame.
#' @param quantile Quantile level(s) in (0, 1).
#' @param ndraw,burn,keep MCMC controls. `ndraw` must be divisible by `keep`.
#' @param prior Optional named list overriding the continuous prior
#'   defaults. Accepted names are `a`, `b`, `b0_mean`, `b0_var`, `mh_sd0`,
#'   `target_acc`, `sigma_init`, `delta_init`, `alpha_delta`, `beta_delta`,
#'   (note `mh_sd0` defaults to 0.5 here, not the binary 0.05: `delta` mixes
#'   over a wider range once `sigma` is live)
#'   plus `tau_init`, `alpha_tau` and `beta_tau`. This is **not** a [prior()] object:
#'   that constructor builds the binary hierarchy and is rejected here.
#' @param standardize Centre and scale **both** `X` and `y`. Defaults to `TRUE`,
#'   which is not the binary default: no proper hierarchy is
#'   exactly response-scale equivariant, so fixing the scale is part of the
#'   specification rather than a preprocessing convenience. Draws are mapped
#'   back to the original units exactly.
#' @param beta_init,beta0_init Starting values.
#' @param q Exponent on sigma in the local-scale prior rate
#'   `sigma^q lambda^2 / 2`. Defaults to **0**, reproducing the published
#'   sigma-free penalty `lambda|beta_j|`. This differs from the adaptive
#'   layer's default of 1 because the two published specifications differ.
#' @param use_boundaries Apply the kernel's numerical clamps.
#' @return An object of class `cbqr`.
#' @seealso [cbqr()] for the common interface, [bbqr_lasso()] for the
#'   binary counterpart.
#' @references
#' Alhamzawi, R., Yu, K. and Benoit, D. F. (2012). Bayesian adaptive lasso
#' quantile regression. *Statistical Modelling*, 12(3), 279--297.
#' \doi{10.1177/1471082X1101200304}
#'
#' Benoit, D. F., Al-Hamzawi, R. and Yu, K. (2013). Bayesian lasso binary
#' quantile regression. *Computational Statistics*, 28(6), 2861--2873.
#' \doi{10.1007/s00180-013-0439-0}
#'
#' Kozumi, H. and Kobayashi, G. (2011). Gibbs sampling methods for Bayesian
#' quantile regression. *Journal of Statistical Computation and Simulation*,
#' 81(11), 1565--1578. \doi{10.1080/00949655.2010.496117}
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(120), x2 = rnorm(120))
#' d$y <- 1 + 2 * d$x1 - d$x2 + rnorm(120)
#' fit <- cbqr_lasso(y ~ x1 + x2, d, ndraw = 400, burn = 100)
#' round(coef(fit), 2)
#' @export
cbqr_lasso <- function(formula, data = NULL, quantile = 0.5,
                       ndraw = 10000, burn = NULL, keep = 1,
                       prior = NULL, q = NULL,
                       standardize = TRUE, use_boundaries = FALSE,
                       beta_init = NULL, beta0_init = 0) {
  .bbqr_check_mcmc(quantile, ndraw, keep)
  burn <- .bbqr_check_burn(burn, ndraw)
  md   <- .cbqr_model_data(formula, data)
  pr   <- .cbqr_prior("lasso", prior)
  q    <- .cbqr_check_q(if (is.null(q)) .cbqr_default_q("lasso") else q)

  std   <- .bbqr_maybe_standardize(md, standardize)
  ystd  <- .cbqr_standardize_y(md$y, standardize)
  binit <- .bbqr_init_beta(beta_init, md$names, md$nvar)
  out_rows <- ndraw / keep

  fn <- .Fortran("qrc_l_mcmc",
    as.integer(md$n), as.integer(md$nvar), as.integer(ndraw), as.integer(keep),
    as.double(ystd$y), as.double(quantile), as.double(std$X),
    as.double(pr$a), as.double(pr$b), as.double(q), as.double(pr$sigma_init),
    as.double(pr$delta_init), as.double(pr$tau_init),
    as.double(pr$mh_sd0), as.double(pr$target_acc),
    as.integer(burn), as.logical(use_boundaries),
    as.double(pr$alpha_delta), as.double(pr$beta_delta),
    as.double(pr$alpha_tau), as.double(pr$beta_tau),
    as.double(pr$b0_mean), as.double(1 / pr$b0_var),
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

  out <- .cbqr_assemble(fn, md, std$scaler, ystd, quantile, "lasso", q,
                        out_rows, ndraw, burn, keep, pr,
                        extra = list(
                          lambdasq = fn$lambda2draw,
                          tau_hyper = fn$tauhyperdraw,
                          delta     = fn$deltadraw,
                          acc       = fn$accdraw))
  out$b0_prior <- c(mean = pr$b0_mean, var = pr$b0_var, prec = 1 / pr$b0_var)
  out
}

## ---------------------------------------------------------------------------
## Unpenalised
## ---------------------------------------------------------------------------
#' Continuous quantile regression, unpenalised
#'
#' Fits the unpenalised hierarchy for an observed response. Unlike its binary
#' counterpart, which holds sigma at the identification anchor and never samples
#' it, this kernel draws sigma every sweep from
#' its Gamma full conditional given the latents, with shape `a + 3n/2`.
#'
#' @param formula A two-sided formula. The intercept is split off and never
#'   penalised, so it must not appear among the penalised columns.
#' @param data Optional data frame.
#' @param quantile Quantile level(s) in (0, 1).
#' @param ndraw,burn,keep MCMC controls. `ndraw` must be divisible by `keep`.
#' @param prior Optional named list overriding the continuous prior
#'   defaults. Accepted names are `a`, `b`, `b0_mean`, `b0_var`, `mh_sd0`,
#'   `target_acc`, `sigma_init`, `delta_init`, `alpha_delta`, `beta_delta`,
#'   (note `mh_sd0` defaults to 0.5 here, not the binary 0.05: `delta` mixes
#'   over a wider range once `sigma` is live)
#'   plus `beta_var`. This is **not** a [prior()] object:
#'   that constructor builds the binary hierarchy and is rejected here.
#' @param standardize Centre and scale **both** `X` and `y`. Defaults to `TRUE`,
#'   which is not the binary default: no proper hierarchy is
#'   exactly response-scale equivariant, so fixing the scale is part of the
#'   specification rather than a preprocessing convenience. Draws are mapped
#'   back to the original units exactly.
#' @param beta_init,beta0_init Starting values.
#' @param use_boundaries Apply the kernel's numerical clamps.
#' @return An object of class `cbqr`.
#' @seealso [cbqr()] for the common interface, [bbqr_none()] for the
#'   binary counterpart.
#' @references
#' Alhamzawi, R., Yu, K. and Benoit, D. F. (2012). Bayesian adaptive lasso
#' quantile regression. *Statistical Modelling*, 12(3), 279--297.
#' \doi{10.1177/1471082X1101200304}
#'
#' Benoit, D. F., Al-Hamzawi, R. and Yu, K. (2013). Bayesian lasso binary
#' quantile regression. *Computational Statistics*, 28(6), 2861--2873.
#' \doi{10.1007/s00180-013-0439-0}
#'
#' Kozumi, H. and Kobayashi, G. (2011). Gibbs sampling methods for Bayesian
#' quantile regression. *Journal of Statistical Computation and Simulation*,
#' 81(11), 1565--1578. \doi{10.1080/00949655.2010.496117}
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(120), x2 = rnorm(120))
#' d$y <- 1 + 2 * d$x1 - d$x2 + rnorm(120)
#' fit <- cbqr_none(y ~ x1 + x2, d, ndraw = 400, burn = 100)
#' round(coef(fit), 2)
#' @export
cbqr_none <- function(formula, data = NULL, quantile = 0.5,
                      ndraw = 10000, burn = NULL, keep = 1,
                      prior = NULL, standardize = TRUE, use_boundaries = FALSE,
                      beta_init = NULL, beta0_init = 0) {
  .bbqr_check_mcmc(quantile, ndraw, keep)
  burn <- .bbqr_check_burn(burn, ndraw)
  md   <- .cbqr_model_data(formula, data)
  pr   <- .cbqr_prior("none", prior)

  std   <- .bbqr_maybe_standardize(md, standardize)
  ystd  <- .cbqr_standardize_y(md$y, standardize)
  binit <- .bbqr_init_beta(beta_init, md$names, md$nvar)
  out_rows <- ndraw / keep

  fn <- .Fortran("qrc_bqr_mcmc",
    as.integer(md$n), as.integer(md$nvar), as.integer(ndraw), as.integer(keep),
    as.double(ystd$y), as.double(quantile), as.double(std$X),
    as.double(pr$beta_var), as.double(pr$a), as.double(pr$b),
    as.double(0), as.double(pr$sigma_init),
    as.logical(use_boundaries),
    as.double(pr$b0_mean), as.double(1 / pr$b0_var),
    betadraw   = double(out_rows * md$nvar),
    beta0draw  = double(out_rows),
    sigmadraw  = double(out_rows),
    beta_init  = binit,
    beta0_init = as.double(beta0_init),
    PACKAGE = "bbqr")

  out <- .cbqr_assemble(fn, md, std$scaler, ystd, quantile, "none", NA_real_,
                        out_rows, ndraw, burn, keep, pr)
  out$b0_prior <- c(mean = pr$b0_mean, var = pr$b0_var, prec = 1 / pr$b0_var)
  out
}

## Resolve a user `clip` argument to the sixteen-vector the continuous kernel
## expects. Accepts a single logical, or a named subset of .CBQR_CLIP_SITES.
.cbqr_clip_vector <- function(clip) {
  full <- stats::setNames(rep(FALSE, 16L), .CBQR_CLIP_SITES)
  if (is.logical(clip) && length(clip) == 1L) {
    full[] <- clip
    return(full)
  }
  if (is.logical(clip) && length(clip) == 16L) {
    if (!is.null(names(clip))) {
      if (!setequal(names(clip), .CBQR_CLIP_SITES))
        stop("'clip' names must be exactly the sixteen continuous clamp sites; ",
             "the binary eighteen-name vector is a DIFFERENT contract and ",
             "would misalign every site.", call. = FALSE)
      return(clip[.CBQR_CLIP_SITES])
    }
    return(stats::setNames(clip, .CBQR_CLIP_SITES))
  }
  if (is.character(clip)) {
    unknown <- setdiff(clip, .CBQR_CLIP_SITES)
    if (length(unknown))
      stop("Unknown clamp site(s): ", paste(sQuote(unknown), collapse = ", "),
           call. = FALSE)
    full[clip] <- TRUE
    return(full)
  }
  stop("'clip' must be a single logical, a logical vector of length 16, or a ",
       "character vector of site names.", call. = FALSE)
}

## ---------------------------------------------------------------------------
## User-facing entry point
## ---------------------------------------------------------------------------
#' Bayesian quantile regression for a continuous response
#'
#' User-facing entry point for the three continuous penalty layers. The binary
#' threshold model is fitted by [bbqr()] instead; the two share no arguments
#' that refer to identification, because a continuous response identifies the
#' scale and a binary one does not.
#'
#' @param formula A two-sided formula. The intercept is split off and never
#'   penalised, so it must not appear among the penalised columns.
#' @param data Optional data frame.
#' @param quantile Quantile level(s) in (0, 1).
#' @param ndraw,burn,keep MCMC controls. `ndraw` must be divisible by `keep`.
#' @param prior Optional named list overriding the continuous prior
#'   defaults. Accepted names are `a`, `b`, `b0_mean`, `b0_var`, `mh_sd0`,
#'   `target_acc`, `sigma_init`, `delta_init`, `alpha_delta`, `beta_delta`,
#'   (note `mh_sd0` defaults to 0.5 here, not the binary 0.05: `delta` mixes
#'   over a wider range once `sigma` is live)
#'   and, by penalty, `omega_init`/`alpha_omega`/`beta_omega`
#'   (`"alasso"`), `tau_init`/`alpha_tau`/`beta_tau` (`"lasso"`) or
#'   `beta_var` (`"none"`). This is **not** a [prior()] object:
#'   that constructor builds the binary hierarchy and is rejected here.
#' @param standardize Centre and scale **both** `X` and `y`. Defaults to `TRUE`,
#'   which is not the binary default: no proper hierarchy is
#'   exactly response-scale equivariant, so fixing the scale is part of the
#'   specification rather than a preprocessing convenience. Draws are mapped
#'   back to the original units exactly.
#' @param beta_init,beta0_init Starting values.
#' @param penalty One of `"alasso"`, `"lasso"`, `"none"`.
#' @param q Exponent on sigma in the local-scale prior. Layer-specific default:
#'   1 for `"alasso"`, 0 for `"lasso"`. Supplying it with `penalty = "none"`
#'   is an error: there is no penalty layer for the scale to enter.
#' @param use_boundaries Apply the kernel's numerical clamps. `"lasso"` and
#'   `"none"` only; the adaptive-lasso kernel takes its clamps through the
#'   per-site `clip` argument instead and ignores this one.
#' @param clip Per-site clamp configuration; `"alasso"` only.
#' @param sampler Latent-draw routine; `"alasso"` only.
#' @return A `cbqr` object, or a `cbqr.list` when `quantile` has length > 1.
#'   The components are:
#'   \describe{
#'     \item{`family`, `penalty`, `q`, `quantile`}{what was fitted.}
#'     \item{`beta`, `beta0`, `sigma`}{the kept draws: a matrix of slopes, and
#'       vectors for the intercept and the ALD inverse-scale. Unlike a binary
#'       fit, `sigma` is always a live parameter here.}
#'     \item{`lambdasq`, `omega`, `delta`}{the shrinkage layer, `"alasso"`
#'       only: the per-coefficient scales, the global scale, and its shape.}
#'     \item{`lambdasq`, `tau_hyper`, `delta`}{the shrinkage layer, `"lasso"`
#'       only: one global scale, its hyperparameter, and the shape. Absent for
#'       `"none"`, which has no shrinkage layer.}
#'     \item{`accept`, `acc`}{Metropolis acceptance rate for `delta`, and its
#'       running trace. `NA` when there is no Metropolis step.}
#'     \item{`ndraw`, `burn`, `keep`, `ndraw_kept`}{the MCMC controls as used.}
#'     \item{`n`, `nvar`, `names`, `terms`, `call`}{the data and the model
#'       frame the fit came from.}
#'     \item{`prior`}{the resolved prior list, including any defaults filled
#'       in.}
#'     \item{`standardized`}{the centring and scaling applied to `X` and `y`,
#'       or `NULL`. Draws are already mapped back, so this is a record rather
#'       than something to undo.}
#'     \item{`b0_prior`}{the intercept prior actually used.}
#'     \item{`clip`, `clip_hits`, `sampler`}{the clamp configuration, how often
#'       each site fired, and the latent-draw routine. `"alasso"` only.}
#'   }
#' @seealso [bbqr()] for the binary threshold model; [cbqr_alasso()],
#'   [cbqr_lasso()] and [cbqr_none()] for the direct forms; [predict.cbqr()],
#'   [summary.cbqr()] and [plot.cbqr()] for the fitted object. [prior()] does
#'   **not** configure these fitters.
#' @references
#' Alhamzawi, R., Yu, K. and Benoit, D. F. (2012). Bayesian adaptive lasso
#' quantile regression. *Statistical Modelling*, 12(3), 279--297.
#' \doi{10.1177/1471082X1101200304}
#'
#' Benoit, D. F., Al-Hamzawi, R. and Yu, K. (2013). Bayesian lasso binary
#' quantile regression. *Computational Statistics*, 28(6), 2861--2873.
#' \doi{10.1007/s00180-013-0439-0}
#'
#' Kozumi, H. and Kobayashi, G. (2011). Gibbs sampling methods for Bayesian
#' quantile regression. *Journal of Statistical Computation and Simulation*,
#' 81(11), 1565--1578. \doi{10.1080/00949655.2010.496117}
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(120), x2 = rnorm(120))
#' d$y <- 1 + 2 * d$x1 - d$x2 + rnorm(120)
#' fit <- cbqr(y ~ x1 + x2, d, ndraw = 400, burn = 100)
#' round(coef(fit), 2)
#' @export
cbqr <- function(formula, data = NULL, quantile = 0.5,
                 penalty = c("alasso", "lasso", "none"),
                 ndraw = 10000, burn = NULL, keep = 1,
                 prior = NULL, q = NULL,
                 standardize = TRUE, use_boundaries = FALSE, clip = NULL,
                 sampler = NULL, beta_init = NULL, beta0_init = 0) {
  penalty <- match.arg(penalty)
  cl <- match.call()

  if (!is.numeric(quantile) || length(quantile) < 1L)
    stop("'quantile' must be a numeric vector of length at least one.",
         call. = FALSE)
  if (!is.null(clip) && penalty != "alasso")
    stop("'clip' is only available for penalty = \"alasso\"; ",
         sQuote(penalty), " takes the all-or-nothing 'use_boundaries'.",
         call. = FALSE)
  if (!is.null(sampler) && penalty != "alasso")
    stop("'sampler' is only available for penalty = \"alasso\".", call. = FALSE)
  if (!is.null(q) && penalty == "none")
    stop("'q' has no effect for penalty = \"none\": there is no penalty layer ",
         "for the scale to enter. Omit it.", call. = FALSE)

  fit_one <- function(qt) {
    out <- switch(penalty,
      alasso = cbqr_alasso(formula, data, quantile = qt, ndraw = ndraw,
                           burn = burn, keep = keep, prior = prior, q = q,
                           standardize = standardize, clip = clip,
                           sampler = sampler,
                           beta_init = beta_init, beta0_init = beta0_init),
      lasso  = cbqr_lasso(formula, data, quantile = qt, ndraw = ndraw,
                          burn = burn, keep = keep, prior = prior, q = q,
                          standardize = standardize,
                          use_boundaries = use_boundaries,
                          beta_init = beta_init, beta0_init = beta0_init),
      none   = cbqr_none(formula, data, quantile = qt, ndraw = ndraw,
                         burn = burn, keep = keep, prior = prior,
                         standardize = standardize,
                         use_boundaries = use_boundaries,
                         beta_init = beta_init, beta0_init = beta0_init))
    out$call <- cl
    out
  }

  if (length(quantile) == 1L) return(fit_one(quantile))
  res <- lapply(quantile, fit_one)
  names(res) <- paste0("tau=", format(quantile))
  class(res) <- c("cbqr.list", "bbqr.list")
  res
}
