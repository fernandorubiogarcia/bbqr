#' Bayesian binary quantile regression
#'
#' Fits a quantile regression model to a binary response by Markov chain Monte
#' Carlo, using the asymmetric Laplace distribution (ALD) and the location-scale
#' mixture representation of Kozumi and Kobayashi (2011). Three penalty layers
#' are available through the `penalty` argument, the identification anchor is
#' exposed through `anchor`, and the prior hierarchy is chosen by `model` on
#' the [prior()] object. The defaults -- `model = "v5"` at `anchor = "sigma1"`
#' -- sample a posterior that is proper unconditionally and identified by
#' construction.
#'
#' @section Derivation version:
#' Each penalty has a published specification and a corrected one, selected by
#' `prior(penalty, model = )`. `"v3"` is the published hierarchy: for the
#' adaptive lasso and the lasso its joint posterior does not exist, because the
#' log-uniform hyperprior on the global shrinkage scale (\eqn{\omega} or
#' \eqn{\tau_h}) is not integrable along the ray where the scale, the
#' \eqn{\lambda_j^2} and the latent \eqn{s_j} vanish together while the
#' likelihood stays bounded below. `"v5"`, the default, replaces every improper
#' component with a proper one -- \eqn{\mathrm{Gamma}(2, 2)} on the global
#' scale, and a \eqn{\tau}-calibrated normal prior on the intercept -- so that
#' the posterior is proper with no condition on the design, on the number of
#' covariates, or on the anchor. `"v4"` is the intermediate step for the
#' adaptive lasso. The fitted object records which hierarchy was sampled in
#' `model`, and whether the chain targets a derived posterior in `derived`.
#' Under `"v4"` and `"v5"` the numerical clamps, argument floors and the
#' historical sampler routine are refused, because the derivations say none of
#' them is part of the target; `"v3"` keeps them reachable so that earlier
#' results can be reproduced. See [prior()] for the hierarchies in full.
#'
#' @section Choosing a penalty:
#' The three penalties share the same latent-variable augmentation and differ
#' only in what sits on top of the slopes:
#' \describe{
#'   \item{`"none"`}{A vague \eqn{N(0, \sigma^2_\beta)} prior on each slope and
#'     no shrinkage, following Benoit and Van den Poel (2012).}
#'   \item{`"lasso"`}{One shrinkage parameter shared by all slopes, following
#'     Benoit, Al-Hamzawi and Yu (2013).}
#'   \item{`"alasso"`}{A separate shrinkage parameter per slope, so that
#'     coefficients are penalised individually rather than uniformly,
#'     following Rubio Garcia (2023).}
#' }
#'
#' @section Identification:
#' A binary threshold model identifies the coefficient vector only up to a
#' positive scale: \eqn{\mathrm{sign}(x'\beta + \varepsilon)} is unchanged if
#' \eqn{\beta} and the error scale are multiplied by the same positive constant.
#' Some anchor must therefore be imposed, and the choice is not innocuous. The
#' available anchors are:
#' \describe{
#'   \item{`"sigma1"`}{Fix the ALD inverse-scale at `sigma_fixed` (default 1)
#'     and never sample it. The default for every penalty. It closes the scale
#'     ray outright, so the posterior summaries of \eqn{\beta} are on a stated
#'     scale, and it is the convention under which the \eqn{\tau}-calibrated
#'     intercept prior of `model = "v5"` is exact. In the simulation study
#'     behind this package its point estimates matched `"free"` while its
#'     credible intervals were markedly narrower, because a sampled
#'     \eqn{\sigma} lets the unidentified scale inflate them.}
#'   \item{`"free"`}{Impose nothing beyond the prior: sample \eqn{\sigma} from
#'     its full conditional and let the prior alone settle the scale. This is
#'     the sampled-\eqn{\sigma} convention of the adaptive-lasso and lasso
#'     derivations, and it is a derived posterior under every `model`; it is
#'     simply not identified by the likelihood, so `||beta||` is a prior
#'     artefact. Not available for `penalty = "none"`: with neither a penalty
#'     nor an anchor the likelihood is flat along the scale ray, so nothing
#'     identifies the model and only the vague prior stops the chain
#'     drifting.}
#'   \item{`"beta1"`}{Hold the first slope at 1 and never sample it. The
#'     first column of the design matrix must be a covariate with a genuinely
#'     non-zero coefficient for this to be sensible.}
#'   \item{`"norm1"`}{Rescale so that \eqn{\|(\beta_0, \beta)\| = 1} after every
#'     sweep.}
#'   \item{`"normslopes"`}{Rescale so that \eqn{\|\beta\| = 1} over the slopes
#'     alone, the convention used in the maximum-score literature.}
#' }
#' The two norm anchors share one projection step and differ only in how the
#' scaling constant is chosen. In both cases \eqn{\beta_0} and the ALD
#' inverse-scale \eqn{\sigma} are rescaled alongside the slopes. At the
#' observed-data level the model is invariant under
#' \eqn{(\beta_0, \beta, \sigma) \mapsto (c\beta_0, c\beta, \sigma/c)}.
#' These projections are identification experiments and are not steps in any
#' of the derivations.
#'
#' Note, however, that this invariance does not extend to the *penalised
#' posterior*: the shrinkage layer is stated on a fixed scale, so repeatedly
#' renormalising the coefficients moves them relative to the penalty. Anchors
#' that constrain the norm of the coefficient vector should therefore be
#' expected to behave differently from `"sigma1"`, which leaves the
#' coefficients where the sampler put them. Under `model = "v4"` or `"v5"`
#' the three projection anchors (`"beta1"`, `"norm1"`, `"normslopes"`) warn
#' and set `derived = FALSE` on the fit: the projection is an exact
#' reparameterisation of the likelihood but not of the prior, so the chain
#' targets no written-down posterior. They remain available because the
#' comparison is the point of exposing the anchor at all.
#'
#' @section Scaling:
#' The prior variances and the penalty are stated on the scale of the design
#' matrix as supplied. If one covariate is measured in units a hundred times
#' larger than another, a shared penalty shrinks them by very different
#' relative amounts, and the recovered direction can be badly wrong. Setting
#' `standardize = TRUE` centres and scales the covariates before sampling and
#' maps the draws back afterwards; the back-transform preserves the linear
#' predictor exactly, but the posterior itself differs because the penalty now
#' acts on comparable coefficients. Unless the covariates are already on a
#' common scale, prefer `standardize = TRUE`.
#'
#' @param formula A model formula. The response must be binary: 0/1, logical,
#'   or a two-level factor (the second level is treated as the success).
#' @param data A data frame in which to interpret `formula`.
#' @param quantile The quantile level(s) to fit, each strictly between 0 and 1.
#'   If more than one is given, the samplers are run once per level and an
#'   object of class `"bbqr.list"` is returned.
#' @param penalty The shrinkage layer: `"alasso"` (adaptive lasso, the
#'   default), `"lasso"`, or `"none"`.
#' @param ndraw Total number of MCMC iterations.
#' @param burn Number of initial iterations to discard. Also the window over
#'   which the Metropolis proposal for the penalty hyperparameter adapts.
#'   Defaults to `ndraw / 10`.
#' @param keep Thinning interval; every `keep`-th draw is retained. `ndraw`
#'   must be divisible by `keep`.
#' @param prior A prior specification from [prior()]. Defaults to
#'   `prior()` for this penalty, the fully proper `model = "v5"`
#'   hierarchy. Pass
#'   `prior(penalty, model = "v3")` to reproduce the published specification.
#' @param anchor The identification anchor; see the Identification section
#'   of [bbqr()]. Defaults to
#'   `"sigma1"` for every penalty. Every anchor is available for every
#'   penalty, with one exception: `penalty = "none"` with `anchor = "free"` is
#'   not identified and is rejected.
#' @param sigma_fixed The value the ALD inverse-scale is held at when
#'   `anchor = "sigma1"`.
#' @param standardize If `TRUE`, centre and scale the covariates before
#'   sampling and map the draws back to the original units afterwards.
#'   **Strongly recommended whenever the covariates are on different scales.**
#'   Note that this is not a mere reparameterization: the prior and the penalty
#'   apply to the standardized coefficients, so the posterior genuinely
#'   changes. That is the point — a single shared penalty is not comparable
#'   across covariates measured in different units, and leaving them unscaled
#'   can pull the estimated direction badly off. See the Scaling section of
#'   [bbqr()].
#' @param use_boundaries If `TRUE`, clamp the latent quantities to a wide
#'   numerical range. Off by default; turn it on only if a chain is producing
#'   non-finite values.
#' @param clip Which numerical clamps to apply, available for
#'   `penalty = "alasso"` only; see [bbqr_alasso()] for the eighteen sites and
#'   the accepted forms, including the `"none"`, `"v7clip"` and `"v4"`
#'   presets. Refused under `model = "v4"` and `"v5"` unless every site is
#'   off.
#' @param sampler Which routine draws the latents `z` and `s`, available for
#'   `penalty = "alasso"` only; `"gig_half"` (default) or `"v4_invgaus"`. See
#'   [bbqr_alasso()]. Only `"gig_half"` is accepted under `model = "v4"` and
#'   `"v5"`.
#'   Supplying it for the other penalties is an error rather than a silent
#'   no-op, because their kernels draw with one fixed routine.
#' @param beta_init,beta0_init Starting values for the slopes and the
#'   intercept. `beta_init` may be named, in which case it is matched to the
#'   columns of the design matrix.
#'
#' @return For a single `quantile`, an object of class `"bbqr"`: a list with
#' components
#' \describe{
#'   \item{`beta`}{Matrix of retained slope draws, one column per covariate.}
#'   \item{`beta0`}{Vector of retained intercept draws.}
#'   \item{`sigma`}{Vector of retained ALD inverse-scale draws. Constant at
#'     `sigma_fixed` under `anchor = "sigma1"`; under the norm anchors it
#'     varies, because the rescaling step absorbs the scale into it.}
#'   \item{`accept`}{Final Metropolis acceptance rate for the penalty
#'     hyperparameter, or `NA` when the sampler has no Metropolis step.}
#'   \item{`penalty`, `anchor`, `quantile`}{The settings used.}
#'   \item{`prior`}{The prior actually applied.}
#'   \item{`model`}{Which hierarchy was sampled: `"v3"`, `"v4"` or `"v5"`.}
#'   \item{`derived`}{`TRUE` when the chain targets a posterior one of the
#'     derivations writes down; `FALSE` when a projection anchor was combined
#'     with `"v4"` or `"v5"`.}
#'   \item{`b0_prior`}{The intercept prior actually used, as `mean`, `var` and
#'     `prec`. Recorded because the `"v5"` reference depends on `quantile`,
#'     so it cannot be reconstructed from the prior object alone.}
#' }
#' Penalty-specific draws are also returned: `lambdasq`, `omega` and `delta`
#' for the adaptive lasso; `lambdasq`, `tauhyper` and `delta` for the lasso.
#' For several quantiles, an object of class `"bbqr.list"` holding one such
#' object per level.
#'
#' @references
#' Benoit, D. F. and Van den Poel, D. (2012). Binary quantile regression: a
#' Bayesian approach based on the asymmetric Laplace distribution.
#' *Journal of Applied Econometrics*, 27(7), 1174--1188.
#' \doi{10.1002/jae.1216}
#'
#' Benoit, D. F., Al-Hamzawi, R. and Yu, K. (2013). Bayesian lasso binary
#' quantile regression. *Computational Statistics*, 28(6), 2861--2873.
#' \doi{10.1007/s00180-013-0439-0}
#'
#' Kozumi, H. and Kobayashi, G. (2011). Gibbs sampling methods for Bayesian
#' quantile regression. *Journal of Statistical Computation and Simulation*,
#' 81(11), 1565--1578. \doi{10.1080/00949655.2010.496117}
#'
#' Rubio Garcia, F. (2023). Bayesian adaptive lasso binary quantile
#' regression with hybrid resampling for classification of imbalanced data.
#' M.S. thesis, Wichita State University, Dept. of Mathematics, Statistics,
#' and Physics.
#' \url{https://soar.wichita.edu/entities/publication/a2f86232-4704-4ec2-b685-751e7b04ec42}
#'
#' @seealso [prior()] for hyperparameters, [summary.bbqr()] for posterior
#'   summaries, [plot.bbqr()] for trace and density plots, and
#'   [predict.bbqr()] for fitted probabilities. [cbqr()] fits the same three
#'   penalties when the response is observed rather than thresholded.
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' X <- matrix(rnorm(n * 4), n, 4, dimnames = list(NULL, paste0("x", 1:4)))
#' y <- as.numeric(X %*% c(1.5, -1, 0, 0) + rnorm(n) > 0)
#' dat <- data.frame(y = y, X)
#'
#' fit <- bbqr(y ~ x1 + x2 + x3 + x4, data = dat, quantile = 0.5,
#'             ndraw = 1000, burn = 200)
#' fit
#' summary(fit)
#'
#' # Compare the three penalties at the same quantile
#' fits <- lapply(c("none", "lasso", "alasso"), function(p)
#'   bbqr(y ~ x1 + x2 + x3 + x4, data = dat, penalty = p,
#'        ndraw = 1000, burn = 200))
#' sapply(fits, coef)
#'
#' # The published adaptive-lasso hierarchy, whose posterior does not exist,
#' # is still available for reproducing earlier results
#' old <- bbqr(y ~ x1 + x2 + x3 + x4, data = dat, ndraw = 1000, burn = 200,
#'             prior = prior("alasso", model = "v3"))
#' c(fit$model, old$model)
#'
#' @export
bbqr <- function(formula, data = NULL, quantile = 0.5,
                 penalty = c("alasso", "lasso", "none"),
                 ndraw = 10000, burn = NULL, keep = 1,
                 prior = NULL,
                 anchor = NULL, sigma_fixed = 1,
                 standardize = FALSE, use_boundaries = FALSE, clip = NULL,
                 sampler = NULL,
                 beta_init = NULL, beta0_init = 0) {

  penalty <- match.arg(penalty)
  cl <- match.call()

  if (!is.numeric(quantile) || length(quantile) < 1L)
    stop("'quantile' must be a numeric vector of length at least one.",
         call. = FALSE)

  ## sigma = 1 for every penalty. The adaptive-lasso and lasso derivations
  ## cover a sampled sigma too ("free"), but the likelihood never identifies
  ## it: with the ray open the norm of beta is a prior artefact and the
  ## intervals widen accordingly, while the tau-calibrated intercept prior of
  ## "v5" is exact only at sigma = 1. "free" stays one argument away.
  if (is.null(anchor)) anchor <- "sigma1"

  ## Only the adaptive-lasso kernel takes per-site switches. Accepting `clip`
  ## and dropping it for the others would silently mislabel an ablation arm.
  if (!is.null(clip) && penalty != "alasso")
    stop("'clip' is only available for penalty = \"alasso\"; ",
         sQuote(penalty), " takes the all-or-nothing 'use_boundaries'.",
         call. = FALSE)

  ## Same reasoning as `clip`: only the adaptive-lasso kernel takes the switch,
  ## and accepting it for the others would mislabel the result.
  if (!is.null(sampler) && penalty != "alasso")
    stop("'sampler' is only available for penalty = \"alasso\"; ",
         sQuote(penalty), " draws its latents with one fixed routine.",
         call. = FALSE)
  if (is.null(sampler)) sampler <- "gig_half"

  fit_one <- function(q) {
    out <- switch(penalty,
      alasso = bbqr_alasso(formula, data, quantile = q, ndraw = ndraw,
                           burn = burn, keep = keep, prior = prior,
                           anchor = anchor, sigma_fixed = sigma_fixed,
                           standardize = standardize,
                           use_boundaries = use_boundaries, clip = clip,
                           sampler = sampler,
                           beta_init = beta_init, beta0_init = beta0_init),
      lasso  = bbqr_lasso(formula, data, quantile = q, ndraw = ndraw,
                          burn = burn, keep = keep, prior = prior,
                          anchor = anchor, sigma_fixed = sigma_fixed,
                          standardize = standardize,
                          use_boundaries = use_boundaries,
                          beta_init = beta_init, beta0_init = beta0_init),
      none   = bbqr_none(formula, data, quantile = q, ndraw = ndraw,
                         burn = burn, keep = keep, prior = prior,
                         anchor = anchor, standardize = standardize,
                         use_boundaries = use_boundaries,
                         beta_init = beta_init, beta0_init = beta0_init))
    out$call <- cl
    out
  }

  if (length(quantile) == 1L) return(fit_one(quantile))

  res <- lapply(quantile, fit_one)
  names(res) <- paste0("tau=", format(quantile, trim = TRUE))
  attr(res, "quantiles") <- quantile
  attr(res, "penalty")   <- penalty
  attr(res, "call")      <- cl
  class(res) <- "bbqr.list"
  res
}
