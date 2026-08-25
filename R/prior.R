#' Prior and hyperparameter settings for the binary fitters
#'
#' Builds the list of prior hyperparameters and Metropolis tuning constants for
#' one of the three penalty layers.
#'
#' This object is accepted by \code{\link[=bbqr]{bbqr()}} and the three direct
#' binary fitters, and by those only. The continuous fitters
#' (\code{\link[=cbqr]{cbqr()}} and friends) take a \strong{plain named list},
#' not a \code{bbqr.prior}: their hierarchy has no anchor, a live
#' \eqn{\sigma}, and the \code{q} exponent, and its reference defaults differ.
#' Passing the result of this function to \code{cbqr()} fails with
#' \dQuote{Unknown prior component(s)}; \code{\link[=cbqr]{cbqr()}} documents
#' the components they do accept. The defaults are the \code{"v5"}
#' hierarchy, in which every prior component is a probability distribution and
#' the posterior is proper unconditionally; \code{model = "v3"} reproduces the
#' specification of the source paper for each penalty instead. In most cases
#' this function only needs to be called to change one or two values.
#'
#' @param penalty Which sampler the prior is for: \code{"alasso"} (adaptive
#'   lasso), \code{"lasso"}, or \code{"none"}.
#' @param model Which prior hierarchy to sample. The three are nested
#'   corrections, each removing one failure of propriety. For the adaptive
#'   lasso they differ only in the prior on \eqn{\omega} and the prior on the
#'   intercept; for the lasso, in the priors on \eqn{\tau_h}, \eqn{\delta} and
#'   the intercept; for the unpenalised model, in the prior on the intercept
#'   alone. \code{"v3"} and \code{"v5"} are defined for every penalty;
#'   \code{"v4"} exists only for \code{"alasso"}, because \eqn{\omega} does.
#'
#'   \describe{
#'     \item{\code{"v5"}}{The default. Every prior component is proper: for
#'       \code{"alasso"}, \code{"v4"} plus \eqn{\beta_0 \sim N(m_{00},
#'       V_{00})}; for \code{"lasso"}, \eqn{\mathrm{Gamma}(2, 2)} priors on
#'       \eqn{\tau_h} and \eqn{\delta} plus the same intercept prior; for
#'       \code{"none"}, the intercept prior alone. The observed likelihood is
#'       bounded by one, so the posterior is proper unconditionally: no
#'       condition on the design, on \eqn{k}, on the anchor, or on \eqn{a}
#'       beyond positivity.}
#'     \item{\code{"v4"}}{\eqn{\omega \sim \mathrm{Gamma}(2, 2)} with a flat
#'       intercept; \code{"alasso"} only. Removes the \eqn{\omega} defect.
#'       Propriety is then conditional: with a free \eqn{\sigma} it requires
#'       \eqn{a > 1} plus an overlap condition on the design, which is why the
#'       default here is \code{a = b = 2} rather than the historical
#'       \code{a = b = 1}. At \code{a = 1} the fit warns.}
#'     \item{\code{"v3"}}{The published specifications, retained as a
#'       baseline. For \code{"alasso"}: \eqn{p(\omega) \propto 1/\omega} with
#'       a flat intercept, the hierarchy Experiments v4 and v6 ran.
#'       \strong{Its joint posterior does not exist.} The prior hierarchy is a
#'       scale family in \eqn{\omega}: on the wedge where \eqn{\omega},
#'       \eqn{\lambda_j^2} and \eqn{s_j} vanish together with
#'       \eqn{\beta_j = O(\sqrt\omega)} the prior measure is
#'       \eqn{d\omega/\omega} while the likelihood is bounded below, so the
#'       normalizing integral diverges. That holds for every prior on
#'       \eqn{\delta} and at either anchor. For \code{"lasso"}: the hierarchy
#'       of Benoit, Al-Hamzawi and Yu (2013), \eqn{p(\tau_h) \propto 1/\tau_h}
#'       and \eqn{p(\delta) \propto 1}, which carries the same defect along
#'       \eqn{\tau_h \to 0}. For \code{"none"}: the flat intercept of Benoit
#'       and Van den Poel (2012), whose posterior is proper if and only if
#'       both outcome classes occur.}
#'   }
#'
#'   \code{"v4"} and \code{"v5"} are defined by their derivations, which state
#'   that no clamp, sampler-argument floor, or
#'   \eqn{\exp(-\kappa_\lambda/\lambda_j^2)} tilt is part of the target. A fit
#'   that names them and enables such a device is refused rather than silently
#'   recorded; see \code{\link{bbqr_alasso}}. \code{"v3"} leaves every switch
#'   reachable, because reproducing the historical build needs them.
#' @param alpha_omega,beta_omega Shape and rate of the \eqn{\mathrm{Gamma}}
#'   prior on the global shrinkage scale \eqn{\omega}, whose full conditional
#'   is \eqn{\mathrm{Gamma}(\alpha_\omega + k\delta, \beta_\omega + \sum_j
#'   \lambda_j^{-2})}. Both zero is \eqn{p(\omega) \propto 1/\omega} exactly,
#'   which is \code{model = "v3"}; \code{"v4"} and \code{"v5"} require both
#'   strictly positive. Defaults are \code{(0, 0)} for \code{"v3"} and
#'   \code{(2, 2)} otherwise.
#' @param b0_mean,b0_var Mean and variance of the \eqn{N(m_{00}, V_{00})} prior
#'   on the intercept. \code{b0_var = Inf} is the flat intercept of
#'   \code{"v3"} and \code{"v4"}, carried into the kernel as a prior precision
#'   of exactly zero rather than as a large finite variance. \code{NA}, the
#'   \code{"v5"} default, means the \eqn{\tau}-calibrated reference, which
#'   cannot be resolved here because \eqn{\tau} is not known until the fit.
#'
#'   That reference is the moment match to the exactly calibrated prior. At the
#'   \eqn{\sigma = 1} anchor with \eqn{x'\beta = 0} the model reduces to
#'   \eqn{\Pr(y = 1) = G_\tau(\beta_0)} with \eqn{G_\tau} the distribution
#'   function of \eqn{-E} for \eqn{E \sim \mathrm{ALD}(0, 1, \tau)}, so by the
#'   probability integral transform \eqn{G_\tau(\beta_0)} is exactly uniform
#'   when \eqn{\beta_0 \sim \mathrm{ALD}(0, 1, 1 - \tau)}. That law is not
#'   conjugate; its first two moments are
#'   \eqn{m_{00} = -\theta = -(1 - 2\tau)/(\tau(1-\tau))} and
#'   \eqn{V_{00} = (1 - 2\tau + 2\tau^2)/(\tau^2(1-\tau)^2)}, giving
#'   \eqn{(0, 8)} at the median.
#'
#'   Note that \eqn{m_{00}} is not zero away from the median. Centring the
#'   intercept prior at the origin puts the prior-predictive mean of
#'   \eqn{\Pr(y = 1)} at 0.73 rather than 0.5 at \eqn{\tau = 0.1}, so it is a
#'   specification error and not a neutral default.
#' @param a,b Shape and rate of the \eqn{\mathrm{Gamma}(a, b)} prior on the ALD
#'   inverse-scale \eqn{\sigma}.  Ignored when the inverse-scale is anchored
#'   (\code{anchor = "sigma1"}, the default) and, \strong{in the binary model},
#'   for \code{penalty = "none"}, where \eqn{\sigma \equiv 1} is the
#'   identification convention.  The continuous unpenalised kernel samples
#'   \eqn{\sigma} and does use \code{a} and \code{b}.  Defaults are
#'   \code{a = b = 2} under \code{"v4"} and \code{"v5"}, the mean-one reference
#'   that also satisfies the \eqn{a > 1} condition a flat intercept needs
#'   against a free \eqn{\sigma}; under \code{"v3"} they follow the source
#'   papers, \code{a = b = 1} for the adaptive lasso and \code{a = b = 0.1} for
#'   the lasso.
#' @param alpha_delta,beta_delta Shape and rate of the
#'   \eqn{\mathrm{Gamma}}{Gamma} prior on the penalty hyperparameter
#'   \eqn{\delta}.
#'
#'   For \code{penalty = "lasso"} the default is \code{(2, 2)} under
#'   \code{"v5"}, which requires \code{beta_delta > 0}: the flat prior is not
#'   integrable at the large-\eqn{\delta} end, so it has no place in a
#'   hierarchy whose every component is proper. Under \code{"v3"} the default
#'   is \code{(1, 0)}, the improper flat prior \eqn{p(\delta) \propto 1} of
#'   Benoit, Al-Hamzawi and Yu (2013): at \code{alpha_delta = 1} and
#'   \code{beta_delta = 0} the prior contribution
#'   \eqn{(\alpha-1)\log\delta - \beta\delta} vanishes and the Metropolis step
#'   targets that paper's Eq. (13) exactly.  There is a single shrinkage
#'   parameter in that model.
#'
#'   For \code{penalty = "alasso"} the default is \code{(2, 2)} under every
#'   \code{model}, and \code{beta_delta} must be strictly positive.  With \eqn{k} adaptive
#'   penalties the marginal posterior of \eqn{\delta} tends to a positive
#'   constant, so a flat prior leaves it non-integrable and the joint posterior
#'   does not exist; any positive rate restores propriety, with no condition on
#'   \eqn{k}.  \eqn{\delta} is the shape of the inverse-Gamma prior on
#'   \eqn{\lambda_j^2} and therefore controls how dispersed the per-coefficient
#'   shrinkage is: both \eqn{\delta \to 0} (where the inverse-Gamma family
#'   degenerates) and \eqn{\delta \to \infty} (where the \eqn{\lambda_j^2}
#'   coincide and the penalty stops being adaptive) are degenerate, so the
#'   prior is chosen to vanish at both ends, which requires
#'   \code{alpha_delta > 1}.  At \code{(2, 2)} the prior has mode \eqn{1/2} and
#'   mean \eqn{1}.
#' @param alpha_tau,beta_tau Shape and rate of the \eqn{\mathrm{Gamma}} prior
#'   on the lasso's global shrinkage hyperparameter \eqn{\tau_h};
#'   \code{penalty = "lasso"} only, where \eqn{\tau_h} plays the role
#'   \eqn{\omega} plays in the adaptive lasso. Its full conditional is
#'   \eqn{\mathrm{Gamma}(\alpha_\tau + \delta, \beta_\tau + \lambda^2)}. Both
#'   zero is \eqn{p(\tau_h) \propto 1/\tau_h}, the hyperprior of Benoit,
#'   Al-Hamzawi and Yu (2013) and the \code{"v3"} default; it leaves the joint
#'   posterior improper along \eqn{\tau_h \to 0} by the same mechanism as
#'   \eqn{\omega}. \code{"v5"} requires both strictly positive and defaults to
#'   \code{(2, 2)}.
#' @param delta_init,omega_init,tau_init Starting values for the penalty
#'   hyperparameters.  \code{omega_init} applies to the adaptive lasso only and
#'   \code{tau_init} to the lasso only.
#' @param mh_sd0 Initial standard deviation of the random-walk Metropolis
#'   proposal for \eqn{\delta}, on the log scale.
#' @param target_acc Target acceptance rate for that Metropolis step.  The
#'   proposal standard deviation is adapted towards this value during burn-in
#'   and then held fixed.
#' @param beta_var Prior variance of each slope under \code{penalty = "none"}.
#'   Benoit and Van den Poel (2012) use 100, a vague prior.
#'
#' @return A list of class \code{"bbqr.prior"}.
#' @seealso \code{\link[=bbqr]{bbqr()}} and the direct binary fitters, which
#'   accept this object; \code{\link[=cbqr]{cbqr()}}, which does not.
#'
#' @references
#' Benoit, D. F. and Van den Poel, D. (2012). Binary quantile regression: a
#' Bayesian approach based on the asymmetric Laplace distribution.
#' \emph{Journal of Applied Econometrics}, 27(7), 1174--1188.
#' \doi{10.1002/jae.1216}
#'
#' Benoit, D. F., Al-Hamzawi, R. and Yu, K. (2013). Bayesian lasso binary
#' quantile regression. \emph{Computational Statistics}, 28(6), 2861--2873.
#' \doi{10.1007/s00180-013-0407-8}
#'
#' @examples
#' # Default adaptive-lasso prior: the fully proper "v5" hierarchy
#' prior("alasso")
#'
#' # The published specification, for reproducing earlier results. Its
#' # posterior does not exist, and the print method says so.
#' prior("alasso", model = "v3")
#'
#' # A vaguer prior on the slopes for the unpenalised sampler
#' prior("none", beta_var = 1000)
#'
#' @export
prior <- function(penalty = c("alasso", "lasso", "none"),
                  model = c("v5", "v4", "v3"),
                  a = NULL, b = NULL,
                  alpha_delta = NULL, beta_delta = NULL,
                  alpha_omega = NULL, beta_omega = NULL,
                  alpha_tau = NULL, beta_tau = NULL,
                  b0_mean = NULL, b0_var = NULL,
                  delta_init = 1, omega_init = 1, tau_init = 1,
                  mh_sd0 = 0.05, target_acc = 0.30,
                  beta_var = 100) {
  penalty <- match.arg(penalty)
  model   <- match.arg(model)

  ## v4 is defined by a proper prior on the adaptive lasso's omega, which the
  ## other two penalties do not have: the lasso's analogous hyperparameter is
  ## tau_h and the unpenalised model has none. So v4 is alasso-only, while v3
  ## (as published) and v5 (every component proper) are defined for all three.
  if (penalty != "alasso" && model == "v4")
    stop("model = \"v4\" is the proper-omega correction, and omega exists only ",
         "for penalty = \"alasso\"; got penalty = ", sQuote(penalty),
         ". Use \"v3\" for the published specification or \"v5\" for the fully ",
         "proper one.", call. = FALSE)

  ## Version presets. Each is a point in the (alpha_omega, beta_omega,
  ## b0_prec) parameter space rather than a code branch, so the three
  ## hierarchies share one sampler and differ only in these numbers.
  ##
  ##   v3  p(omega) prop 1/omega, flat intercept.   IMPROPER -- retained only
  ##       as the baseline that Experiments v4 and v6 actually ran.
  ##   v4  omega ~ Gamma(2, 2), flat intercept.     Proper iff a > 1 with a
  ##       free sigma, hence a = b = 2 rather than the historical a = b = 1.
  ##   v5  v4 plus beta0 ~ N(m00, V00).             Proper unconditionally.
  ##
  ## `b0_var = NA_real_` is the sentinel for "the tau-calibrated reference",
  ## which cannot be resolved here because tau is not known until the fit; see
  ## .bbqr_b0_reference().
  ## Under "v5" the intercept prior is proper, so no condition on `a` is
  ## required and the published values would be admissible; they are still
  ## moved to the mean-one reference so that one reading of the specification
  ## covers sigma, omega, tau_h and delta together.
  if (is.null(a))
    a <- switch(penalty, alasso = if (model == "v3") 1 else 2,
                lasso = if (model == "v3") 0.1 else 2, none = 1)
  if (is.null(b))
    b <- switch(penalty, alasso = if (model == "v3") 1 else 2,
                lasso = if (model == "v3") 0.1 else 2, none = 1)
  ## The lasso keeps the flat prior of Benoit, Al-Hamzawi and Yu (2013) under
  ## model = "v3", and takes a proper Gamma(2, 2) under "v5": an improper
  ## component anywhere in the hierarchy leaves the joint posterior undefined,
  ## and "v5" means every component is a probability distribution.
  if (is.null(alpha_delta))
    alpha_delta <- switch(penalty, alasso = 2,
                          lasso = if (model == "v5") 2 else 1, none = 1)
  if (is.null(beta_delta))
    beta_delta <- switch(penalty, alasso = 2,
                         lasso = if (model == "v5") 2 else 0, none = 1)
  if (is.null(alpha_omega)) alpha_omega <- if (model == "v3") 0 else 2
  if (is.null(beta_omega))  beta_omega  <- if (model == "v3") 0 else 2
  ## The lasso's tau_h plays omega's role: its published conditional is
  ## Gamma(delta, lambda^2) with no additive constants, the signature of
  ## p(tau_h) prop 1/tau_h, and it carries the same non-integrable zero-scale
  ## ray. (0, 0) reproduces the published sampler exactly.
  if (is.null(alpha_tau)) alpha_tau <- if (model == "v5") 2 else 0
  if (is.null(beta_tau))  beta_tau  <- if (model == "v5") 2 else 0
  if (is.null(b0_mean)) b0_mean <- if (model == "v5") NA_real_ else 0
  if (is.null(b0_var))  b0_var  <- if (model == "v5") NA_real_ else Inf

  ## Propriety. With k adaptive penalties the marginal posterior of delta tends
  ## to a positive constant, so a flat prior leaves it non-integrable and the
  ## joint posterior does not exist. Any Gamma prior with a positive rate is
  ## enough; no condition on k is required. The single-penalty lasso keeps the
  ## flat prior of Benoit, Al-Hamzawi and Yu (2013), where there is one
  ## shrinkage parameter and nothing for a drifting delta to homogenise.
  if (penalty == "alasso" && !(is.numeric(beta_delta) &&
                               length(beta_delta) == 1L &&
                               is.finite(beta_delta) && beta_delta > 0))
    stop("'beta_delta' must be a single positive number for penalty = \"alasso\": ",
         "the adaptive-lasso posterior is improper under a flat prior on delta. ",
         "The default is beta_delta = 2.", call. = FALSE)

  chk <- function(x, nm) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
      stop(sprintf("'%s' must be a single positive finite number.", nm),
           call. = FALSE)
    as.double(x)
  }
  ## beta_delta alone may be exactly 0, which together with alpha_delta = 1
  ## gives the improper flat prior p(delta) ~ 1 used in the 2013 paper.
  ## alpha_omega and beta_omega may likewise be exactly 0, which together give
  ## V3's p(omega) prop 1/omega.
  chk0 <- function(x, nm) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0)
      stop(sprintf("'%s' must be a single non-negative finite number.", nm),
           call. = FALSE)
    as.double(x)
  }
  ## b0_mean is unrestricted in sign; b0_var may be Inf (the flat intercept of
  ## V3 and V4) or NA (resolve to the tau-calibrated reference at fit time).
  chk_b0m <- function(x) {
    if (!is.numeric(x) || length(x) != 1L || (!is.na(x) && !is.finite(x)))
      stop("'b0_mean' must be a single finite number, or NA for the ",
           "tau-calibrated reference.", call. = FALSE)
    as.double(x)
  }
  chk_b0v <- function(x) {
    if (!is.numeric(x) || length(x) != 1L || (!is.na(x) && x <= 0))
      stop("'b0_var' must be a single positive number, Inf for a flat ",
           "intercept, or NA for the tau-calibrated reference.", call. = FALSE)
    as.double(x)
  }
  out <- list(
    penalty     = penalty,
    model       = model,
    a           = chk(a, "a"),
    b           = chk(b, "b"),
    alpha_delta = chk(alpha_delta, "alpha_delta"),
    beta_delta  = chk0(beta_delta, "beta_delta"),
    alpha_omega = chk0(alpha_omega, "alpha_omega"),
    beta_omega  = chk0(beta_omega, "beta_omega"),
    alpha_tau   = chk0(alpha_tau, "alpha_tau"),
    beta_tau    = chk0(beta_tau, "beta_tau"),
    b0_mean     = chk_b0m(b0_mean),
    b0_var      = chk_b0v(b0_var),
    delta_init  = chk(delta_init, "delta_init"),
    omega_init  = chk(omega_init, "omega_init"),
    tau_init    = chk(tau_init, "tau_init"),
    mh_sd0      = chk(mh_sd0, "mh_sd0"),
    target_acc  = chk(target_acc, "target_acc"),
    beta_var    = chk(beta_var, "beta_var"))

  if (out$target_acc >= 1)
    stop("'target_acc' must lie strictly between 0 and 1.", call. = FALSE)

  ## Propriety of the omega block. With p(omega) prop 1/omega the joint
  ## posterior does not exist, whatever the prior on delta and whichever
  ## anchor is used: on the wedge where omega, lambda_j^2 and s_j vanish
  ## together with beta_j = O(sqrt(omega)) the prior measure is d omega /
  ## omega and the likelihood is bounded below. Both hyperparameters must be
  ## positive together -- alpha_omega alone controls the origin and
  ## beta_omega alone the upper tail.
  if (penalty == "lasso" && model == "v5" &&
      !(out$alpha_tau > 0 && out$beta_tau > 0))
    stop("'alpha_tau' and 'beta_tau' must both be strictly positive for ",
         "model = \"v5\": the published p(tau_h) prop 1/tau_h leaves the joint ",
         "posterior improper along tau_h -> 0. Use model = \"v3\" to reproduce ",
         "that specification deliberately.", call. = FALSE)
  if (penalty == "lasso" && model == "v3" &&
      (out$alpha_tau > 0 || out$beta_tau > 0))
    stop("model = \"v3\" is the published lasso, defined by ",
         "p(tau_h) prop 1/tau_h, i.e. alpha_tau = beta_tau = 0.", call. = FALSE)
  ## The lasso's delta is the one place a flat prior is admissible at all --
  ## under "v3", as published -- and "v5" means every component is proper, so
  ## a v5 prior carrying beta_delta = 0 would be a mislabelled arm.
  if (penalty == "lasso" && model == "v5" && !(out$beta_delta > 0))
    stop("'beta_delta' must be strictly positive for model = \"v5\": the flat ",
         "prior p(delta) prop 1 of the published lasso is not integrable at the ",
         "large-delta end, so the joint posterior is improper. Use model = ",
         "\"v3\" to reproduce that specification deliberately.", call. = FALSE)
  if (penalty == "alasso" && model != "v3" &&
      !(out$alpha_omega > 0 && out$beta_omega > 0))
    stop("'alpha_omega' and 'beta_omega' must both be strictly positive for ",
         "model = ", sQuote(model), ": a log-uniform prior on omega leaves ",
         "the joint posterior improper. Use model = \"v3\" to reproduce that ",
         "specification deliberately.", call. = FALSE)
  if (penalty == "alasso" && model == "v3" &&
      (out$alpha_omega > 0 || out$beta_omega > 0))
    stop("model = \"v3\" is defined by p(omega) prop 1/omega, i.e. ",
         "alpha_omega = beta_omega = 0. Use model = \"v4\" for a proper ",
         "prior on omega.", call. = FALSE)

  ## The intercept prior is what separates V5 from V4, so a mismatch between
  ## the two is a mislabelled arm rather than an unusual setting.
  if (model == "v5" && identical(out$b0_var, Inf))
    stop("model = \"v5\" is defined by a proper intercept prior, so ",
         "'b0_var' cannot be Inf. Use model = \"v4\" for a flat intercept.",
         call. = FALSE)
  if (model != "v5" && !identical(out$b0_var, Inf))
    stop("model = ", sQuote(model), " has a flat intercept, so 'b0_var' must ",
         "be Inf. Use model = \"v5\" for a proper intercept prior.",
         call. = FALSE)

  class(out) <- "bbqr.prior"
  out
}

#' The tau-calibrated reference intercept prior
#'
#' At the `sigma = 1` anchor and with `x'beta = 0`, the model reduces to
#' `Pr(y = 1) = G_tau(beta0)` with `G_tau` the distribution function of `-E`
#' for `E ~ ALD(0, 1, tau)`. By the probability integral transform,
#' `G_tau(beta0)` is exactly `Uniform(0, 1)` when `beta0 ~ ALD(0, 1, 1 - tau)`,
#' so the reference normal prior is that law's moment match:
#'
#'   `mean = -theta = -(1 - 2 tau) / (tau (1 - tau))`
#'   `var  = (1 - 2 tau + 2 tau^2) / (tau^2 (1 - tau)^2)`
#'
#' The mean is not zero away from the median. Centring the intercept prior at
#' the origin puts the prior-predictive mean of `Pr(y = 1)` at 0.73 rather than
#' 0.5 at `tau = 0.1`, so it is a specification error rather than a neutral
#' default. The variance is within 5% of the value that minimises the
#' Kolmogorov-Smirnov distance to `Uniform(0, 1)` over scalings of itself, so
#' no optimisation is performed: the closed form is used as it stands.
#'
#' @param quantile The quantile level `tau`, in `(0, 1)`.
#' @return A list with elements `mean` and `var`.
#' @keywords internal
.bbqr_b0_reference <- function(quantile) {
  tau <- as.double(quantile)
  if (!is.finite(tau) || tau <= 0 || tau >= 1)
    stop("'quantile' must lie strictly between 0 and 1.", call. = FALSE)
  list(mean = -(1 - 2 * tau) / (tau * (1 - tau)),
       var  = (1 - 2 * tau + 2 * tau^2) / (tau^2 * (1 - tau)^2))
}

#' @export
print.bbqr.prior <- function(x, ...) {
  cat("bbqr prior for penalty =", sQuote(x$penalty),
      ", model =", sQuote(x$model))
  cat(switch(paste(x$penalty, x$model),
    "alasso v3" = "  [IMPROPER: p(omega) prop 1/omega]",
    "alasso v4" = "  [proper omega, flat intercept]",
    "alasso v5" = "  [proper omega and intercept]",
    "lasso v3"  = "  [IMPROPER: p(tau_h) prop 1/tau_h, flat delta]",
    "lasso v5"  = "  [proper tau_h, delta and intercept]",
    "none v3"   = "  [flat intercept: proper iff both classes observed]",
    "none v5"   = "  [proper intercept]",
    ""))
  cat("\n\n")
  keep <- switch(x$penalty,
    alasso = c("a", "b", "alpha_delta", "beta_delta",
               "alpha_omega", "beta_omega", "b0_mean", "b0_var",
               "delta_init", "omega_init", "mh_sd0", "target_acc"),
    lasso  = c("a", "b", "alpha_delta", "beta_delta",
               "alpha_tau", "beta_tau", "b0_mean", "b0_var",
               "delta_init", "tau_init", "mh_sd0", "target_acc"),
    none   = c("beta_var", "b0_mean", "b0_var"))
  for (nm in keep) {
    v <- x[[nm]]
    lab <- if (is.na(v)) "<tau-calibrated>" else sprintf("%g", v)
    cat(sprintf("  %-12s %s\n", nm, lab))
  }
  invisible(x)
}
