# bbqr 0.3.0

## Continuous-response quantile regression

* **The continuous adaptive lasso's `omega` hyperprior now defaults to
  `Gamma(0.1, 0.1)`, not `Gamma(2, 2)`.** `Gamma(a, a)` has mean 1 for every
  `a`, so this holds the prior mean fixed and only widens it -- but the shape
  also fixes the density at the origin, and that is the operative change:
  shape > 1 sends the density to **zero** as `omega -> 0`, shape < 1 sends it
  to infinity. Since small `omega` means strong shrinkage, the old default
  forbade exactly the region the data want. Measured in Experiment C1's
  6 x 2 factorial (3,200 paired scenarios per cell, n = 1000): the posterior
  places ~23% of `omega`'s mass below `1e-3` and `Gamma(2, 2)` permitted
  **none** of it. Loosening moves RMSE `0.0895 -> 0.0843` (-5.8%), coverage
  `0.8008 -> 0.8234` and Winkler `0.7074 -> 0.6664`; it wins at every `tau`
  and never trades one quantile against another.

  **Propriety is unaffected.** Proposition 1C requires only that the
  hyperparameters be strictly positive -- there is no lower bound, and its
  normalising bound does not involve `alpha_omega` or `beta_omega` at all.
  `Gamma(0.1, 0.1)` is exactly as proper as `Gamma(2, 2)`.

  **Scope.** This changes `cbqr_alasso()` only. The continuous *lasso* layer's
  `tau` hyperprior and the *binary* hierarchy in `prior()` are unchanged: the
  factorial covered the continuous adaptive lasso alone, and the binary
  sensitivity study reached a different conclusion about what the binding
  costs. Results from earlier versions will differ; pass
  `prior = list(alpha_omega = 2, beta_omega = 2)` to reproduce them.

* New `cbqr()`, with `cbqr_alasso()`, `cbqr_lasso()` and `cbqr_none()`, fitting
  the three penalty layers to an **observed continuous** response, against
  three new Fortran kernels. The binary path is unchanged.

* There is **no `anchor` and no `sigma_fixed`**. A continuous response
  identifies the scale, so `sigma` is sampled every sweep and the anchors are
  misspecifications rather than competing conventions. The unpenalised kernel
  therefore *gains* a `sigma` draw, which its binary counterpart does not have
  at all.

* New `q` argument: the exponent on `sigma` in the local-scale prior rate,
  `sigma^q / (2 lambda_j^2)`. **The layer defaults differ, and this is
  inherited rather than chosen** -- the published lasso penalty is
  `lambda|beta_j|` (`q = 0`) and the published adaptive penalty is
  `(sqrt(sigma)/lambda_j)|beta_j|` (`q = 1`). In the binary model `sigma` is
  anchored, so the two give the same sampler and the difference is invisible;
  it is not invisible here. `q = 2` is the response-scale equivariant choice,
  and draws `sigma` by slice sampling rather than conjugately.

* `standardize` defaults to `TRUE` and now scales **`y` as well as `X`**. No
  proper hierarchy is exactly response-scale equivariant, so fixing the scale
  is part of the specification rather than a preprocessing convenience. Draws are mapped back to the original units exactly.

* `predict()` on a `cbqr` fit returns the fitted conditional quantile. Note
  that the binary `predict()` returns `P(y = 1 | x)` from the ALD distribution
  function, which on a continuous fit saturates at 1 -- `cbqr` objects get
  their own method for exactly this reason.

* `plot()` accepts `sigma` as a plottable parameter, which the binary model has
  nothing to show for.

## The samplers are now tested against the posterior they claim to target

* Geweke's joint-distribution test (Geweke 2004, *JASA* **99**, 799-804) is
  added for the unpenalised binary kernel and for all three continuous
  kernels. Every other test in the suite checks a *consequence* of the sampler
  being right -- constraints hold, coefficients are recovered, draws are
  finite. A kernel with a wrong full conditional passes all of them and still
  converges to the wrong distribution. This one draws the joint
  `p(theta, y)` twice, once ancestrally from the prior and once by a two-block
  Gibbs chain that uses the package's own sweep as its `theta` transition, and
  the two agree only if that sweep leaves the posterior invariant.

* Largest `z` over eleven to fifteen moments at `tau` = 0.3, 0.5, 0.75:
  `bbqr_none` 1.36 / 1.25 / 1.18; `cbqr_none` 2.05 / 2.78 / 2.37;
  `cbqr_lasso` 1.66 / 1.67 / 1.61; `cbqr_alasso` 1.42 / 1.59 / 2.62. The last
  is the one the V5 correction lives in.

* The test is shown to have power rather than merely to pass: a prior shape
  perturbed from 4 to 5 reads 19.6 to 27.4, a quantile off by 0.05 reads 24.7,
  against controls of 1.8 to 2.9.

* The **binary** penalised kernels are still not covered. Their sweeps carry
  `lambda^2`, `s`, `omega` and `delta` between sweeps with no argument to
  initialise them, so a run cannot continue a chain. `cbqr` takes
  `sigma_init`, `delta_init` and `omega_init`/`tau_init` as prior components,
  which is why the continuous penalised kernels can be tested and their binary
  twins cannot.

* Both test files are `skip_on_cran()`. `benchmarks/geweke_joint.R` and
  `benchmarks/geweke_joint_continuous.R` are the full-power versions.

## One package, one description

* The release decision is taken: continuous and binary stay in **one** package
  rather than splitting into siblings, because they share `prior()`, the
  hierarchies, the derivation and the kernels' structure. `DESCRIPTION`'s
  `Title` drops "Binary" and its `Description` now covers both response types;
  `?bbqr-package` gains a *Two response families* section; the README opens
  with the two entry points and gains a continuous section; and the vignette
  gains *When the response is observed*. Four of the package's eight exported
  fitters used to appear nowhere outside `NEWS.md`.

* An audit of every help page against that change found sixteen statements
  that were wrong rather than merely thin, and they are corrected here. Among
  them: `plot.cbqr()` had no `alue`; the continuous kernels were described
  as having "no latent draw" when they still draw `z` and `s` every sweep (it
  is the latent *response* that goes); `prior()` was said to serve both
  families when `cbqr()` rejects it; and the continuous adaptive lasso at its
  default `q = 1` reproduces Alhamzawi, Yu and Benoit (2012), not the 2023
  thesis. The four `cbqr` pages gained examples, references and cross-links,
  so `R CMD check` now exercises the continuous code it previously skipped.

# bbqr 0.2.0

## The default hierarchy is now the one whose posterior exists

* `prior()` defaults to `model = "v5"` for every penalty. The published
  adaptive-lasso and lasso hierarchies (`model = "v3"`) have no posterior:
  the log-uniform hyperprior on the global shrinkage scale (`omega`, or
  `tau_h` for the lasso) is not integrable along the ray where the scale,
  the `lambda_j^2` and the latent `s_j` vanish together, for every prior on
  `delta` and at either anchor. `"v5"` replaces every improper component with
  a proper one -- `Gamma(2, 2)` on the global scale and, for the lasso, on
  `delta`; a `tau`-calibrated normal prior on the intercept -- so the
  posterior is proper with no condition on the design, on `k`, or on the
  anchor. `"v3"` remains available for reproducing earlier results and is
  labelled `IMPROPER` wherever it is printed.

* The default anchor is `"sigma1"` for every penalty (previously `"free"` for
  the adaptive lasso and the lasso). Holding `sigma` at 1 closes the scale
  ray the likelihood never identifies, gives the coefficients a stated scale,
  keeps the credible intervals from inheriting the width of an unidentified
  direction, and is the convention under which the `"v5"` intercept prior is
  exact. `"free"` -- `sigma` sampled from its full conditional, as the
  derivations also allow -- is one argument away and is a derived posterior
  under every `model`. `bbqr_alasso()` and `bbqr_lasso()` changed
  accordingly.

* A `"v5"` lasso prior now requires `beta_delta > 0`, as the adaptive lasso
  already did: a flat prior on `delta` is admissible only under `"v3"`, where
  it is the published specification.

## Provenance on the fitted object

* `print()` and `summary()` report which hierarchy a fit sampled and flag a
  chain that does not target a derived posterior (a projection anchor
  combined with `"v4"` or `"v5"`). Objects saved by earlier builds carry no
  `model` and print no hierarchy line.

## Documentation

* `prior()` documents `model` for all three penalties and the lasso's
  `alpha_tau`/`beta_tau`; `bbqr_alasso()` documents `kappa_lam`; `bbqr()`
  gains a *Derivation version* section and documents the `model`, `derived`
  and `b0_prior` components of the fit. The clip-site count is eighteen,
  not fifteen.

* The README and the vignette describe the two hierarchies and the new
  defaults.

* The help pages now render their code spans, emphasis and cross-links
  (roxygen markdown was never switched on, so these appeared as literal
  backticks, asterisks and brackets), and the Kozumi and Kobayashi (2011)
  DOI is corrected.

## Benchmarks

* The drivers that predate the hierarchy split (`Experiment_v4.R` through
  `Experiment_v9.R`, the real-data and churn benchmarks, and
  `diagnose_hyperparams.R`) now pass `model = "v3"` explicitly, so re-running
  them reproduces what they reported rather than picking up the new default.

# bbqr 0.1.0

* First working version: adaptive-lasso, lasso and unpenalised binary
  quantile regression with five identification anchors behind one interface.
  Neither this nor 0.2.0 was released to CRAN; 0.3.0 is the first submission.
