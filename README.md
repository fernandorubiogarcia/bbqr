# bbqr

<!-- badges: start -->
[![License: GPL (>= 2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%202%29-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![R >= 4.2](https://img.shields.io/badge/R-%3E%3D%204.2-blue.svg)](https://cran.r-project.org/)
<!-- badges: end -->

Bayesian quantile regression with three penalty layers behind a single
interface, for two kinds of response:

| response | entry point | scale |
|---|---|---|
| **binary**, observed 0/1 through a threshold | `bbqr()` | not identified; five `anchor`s fix it |
| **continuous**, observed directly | `cbqr()` | identified; `sigma` is sampled |

In the binary family each penalty also ships in two prior hierarchies — the one
as published and the one whose posterior exists. The continuous family ships
the second only.

Most of what follows is written for the binary case, which is the harder one
and the reason the package exists. [Continuous responses](#continuous-responses)
covers what changes.

Quantile regression on a binary outcome has a property that continuous quantile
regression does not: the model

$$y_i = \mathbb{1}\\{\beta_0 + x_i'\beta + \varepsilon_i > 0\\}$$

is unchanged if $\beta_0$, $\beta$ and the error scale are all multiplied by the
same positive constant. The data identify the *direction* of $\beta$, never its
length, so every method has to pin the scale down somehow.

Existing packages each make that choice internally and do not expose it. `bbqr`
makes it an argument, so the choice can be compared rather than assumed. It
also separates that question — *identification*, which only a constraint can
settle — from a second one the published hierarchies get wrong: *propriety*,
which only the prior can settle.

## Installation

```r
# install.packages("remotes")
remotes::install_github("fernandorubiogarcia/bbqr")
```

Building from source requires a Fortran compiler (Rtools on Windows,
`gfortran` elsewhere).

## Usage

```r
library(bbqr)

set.seed(1)
n <- 300
X <- matrix(rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
y <- as.numeric(X %*% c(1.5, -1, 0.8, 0, 0) + rnorm(n) > 0)
dat <- data.frame(y = y, X)

fit <- bbqr(y ~ x1 + x2 + x3 + x4 + x5, data = dat,
            quantile = 0.5, penalty = "alasso",
            ndraw = 10000, burn = 2000)

fit
#> Bayesian binary quantile regression
#>   Quantile        0.5
#>   Penalty         adaptive lasso
#>   Hierarchy       model = "v5"
#>   Identification  sigma = 1
#>   ...
summary(fit)
```

`summary()` reports posterior means, credible intervals, and flags coefficients
whose interval excludes zero. `coef()`, `predict()` and `plot()` behave as
expected; `quantile` accepts a vector to fit several levels at once.

The defaults — `model = "v5"` on the prior and `anchor = "sigma1"` — sample a
posterior that is proper unconditionally and identified by construction. Both
are explained below.

## Penalties

The three penalties share the latent-variable augmentation of Kozumi and
Kobayashi (2011) and differ only in what sits above the slopes.

| `penalty` | Shrinkage | Implements |
|---|---|---|
| `"none"` | vague normal prior, no shrinkage | Benoit and Van den Poel (2012) |
| `"lasso"` | one parameter shared by all slopes | Benoit, Al-Hamzawi and Yu (2013) |
| `"alasso"` | one parameter per slope | Rubio Garcia (2023) |

## Two hierarchies per penalty

Each penalty ships in two forms, chosen by `model` on the `prior()` object.

**`model = "v3"` is the hierarchy as published.** For the adaptive lasso and
the lasso its joint posterior does not exist. The global shrinkage scale
($\omega$ for the adaptive lasso, $\tau_h$ for the lasso) carries the
log-uniform hyperprior $p(\omega) \propto 1/\omega$, and along the ray where
$\omega$, the $\lambda_j^2$ and the latent $s_j$ vanish together with
$\beta_j = O(\sqrt{\omega})$, the prior measure is $d\omega/\omega$ while the
likelihood stays bounded below. The normalising integral diverges, for every
prior on $\delta$ and at either anchor. A Gibbs sampler on such a target runs
without complaint; what it produces is a function of chain length and seed,
not of the data and prior. For `penalty = "none"` the published flat intercept
is proper if and only if both outcome classes occur.

**`model = "v5"`, the default, replaces every improper component with a proper
one:** $\mathrm{Gamma}(2, 2)$ on the global shrinkage scale (and, for the
lasso, on $\delta$), and a $\tau$-calibrated normal prior on the intercept —
the moment match to the prior under which $\Pr(y = 1)$ is exactly uniform at
$x'\beta = 0$, so that it is centred away from zero when $\tau \neq 0.5$. The
observed likelihood is bounded by one, so the posterior is proper with no
condition on the design, on the number of covariates, or on the anchor. In the
simulation study behind the package the corrected hierarchy beat the published
one on every accuracy, calibration and mixing metric, paired scenario by
scenario.

`"v4"` is the intermediate step for the adaptive lasso — proper $\omega$, flat
intercept — and is proper only conditionally.

```r
# the published hierarchy, for reproducing earlier results
old <- bbqr(y ~ x1 + x2 + x3 + x4 + x5, data = dat,
            prior = prior("alasso", model = "v3"))
prior("alasso", model = "v3")
#> bbqr prior for penalty = 'alasso' , model = 'v3'  [IMPROPER: p(omega) prop 1/omega]
```

Under `"v4"` and `"v5"` the package refuses every device its derivations
exclude: the numerical clamps, the sampler-argument floors, the historical
inverse-Gaussian routine for the latents, and the $\exp(-\kappa/\lambda_j^2)$
tilt. Those devices exist because the published hierarchy needed them — a
numerical guard that fires often is a measurement about the model, not about
the arithmetic — and `"v3"` keeps them reachable so that earlier runs can be
reproduced.

`prior("alasso")` uses a proper $\mathrm{Gamma}(2, 2)$ prior on $\delta$ under
every `model`, and requires `beta_delta > 0`. With $k$ adaptive penalties the
marginal posterior of $\delta$ tends to a positive constant, so the flat prior
leaves it non-integrable; any positive rate restores propriety. Since $\delta$
is the shape of the inverse-Gamma prior on $\lambda_j^2$, it controls how
dispersed the per-coefficient shrinkage is, and both $\delta \to 0$ and
$\delta \to \infty$ are degenerate — hence a prior that vanishes at both ends.

## Identification anchors

| `anchor` | Restriction |
|---|---|
| `"sigma1"` | hold the ALD inverse-scale at 1 *(default)* |
| `"free"` | impose nothing beyond the prior; sample $\sigma$ |
| `"beta1"` | hold the first slope at 1 |
| `"norm1"` | rescale so $\lVert(\beta_0, \beta)\rVert = 1$ |
| `"normslopes"` | rescale so $\lVert\beta\rVert = 1$, slopes only |

All five are available for every penalty except `penalty = "none"` with
`anchor = "free"`, which is refused: with neither a penalty nor an anchor the
likelihood is flat along the scale ray and nothing identifies the model.

`"sigma1"` is the default because it closes the scale ray outright. With the
ray closed, four hierarchy arms in the simulation study agreed on raw RMSE to
within 4%; with it open (`"free"`) they spread over 4.5×, and every bit of that
spread was impropriety rather than non-identification. `"free"` has equally
good point estimates but credible intervals about half again as wide at the
median, because a sampled $\sigma$ lets the unidentified scale inflate them.
It is the sampled-$\sigma$ convention of the adaptive-lasso and lasso
derivations, and it stays one argument away. **Identification is fixed by the
anchor, propriety by the prior. Both are needed; neither substitutes.**

```r
# compare anchors on the same data
sapply(c("sigma1", "free", "beta1", "norm1", "normslopes"), function(a)
  coef(bbqr(y ~ x1 + x2 + x3 + x4 + x5, dat, anchor = a,
            ndraw = 4000, burn = 1000)))
```

The three projection anchors rescale $\beta_0$ and the ALD inverse-scale
$\sigma$ alongside the slopes after every sweep. At the observed-data level the
model is invariant under $(\beta_0, \beta, \sigma) \mapsto (c\beta_0, c\beta,
\sigma/c)$, but that invariance does **not** extend to the penalised posterior,
because the shrinkage layer is stated on a fixed scale: the projection moves
the coefficients relative to the penalty every sweep, and the chain targets no
written-down model. Under `"v4"` and `"v5"` such a fit therefore warns and is
recorded with `derived = FALSE`. The projections are identification
experiments, not steps in any derivation; they are offered because the
comparison is the point of exposing the anchor. In practice they are harmless
when the prior on $\beta$ is vague and damaging when it is informative on the
scale the coefficients occupy. `vignette("bbqr")` works through this.

## Two things to know

**Scale your covariates.** The prior and the penalty are stated on the scale of
the design matrix as supplied, so covariates in very different units are shrunk
by very different relative amounts. Use `standardize = TRUE` unless they already
share a scale. This is a requirement rather than a convenience under `"v5"`:
the only prior on $\omega$ that would make the hierarchy equivariant in the
predictor scale is the improper one.

**Compare directions, not magnitudes.** Coefficients are identified only up to a
positive scale, so their magnitudes are interpretable only relative to the anchor
you chose. Ratios and sign patterns carry across anchors; levels do not.

## Continuous responses

When the response is observed rather than thresholded, use `cbqr()`.
The three penalties carry over, and so does the corrected hierarchy. The
`prior()` constructor and the `model` switch do **not**: continuous fits
take a plain named list, and ship one hierarchy rather than two, because
there is no published continuous specification being reproduced alongside
it. What else changes is forced by the response being observed.

```r
fit <- cbqr(y ~ x1 + x2 + x3, data = d, quantile = 0.5, penalty = "alasso")
```

**There is no `anchor`, and no `sigma_fixed`.** A continuous response
identifies the scale, so `sigma` is sampled every sweep and the anchors would
be misspecifications rather than competing conventions. The unpenalised kernel
*gains* a `sigma` draw, which its binary counterpart does not have at all.

**A new `q`, and its default differs by penalty.** `q` is the exponent on
`sigma` in the local-scale prior rate, $\sigma^q / (2\lambda_j^2)$. The
defaults are inherited rather than chosen: the published lasso penalty is
$\lambda|\beta_j|$, giving `q = 0`, and the published adaptive penalty is
$(\sqrt{\sigma}/\lambda_j)|\beta_j|$, giving `q = 1`. With `sigma` anchored the
two produce the same sampler, so the difference never surfaced in the binary
model; here it does. `q = 2` is the response-scale equivariant choice and draws
`sigma` by slice sampling rather than conjugately.

**`standardize` defaults to `TRUE` and scales `y` as well as `X`.** No proper
hierarchy is exactly response-scale equivariant, so fixing the scale is part of
the specification rather than a preprocessing convenience. Draws are mapped back
to the original units exactly.

**`predict()` returns the fitted conditional quantile**, not a probability — a
`cbqr` fit gets its own method because the binary `predict()` returns
$P(y = 1 \mid x)$ from the ALD distribution function, which saturates at 1 on a
continuous response.

## Relationship to bayesQR and Brq

`bayesQR` and `Brq` both fit binary quantile regression, and both fix the
identification convention internally. `bbqr` differs in exposing it, in offering
the three penalties through one interface, in shipping a hierarchy whose
posterior exists alongside the published one, and in sampler design: the
unpenalised kernel uses a Gibbs step via the Kozumi–Kobayashi augmentation where
Benoit and Van den Poel (2012) use Metropolis–Hastings. It targets the same
posterior by a different algorithm.

The Fortran RNG bridge (`src/wrapper.c`) and the package layout derive from
`bayesQR`, which is GPL (>= 2) and is credited in `Authors@R`. The three MCMC
kernels
were written for this package; see `inst/COPYRIGHTS` for the full breakdown.

## References

Rubio Garcia, F. (2023). *Bayesian adaptive lasso binary quantile regression
with hybrid resampling for classification of imbalanced data.* M.S. thesis,
Wichita State University, Dept. of Mathematics, Statistics, and Physics.
[soar.wichita.edu](https://soar.wichita.edu/entities/publication/a2f86232-4704-4ec2-b685-751e7b04ec42)

Benoit, D. F. and Van den Poel, D. (2012). Binary quantile regression: a
Bayesian approach based on the asymmetric Laplace distribution. *Journal of
Applied Econometrics*, 27(7), 1174–1188.
[doi:10.1002/jae.1216](https://doi.org/10.1002/jae.1216)

Benoit, D. F., Al-Hamzawi, R. and Yu, K. (2013). Bayesian lasso binary quantile
regression. *Computational Statistics*, 28(6), 2861–2873.
[doi:10.1007/s00180-013-0407-8](https://doi.org/10.1007/s00180-013-0407-8)

Kozumi, H. and Kobayashi, G. (2011). Gibbs sampling methods for Bayesian
quantile regression. *Journal of Statistical Computation and Simulation*,
81(11), 1565–1578.
[doi:10.1080/00949655.2010.496117](https://doi.org/10.1080/00949655.2010.496117)

## License

GPL (>= 2). See `inst/COPYRIGHTS` for third-party attribution.
