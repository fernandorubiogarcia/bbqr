#' @keywords internal
#' @aliases bbqr-package
#'
#' @section Two response families:
#' The package fits the same three penalty layers to two kinds of response,
#' through two parallel entry points.
#'
#' \describe{
#'   \item{[bbqr()]}{A **binary** response, observed only as 0/1 through a
#'     threshold. The latent response is drawn each sweep. Because the
#'     threshold identifies the coefficient vector only up to a positive
#'     scale, these fits take an `anchor` that says how the scale is fixed.
#'     `bbqr_alasso()`, `bbqr_lasso()` and `bbqr_none()` are the direct forms.}
#'   \item{[cbqr()]}{An **observed continuous** response. There is no latent
#'     *response* draw and no `anchor`: the data identify the scale, so `sigma`
#'     is sampled every sweep. The location-scale latents `z` and `s` are still
#'     drawn. These fits instead take `q`, the exponent on
#'     `sigma` in the local-scale prior rate, which the binary model cannot
#'     see because its `sigma` is held fixed. `cbqr_alasso()`, `cbqr_lasso()`
#'     and `cbqr_none()` are the direct forms.}
#' }
#'
#' The two families do **not** share a prior constructor. [prior()] builds a
#' `bbqr.prior` object for the binary fitters; the continuous fitters take a
#' plain named list, because their hierarchy has different components (no
#' anchor, a live `sigma`, and the `q` exponent) and different reference
#' defaults. Passing a [prior()] object to `cbqr()` is an error, and says so.
#'
#' @section Method sources:
#' The binary samplers implement Rubio Garcia (2023) for the adaptive lasso,
#' and Benoit, Al-Hamzawi and Yu (2013) and Benoit and Van den Poel (2012) for
#' the lasso and unpenalised layers. On the continuous side the published
#' specification each layer reproduces is a different one: at its default
#' `q = 1` the adaptive lasso is the penalty of Alhamzawi, Yu and Benoit
#' (2012), the lasso at `q = 0` is Benoit, Al-Hamzawi and Yu (2013) carried
#' over to an observed response, and the unpenalised layer is the sampler of
#' Kozumi and Kobayashi (2011). See `citation("bbqr")`.
#'
#' **In the binary family**, each of the three penalties is shipped in two
#' forms, selected by `prior(penalty, model = )`: `"v3"`, the hierarchy as
#' published, and `"v5"`, the default, in which every improper prior component
#' has been replaced by a proper one so that the posterior exists
#' unconditionally. The published adaptive-lasso and lasso hierarchies have no
#' posterior, for the reason given in [prior()]; they are kept so that earlier
#' results can be reproduced, and the print method labels them. The continuous
#' family ships one hierarchy, the proper one: there is no `model` argument and
#' no published continuous specification being reproduced alongside it.
#'
#' @section Getting started:
#' [bbqr()] is the entry point for a 0/1 response; see its help page for the
#' three penalties, the identification anchors and the derivation versions.
#' [cbqr()] is the entry point for a continuous one. [prior()] sets
#' hyperparameters and selects the hierarchy for the binary fitters only; the
#' continuous fitters take a plain named list, documented at [cbqr()].
#'
#' @section Relationship to bayesQR:
#' The Fortran RNG wrapper (`src/wrapper.c`) and the general layout of the
#' package follow Benoit's \pkg{bayesQR}, which is GPL (>= 2) licensed and is
#' acknowledged in the `Authors@R` field of `DESCRIPTION`. The six MCMC
#' kernels shipped here -- three binary, three continuous -- were written for
#' this package.
#'
#' @importFrom stats coef density model.frame model.matrix model.response
#'   quantile sd setNames delete.response
#' @importFrom graphics abline par
#' @importFrom utils tail
#' @useDynLib bbqr, .registration = TRUE
"_PACKAGE"
