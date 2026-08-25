## Methods for class `cbqr`.
##
## `cbqr` objects inherit from `bbqr` so that generics with no response-specific
## content (currently only `coef`) are shared. Everything else MUST be defined
## here, because the binary methods are wrong for a continuous fit in two
## different ways and one of them is silent:
##
##   * print/summary call .bbqr_anchor_label(x$anchor) on an object that has no
##     anchor, and error with "EXPR must be a length 1 vector". Loud, harmless.
##   * predict.bbqr's "response" type returns P(y = 1 | x) computed from the ALD
##     CDF. On a continuous fit with a linear predictor far from zero that
##     saturates at 1 and returns a vector of ones WITHOUT WARNING. This is the
##     dangerous one and is the reason predict.cbqr exists rather than being
##     left to inherit.

#' @export
print.cbqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Bayesian quantile regression (continuous response)\n\n")
  if (!is.null(x$call)) { cat("Call:\n"); print(x$call); cat("\n") }
  cat(sprintf("  Quantile        %g\n", x$quantile))
  cat(sprintf("  Penalty         %s\n", .cbqr_penalty_label(x$penalty)))
  if (!is.na(x$q))
    cat(sprintf("  Scale exponent  q = %g%s\n", x$q, .cbqr_q_note(x)))
  cat(sprintf("  Observations    %d\n", x$n))
  cat(sprintf("  Draws           %d kept of %d (burn %d, thin %d)\n",
              x$ndraw_kept, x$ndraw, x$burn, x$keep))
  cat(sprintf("  Standardized    X %s, y %s\n",
              if (isTRUE(x$standardized$x)) "yes" else "no",
              if (isTRUE(x$standardized$y)) "yes" else "no"))
  cat(sprintf("  sigma (inverse scale)  mean %.4g\n", mean(x$sigma)))
  cat("\nPosterior mean coefficients:\n")
  print(round(stats::coef(x), digits))
  invisible(x)
}

.cbqr_penalty_label <- function(p)
  switch(p, none = "none", lasso = "lasso (shared lambda^2)",
         alasso = "adaptive lasso (lambda_j^2, omega/delta sampled)", p)

## Flag the two published exponents so a non-default q is never mistaken for
## the reference specification.
.cbqr_q_note <- function(x) {
  d <- .cbqr_default_q(x$penalty)
  if (isTRUE(all.equal(x$q, d))) " (layer default)" else
    sprintf(" (NON-DEFAULT; the %s layer's published value is q = %g)",
            x$penalty, d)
}

#' Summarise a fitted continuous quantile regression
#'
#' Posterior means, standard deviations and equal-tailed credible intervals for
#' every coefficient, plus the ALD inverse-scale. A coefficient whose interval
#' excludes zero is flagged.
#'
#' Unlike the binary summary, there is no identification caveat: a continuous
#' response identifies the scale, so the coefficients are on the scale of `y`
#' and are directly interpretable.
#'
#' @param object An object of class `"cbqr"`.
#' @param level Credible level for the intervals.
#' @param ... Ignored.
#' @return An object of class `"cbqr.summary"`.
#' @export
summary.cbqr <- function(object, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("'level' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  draws <- cbind(`(Intercept)` = object$beta0, object$beta)

  qs <- t(apply(draws, 2L, stats::quantile, probs = probs))
  tab <- data.frame(
    mean     = colMeans(draws),
    sd       = apply(draws, 2L, stats::sd),
    lower    = qs[, 1L],
    upper    = qs[, 2L],
    selected = !(qs[, 1L] <= 0 & qs[, 2L] >= 0),
    row.names = colnames(draws))

  sq <- stats::quantile(object$sigma, probs = probs)
  sigma_row <- data.frame(mean = mean(object$sigma), sd = stats::sd(object$sigma),
                          lower = sq[1L], upper = sq[2L], row.names = "sigma")

  out <- list(coefficients = tab, sigma = sigma_row, level = level,
              quantile = object$quantile, penalty = object$penalty, q = object$q,
              n = object$n, ndraw = object$ndraw, burn = object$burn,
              keep = object$keep, ndraw_kept = object$ndraw_kept,
              accept = object$accept, standardized = object$standardized,
              call = object$call)
  class(out) <- "cbqr.summary"
  out
}

#' @export
print.cbqr.summary <- function(x, digits = max(3L, getOption("digits") - 3L),
                               ...) {
  cat("Bayesian quantile regression (continuous response)\n\n")
  if (!is.null(x$call)) { cat("Call:\n"); print(x$call); cat("\n") }
  cat(sprintf("  Quantile        %g\n", x$quantile))
  cat(sprintf("  Penalty         %s\n", .cbqr_penalty_label(x$penalty)))
  if (!is.na(x$q)) cat(sprintf("  Scale exponent  q = %g\n", x$q))
  cat(sprintf("  Observations    %d\n", x$n))
  cat(sprintf("  Draws           %d kept of %d (burn %d, thin %d)\n",
              x$ndraw_kept, x$ndraw, x$burn, x$keep))
  if (!is.null(x$accept) && !is.na(x$accept))
    cat(sprintf("  MH acceptance   %.1f%%\n", 100 * x$accept))

  cat(sprintf("\nPosterior summary (%.0f%% credible intervals):\n",
              100 * x$level))
  tab <- x$coefficients
  disp <- data.frame(
    Mean  = round(tab$mean, digits),
    SD    = round(tab$sd, digits),
    Lower = round(tab$lower, digits),
    Upper = round(tab$upper, digits),
    `  `  = ifelse(tab$selected, "*", " "),
    row.names = rownames(tab), check.names = FALSE)
  print(disp)
  cat("\n* interval excludes zero\n")

  cat(sprintf("\nALD inverse scale: mean %.4g, %.0f%% CI [%.4g, %.4g]\n",
              x$sigma$mean, 100 * x$level, x$sigma$lower, x$sigma$upper))
  if (isTRUE(x$standardized$y))
    cat("Coefficients and sigma are reported on the ORIGINAL y scale.\n")
  invisible(x)
}

#' Predictions from a continuous quantile regression
#'
#' Returns the fitted conditional quantile of the response, which is the
#' linear predictor. That is a different quantity from what [predict.bbqr()]
#' returns -- `P(y = 1 | x)` from the ALD distribution function -- which on a
#' continuous fit saturates at 1 and is silently meaningless.
#'
#' @param object An object of class `"cbqr"`.
#' @param newdata A data frame of predictors. Required: the fitted object does
#'   not store the design matrix.
#' @param interval Return posterior credible intervals for the fitted quantile
#'   as well as the point estimate. Uses the full retained draws, so the cost is
#'   `nrow(newdata) * ndraw_kept` doubles.
#' @param level Credible level when `interval = TRUE`.
#' @param ... Ignored.
#' @return A numeric vector, or a matrix with columns `fit`, `lower`, `upper`
#'   when `interval = TRUE`.
#' @export
predict.cbqr <- function(object, newdata = NULL, interval = FALSE,
                         level = 0.95, ...) {
  if (is.null(newdata))
    stop("'newdata' must be supplied; the fitted object does not store the ",
         "design matrix.", call. = FALSE)
  tt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(tt, newdata, xlev = NULL)
  X  <- stats::model.matrix(tt, mf)
  if (!"(Intercept)" %in% colnames(X)) X <- cbind(`(Intercept)` = 1, X)
  cf <- stats::coef(object)
  miss <- setdiff(names(cf), colnames(X))
  if (length(miss))
    stop("'newdata' is missing required columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  Xm <- X[, names(cf), drop = FALSE]

  if (!isTRUE(interval)) return(as.vector(Xm %*% cf))

  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("'level' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  D <- cbind(object$beta0, object$beta)          # draws x (1 + nvar)
  eta <- Xm %*% t(D)                             # obs x draws
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  qs <- t(apply(eta, 1L, stats::quantile, probs = probs))
  out <- cbind(fit = rowMeans(eta), lower = qs[, 1L], upper = qs[, 2L])
  rownames(out) <- rownames(newdata)
  out
}

#' Trace and density plots for a continuous quantile regression
#'
#' Same interface as [plot.bbqr()], with one addition: `sigma` is a plottable
#' parameter here. In the binary model it is constant under the default
#' `anchor = "sigma1"`, so there is usually nothing to look at; here it is
#' always a live parameter.
#'
#' @param x An object of class `"cbqr"`.
#' @param which Coefficient names or indices; defaults to all, plus `sigma`.
#' @param type `"both"`, `"trace"` or `"density"`.
#' @param ... Passed to the underlying plot calls.
#' @return `x`, invisibly. Called for the side effect of drawing.
#' @seealso [plot.bbqr()] for the binary counterpart.
#' @export
plot.cbqr <- function(x, which = NULL, type = c("both", "trace", "density"),
                      ...) {
  type <- match.arg(type)
  draws <- cbind(`(Intercept)` = x$beta0, x$beta, sigma = x$sigma)
  if (is.null(which)) which <- colnames(draws)
  if (is.numeric(which)) which <- colnames(draws)[which]
  bad <- setdiff(which, colnames(draws))
  if (length(bad))
    stop("Unknown parameter(s): ", paste(bad, collapse = ", "), call. = FALSE)

  ncol_pan <- if (type == "both") 2L else 1L
  op <- graphics::par(mfrow = c(length(which), ncol_pan),
                      mar = c(4, 4, 2, 1))
  on.exit(graphics::par(op), add = TRUE)

  for (nm in which) {
    v <- draws[, nm]
    if (type %in% c("trace", "both"))
      graphics::plot(v, type = "l", xlab = "Iteration", ylab = nm,
                     main = paste("Trace:", nm), ...)
    if (type %in% c("density", "both")) {
      d <- stats::density(v)
      graphics::plot(d, main = paste("Density:", nm), xlab = nm, ...)
      ## Zero is a meaningful reference for a coefficient and not for a scale.
      if (nm != "sigma")
        graphics::abline(v = 0, lty = 2, col = "grey50")
    }
  }
  invisible(x)
}

#' @export
print.cbqr.list <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Bayesian quantile regression (continuous response)\n")
  cat(sprintf("%d quantiles: %s\n\n", length(x),
              paste(vapply(x, function(f) format(f$quantile), character(1)),
                    collapse = ", ")))
  print(round(stats::coef(x), digits))
  invisible(x)
}

#' @export
summary.cbqr.list <- function(object, level = 0.95, ...)
  lapply(object, summary, level = level)
