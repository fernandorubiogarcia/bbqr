## One line naming the hierarchy a fit sampled. Objects from builds before
## the version split carry no `model`; they were the V3 hierarchy, and saying
## nothing is better than guessing on their behalf.
.bbqr_print_hierarchy <- function(x) {
  if (is.null(x$model)) return(invisible(NULL))
  cat(sprintf("  Hierarchy       model = \"%s\"%s\n", x$model,
              if (isTRUE(x$derived)) "" else "  (not a derived posterior)"))
  invisible(NULL)
}

#' @export
print.bbqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Bayesian binary quantile regression\n\n")
  if (!is.null(x$call)) {
    cat("Call:\n"); print(x$call); cat("\n")
  }
  cat(sprintf("  Quantile        %g\n", x$quantile))
  cat(sprintf("  Penalty         %s\n", .bbqr_penalty_label(x$penalty)))
  .bbqr_print_hierarchy(x)
  cat(sprintf("  Identification  %s\n", .bbqr_anchor_label(x$anchor)))
  cat(sprintf("  Observations    %d\n", x$n))
  cat(sprintf("  Draws           %d kept of %d (burn %d, thin %d)\n",
              x$ndraw_kept, x$ndraw, x$burn, x$keep))
  cat("\nPosterior mean coefficients:\n")
  print(round(stats::coef(x), digits))
  invisible(x)
}

#' Extract posterior mean coefficients
#'
#' @param object An object of class `"bbqr"`, `"cbqr"`, `"bbqr.list"` or
#'   `"cbqr.list"`. Continuous fits inherit from `"bbqr"`, so this method
#'   serves them too.
#' @param ... Ignored.
#' @return A named numeric vector of posterior means, intercept first. For a
#'   `"bbqr.list"`, a matrix with one column per quantile.
#' @export
coef.bbqr <- function(object, ...) {
  stats::setNames(c(mean(object$beta0), colMeans(object$beta)),
                  c("(Intercept)", object$names))
}

#' @rdname coef.bbqr
#' @export
coef.bbqr.list <- function(object, ...) {
  m <- vapply(object, stats::coef, numeric(length(object[[1L]]$names) + 1L))
  colnames(m) <- names(object)
  m
}

#' Summarise a fitted binary quantile regression
#'
#' Posterior means, standard deviations and equal-tailed credible intervals for
#' every coefficient. A coefficient whose interval excludes zero is flagged as
#' selected, which is the variable-selection rule used with these samplers.
#'
#' @param object An object of class `"bbqr"`.
#' @param level Credible level for the intervals.
#' @param ... Ignored.
#' @return An object of class `"bbqr.summary"`.
#' @examples
#' set.seed(1)
#' n <- 150
#' X <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
#' y <- as.numeric(X %*% c(1, -1, 0) + rnorm(n) > 0)
#' fit <- bbqr(y ~ x1 + x2 + x3, data = data.frame(y = y, X),
#'             ndraw = 500, burn = 100)
#' summary(fit)
#' @export
summary.bbqr <- function(object, level = 0.95, ...) {
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

  out <- list(coefficients = tab, level = level, quantile = object$quantile,
              penalty = object$penalty, anchor = object$anchor,
              n = object$n, ndraw = object$ndraw, burn = object$burn,
              keep = object$keep, ndraw_kept = object$ndraw_kept,
              accept = object$accept, model = object$model,
              derived = object$derived, call = object$call)
  class(out) <- "bbqr.summary"
  out
}

#' @export
print.bbqr.summary <- function(x, digits = max(3L, getOption("digits") - 3L),
                               ...) {
  cat("Bayesian binary quantile regression\n\n")
  if (!is.null(x$call)) { cat("Call:\n"); print(x$call); cat("\n") }
  cat(sprintf("  Quantile        %g\n", x$quantile))
  cat(sprintf("  Penalty         %s\n", .bbqr_penalty_label(x$penalty)))
  .bbqr_print_hierarchy(x)
  cat(sprintf("  Identification  %s\n", .bbqr_anchor_label(x$anchor)))
  cat(sprintf("  Observations    %d\n", x$n))
  cat(sprintf("  Draws           %d kept of %d (burn %d, thin %d)\n",
              x$ndraw_kept, x$ndraw, x$burn, x$keep))
  if (!is.na(x$accept))
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
  cat("Note: coefficients are identified only up to a positive scale;\n")
  cat("      see the Identification section of ?bbqr.\n")
  invisible(x)
}

#' @export
print.bbqr.list <- function(x, digits = max(3L, getOption("digits") - 3L),
                            ...) {
  cat("Bayesian binary quantile regression at",
      length(x), "quantiles\n\n")
  if (!is.null(attr(x, "call"))) { print(attr(x, "call")); cat("\n") }
  cat(sprintf("  Penalty  %s\n", .bbqr_penalty_label(attr(x, "penalty"))))
  cat("\nPosterior mean coefficients:\n")
  print(round(stats::coef(x), digits))
  invisible(x)
}

#' @export
summary.bbqr.list <- function(object, level = 0.95, ...) {
  lapply(object, summary, level = level)
}

#' Predicted probabilities and indices from a binary quantile regression
#'
#' @param object An object of class `"bbqr"`.
#' @param newdata A data frame of covariates. Required: fitted objects do not
#'   retain the design matrix, so that they stay small enough to hold thousands
#'   of fits in memory during a simulation study.
#' @param type `"response"` returns \eqn{P(y = 1 \mid x)} implied by the ALD at
#'   the fitted quantile; `"link"` returns the latent index
#'   \eqn{\beta_0 + x'\beta}; `"class"` returns the 0/1 label given by the sign
#'   of that index.
#' @param ... Ignored.
#' @return A numeric vector, one entry per row of `newdata`.
#' @details
#' Predictions use the posterior mean coefficients. Because the coefficient
#' vector is identified only up to a positive scale, `type = "response"`
#' depends on the identification anchor, whereas `type = "class"` does not.
#' @examples
#' set.seed(1)
#' n <- 150
#' X <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
#' y <- as.numeric(X %*% c(1, -1, 0) + rnorm(n) > 0)
#' d <- data.frame(y = y, X)
#' fit <- bbqr(y ~ x1 + x2 + x3, data = d, ndraw = 500, burn = 100)
#' head(predict(fit, newdata = d, type = "response"))
#' mean(predict(fit, newdata = d, type = "class") == y)
#' @export
predict.bbqr <- function(object, newdata = NULL,
                         type = c("response", "link", "class"), ...) {
  type <- match.arg(type)
  cf <- stats::coef(object)

  if (is.null(newdata)) {
    stop("'newdata' must be supplied; the fitted object does not store the ",
         "design matrix.", call. = FALSE)
  }
  tt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(tt, newdata, xlev = NULL)
  X  <- stats::model.matrix(tt, mf)
  if (!"(Intercept)" %in% colnames(X)) X <- cbind(`(Intercept)` = 1, X)
  miss <- setdiff(names(cf), colnames(X))
  if (length(miss))
    stop("'newdata' is missing required columns: ",
         paste(miss, collapse = ", "), call. = FALSE)

  eta <- as.vector(X[, names(cf), drop = FALSE] %*% cf)
  if (type == "link")  return(eta)
  if (type == "class") return(as.integer(eta > 0))

  ## P(y = 1 | x) = 1 - F_ALD(-eta) with scale sigma and skewness tau.
  tau   <- object$quantile
  sigma <- mean(object$sigma)
  ifelse(eta > 0,
         1 - tau * exp(-(1 - tau) * eta / sigma),
         (1 - tau) * exp(tau * eta / sigma))
}

#' Trace and density plots for a fitted binary quantile regression
#'
#' @param x An object of class `"bbqr"`.
#' @param which Which coefficients to plot, by name or index. Index 1 is the
#'   intercept. Defaults to all.
#' @param type `"trace"`, `"density"`, or `"both"` (the default), which draws
#'   the two side by side.
#' @param ... Passed to the underlying plotting calls.
#' @return `x`, invisibly. Called for its side effect.
#' @examples
#' set.seed(1)
#' n <- 150
#' X <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
#' y <- as.numeric(X %*% c(1, -1, 0) + rnorm(n) > 0)
#' fit <- bbqr(y ~ x1 + x2 + x3, data = data.frame(y = y, X),
#'             ndraw = 500, burn = 100)
#' op <- graphics::par(no.readonly = TRUE)
#' plot(fit, which = "x1")
#' graphics::par(op)
#' @export
plot.bbqr <- function(x, which = NULL,
                      type = c("both", "trace", "density"), ...) {
  type <- match.arg(type)
  draws <- cbind(`(Intercept)` = x$beta0, x$beta)
  if (is.null(which)) which <- colnames(draws)
  if (is.numeric(which)) which <- colnames(draws)[which]
  bad <- setdiff(which, colnames(draws))
  if (length(bad))
    stop("Unknown coefficient(s): ", paste(bad, collapse = ", "),
         call. = FALSE)

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
      graphics::abline(v = 0, lty = 2, col = "grey50")
    }
  }
  invisible(x)
}
