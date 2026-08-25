## Small, fast simulated data set shared by the tests.
sim_binary <- function(n = 200, beta = c(1.5, -1, 0), seed = 1) {
  set.seed(seed)
  p <- length(beta)
  X <- matrix(stats::rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))
  y <- as.numeric(as.vector(X %*% beta) + stats::rnorm(n) > 0)
  list(data = data.frame(y = y, X),
       formula = stats::as.formula(
         paste("y ~", paste(paste0("x", seq_len(p)), collapse = " + "))),
       beta = beta)
}

## The published (V3) hierarchy, for tests that exercise devices the
## derivation-faithful default refuses: the clamps, the argument floors, the
## historical z/s routine and the kappa tilt. Under model = "v5" (the package
## default) every one of those is an error, so a test that turns one on must
## say which hierarchy it is testing.
v3_prior <- function(penalty = "alasso", ...) prior(penalty, model = "v3", ...)

## Cosine similarity between an estimated slope vector and the truth. The
## coefficient vector is identified only up to positive scale, so direction is
## the only thing worth comparing.
cos_sim <- function(est, truth) {
  sum(est * truth) / sqrt(sum(est^2) * sum(truth^2))
}
