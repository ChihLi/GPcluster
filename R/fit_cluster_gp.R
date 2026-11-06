# Fit one cluster with a consistent covariance model in both execution paths.
fit_cluster_gp <- function(X, Y, nugget, verbose) {
  estimate <- is.null(nugget)
  # A positive start enables mlegp estimation even without replicated inputs.
  # Y is on the globally standardized response scale used by GPcluster_fit.
  start <- if (estimate) max(0.01 * stats::var(Y), 1e-6) else nugget
  fit <- mlegp::mlegp(X, Y, nugget = start,
                      nugget.known = if (estimate) 0 else 1,
                      verbose = as.integer(verbose))
  if (!inherits(fit, "gp")) {
    stop("mlegp failed to fit a cluster Gaussian process.", call. = FALSE)
  }
  fit
}
