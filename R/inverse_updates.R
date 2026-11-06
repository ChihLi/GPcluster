# Schur-complement updates on the covariance inverse, in observation order.
inverse_delete <- function(Q, index) {
  keep <- seq_len(nrow(Q)) != index
  pivot <- Q[index, index]
  if (!is.finite(pivot) || pivot <= 0)
    stop("Invalid covariance inverse pivot during deletion.", call. = FALSE)
  out <- Q[keep, keep, drop = FALSE] -
    tcrossprod(Q[keep, index]) / pivot
  (out + t(out)) / 2
}

inverse_add <- function(Q, covariance, variance, index, covariance_matrix = NULL) {
  v <- as.vector(Q %*% covariance)
  schur <- variance - sum(covariance * v)
  # Near-singular additions amplify roundoff in the earlier inverse update.
  # The optional full matrix is evaluated lazily only for this fallback.
  if ((!is.finite(schur) || schur <= 1e-5 * variance) &&
      !is.null(covariance_matrix))
    return(chol2inv(chol(covariance_matrix)))
  if (!is.finite(schur) || schur <= 0)
    stop("Nonpositive conditional variance during covariance update; increase the nugget.", call. = FALSE)
  out <- matrix(0, nrow(Q) + 1L, nrow(Q) + 1L)
  keep <- seq_len(nrow(out)) != index
  out[keep, keep] <- Q + tcrossprod(v) / schur
  out[index, keep] <- out[keep, index] <- -v / schur
  out[index, index] <- 1 / schur
  (out + t(out)) / 2
}
