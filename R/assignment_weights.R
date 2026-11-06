# Recompute weights for every observation; zero-weight gates cannot receive mass.
assignment_weights <- function(density, gate, log = FALSE) {
  if (any(!is.finite(gate)) || any(gate < 0))
    stop("Invalid cluster probabilities.", call. = FALSE)
  if (!log) density <- log(density)
  weights <- numeric(length(gate))
  infinite <- is.infinite(density) & density > 0 & gate > 0
  if (any(infinite)) {
    weights[infinite] <- gate[infinite]
  } else {
    valid <- is.finite(density) & gate > 0
    if (any(valid)) {
      logweights <- density[valid] + log(gate[valid])
      weights[valid] <- exp(logweights - max(logweights))
    }
  }
  weights
}
