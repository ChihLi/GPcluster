F.mixed <- function(x, w, u, s) {
  cdf <- numeric(length(w))
  positive <- s > 0
  cdf[positive] <- pnorm(x, mean = u[positive], sd = s[positive])
  cdf[!positive] <- as.numeric(x >= u[!positive])
  sum(w * cdf)
}

# Component quantiles bracket the mixture quantile. Bisection also handles
# jumps in the CDF from zero-variance components (point masses).
F_inv <- function(p, w, u, s, br = NULL) {
  if (length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1 ||
      !length(w) || length(w) != length(u) || length(w) != length(s) ||
      any(!is.finite(c(w, u, s))) || any(w < 0) || any(s < 0) ||
      !is.finite(sum(w)) || sum(w) <= 0) {
    stop("Invalid mixture quantile inputs.", call. = FALSE)
  }
  keep <- w > 0
  w <- w[keep] / sum(w)
  u <- u[keep]
  s <- s[keep]
  q <- u + s * stats::qnorm(p)
  lo <- min(q)
  hi <- max(q)
  if (F.mixed(lo, w, u, s) >= p) return(lo)
  for (iteration in seq_len(200)) {
    mid <- lo / 2 + hi / 2
    if (mid == lo || mid == hi) break
    if (F.mixed(mid, w, u, s) >= p) hi <- mid else lo <- mid
  }
  hi
}
