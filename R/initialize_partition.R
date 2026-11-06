initialize_partition <- function(X, K, minimum, maximum, random = FALSE) {
  n <- nrow(X)
  if (K == 1) return(list(seq_len(n)))
  if (random) {
    labels <- sample(rep(seq_len(K), length.out = n))
  } else {
    labels <- stats::kmeans(X, centers = K, iter.max = 100, nstart = 100)$cluster
  }
  sizes <- tabulate(labels, nbins = K)
  centers <- t(matrix(vapply(seq_len(K), function(k)
    colMeans(X[labels == k, , drop = FALSE]), numeric(ncol(X))),
    nrow = ncol(X), ncol = K))
  # Each move reduces an excess or deficit, so at most n moves are needed.
  while (any(sizes < minimum) || any(sizes > maximum)) {
    destination <- if (any(sizes < minimum)) which(sizes < minimum)[1] else
      which(sizes < maximum)[1]
    donors <- if (any(sizes > maximum)) which(sizes > maximum) else
      which(sizes > minimum)
    candidates <- which(labels %in% donors)
    distances <- rowSums(sweep(X[candidates, , drop = FALSE], 2,
                              centers[destination, ], "-")^2)
    index <- candidates[which.min(distances)]
    donor <- labels[index]
    labels[index] <- destination
    sizes[donor] <- sizes[donor] - 1L
    sizes[destination] <- sizes[destination] + 1L
  }
  lapply(seq_len(K), function(k) which(labels == k))
}
