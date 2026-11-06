library(GPcluster)
local({
  cl <- parallel::makePSOCKcluster(2)
  on.exit({parallel::stopCluster(cl); foreach::registerDoSEQ()}, add = TRUE)
  parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
  doParallel::registerDoParallel(cl)
  set.seed(15)
  X <- matrix(seq(0, 1, length.out = 30))
  Y <- sin(6 * X[, 1]) + rnorm(30, sd = .15)
  for (K in c(1, 2)) {
    set.seed(8)
    serial <- GPcluster_fit(X, Y, K = K, nugget = NULL, iter_max = 1)
    set.seed(8)
    workers <- GPcluster_fit(X, Y, K = K, nugget = NULL, iter_max = 1, parallel = TRUE)
    stopifnot(isTRUE(all.equal(serial$LOOCV, workers$LOOCV)),
              isTRUE(all.equal(predict(serial, X[1:3, , drop = FALSE]),
                               predict(workers, X[1:3, , drop = FALSE], parallel = TRUE))))
  }
})
cat('Two-worker PSOCK checks passed\n')
