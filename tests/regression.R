library(GPcluster)
set.seed(1)
X <- matrix(seq(-2, 2, length.out = 30), ncol = 1)
Y <- sin(X[, 1]) + 0.2 * cos(4 * X[, 1])
for (iterations in c(0, 1)) {
    fit <- GPcluster_fit(X, Y, K = 2, iter_max = iterations,
                         save.loocvpred = TRUE, nstart = 2)
    stopifnot(length(fit$ypred.loocv) == length(Y),
              isTRUE(all.equal(sqrt(mean((Y - fit$ypred.loocv)^2)),
                               min(fit$LOOCV))))
}
pred <- predict(fit, X[c(5, 15, 25), , drop = FALSE])
stopifnot(all(is.finite(unlist(pred))), all(pred$LCL <= pred$UCL))
# Exercise foreach's result simplification for a single predictive location.
foreach::registerDoSEQ()
one <- X[5, , drop = FALSE]
stopifnot(isTRUE(all.equal(predict(fit, one, parallel = TRUE),
                           predict(fit, one), tolerance = 1e-6)))
for (gate in c("lda", "qda")) {
    set.seed(1)
    fit <- GPcluster_fit(X, Y, K = 2, iter_max = 1, gfun = gate)
    stopifnot(all(is.finite(unlist(predict(fit, one)))))
}

# Candidate selection must return the actual winning fit, without refitting it.
set.seed(42)
selected <- GPcluster_fit(X, Y, K = c(3, 2, 3), iter_max = 1,
                          nstart = 2, save.loocvpred = TRUE)
set.seed(42)
fits <- lapply(2:3, function(k) GPcluster_fit(X, Y, K = k, iter_max = 1,
                                            nstart = 2, save.loocvpred = TRUE))
scores <- vapply(fits, function(f) min(f$LOOCV), numeric(1))
winner <- fits[[which.min(scores)]]
stopifnot(identical(selected$K_selection$K, c(2, 3)),
          isTRUE(all.equal(selected$K_selection$LOOCV, scores)),
          sum(selected$K_selection$selected) == 1,
          selected$K == winner$K,
          identical(selected$P, winner$P),
          identical(selected$LOOCV, winner$LOOCV),
          identical(selected$ypred.loocv, winner$ypred.loocv),
          isTRUE(all.equal(predict(selected, one), predict(winner, one))),
          is.null(winner$K_selection))
for (bad in list(numeric(0), NA, Inf, 2.5, c(2, 0), "2")) {
    stopifnot(inherits(try(GPcluster_fit(X, Y, K = bad), silent = TRUE),
                       "try-error"))
}
stopifnot(inherits(try(GPcluster_fit(X, Y, K = 2:3, n_max = 10),
                      silent = TRUE), "try-error"))

# Nugget estimation must work without replicated inputs, in one or two dimensions.
set.seed(15)
xnoise <- matrix(seq(0, 1, length.out = 30), ncol = 1)
ynoise <- sin(6 * xnoise[, 1]) + rnorm(30, sd = 0.15)
fit_one <- getFromNamespace("fit_cluster_gp", "GPcluster")
for (design in list(xnoise, cbind(xnoise, cos(4 * xnoise)))) {
    estimated <- fit_one(design, ynoise, NULL, FALSE)
    initial <- max(0.01 * var(ynoise), 1e-6)
    direct <- mlegp::mlegp(design, ynoise, nugget = initial,
                           nugget.known = 0, verbose = 0)
    stopifnot(is.finite(estimated$nugget), estimated$nugget > 0,
              !isTRUE(all.equal(estimated$nugget, initial)),
              isTRUE(all.equal(estimated$nugget, direct$nugget)))
}
for (nug in list(NULL, 1e-4)) {
    set.seed(8)
    serial <- GPcluster_fit(xnoise, ynoise, K = 2, nugget = nug,
                            iter_max = 1, save.loocvpred = TRUE)
    stopifnot(all(is.finite(serial$LOOCV)),
              all(is.finite(unlist(predict(serial, xnoise[1:3, , drop = FALSE])))))
    # Sequential foreach exercises the same worker body without opening sockets.
    set.seed(8)
    parfit <- suppressWarnings(GPcluster_fit(xnoise, ynoise, K = 2, nugget = nug,
                                            iter_max = 1, parallel = TRUE,
                                            save.loocvpred = TRUE))
    stopifnot(isTRUE(all.equal(serial$LOOCV, parfit$LOOCV)),
              isTRUE(all.equal(serial$GP.model, parfit$GP.model)))
}
for (bad in list(NA, Inf, -1, c(0.1, 0.2), "estimate")) {
    stopifnot(inherits(try(GPcluster_fit(X, Y, K = 2, nugget = bad),
                          silent = TRUE), "try-error"))
}

# Quantiles far outside the old fixed bracket and degenerate components.
quantile_mix <- getFromNamespace("F_inv", "GPcluster")
cdf_mix <- getFromNamespace("F.mixed", "GPcluster")
for (p in c(0.001, 0.025, 0.5, 0.975, 0.999)) {
    stopifnot(isTRUE(all.equal(quantile_mix(p, 1, 100, 20),
                               qnorm(p, 100, 20))))
    q <- quantile_mix(p, c(.3, .7), c(-100, 200), c(2, 30))
    stopifnot(abs(cdf_mix(q, c(.3, .7), c(-100, 200), c(2, 30)) - p) < 1e-10)
}
stopifnot(quantile_mix(.5, c(.8, .2), c(50, 100), c(0, 0)) == 50,
          quantile_mix(.9, c(.8, .2), c(50, 100), c(0, 0)) == 100,
          quantile_mix(.5, c(1, 0), c(42, -1e6), c(0, 1)) == 42)
weights <- getFromNamespace("assignment_weights", "GPcluster")
stopifnot(identical(weights(c(Inf, 2), c(0, 1)), c(0, 1)),
          identical(weights(c(NaN, 2), c(.5, .5)), c(0, 1)),
          identical(weights(c(Inf, 2), c(.5, .5)), c(.5, 0)),
          identical(weights(c(0, 0), c(.5, .5)), c(0, 0)))
expect_error <- function(expr, pattern) {
    error <- tryCatch({force(expr); NULL}, error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
expect_error(GPcluster_fit(X, rep(1, nrow(X)), K = 2), "nonconstant")
expect_error(GPcluster_fit(X, Y[-1], K = 2), "one value")
expect_error(GPcluster_fit(X, Y, K = 2, iter_max = -1), "iter_max")
expect_error(GPcluster_fit(X, Y, K = 2, nstart = 0), "nstart")
expect_error(GPcluster_fit(X, Y, K = 2, parallel = NA), "parallel")
expect_error(GPcluster_fit(X, Y, K = 11), "infeasible")
expect_error(GPcluster_fit(X, Y, K = 2, n_max = 1), "infeasible")

expect_error(predict(selected, matrix(1, 1, 2)), "same number")
expect_error(predict(selected, NA_real_), "finite")
expect_error(predict(selected, one, conf.level = 1), "conf.level")

# Inverse updates agree with full inversion, including ill-conditioned covariances.
delete_inverse <- getFromNamespace("inverse_delete", "GPcluster")
add_inverse <- getFromNamespace("inverse_add", "GPcluster")
for (nugget in c(1e-6, .1)) {
    xx <- c(0, 1e-4, .2, .4, .7, 1)
    V <- exp(-outer(xx, xx, "-")^2 / .1) + diag(nugget, length(xx))
    for (index in seq_along(xx)) {
        keep <- seq_along(xx) != index
        Q <- delete_inverse(solve(V), index)
        stopifnot(isTRUE(all.equal(Q, solve(V[keep, keep]), tolerance = 1e-5)))
        restored <- add_inverse(Q, V[keep, index], V[index, index], index, covariance_matrix = V)
        stopifnot(isTRUE(all.equal(restored, solve(V), tolerance = 1e-5)))
    }
}
# Both ordinary densities underflow, but their relative log weights remain usable.
logweights <- weights(dnorm(c(100, 101), log = TRUE), c(.5, .5), log = TRUE)
stopifnot(logweights[1] == 1, logweights[2] > 0, logweights[2] < 1)
initialize <- getFromNamespace("initialize_partition", "GPcluster")
set.seed(1)
partition <- initialize(matrix(c(seq(0, .1, length.out = 29), 100)), 2, 3, 20)
stopifnot(all(lengths(partition) >= 3), all(lengths(partition) <= 20),
          identical(sort(unlist(partition)), seq_len(30)))
# A single GP is a valid candidate with all supported gate choices.
for (gate in c("logistic", "lda", "qda")) {
    single <- GPcluster_fit(X, Y, K = 1, gfun = gate, iter_max = 0,
                             save.loocvpred = TRUE)
    stopifnot(single$K == 1, is.null(single$gfun.model),
              all(is.finite(unlist(predict(single, one)))),
              isTRUE(all.equal(predict(single, one),
                               predict(single, one, parallel = TRUE))))
}
set.seed(10)
baseline_selection <- GPcluster_fit(X, Y, K = 1:3, iter_max = 1)
stopifnot(baseline_selection$K == baseline_selection$K_selection$K[
    which.min(baseline_selection$K_selection$LOOCV)])

# Verify inverse updates after actual successive E-step cluster transfers.
xe <- matrix(seq(0, 1, length.out = 30))
set.seed(4)
ye <- sin(8 * xe[, 1]) + rnorm(30, sd = .2)
partition <- list(seq(1, 30, 2), seq(2, 30, 2))
mstep <- getFromNamespace("GPcluster_Mstep", "GPcluster")
estep <- getFromNamespace("GPcluster_Estep", "GPcluster")
models <- mstep(xe, ye, partition, 2, .1, "logistic", FALSE, FALSE)
set.seed(5)
moved <- estep(xe, ye, partition, 2, models, "logistic", 3, 25, FALSE)
stopifnot(!identical(moved$P, partition))
for (k in 1:2) {
    gp <- moved$M.models$GP.model[[k]]
    V <- mlegp::calcVarMatrix(gp$X, gp$beta, gp$a, gp$nugget, gp$sig2, 0, nrow(gp$X))
    stopifnot(isTRUE(all.equal(gp$invVarMatrix, solve(V), tolerance = 1e-6)),
              identical(as.numeric(gp$Z), ye[moved$P[[k]]]))
}
