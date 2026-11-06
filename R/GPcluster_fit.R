#' Fit a Clustered Gaussian Process Model
#'
#' @description The function performs fitting procedure for a clustered Gaussian process model.
#'
#' @param X a finite numeric design matrix with dimension \code{n} by \code{d},
#'   a numeric data frame, or a vector for one-dimensional inputs.
#' @param Y a finite, nonconstant numeric response vector of size \code{n}.
#' @param K a positive integer specifying the number of clusters, or a
#'   vector of candidate integers. For multiple candidates, the model with the
#'   smallest LOOCV root mean squared error is returned.
#' @param n_max a positive integer specifying the maximum number of observations in each cluster. The default is 200.
#' @param nugget a nonnegative scalar fixing the nugget variance on the globally
#'   standardized response scale. The default is 1e-6. Set to \code{NULL}
#'   to estimate a separate nugget in each cluster jointly with its correlation
#'   parameters using \pkg{mlegp}. The starting value is one percent of the
#'   within-cluster response variance, with a minimum starting value of 1e-6;
#'   this minimum is not a lower bound on the estimate.
#' @param gfun a character string specifying the latent class distribution of the cluster assignment. It can be "logistic" for multinomial logistic regression, "lda" for linear discriminant analysis, and "qda" for quadratic discriminant analysis.
#' @param parallel logical. If \code{TRUE}, use a parallel backend registered with \pkg{foreach}.
#' @param iter_max a nonnegative integer specifying the maximum numbers of iterations in the stochastic EM algorithm. The default is 10.
#' @param save.loocvpred logical.  If \code{TRUE}, the LOOCV predictions will be saved and returned.
#' @param verbose logical. If \code{TRUE}, additional diagnostics are printed.
#' @param nstart a positive integer specifying how many random initializations are performed. The default is 1.
#'
#' @details A clustered Gaussian process model segments the input data into multiple clusters, in each of which a Gaussian process is performed.
#' The unknown parameters are estimated by a stochastic EM algorithm, where \code{iter_max} can assign the number of iterations in the algorithm. The function also returns the LOOCV error in each iteration, which can assess the prediction performance. The iteration with minimum LOOCV will be selected as the final assignment. More details can be seen in Sung et al. (2023), \doi{10.5705/ss.202020.0456}.
#'
#' Candidate values are sorted and deduplicated before fitting. Exact ties select
#' the smaller number of clusters. Each candidate uses the supplied fitting
#' options, including \code{nstart}. Selection uses the existing LOOCV diagnostic
#' (with fitted partitions and model parameters), not a full refit of the
#' clustering procedure for every held-out observation. Set a seed before fitting
#' for reproducibility. Candidate fitting errors stop selection and identify the
#' failing value. In a candidate search, every value must satisfy
#' \code{K * n_max >= length(Y)}; candidates are not silently increased.
#'
#' The existing minimum cluster size is \code{3 * ncol(X)}. Infeasible size
#' constraints are rejected before fitting. Initial partitions are repaired by
#' transferring observations from oversized or sufficiently large clusters,
#' choosing transfers by distance to the destination center. A single cluster
#' bypasses the classification model and provides a standard GP baseline.
#'
#' @return
#' \item{X}{data \code{X}.}
#' \item{Y}{data \code{Y}.}
#' \item{Z}{cluster assignments.}
#' \item{P}{the set of indices of the observations in each cluster.}
#' \item{K}{the number of clusters in the returned model.}
#' \item{LOOCV}{LOOCV root mean squared errors for the winning initialization
#' of the returned model, including the initial fit.}
#' \item{K_selection}{For a vector of candidates, a data frame containing each
#' candidate \code{K}, its minimum \code{LOOCV}, and a logical \code{selected}
#' indicator. Absent for a scalar \code{K}.}
#' \item{ypred.loocv}{If \code{save.loocvpred} is \code{TRUE}, then it returns the LOOCV predictions, otherwise it returns \code{NULL}.}
#' \item{nugget}{\code{nugget}.}
#' \item{gfun}{\code{gfun}.}
#' \item{n_max}{\code{n_max}.}
#' \item{GP.model}{a list of fitted Gaussian process models.}
#' \item{gfun.model}{fitted latent class distribution.}
#'
#' @seealso \code{\link{predict.GPcluster}} for prediction of the clustered Gaussian process model.
#' @author Chih-Li Sung <iamdfchile@gmail.com>
#'
#' @import nnet
#' @import mlegp
#' @import foreach
#' @importFrom stats dnorm kmeans pnorm predict rmultinom sd uniroot var
#' @examples
#' set.seed(1)
#' X <- matrix(seq(-2, 2, length.out = 30), ncol = 1)
#' Y <- sin(X[, 1]) + 0.2 * cos(4 * X[, 1])
#' fit <- GPcluster_fit(X, Y, K = 2, iter_max = 1)
#' pred <- predict(fit, X[c(5, 15, 25), , drop = FALSE])
#' pred$yhat
#' selected <- GPcluster_fit(X, Y, K = 2:3, iter_max = 1)
#' selected$K
#' selected$K_selection
#' @export
#'


GPcluster_fit <- function(X, Y, K = 10, n_max = 200, nugget = 1e-6, gfun = c("logistic", "lda", "qda")[1], parallel = FALSE, iter_max = 10, save.loocvpred = FALSE, verbose = FALSE, nstart = 1){

  if (is.data.frame(X)) X <- as.matrix(X)
  if (is.numeric(X) && is.null(dim(X))) X <- matrix(X, ncol = 1)
  if (!is.matrix(X) || !is.numeric(X) || !nrow(X) || !ncol(X) ||
      any(!is.finite(X))) stop("X must be a nonempty finite numeric matrix or vector.", call. = FALSE)
  if (!is.numeric(Y) || !is.null(dim(Y)) || length(Y) != nrow(X) ||
      any(!is.finite(Y)) || length(Y) < 2L || !is.finite(sd(Y)) || sd(Y) <= 0)
    stop("Y must be a finite nonconstant numeric vector with one value per row of X.", call. = FALSE)
  for (name in c("n_max", "iter_max", "nstart")) {
    value <- get(name)
    lower <- if (name == "iter_max") 0 else 1
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < lower || value != floor(value))
      stop(name, " must be an integer >= ", lower, ".", call. = FALSE)
  }
  for (name in c("parallel", "verbose", "save.loocvpred")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  gfun <- match.arg(gfun, c("logistic", "lda", "qda"))
  if (!is.null(nugget) &&
      (!is.numeric(nugget) || length(nugget) != 1L ||
       !is.finite(nugget) || nugget < 0)) {
    stop("nugget must be NULL or a finite nonnegative scalar.", call. = FALSE)
  }
  if (!is.numeric(K) || !length(K) || anyNA(K) ||
      any(!is.finite(K)) || any(K < 1 | K != floor(K))) {
    stop("K must contain positive integers.", call. = FALSE)
  }
  if (length(K) > 1L) {
    candidates <- sort(unique(K))
    if (!is.numeric(n_max) || length(n_max) != 1L ||
        !is.finite(n_max) || n_max < 1 || n_max != floor(n_max)) {
      stop("n_max must be a positive integer.", call. = FALSE)
    }
    if (any(candidates * n_max < length(Y))) {
      stop("Every candidate K must satisfy K * n_max >= length(Y).",
           call. = FALSE)
    }
    scores <- numeric(length(candidates))
    best <- NULL
    best.score <- Inf
    for (j in seq_along(candidates)) {
      candidate <- candidates[j]
      fit <- tryCatch(
        GPcluster_fit(X, Y, K = candidate, n_max = n_max, nugget = nugget,
                      gfun = gfun, parallel = parallel, iter_max = iter_max,
                      save.loocvpred = save.loocvpred, verbose = verbose,
                      nstart = nstart),
        error = function(e) stop("Fitting K = ", candidate, " failed: ",
                                 conditionMessage(e), call. = FALSE)
      )
      score <- min(fit$LOOCV)
      if (length(score) != 1L || !is.finite(score)) {
        stop("Non-finite LOOCV score for K = ", candidate, ".", call. = FALSE)
      }
      scores[j] <- score
      if (score < best.score) {
        best <- fit
        best.score <- score
      }
    }
    best$K_selection <- data.frame(K = candidates, LOOCV = scores,
                                  selected = candidates == best$K)
    return(best)
  }

  ### setting
  if(parallel){
    if (!requireNamespace("foreach", quietly = TRUE)) {
      # EXCLUDE COVERAGE START
      stop("foreach package required for parallel GPcluster_fit operation",
           call. = FALSE)
      # EXCLUDE COVERAGE END
    }
    if (foreach::getDoParWorkers() == 1) {
      # EXCLUDE COVERAGE START
      warning("No parallel backend registered", call. = TRUE)
      # EXCLUDE COVERAGE END
    }
  }

  ### 0. initialization
  # if(!is.null(subsample_num)){ # under development
  #   sample.index <- sort(sample(1:length(Y), subsample_num))
  #   X <- X[sample.index,]
  #   Y <- Y[sample.index]
  # }
  n <- length(Y)
  if(n_max * K < n){
    warning("Because n_max * K < n, K = ceiling(n/n_max).")
    K <- ceiling(n/n_max)
  }
  if (is.null(ncol(X))) X <- matrix(X, ncol = 1)

  Y.ori <- Y
  n_min <- ncol(X) + 2 + (ncol(X)-1)*2
  if (n_max < n_min || K > floor(n / n_min))
    stop("Cluster sizes are infeasible: require n_max >= ", n_min,
         " and at least ", n_min, " observations per cluster.", call. = FALSE)
  if (nrow(unique(X)) < K)
    stop("K exceeds the number of distinct input rows.", call. = FALSE)
  P <- vector("list", K) # index of membership
  Y <- scale(Y)

  if(!save.loocvpred) ypred.loocv <- NULL

  M.ls <- P.ls <- LOOCV.ls <- vector("list", nstart)
  ypred.ls <- vector("list", nstart)
  dev.vt <- rep(NA, nstart)
  for(ii in 1:nstart){
    P <- initialize_partition(X, K, n_min, n_max, random = ii > 1)
    M.models <- GPcluster_Mstep(X, Y, P, K, nugget, gfun, parallel, verbose)
    yhat.m <- matrix(0, nrow = nrow(X), ncol = K)

    if(parallel){
      yhat.m <- foreach::foreach(k = seq(K), .combine = cbind, .packages = "mlegp")%dopar%{
        yhat <- predict(M.models$GP.model[[k]], newData = X, se.fit = FALSE)
        yhat[P[[k]]] <- M.models$GP.model[[k]]$cv[,1]
        return(yhat)
      }
    }else{
      for(k in 1:K){
        yhat.m[,k] <- predict(M.models$GP.model[[k]], newData = X, se.fit = FALSE)
        yhat.m[P[[k]],k] <- M.models$GP.model[[k]]$cv[,1]
      }
    }
    yhat.m <- matrix(yhat.m, nrow = nrow(X), ncol = K)
    mu.pred <- apply(yhat.m * M.models$gfun.prob, 1, sum)
    dev <- sqrt(mean((mu.pred - Y)^2)) * c(sd(Y.ori))
    dev.new <- dev + 1

    iter <- counter <- 0
    P.final <- P
    M.models.final <- M.models
    dev.final <- dev
    if(save.loocvpred) ypred.loocv <- mu.pred * c(sd(Y.ori)) + mean(Y.ori)
    loocv.error <- rep(0,iter_max+1) ## save loocv
    loocv.error[1] <- dev
    while(iter < iter_max){
      #if(dev > dev.new & dev - dev.new < conv_threshold) counter <- counter + 1 else counter <- 0 # under development
      dev <- dev.new
      iter <- iter + 1

      ### 1. E-step
      P <- GPcluster_Estep(X, Y, P, K, M.models, gfun, n_min, n_max, parallel)$P

      ### 2. M-step
      M.models <- GPcluster_Mstep(X, Y, P, K, nugget, gfun, parallel, verbose)

      if(parallel){
        yhat.m <- foreach::foreach(k = seq(K), .combine = cbind, .packages = "mlegp")%dopar%{
          yhat <- predict(M.models$GP.model[[k]], newData = X, se.fit = FALSE)
          yhat[P[[k]]] <- M.models$GP.model[[k]]$cv[,1]
          return(yhat)
        }
      }else{
        for(k in 1:K){
          yhat.m[,k] <- predict(M.models$GP.model[[k]], newData = X, se.fit = FALSE)
          yhat.m[P[[k]],k] <- M.models$GP.model[[k]]$cv[,1]
        }
      }

      yhat.m <- matrix(yhat.m, nrow = nrow(X), ncol = K)
    mu.pred <- apply(yhat.m * M.models$gfun.prob, 1, sum)
      dev.new <- sqrt(mean((mu.pred - Y)^2)) * c(sd(Y.ori))
      loocv.error[iter+1] <- dev.new

      if(dev.new < dev.final){
        P.final <- P
        M.models.final <- M.models
        dev.final <- dev.new
        if(save.loocvpred) ypred.loocv <- mu.pred * c(sd(Y.ori)) + mean(Y.ori)
      }
      if(verbose) print(dev.final)
    }
    if(save.loocvpred) ypred.ls[[ii]] <- ypred.loocv
    dev.vt[ii] <- dev.final
    P.ls[[ii]] <- P.final
    M.ls[[ii]] <- M.models.final
    LOOCV.ls[[ii]] <- loocv.error
  }
  if(verbose) print(min(dev.vt))
  if(save.loocvpred) ypred.loocv <- ypred.ls[[which.min(dev.vt)]]
  M.models <- M.ls[[which.min(dev.vt)]]
  P <- P.ls[[which.min(dev.vt)]]
  LOOCV <- LOOCV.ls[[which.min(dev.vt)]]

  Z <- rep(0, n)
  for (k in 1:K) Z[P[[k]]] <- k

  out <- list(X = X, Y = Y.ori, Z = Z, P = P, K = K, LOOCV = LOOCV, ypred.loocv = ypred.loocv, nugget = nugget, gfun = gfun, n_max = n_max, GP.model = M.models$GP.model, gfun.model = M.models$gfun.model)
  class(out) <- "GPcluster"
  return(out)
}
