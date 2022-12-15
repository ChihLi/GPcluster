#' Fit a Clustered Gaussian Process Model
#'
#' @description The function performs fitting procedure for a clustered Gaussian process model.
#'
#' @param X a design matrix with dimension \code{n} by \code{d}.
#' @param Y a response vector of size \code{n}.
#' @param K a positive integer specifying number of clusters.
#' @param n_max a positive integer specifying the maximum number of observations in each cluster. The default is 200.
#' @param nugget a positive value specifying the nugget for fitting a Gaussian process model. The default is 1e-6.
#' @param gfun a character string specifying the latent class distribution of the cluster assignment. It can be "logistic" for multinomial logistric regression, "lda" for linear discriminant analysis, and "qda" for quadratic discriminant analysis.
#' @param parallel logical. If \code{TRUE}, apply function in parallel in \code{ldply} using parallel backend provided by foreach.
#' @param iter_max a positive integer specifying the maximum numbers of iterations in the stochastic EM algorithm. The default is 10.
#' @param save.loocvpred logical.  If \code{TRUE}, the LOOCV predictions will be saved and returned.
#' @param verbose logical. If \code{TRUE}, additional diagnostics are printed.
#' @param nstart a positive integer specifying how many random initializations are performed. The default is 1.
#'
#' @details A clustered Gaussian process model segments the input data into multiple clusters, in each of which a Gaussian process is performed.
#' The unknown parameters are estimated by a stochastic EM algorithm, where \code{iter_max} can assign the number of iterations in the algorithm. The function also returns the LOOCV error in each iteration, which can assess the prediction performance. The iteration with minimum LOOCV will be selected as the final assignment. More details can be seen in Sung et al. (2019).
#'
#' @return
#' \item{X}{data \code{X}.}
#' \item{Y}{data \code{Y}.}
#' \item{Z}{cluster assignments.}
#' \item{P}{the set of indices of the observations in each cluster.}
#' \item{K}{\code{K}.}
#' \item{LOOCV}{Leave-one-out cross-validation of each iteration.}
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
#' @examples \dontrun{
#'
#' library(GPcluster)
#' ##### Gramacy and Lee (2009)
#' data.df <- read.csv(system.file("extdata", "GL2009.csv", package = "GPcluster"))
#' X <- data.df$X
#' Y <- data.df$Y
#' X.test <- matrix(sort(c(X, seq(0,20,0.1))), ncol = 1)
#'
#' set.seed(1)
#' fit.object <- GPcluster_fit(matrix(X, ncol = 1), Y, K = 2, iter_max = 100)
#' pred.out <- predict(fit.object, X.test, conf.level = 0.95)
#' yhat <- pred.out$yhat
#' plot(X, Y, ylim = c(-1.5,1.5), xlim = c(0,20), xlab = "x", ylab = "y", type = "n")
#' polygon(c(rev(X.test[,1]), X.test[,1]), c(rev(pred.out$UCL), pred.out$LCL), col = 'grey80', border = NA)
#'
#' lines(X.test[,1], yhat, col = 4, lwd = 2, lty = 2)
#' curve(sin(pi*x/5) + 0.2*cos(4*pi*x/5), from = 0, to = 10, add = TRUE, lwd = 1)
#' curve(x/10 - 1, from = 10, to = 20, add = TRUE, lwd = 1)
#' points(X, Y, cex = 1.5)
#' for(i in 1:2) points(X[fit.object$P[[i]]], Y[fit.object$P[[i]]], col = c(2,3)[i], pch = 18, cex = 1.5)
#'
#'
#' #####  Montagna and Tokdar (2016)
#' X <- seq(-2,2,length.out = 15)
#' Y <- sin(X) + 2*exp(-30*X^2)
#' X.test <- seq(-2,2,length.out = 500)
#' set.seed(1)
#'
#' K <- 3
#' fit.object <- GPcluster_fit(matrix(X, ncol = 1), Y, K = K, iter_max = 100)
#'
#' pred.out <- predict(fit.object, X.test, conf.level = 0.95)
#' yhat <- pred.out$yhat
#' plot(X, Y, ylim = c(-1.2,2), xlim = c(-2,2), type = "n")
#' polygon(c(rev(X.test), X.test), c(rev(pred.out$UCL), pred.out$LCL), col = 'grey80', border = NA)
#' curve(sin(x) + 2*exp(-30*x^2), from = -2, to = 2, add = TRUE, lwd = 1)
#' lines(X.test, yhat, col = 4, lwd = 2, lty = 2)
#' points(X, Y, cex = 1.5)
#' for(i in 1:K) points(X[fit.object$P[[i]]], Y[fit.object$P[[i]]], col = i+1, pch = 18, cex = 1.5)
#'
#' #####             Testing function: borehole function                        #####
#' #####  Thanks to Sonja Surjanovic and Derek Bingham, Simon Fraser University #####
#' borehole <- function(xx)
#' {
#'   rw <- 0.05  + xx[1] * 0.1
#'   r  <- 100   + xx[2] * 49900
#'   Tu <- 63070 + xx[3] * 52530
#'   Hu <- 990   + xx[4] * 120
#'   Tl <- 63.1  + xx[5] * 52.9
#'   Hl <- 700   + xx[6] * 120
#'   L  <- 1120  + xx[7] * 560
#'   Kw <- 9855  + xx[8] * 2190
#'
#'   frac1 <- 2 * pi * Tu * (Hu-Hl)
#'
#'   frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
#'   frac2b <- Tu / Tl
#'   frac2 <- log(r/rw) * (1+frac2a+frac2b)
#'
#'   y <- frac1 / frac2
#'   return(y)
#' }
#'
#' #####   Training data and testing data   #####
#' set.seed(1)
#' n <- 1000; n_new <- 10000; d <- 8
#' X.train <- matrix(runif(d*n), ncol = d)
#' Y.train <- apply(X.train, 1, borehole)
#' X.test <- matrix(runif(d*n_new), ncol = d)
#' Y.test <- apply(X.test, 1, borehole)
#'
#' #####   Fitting    #####
#' fit.object <- GPcluster_fit(X.train, Y.train, K = 5, iter_max = 5)
#'
#' #####   Prediction   ######
#' Y.pred <- predict(fit.object, X.test)$yhat
#' print(sqrt(mean((Y.test - Y.pred)^2)))
#'
#' #####   borehole function: parallel computing   #####
#' library(doParallel)
#' library(snowfall)
#' sfInit(parallel = TRUE, cpus = 5L) # request 5 CPUs
#' cl <- sfGetCluster()
#' registerDoParallel(cl)
#'
#'#' #####   Fitting    #####
#' fit.object <- GPcluster_fit(X.train, Y.train, K = 5, iter_max = 5, parallel = TRUE)
#'
#' #####   Prediction   ######
#' Y.pred <- predict(fit.object, X.test, parallel = TRUE)$yhat
#' print(sqrt(mean((Y.test - Y.pred)^2)))
#' sfStop()
#' }
#' @export
#'


GPcluster_fit <- function(X, Y, K = 10, n_max = 200, nugget = 1e-6, gfun = c("logistic", "lda", "qda")[1], parallel = FALSE, iter_max = 10, save.loocvpred = FALSE, verbose = FALSE, nstart = 1){

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
  P <- vector("list", K) # index of membership
  Y <- scale(Y)

  if(!save.loocvpred) ypred.loocv <- NULL

  M.ls <- P.ls <- LOOCV.ls <- vector("list", nstart)
  dev.vt <- rep(NA, nstart)
  for(ii in 1:nstart){
    if(ii == 1){
      kmeans.out <- kmeans(X, centers = K, iter.max = 100, nstart = 100)
      for(k in 1:K)  P[[k]] <- as.numeric(which(kmeans.out$cluster == k))
    }else{
      Z <- c(rep(1:K, each = floor(n/K)), rep(K, each = n%%K))
      Z <- Z[sample(1:n)]
      for(k in 1:K)  P[[k]] <- as.numeric(which(Z == k))
    }

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
    mu.pred <- apply(yhat.m * M.models$gfun.prob, 1, sum)
    dev <- sqrt(mean((mu.pred - Y)^2)) * c(sd(Y.ori))
    dev.new <- dev + 1

    iter <- counter <- 0
    P.final <- P
    M.models.final <- M.models
    dev.final <- dev
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

      mu.pred <- apply(yhat.m * M.models$gfun.prob, 1, sum)
      dev.new <- sqrt(mean((mu.pred - Y)^2)) * c(sd(Y.ori))
      loocv.error[iter+1] <- dev.new

      if(dev.new < dev.final){
        P.final <- P
        M.models.final <- M.models
        dev.final <- dev.new
        if(save.loocvpred) ypred.loocv <- mu.pred * c(sd(Y.ori)) + mean(Y.ori)
      }
      print(dev.final)
    }
    dev.vt[ii] <- dev.final
    P.ls[[ii]] <- P.final
    M.ls[[ii]] <- M.models.final
    LOOCV.ls[[ii]] <- loocv.error
  }
  print(min(dev.vt))
  M.models <- M.ls[[which.min(dev.vt)]]
  P <- P.ls[[which.min(dev.vt)]]
  LOOCV <- LOOCV.ls[[which.min(dev.vt)]]

  Z <- rep(0, n)
  for (k in 1:K) Z[P[[k]]] <- k

  out <- list(X = X, Y = Y.ori, Z = Z, P = P, K = K, LOOCV = LOOCV, ypred.loocv = ypred.loocv, nugget = nugget, gfun = gfun, n_max = n_max, GP.model = M.models$GP.model, gfun.model = M.models$gfun.model)
  class(out) <- "GPcluster"
  return(out)
}
