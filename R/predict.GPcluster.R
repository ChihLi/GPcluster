#' Prediction of a Clustered Gaussian Process Model
#'
#' @description The function computes the predicted responses.
#'
#' @param object a class GPcluster object estimated by \code{GPcluster_fit}.
#' @param xnew a testing matrix with dimension \code{n_new} by \code{d} in which each row corresponds to a predictive location.
#' @param conf.level a value specifying confidence level of the confidence interval. The default is 0.95.
#' @param conf.out logical. If \code{TRUE}, predictive confidence intervals will be returned.
#' @param parallel logical. If \code{TRUE}, apply function in parallel in \code{ldply} using parallel backend provided by foreach.
#' @return
#' \item{yhat}{a vector displaying predicted responses at locations \code{xnew}.}
#' \item{s2}{a vector displaying predicted variances at locations \code{xnew}.}
#' \item{LCL}{a vector with length \code{n_new} displaying lower bound of predictive confidence intervals at locations \code{xnew}.}
#' \item{UCL}{a vector with length \code{n_new} displaying upper bound of predictive confidence intervals at locations \code{xnew}.}
#'
#' @seealso \code{\link{GPcluster_fit}} for fitting a clustered Gaussian process model.
#' @author Chih-Li Sung <iamdfchile@gmail.com>
#'
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
#' @export predict.GPcluster
#' @export


predict.GPcluster <- function(object, xnew, conf.level = 0.95, conf.out = TRUE, parallel = FALSE){

  if (is.null(ncol(xnew))) xnew <- matrix(xnew, ncol = 1)

  X <- object$X
  Y <- object$Y
  K <- object$K
  Z <- object$Z
  K <- object$K
  gfun <- object$gfun
  GP.model <- object$GP.model
  gfun.model <- object$gfun.model

  if(parallel){
    out.m <- foreach::foreach(k = seq(K), .combine = cbind, .packages = "mlegp")%dopar%{
      pred.out <- predict(GP.model[[k]], newData = xnew, se.fit = TRUE)
      return(c(pred.out$se.fit, pred.out$fit))
    }
    RMSE.m <- out.m[1:nrow(xnew),]
    yhat.m <- out.m[(nrow(xnew) + 1):(2 * nrow(xnew)),]
  }else{
    RMSE.m <- yhat.m <- matrix(0, nrow = nrow(xnew), ncol = K)
    for(k in 1:K){
      pred.out <- fastpredict.gp(GP.model[[k]], newData = xnew)
      RMSE.m[,k] <- pred.out$se.fit
      yhat.m[,k] <- pred.out$fit
    }
  }


  xnew <- data.frame(xnew)
  colnames(xnew) <- paste0("X", 1:ncol(xnew))
  if(gfun == "logistic"){
    # colnames(xnew) <- gfun.model$coefnames[-1]
    gfun.prob <- predict(gfun.model, newdata = xnew, type = "prob")
    if(K == 2){
      gfun.prob <- matrix(gfun.prob, ncol = 1)
      gfun.prob <- cbind(1 - apply(gfun.prob, 1, sum), gfun.prob)
    }
  }else if(gfun == "lda"){
    gfun.prob <- predict(gfun.model, newdata = xnew)$posterior
  }else if(gfun == "qda"){
    gfun.prob <- predict(gfun.model, newdata = xnew)$posterior
    # }else if(gfun == "bagging"){    # under development
    #   # colnames(xnew) <- attr(gfun.model$terms, "term.labels")
    #   gfun.prob <- predict(gfun.model, newdata = xnew, type = "prob")
    # }else if(gfun == "forest"){
    #   # colnames(xnew) <- attr(gfun.model$terms, "term.labels")
    #   gfun.prob <- predict(gfun.model, newdata = xnew, type = "prob")
    # }else if(gfun == "gp"){
    #   gfun.prob <- predict(gfun.model, newdata = xnew, type = "probabilities")
    # }else if(gfun == "adaboost"){
    #   gfun.prob <- predict(gfun.model, newdata = xnew)$prob
    # }else if(gfun == "constant"){
    #   gfun.prob <- matrix(1/K, ncol = K, nrow = nrow(xnew))
    # }else if(gfun == "forest.SRC"){
    #   gfun.prob <- predict(gfun.model, newdata = xnew)$predicted
  }

  if(nrow(xnew) == 1) gfun.prob <- matrix(gfun.prob, nrow = 1)


  if(parallel){
    ytimesp <- yhat.m * gfun.prob
    mu.pred <- foreach(i = seq(nrow(xnew)), .combine = c)%dopar%{
      return(sum(ytimesp[i,]))
    }
  }else{
    mu.pred <- apply(yhat.m * gfun.prob, 1, sum)
  }

  var.pred <- apply(RMSE.m^2 * gfun.prob, 1, sum) + apply(yhat.m^2 * gfun.prob, 1, sum)
  var.pred <- var.pred - mu.pred^2
  mu.pred <- mu.pred * c(sd(Y)) + mean(Y)
  var.pred <- var.pred * c(var(Y))

  if(conf.out){
    if(parallel){
      quantile.pred <- foreach(i = seq(nrow(xnew)), .combine = rbind, .export = c("F.mixed", "F_inv"))%dopar%{
        q1 <- F_inv((1-conf.level)/2, gfun.prob[i,], yhat.m[i,], RMSE.m[i,])
        q2 <- F_inv(1-(1-conf.level)/2, gfun.prob[i,], yhat.m[i,], RMSE.m[i,])
        return(c(q1,q2))}
      rownames(quantile.pred) <- NULL
      quantile.pred <- t(quantile.pred)
    }else{
      quantile.pred <- apply(cbind(gfun.prob, yhat.m, RMSE.m), 1, function(x.vt) {
        q1 <- F_inv((1-conf.level)/2, x.vt[1:K], x.vt[(K+1):(2*K)], x.vt[(2*K+1):(3*K)])
        q2 <- F_inv(1-(1-conf.level)/2, x.vt[1:K], x.vt[(K+1):(2*K)], x.vt[(2*K+1):(3*K)])
        return(c(q1,q2))})
    }
    quantile.pred <- quantile.pred * c(sd(Y)) + mean(Y)
    return(list(yhat = mu.pred, s2=var.pred, LCL = quantile.pred[1,], UCL = quantile.pred[2,]))
  }else{
    return(list(yhat= mu.pred, s2=var.pred))
  }


}

