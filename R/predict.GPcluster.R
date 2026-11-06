#' Prediction of a Clustered Gaussian Process Model
#'
#' @description The function computes the predicted responses.
#'
#' @param object a class GPcluster object estimated by \code{GPcluster_fit}.
#' @param xnew a testing matrix with dimension \code{n_new} by \code{d} in which each row corresponds to a predictive location.
#' @param conf.level a value specifying confidence level of the confidence interval. The default is 0.95.
#' @param conf.out logical. If \code{TRUE}, predictive confidence intervals will be returned.
#' @param parallel logical. If \code{TRUE}, use a parallel backend registered with \pkg{foreach}.
#' @param ... additional arguments, currently unused.
#' @return
#' \item{yhat}{a vector displaying predicted responses at locations \code{xnew}.}
#' \item{s2}{a vector displaying predicted variances at locations \code{xnew}.}
#' \item{LCL}{a vector with length \code{n_new} displaying lower bound of predictive confidence intervals at locations \code{xnew}.}
#' \item{UCL}{a vector with length \code{n_new} displaying upper bound of predictive confidence intervals at locations \code{xnew}.}
#'
#' @seealso \code{\link{GPcluster_fit}} for fitting a clustered Gaussian process model.
#' @author Chih-Li Sung <iamdfchile@gmail.com>
#'
#' @examples
#' set.seed(1)
#' X <- matrix(seq(-2, 2, length.out = 30), ncol = 1)
#' Y <- sin(X[, 1]) + 0.2 * cos(4 * X[, 1])
#' fit <- GPcluster_fit(X, Y, K = 2, iter_max = 1)
#' pred <- predict(fit, X[c(5, 15, 25), , drop = FALSE])
#' pred$yhat
#' @export predict.GPcluster
#' @export


predict.GPcluster <- function(object, xnew, conf.level = 0.95, conf.out = TRUE, parallel = FALSE, ...){

  i <- NULL # foreach iteration variable

  if (is.data.frame(xnew)) xnew <- as.matrix(xnew)
  if (is.numeric(xnew) && is.null(dim(xnew))) xnew <- matrix(xnew, ncol = 1)
  if (!is.matrix(xnew) || !is.numeric(xnew) || !nrow(xnew) ||
      ncol(xnew) != ncol(object$X) || any(!is.finite(xnew)))
    stop("xnew must contain finite numeric rows with the same number of columns as X.", call. = FALSE)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be strictly between zero and one.", call. = FALSE)
  for (flag in list(conf.out, parallel)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag))
      stop("conf.out and parallel must be TRUE or FALSE.", call. = FALSE)
  }

  X <- object$X
  Y <- object$Y
  K <- object$K
  Z <- object$Z
  K <- object$K
  gfun <- object$gfun
  GP.model <- object$GP.model
  gfun.model <- object$gfun.model

  if(parallel){
    out.m <- foreach::foreach(k = seq(K), .combine = cbind, .packages = "mlegp",
                              .export = "fastpredict.gp")%dopar%{
      pred.out <- fastpredict.gp(GP.model[[k]], newData = xnew)
      return(c(pred.out$se.fit, pred.out$fit))
    }
    out.m <- matrix(out.m, nrow = 2 * nrow(xnew), ncol = K)
    RMSE.m <- out.m[1:nrow(xnew), , drop = FALSE]
    yhat.m <- out.m[(nrow(xnew) + 1):(2 * nrow(xnew)), , drop = FALSE]
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
  if (K == 1) {
    gfun.prob <- matrix(1, nrow(xnew), 1)
  }else if(gfun == "logistic"){
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
      quantile.pred <- t(matrix(quantile.pred, ncol = 2))
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

