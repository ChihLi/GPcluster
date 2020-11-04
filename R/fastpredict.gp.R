fastpredict.gp <- function (object, newData) 
{
  if (is.null(dim(newData))) {
    r = calcCorOneObs(object$X, object$beta, object$a, newData)
    rVinv = r %*% object$invVarMatrix
    ans = predictMu(object, newData) + object$sig2 * rVinv %*% (object$Z - object$mu)
    se = object$sig2 + object$nugget - object$sig2 * rVinv %*% t(r) * object$sig2
    if(se < 0) se = 0
    se = sqrt(se)
    return(list(fit = ans, se.fit = se))
  }
  ans = matrix(0, dim(newData)[1])
  se = matrix(0, dim(newData)[1])
  
  for (i in 1:dim(newData)[1]) {
    r = calcCorOneObs(object$X, object$beta, object$a, newData[i,])
    rVinv = r %*% object$invVarMatrix
    ans[i] = predictMu(object, newData[i,]) + object$sig2 * rVinv %*% (object$Z - object$mu)
    se[i] = object$sig2 + object$nugget - object$sig2 * rVinv %*% t(r) * object$sig2
  }
  se[se<0] = 0
  se = sqrt(se)
  
  return(list(fit = ans, se.fit = se))
}