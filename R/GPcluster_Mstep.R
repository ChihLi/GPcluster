GPcluster_Mstep <- function(X, Y, P, K, nugget, gfun, parallel, verbose){
  
  Z <- rep(0, length(Y))
  for (k in 1:K) Z[P[[k]]] <- k

  if(parallel){
    GP.model <- foreach::foreach(k = seq(K), .packages = c("laGP", "mlegp"), .combine = c, .verbose = verbose)%dopar%{
      if(is.null(nugget)){ #mle nugget
        d <- darg(NULL, X[P[[k]], , drop = FALSE])
        g <- garg(list(mle=TRUE), Y[P[[k]]])
        gpi <- newGP(X[P[[k]], , drop = FALSE], Y[P[[k]]], d=rep(d$start, ncol(X)), g=g$start, dK=TRUE)
        mle.obj <- mleGP(gpi, "g", g$min, g$max)
        nugget.mle <- mle.obj$g
        deleteGP(gpi)
        list(mlegp(X[P[[k]], ], Y[P[[k]]], nugget = nugget.mle, nugget.known = 1, verbose = verbose))
      }else{
        list(mlegp(X[P[[k]], ], Y[P[[k]]], nugget = nugget, nugget.known = 1, verbose = verbose))
      }
    }
  }else{
    GP.model <- vector("list", K) 
    for(k in 1:K) {
      if(is.null(nugget)){ #mle nugget
        d <- darg(NULL, X[P[[k]], , drop = FALSE])
        g <- garg(list(mle=TRUE), Y[P[[k]]])
        gpi <- newGP(X[P[[k]], , drop = FALSE], Y[P[[k]]], d=rep(d$start, ncol(X)), g=g$start, dK=TRUE)
        mle.obj <- mleGP(gpi, "g", g$min, g$max)
        nugget.mle <- mle.obj$g
        deleteGP(gpi)
        GP.model[[k]] <- mlegp(X[P[[k]], ], Y[P[[k]]], nugget = nugget.mle, nugget.known = 1, verbose = verbose)
      }else{
        GP.model[[k]] <- mlegp(X[P[[k]], ], Y[P[[k]]], nugget = nugget, nugget.known = 1, verbose = verbose)
        if(length(P[[k]]) <= 3 & ncol(X) == 1) {
          if(GP.model[[k]]$beta/var(GP.model[[k]]$X) > 1000){
            GP.model[[k]]$beta <- 100  * var(GP.model[[k]]$X)
            GP.model[[k]]$invVarMatrix <- solve(calcVarMatrix(GP.model[[k]]$X, GP.model[[k]]$beta, GP.model[[k]]$a,
                                                              GP.model[[k]]$nugget, GP.model[[k]]$sig2, 
                                                              0, GP.model[[k]]$numObs) + diag(nugget, nrow(GP.model[[k]]$X)))
            
          }
        }
      }
    } 
  }
  
  gfun.model <- gfun.fit(X, Z, K, gfun)
  
  if(gfun == "logistic"){
    gfun.prob <- gfun.model$fitted.values
    if(K == 2) gfun.prob <- cbind(1 - apply(gfun.prob, 1, sum), gfun.prob)
  }else if(gfun == "lda"){
    gfun.prob <- predict(gfun.model)$posterior
  }else if(gfun == "qda"){
    gfun.prob <- predict(gfun.model)$posterior
  }else if(gfun == "bagging"){
    gfun.prob <- gfun.model$votes
  }else if(gfun == "forest"){
    gfun.prob <- gfun.model$votes
  }else if(gfun == "gp"){
    X.tmp <- X
    colnames(X.tmp) <- paste0("X", 1:ncol(X))
    gfun.prob <- predict(gfun.model, X.tmp, type = "probabilities")
    rm("X.tmp")
  }else if(gfun == "adaboost"){
    X.tmp <- X
    colnames(X.tmp) <- paste0("X", 1:ncol(X))
    gfun.prob <- predict(gfun.model, X.tmp)$prob
    rm("X.tmp")
  }else if(gfun == "constant"){
    gfun.prob <- matrix(1/K, ncol = K, nrow = nrow(X))
  }else if(gfun == "forest.SRC"){
    ### random forest ###
    gfun.prob <- gfun.model$predicted
  }
  
  prob.vt <- rep(0, length(Y))
  for(k in 1:K) prob.vt[P[[k]]] <- gfun.prob[P[[k]],k]
  deviance <- -2 * sum(log(prob.vt))
  
  return(list("GP.model" = GP.model, "gfun.model" = gfun.model, 
              "gfun.prob" = gfun.prob, "deviance" = deviance))
}