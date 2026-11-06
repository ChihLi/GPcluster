GPcluster_Mstep <- function(X, Y, P, K, nugget, gfun, parallel, verbose){
  
  Z <- rep(0, length(Y))
  for (k in 1:K) Z[P[[k]]] <- k

  if(parallel){
    GP.model <- foreach::foreach(k = seq_len(K), .packages = "mlegp",
                                .export = "fit_cluster_gp",
                                .verbose = verbose) %dopar% {
      fit_cluster_gp(X[P[[k]], , drop = FALSE], Y[P[[k]]], nugget, verbose)
    }
  }else{
    GP.model <- lapply(seq_len(K), function(k) {
      fit_cluster_gp(X[P[[k]], , drop = FALSE], Y[P[[k]]], nugget, verbose)
    })
  }

  gfun.model <- if (K == 1) NULL else gfun.fit(X, Z, K, gfun, verbose)
  
  if (K == 1) {
    gfun.prob <- matrix(1, nrow(X), 1)
  }else if(gfun == "logistic"){
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