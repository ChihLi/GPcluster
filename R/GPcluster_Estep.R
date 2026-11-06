GPcluster_Estep <- function(X, Y, P, K, M.models, gfun, n_min, n_max, parallel){
  
  GP.model <- M.models$GP.model
  gfun.prob <- M.models$gfun.prob
  
  ### stochastic hard assignment
  for(i in 1:nrow(X)){
    # print(i)
    x <- X[i,]
    y <- Y[i]
    ori.index <- which(sapply(P, function(x) any(x == i)))
    density.vt <- rep(-Inf, K)
    
      for(k in which(gfun.prob[i,]>0)){
        if(ori.index == k){
          rm.index <- which(i == P[[k]])
          mu_hat <- GP.model[[k]]$mu[1]
          Q_r <- GP.model[[k]]$invVarMatrix[rm.index,]
          Y_rm <- GP.model[[k]]$Z[-rm.index,1]
          mu.new <- mu_hat - sum(Q_r[-rm.index] * (Y_rm - mu_hat)) / Q_r[rm.index]
          if(Q_r[rm.index] < 0) {
            sigma.new <- 0
          }else {
            sigma.new <- sqrt(1/Q_r[rm.index])
          }
          density.vt[k] <- dnorm(y, mean = mu.new, sd = sigma.new, log = TRUE)
        }else{
          pred.out <- fastpredict.gp(GP.model[[k]], newData = x)
          mu.new <- pred.out$fit
          sigma.new <- pred.out$se.fit
          density.vt[k] <- dnorm(y, mean = mu.new, sd = sigma.new, log = TRUE)
        }
      }
    
    gfun.p <- gfun.prob[i,]
    assign.prob <- assignment_weights(density.vt, gfun.p, log = TRUE)

    n.group <- sapply(P,length)
    if(any(n.group >= n_max)){
      assign.prob[n.group >= n_max & seq_len(K) != ori.index] <- 0
    }
    if(n.group[ori.index] <= n_min | all(assign.prob == 0)){
      assign.prob[ori.index] <- 1
      assign.prob[-ori.index] <- 0
    }
    assign.prob <- assign.prob/sum(assign.prob)
    
    hard.assign <- rmultinom(1, size = 1, prob = assign.prob)
    new.index <- which(hard.assign[,1] == 1)
    
    if(new.index != ori.index) {
      ########## update GP.model and P #####
      
      ### 1.  delete i in ori.index ###
      keep.index <- P[[ori.index]] != i
      # update GP.model
      GP.model[[ori.index]]$Z <- matrix(GP.model[[ori.index]]$Z[keep.index,1], ncol = 1)
      GP.model[[ori.index]]$mu <- matrix(GP.model[[ori.index]]$mu[keep.index,1], ncol = 1)
      GP.model[[ori.index]]$X <- GP.model[[ori.index]]$X[keep.index, , drop= FALSE]
      GP.model[[ori.index]]$invVarMatrix <- inverse_delete(
        GP.model[[ori.index]]$invVarMatrix, which(!keep.index))

      # update P
      P[[ori.index]] <- P[[ori.index]][keep.index]
      
      ### 2.  add i in new index ###
      # update P
      P[[new.index]] <- sort(c(P[[new.index]], i))
      # update GP.model
      add.index <- P[[new.index]] == i
      newY.tmp <- GP.model[[new.index]]$Z[,1]
      newY.tmp <- append(newY.tmp, y, after = which(add.index) - 1)
      GP.model[[new.index]]$Z <- matrix(newY.tmp, ncol = 1)
      GP.model[[new.index]]$mu <- matrix(rep(GP.model[[new.index]]$mu[1,1], length(P[[new.index]])), ncol = 1)
      newX.tmp <- matrix(0, nrow = nrow(GP.model[[new.index]]$X) + 1, ncol = ncol(GP.model[[new.index]]$X))
      newX.tmp[!add.index,] <- GP.model[[new.index]]$X
      newX.tmp[add.index,] <- x
      GP.model[[new.index]]$X <- newX.tmp
      gp <- GP.model[[new.index]]
      r <- as.numeric(mlegp::calcCorOneObs(gp$X[!add.index, , drop = FALSE],
                                           gp$beta, gp$a, x)) * gp$sig2
      GP.model[[new.index]]$invVarMatrix <- inverse_add(
        gp$invVarMatrix, r, gp$sig2 + gp$nugget, which(add.index),
        covariance_matrix = mlegp::calcVarMatrix(gp$X, gp$beta, gp$a,
                                                gp$nugget, gp$sig2, 0, nrow(gp$X)))

    }
  }
  
  M.models$GP.model <- GP.model
  return(list(P = P, M.models = M.models))
  
}