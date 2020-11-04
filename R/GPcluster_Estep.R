GPcluster_Estep <- function(X, Y, P, K, M.models, gfun, n_min, n_max, parallel){
  
  GP.model <- M.models$GP.model
  gfun.prob <- M.models$gfun.prob
  
  ### stochastic hard assignment
  for(i in 1:nrow(X)){
    # print(i)
    x <- X[i,]
    y <- Y[i]
    ori.index <- which(sapply(P, function(x) any(x == i)))
    density.vt <- rep(0, K)
    
    #if(parallel){   # parallel is slower
    if(FALSE){
      density.vt <- foreach::foreach(k = seq(K), .combine = c, .packages = "mlegp")%dopar%{
        if(ori.index == k){
          rm.index <- which(i == P[[k]])
          mu_hat <- GP.model[[k]]$mu[1]
          Q_r <- GP.model[[k]]$invVarMatrix[rm.index,]
          Y_rm <- GP.model[[k]]$Z[-rm.index,1]
          mu.new <- mu_hat - sum(Q_r[-rm.index] * (Y_rm - mu_hat)) / Q_r[rm.index]
          sigma.new <- sqrt(1/Q_r[rm.index])
          out <- dnorm(y, mean = mu.new, sd = sigma.new)
        }else{
          pred.out <- predict(GP.model[[k]], newData = x, se.fit = TRUE)
          mu.new <- pred.out$fit
          sigma.new <- pred.out$se.fit
          out <- dnorm(y, mean = mu.new, sd = sigma.new)
        }
        return(out)
      }
    }else{
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
          density.vt[k] <- dnorm(y, mean = mu.new, sd = sigma.new)
        }else{
          pred.out <- fastpredict.gp(GP.model[[k]], newData = x)
          mu.new <- pred.out$fit
          sigma.new <- pred.out$se.fit
          density.vt[k] <- dnorm(y, mean = mu.new, sd = sigma.new)
        }
      }
    }
    
    gfun.p <- gfun.prob[i,]
    if(any(!is.finite(density.vt))){
      if(any(gfun.p[!is.finite(density.vt)] > 0)) {
        assign.prob <- rep(0, length(density.vt))
        assign.prob[!is.finite(density.vt)] <- gfun.p[!is.finite(density.vt)]
      }else{
        density.vt[!is.finite(density.vt)] <- 0 
      }
    }else{
      assign.prob <- density.vt * gfun.p
    }
    
    n.group <- sapply(P,length)
    if(any(n.group >= n_max)){
      assign.prob[n.group >= n_max] <- 0
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
      if(TRUE){
        # use Woodbury formula: see https://math.stackexchange.com/questions/1248220/find-the-inverse-of-a-submatrix-of-a-given-matrix
        r <- mlegp::calcCorOneObs(GP.model[[ori.index]]$X, GP.model[[ori.index]]$beta, 
                                  GP.model[[ori.index]]$a, x)
        U <- matrix(0, ncol = 2, nrow = length(P[[ori.index]]))
        U[!keep.index, 1] <- 1
        U[keep.index, 2] <- r
        V <- matrix(0, ncol = length(P[[ori.index]]), nrow = 2)
        V[1, keep.index] <- r
        V[2, !keep.index] <- 1
        R.inv <- GP.model[[ori.index]]$invVarMatrix * GP.model[[ori.index]]$sig2
        R.inv <- R.inv + R.inv %*% U %*% solve(diag(1,2) - V %*% R.inv %*% U , V) %*% R.inv
        GP.model[[ori.index]]$invVarMatrix <- R.inv[keep.index, keep.index]/GP.model[[ori.index]]$sig2
      }else{
        GP.model[[ori.index]]$invVarMatrix <- solve(calcVarMatrix(GP.model[[ori.index]]$X, GP.model[[ori.index]]$beta, GP.model[[ori.index]]$a,GP.model[[ori.index]]$nugget, GP.model[[ori.index]]$sig2, 0, dim(GP.model[[ori.index]]$X)[1]))
      }
      
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
      R.inv <- GP.model[[new.index]]$invVarMatrix * GP.model[[new.index]]$sig2 
      r <- mlegp::calcCorOneObs(GP.model[[new.index]]$X, GP.model[[new.index]]$beta, 
                                GP.model[[new.index]]$a, x)
      U <- r
      U <- U[,!add.index, drop = FALSE]
      s <- 1/(1 + GP.model[[new.index]]$nugget/GP.model[[new.index]]$sig2 - U %*% R.inv %*% t(U))
      newRinv.tmp <- matrix(0, ncol(R.inv) + 1, ncol(R.inv) + 1)
      newRinv.tmp[add.index, add.index] <- s[1]
      newRinv.tmp[add.index, !add.index] <- newRinv.tmp[!add.index, add.index] <- - U %*% R.inv * s[1]
      newRinv.tmp[!add.index, !add.index] <- R.inv + R.inv %*% t(U) %*% U %*% R.inv * s[1]
      GP.model[[new.index]]$invVarMatrix <- newRinv.tmp / GP.model[[new.index]]$sig2
    }
  }
  
  return(list(P = P, M.models = M.models))
  
}