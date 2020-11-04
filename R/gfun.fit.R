gfun.fit <- function(X, Z, K, gfun){

  fit.data <- data.frame(X, "Z" = factor(Z))
  colnames(fit.data)[1:(ncol(fit.data)-1)] <- paste0("X", 1:(ncol(fit.data)-1))

  if(gfun == "logistic"){
    ### multinomial logistic regression model ###
    if(nrow(fit.data) <= 10000){ # fix nnet weight issue
      out <- nnet::multinom(Z ~ ., data = fit.data, type = "prob")
    }else{
      out <- nnet::multinom(Z ~ ., data = fit.data, type = "prob", MaxNWts = nrow(fit.data))
    }
  }else if(gfun == "lda"){
    ### LDA ###
    out <- MASS::lda(Z ~ ., data = fit.data)
  }else if(gfun == "qda"){
    ### QDA ###
    out <- MASS::qda(Z ~ ., data = fit.data)
  # }else if(gfun == "bagging"){ # under development
  #   ### bagging ###
  #   out <- randomForest(Z ~ ., data = fit.data, mtry = ncol(fit.data)-1)
  # }else if(gfun == "forest"){
  #   ### random forest ###
  #   out <- randomForest(Z ~ ., data = fit.data)
  #   #out$deviance <- 0
  # }else if(gfun == "gp"){
  #   ### binary gaussian process ###
  #   out <- gausspr(Z ~ ., data = fit.data)
  #   #out$deviance <- 0
  # }else if(gfun == "adaboost"){
  #   ### adaboost ###
  #   out <- adaboost(Z ~., data = fit.data, 10)
  # }else if(gfun == "constant"){
  #   out <- NULL
  # }else if(gfun == "forest.SRC"){
  #   ### random forest ###
  #   out <- rfsrc(Z ~ ., data = fit.data, sampsize = 1e4, ntree = 100)
  }
  return(out)
}
