library(dlm)
library(numDeriv)

generate_dlm_data <- function(n_subjects,n_obs,m0=0,C0=1,V=1,W=1,GG=1,FF=1,intercept=0) {
  Theta <- matrix(NA,n_subjects,n_obs)
  Y <- matrix(NA,n_subjects,n_obs)
  
  Theta0 <- m0 + sqrt(C0)*rnorm(n_subjects)
  #cat("Generated m0:", mean(Theta0),"\n")
  #cat("Generated C0:", var(Theta0), "\n")
  
  for (i in 1:ncol(Y)) {
    if (i==1) {Theta[,i] <- GG*Theta0 + intercept + sqrt(W)*rnorm(n_subjects)}
    else {Theta[,i] <- GG*Theta[,i-1] + intercept + sqrt(W)*rnorm(n_subjects)}
    
    Y[,i] <- FF*Theta[,i] + sqrt(V)*rnorm(n_subjects)
  }
  return(list(Theta=Theta,Y=Y))
}

model_to_estimate <- function(x,GG) {
  names(x) <- NULL
  dlm(m0=c(x[1],1),
      C0=diag(c(exp(x[2]),0)),
      V=exp(x[3]), W=diag(c(exp(x[4]),0)),
      GG=matrix(c(ifelse(is.null(GG),x[6],GG),0,x[5],1),2,2),
      FF=matrix(c(1,0),ncol=2))
}

model_to_estimate_weighted <- function(x,GG, weight) {
  names(x) <- NULL
  
  X=cbind(rep(exp(x[4]),length(weight)) * weight,
          rep(x[5],length(weight)) * weight)
  
  dlm(m0=c(x[1],1),
      C0=diag(c(exp(x[2]),0)),
      V=exp(x[3]), 
      W=diag(c(exp(x[4]),0)),
      JW=diag(c(1,0)),
      X=X,
      JGG=matrix(c(0,0,2,0),2,2),
      GG=matrix(c(1,0,x[5],1),2,2),
      FF=matrix(c(1,0),ncol=2))
}

logLik_dlmTrend <- function(parm, GG, Y, weight) {
  if (is.null(weight)) {
    mod <- model_to_estimate(parm, GG)
  } else {
    mod <- model_to_estimate_weighted(parm, GG, weight)
  }
  # tmp <- La.svd(mod$V, nu = 0)
  # Dv <- sqrt(tmp$d)
  # if (any(Dv < .Machine$double.eps^0.3)) {
  #   stop("Degenerated V")
  # }
  
  LL<-0
  for (i in 1:nrow(Y)) {
    suppressWarnings(LL <- LL + dlmLL(y = Y[i,], mod = mod, debug = FALSE))
  }
  return(LL)
}

#' Estimate dlm model with trend using Method of Moments estimator
#' @param Y matrix; Input data
#' @param GG double; Transistion parameter (set to NULL to estimate)
#' @param mle logical; Should the mme estimate be used to initialize an maximum likelihood estimation?
dlmTrend_mme <- function(Y,GG=NULL,mle=FALSE)  {
  
  #stop("Weighting is not properly implemented yet.")
  
  #if(!is.null(GG)) {stop("Not implemented because of the weighting")}
  #if (is.null(weight)) { weight <- rep(1,ncol(Y)) }
  
  M <- colMeans(Y,na.rm = TRUE)
  Q <- cov(Y,use="pairwise.complete.obs")
  a <- diag(Q[-1,,drop=FALSE])/diag(Q[-nrow(Q),,drop=FALSE])
  if (!is.null(GG)) {a<-rep(GG,length(a))}
  b <- (M[-1] - a*M[-length(M)])
  
  D <- cov(Y[,-1]-sweep(sweep(Y[,-ncol(Y)],2,a,"*"),2,b,"-"),use="pairwise.complete.obs") # Misschien 1/(N-2) ipv 1/(N-1)
  Dt <- diag(D)
  Dt_tplus <- diag(D[-1,,drop=FALSE]) 
  rho <- Dt_tplus / sqrt(Dt[-length(Dt)]*Dt[-1])
  
  D_hat <- mean(Dt,na.rm=TRUE)
  a_hat <- mean(a,na.rm=TRUE)
  b_hat <- mean(b,na.rm=TRUE)
  rho_hat <- mean(rho,na.rm=TRUE)
  
  V <- -rho_hat * D_hat/a_hat
  W <- D_hat - (1+a_hat^2)*V
  C <- (Q[1,1] - V - W)/(a_hat^2)
  m <- (M[1] - b_hat)/a_hat
  
  if (is.null(GG)) {
    init <- c(m,log(C),log(V),log(W),b_hat,a_hat)
  } else {
    init <- c(m,log(C),log(V),log(W),b_hat)
  }
  
  if (mle) {
    opt <- optim(init, logLik_dlmTrend, method = "L-BFGS-B", GG=GG, Y=Y, weight=NULL)
    init <- opt$par
  }
  
  out <- list()
  out$hessian <- hessian(logLik_dlmTrend,x = opt$par,GG=GG, Y=Y, weight=NULL)
  out$loglikelihood <- logLik_dlmTrend(init,GG,Y,weight=NULL)        
  out$model <- model_to_estimate_weighted(opt$par, weight=NULL)
  out$opt$par <- init
  
  class(out) <- "dlmTrend"
  out
}

# data <- generate_dlm_data(n_subjects=500,n_obs=10,m0=35,C0=45,V=5,W=10,GG=1,FF=1,intercept=10)
# 
# dlmTrend_mme(data$Y,1)
# dlmTrend_mme_mle(data$Y,1)
# 
# print(dlmTrend_mme(data$Y))
# print(dlmTrend_mme_mle(data$Y))

#' Estimate dlm model with trend using Maximum Likelhood estimator
#' @param Y matrix; Input data
#' @param GG double; Transistion parameter (set to NULL to estimate)
dlmTrend_mle <- function(Y, GG=NULL, init=NULL, weight=NULL) {
  
  if (!is.null(weight) & ncol(Y)!=length(weight)) stop("Incorrect number of weights")
  
  if (is.null(GG)) {
    if (is.null(init)) {init <- c(runif(5)*4 - 2,1) }
    opt <- optim(init, logLik_dlmTrend, method = "L-BFGS-B", GG=GG, Y=Y,weight=weight)
  } else {
    if (is.null(init)) {init <- runif(5)*4 -2 }
    opt <- optim(init, logLik_dlmTrend, method = "L-BFGS-B", GG=GG, Y=Y,weight=weight)
  }
  
  out <- list()
  out$hessian <- hessian(logLik_dlmTrend,x = opt$par,GG=GG, Y=Y, weight=weight)
  out$loglikelihood <- opt$value     
  out$model <- model_to_estimate_weighted(opt$par, weight=weight)
  out$opt <- opt
  
  class(out) <- "dlmTrend"
  out
}

dlmTrend_mle_max <- function(Y, GG=NULL, init=NULL, weight=NULL, restart=1) {
  
  models <- map(1:restart, 
                function(x) { dlmTrend_mle(Y, GG, init, weight) }
  )

  best_model <- which.max(sapply(models, function(x) x$loglikelihood))
  out <- model[[best_model]]
  out$models <- models
  class(out) <- "dlmTrend"
  out
}

summary.dlmTrend <- function(object,...) {
  print(object$model)
}

logLik.dlmTrend <- function(object,...) {
  -object$loglikelihood
}

glance.dlmTrend <- function(object,...) {
  suppressWarnings(tibble(
    V = c(V(object$model)),
    W = W(object$model)[1,1],
    Trend = GG(object$model)[1,2],
    G = GG(object$model)[1,1],
    m0 = m0(object$model)[1],
    C0 = C0(object$model)[1,1],
    LogLikelihood = -object$loglikelihood
  ))
}

tidy.dlmTrend <- function(object,...,conf=TRUE) {
  par_un <- object$opt$par
  
  if (conf) {
    H <- object$hessian
    R <- chol(H,pivot=TRUE)
    V <- chol2inv(R)
    rownames(V) <- colnames(V) <- colnames(H)
    
    SIMS <- 1000
    len <- length(par_un)
    unconstrained <- par_un + t(chol(V)) %*% 
      matrix(rnorm(SIMS * len), nrow = len, ncol = SIMS)
    theta_sims <- t(apply(unconstrained, 2, FUN = function(upars) {
      c(upars[1],exp(upars[2]),exp(upars[3]),exp(upars[4]),upars[5])
    }))
   
    conf_int <- apply(theta_sims,2,function(x) {quantile(x,c(0.025,0.975))})
    se <- apply(theta_sims, 2, sd)
  }
  if (conf) {
    tibble(
      term = c("V","W","Trend","m0","C0","G"),
      estimate = suppressWarnings(c(V(object$model), W(object$model)[1,1], GG(object$model)[1,2],
                 m0(object$model)[1],C0(object$model)[1,1], GG(object$model)[1,1])),
      std.error = c(se[c(3,4,5,1,2)],NA),
      low=c(conf_int[1,c(3,4,5,1,2)],NA),
      high=c(conf_int[2,c(3,4,5,1,2)],NA)
    )
  } else {
    tibble(
      term = c("V","W","Trend","m0","C0","G"),
      estimate = suppressWarnings(c(V(object$model), W(object$model)[1,1], GG(object$model)[1,2],
                                    m0(object$model)[1],C0(object$model)[1,1], GG(object$model)[1,1]))
    )
  }
}

cnfsims <- function(object,...) {
  par_un <- object$opt$par
  H <- object$hessian
  R <- chol(H,pivot=TRUE)
  V <- chol2inv(R)
  rownames(V) <- colnames(V) <- colnames(H)
  
  SIMS <- 1000
  len <- length(par_un)
  unconstrained <- par_un + t(chol(V)) %*% 
    matrix(rnorm(SIMS * len), nrow = len, ncol = SIMS)
  theta_sims <- t(apply(unconstrained, 2, FUN = function(upars) {
    c(upars[1],exp(upars[2]),exp(upars[3]),exp(upars[4]),upars[5])
  }))
  
  colnames(theta_sims) <- c("V","W","Trend","m0","C0","G")
  
  tibble(
    term = c("V","W","Trend","m0","C0","G"),
    estimate = suppressWarnings(c(V(object$model), W(object$model)[1,1], GG(object$model)[1,2],
                                  m0(object$model)[1],C0(object$model)[1,1], GG(object$model)[1,1])),
    std.error = c(se[c(3,4,5,1,2)],NA),
    low=c(conf_int[1,c(3,4,5,1,2)],NA),
    high=c(conf_int[2,c(3,4,5,1,2)],NA)
  )
}
