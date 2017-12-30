## Inicializando algoritmo EM

In_EM <- function(X, K, y){
  library(mclust)
  beta <- array()
  X_lag <- matrix(nr = nrow(X), nc = ncol(X))
  for(i in 1:K){
    lin_mod <- lm(X[,i]~y)  
    beta[i] <- lin_mod$coefficients[2]
    X_lag[,i] <- X[,i] - beta[i]*y
  }
  betas <- rep(list(beta), K)

  fit <- Mclust(X_lag, G=K)
  alpha <- fit$parameters$pro
  mus <- lapply(apply(fit$parameters$mean, 2, list), unlist)
  variances <- list()
  vars <- fit$parameters$variance$sigma
  ini <- 1
  end <- K*K
  for(i in 1:K){
    variances[[i]] <- matrix(vars[ini:end], nrow=K)
    ini <- end+1
    end <- end+K*K
  }
  
  #   betas <- list(c(1,1), c(0.5,0.5))
  #   mus = list(c(20,30), c(25,35))
  #   betas = list(c(1,1), c(0.5,0.5))
  #   variances <- list(1000*diag(2), 1000*diag(2))
  #   alpha <- c(0.5,0.5)
  
  inicialization <- list(mu = mus, sig = variances, alpha = alpha, beta = betas)
  return(inicialization)
}


## E-Step : Calculating probabilities of each cluster

E.step <- function(X, y, phi, K) {
  prob <- matrix(nrow=nrow(X), ncol=K)
  for(grupo in 1:K){
    for(i in 1:nrow(X)){
      prob[i,grupo] <- dmvnorm(X[i,],
                               phi$mu[[grupo]]+phi$beta[[grupo]]*y[i],
                               phi$sig[[grupo]])+10^(-10)
    }
  }
  return(prob/rowSums(prob))
}


## M-Step

M.step <- function(X, y, phi, probs, K) {
  
  x_beta <- list(cbind(phi$beta[[1]][1]*y, 
                       phi$beta[[1]][2]*y),
                 cbind(phi$beta[[2]][1]*y, 
                       phi$beta[[2]][2]*y))

  mus <- list(
    colSums((X-x_beta[[1]])*probs[,1])/
      colSums(probs)[1], 
    colSums((X-x_beta[[2]])*probs[,2])/
      colSums(probs)[2]        
  )
  
  betas <- list(
    colSums(sweep(X, 2, phi$mu[[1]])*y*as.vector(matrix(rep(probs[,1],2), nc=2, byrow = F)))/
      sum(y^(2)*probs[,1]),
    colSums(sweep(X, 2, phi$mu[[2]])*y*as.vector(matrix(rep(probs[,2],2), nc=2, byrow = F)))/
      sum(y^(2)*probs[,2]))
  
  covs <- lapply(1:K, function(i) cov.wt((X-x_beta[[i]]), probs[,i]))
  phi$sig <- lapply(covs, "[[", "cov")
  # mat1 <- matrix(0,nr=2,nc=2)
  # mat2 <- matrix(0,nr=2,nc=2)
  # for(i in 1:nrow(X)){
  #   m1 <- as.matrix(X[i,]-phi$mu[[1]]-phi$beta[[1]]*y[i])
  #   m2 <- as.matrix(X[i,]-phi$mu[[2]]-phi$beta[[2]]*y[i])
  #   mat1 <- mat1 + probs[i,1]*t(m1)%*%m1
  #   mat2 <- mat2 + probs[i,2]*t(m2)%*%m2
  # }
  # mat1 = mat1/sum(probs[,1])
  # mat2 = mat2/sum(probs[,2])
  # phi$sig <- list(mat1, mat2)

  phi$alpha <- colMeans(probs)
  
  phi$beta <- betas
  phi$mu <- mus

  return(phi)
}


log.like <- function(X, y, phi, K) {
  prob <- matrix(nrow=nrow(X), ncol=K)
  for(grupo in 1:K){
    for(i in 1:nrow(X)){
      prob[i,grupo] <- phi$alpha[grupo]*dmvnorm(X[i,],
                                                phi$mu[[grupo]]+phi$beta[[grupo]]*y[i],
                                                phi$sig[[grupo]])
    }
  }
  sum(log(rowSums(prob)))
}


run.em <- function(X, y, K, max_iter=100) {
  #learning_rate = 1
  phi <- In_EM(X,K = K, y=y)
  
  for(i in 1:max_iter) {
    oldphi <- phi
    probs <- E.step(X, y, phi, K)
    #probs <- t(apply(probs, 1, function(x) 0.5 + (x - 0.5) * (1-learning_rate) / (1+learning_rate)))
    #probs <- probs/rowSums(probs)
    phi <- M.step(X, y, phi, probs, K)
    #learning_rate = max(learning_rate - 0.01,0)
    if((log.like(X, y, phi, K) - log.like(X, y, oldphi, K)) < 0.01)
      break
  }
  return(phi)
}




require(MASS)
require(mvtnorm)

## Probability of each cluster for each point
E.step.old <- function(X, phi, N) {
  h <-
    with(phi, do.call(cbind,
                      lapply(1:N, function(i)
                        dmvnorm(X, mu[[i]], sig[[i]]))))
  h/rowSums(h)                       #Normalize
}

## Given the probability of each cluster for each point, we try to
## find the values of mu, sigma, and alpha that maximize the likelihood
M.step.old <- function(X, h, N) {
  covs <- lapply(1:N, function(i) cov.wt(X, h[,i]))
  mu <- lapply(covs, "[[", "center")
  sig <- lapply(covs, "[[", "cov")
  alpha <- colMeans(h)
  list(mu = mu, sig = sig, alpha = alpha)
}

log.like.old <- function(X, phi, N) {
  probs <- 
    with(phi, do.call(cbind,
                      lapply(1:N, function(i)
                        alpha[i] * dmvnorm(X, mu[[i]], sig[[i]]))))
  sum(log(rowSums(probs)))
}

run.em.old <- function(X, K, max_iter=100) {
  N <- nrow(X)
  covs <- replicate(N, list(cov.wt(X[sample(nrow(X), 30),])))
  mu <- lapply(covs, "[[", "center")
  sig <- lapply(covs, "[[", "cov")
  alpha <- rep(1/N, N)
  phi <<- list(mu = mu, sig = sig, alpha = alpha, beta = list(c(0,0), c(0,0)))
  
  for(i in 1:max_iter) {
    oldphi <- phi
    probs <- E.step.old(X, phi, K)
    phi <- M.step.old(X, probs, K)
    if((log.like.old(X, phi, K) - log.like.old(X, oldphi, K)) < 0.01)
      break
  }
  return(phi)
}


clusters_perf <- function(real, probs){
  cl <- apply(probs, 1, which.max)
  pair <- comembership_table(real, cl)
  return(c(
    (pair[[1]][1]+pair[[4]][1])/sum(unlist(pair)),
    adjustedRandIndex(cl, real)))
}


euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

dist_parameters <- function(est, sim){
  return(min(
    euc_dist(est[[1]], sim[[1]]) + euc_dist(est[[2]], sim[[2]]),
    euc_dist(est[[1]], sim[[2]]) + euc_dist(est[[2]], sim[[1]])))
}

simulations <- function(seed, n=100){
  library(MASS)
  set.seed(seed)
  x = rnorm(n, mean=runif(1,0,10), sd = runif(1,0,5))
  
  set.seed(seed*100)
  betas1 = c(runif(1,2,5), runif(1,2,5))
  betas2 = c(runif(1,-5,-2), runif(1,-5,-2))
  mean1 = c(0,0)
  mean2 = c(10,0)
  
  set.seed(runif(1,0,1000))
  s1 <- runif(1,1,10)
  s2 <- runif(1,1,10)
  sigma1 <- matrix(c(s1,0,0,s1), nrow = 2)
  sigma2 <- matrix(c(s2,0,0,s2), nrow = 2)
  
  df = matrix(nrow=n, ncol=3)
  cont <- 1
  set.seed(seed)
  for(i in 1:(n/2)){
    df[cont,] = cbind(matrix(mvrnorm(1,
                                c(mean1[1]+betas1[1]*x[cont],
                                  mean1[2]+betas1[2]*x[cont]),
                                sigma1), nr=1), 1)
    cont<-cont+1
    df[cont,] = cbind(matrix(mvrnorm(1,
                                c(mean2[1]+betas2[1]*x[cont],
                                  mean2[2]+betas2[2]*x[cont]),
                                sigma2), nr=1), 2)
    cont<-cont+1
  }
  df <- data.frame(df)
  output <- list(x, mean1, mean2, s1, s2, betas1, betas2, df)
  names(output) <- c("x", "mean1", "mean2", "s1", "s2",
                     "betas1", "betas2", "df") 
  return(output)
}

perf_simulations <- function(n=200, K=2){
  output <- matrix(nrow = n, ncol = 28)
  for(i in 1:n){
    dataset <- simulations(seed=i)
    X <- dataset$df[,-ncol(dataset$df)]
    resp <- dataset$df[,ncol(dataset$df)]
    y <- dataset$x
    
    proposal <- run.em(X, y, K)
    current <- run.em.old(X, K)
    alldata <- run.em.old(cbind(X,y), K)
    partial <- In_EM(X,K,y)
    
    proposal.perf <- clusters_perf(resp, E.step(X, y, proposal, K))
    current.perf <- clusters_perf(resp, E.step.old(X, current, K))
    all.perf <- clusters_perf(resp, E.step.old(cbind(X,y), alldata, K))
    partial.perf <- clusters_perf(resp, E.step(X, y, partial, K))
    
    betas_sim <- list(dataset$betas1, dataset$betas2)
    distBeta.prop <- dist_parameters(proposal$beta, betas_sim)
    distBeta.curr <- dist_parameters(rep(list(c(0,0)), 2), betas_sim)
    distBeta.alll <- dist_parameters(rep(list(c(0,0)), 2), betas_sim)
    distBeta.part <- dist_parameters(partial$beta, betas_sim)
    
    mus_sim <- list(dataset$mean1, dataset$mean2)
    distMu.prop <- dist_parameters(proposal$mu, mus_sim)
    distMu.curr <- dist_parameters(current$mu, mus_sim)
    distMu.alll <- dist_parameters(lapply(alldata$mu, 
                                          function(x) x[1:2]), mus_sim)
    distMu.part <- dist_parameters(partial$mu, mus_sim)
    
    output[i,] <- c(mean(y),
                    sd(y),
                    dataset$mean1[1],
                    dataset$mean1[2],
                    dataset$mean2[1],
                    dataset$mean2[2],
                    dataset$s1,
                    dataset$s2,
                    dataset$betas1[1],
                    dataset$betas1[2],
                    dataset$betas2[1],
                    dataset$betas2[2],
                    proposal.perf[1],
                    proposal.perf[2],
                    current.perf[1],
                    current.perf[2],
                    all.perf[1],
                    all.perf[2],
                    partial.perf[1],
                    partial.perf[2],
                    distBeta.prop,
                    distBeta.curr,
                    distBeta.alll,
                    distBeta.part,
                    distMu.prop,
                    distMu.curr,
                    distMu.alll,
                    distMu.part)
  }
  colnames(output) <- c('mean_co',
                     'sd_co',
                     'mean1_1',
                     'mean1_2',
                     'mean2_1',
                     'mean2_2',
                     'sigma1',
                     'sigma2',
                     'beta1_1',
                     'beta1_2',
                     'beta2_1',
                     'beta2_2',
                     'jacc_cocem',
                     'rand_cocem',
                     'jacc_cem',
                     'rand_cem',
                     'jacc_all',
                     'rand_all',
                     'jacc_mix',
                     'rand_mix',
                     'beta_dist_cocem',
                     'beta_dist_cem',
                     'beta_dist_all',
                     'beta_dist_mix',
                     'mu_dist_cocem',
                     'mu_dist_cem',
                     'mu_dist_all',
                     'mu_dist_mix')
  return(output)
}

res <- perf_simulations()

jaccs <- res[,c(14,16,18, 20)]
jaccs <- data.frame(jaccs)
jaccs <- apply(jaccs, c(1,2), as.numeric)

plot(jaccs[,2], jaccs[,1], xlab='CEM', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index')
abline(a=0,b=1)
plot(jaccs[,3], jaccs[,1], xlab='All', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index')
abline(a=0,b=1)
plot(jaccs[,4], jaccs[,1], xlab='2step', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index')
abline(a=0,b=1)

df <- data.frame(res)
df[,'dist_betas'] <- sqrt(((df$beta1_1-df$beta2_1)^2+
                            (df$beta1_2-df$beta2_2)^2)/2)
df[,'dist_means'] <- sqrt(((df$mean1_1-df$mean2_1)^2+
                             (df$mean1_2-df$mean2_2)^2)/2)

df[,'class'] <- as.factor(ifelse(df$rand_cocem < df$rand_cem, 1, 0))

library(randomForest)
fit <- randomForest(as.formula('class ~ 
                               dist_betas+dist_means+mean_co+sd_co'),
             data=df, ntree=200, importance=T)

fit$importance

diag = function(obs, prev, label = '', plot=TRUE){
  library(ROCR)
  pred = prediction(prev, obs)
  perf <- performance(pred,"tpr","fpr")
  if(plot==TRUE){
    plot(perf, main=paste0(label, '  (ROC Curve)'), col='red', lwd=2)
    abline(a=0,b=1, lwd=2)
  }
  KS <- max(attr(perf,'y.values')[[1]]-attr(perf,'x.values')[[1]])
  perf = performance(pred,"auc")
  GINI<-2*attr(perf,'y.values')[[1]]-1
  ROC<-attr(perf,'y.values')
  saida = list(GINI, unlist(ROC), KS)
  names(saida) <- c('GINI','ROC', 'KS')
  return(saida)
}

pred <- predict(fit, newdata=df, type='prob')[,2]
diag(df$class, pred)


fit <- glm(as.formula('class ~ 
                               dist_betas+dist_means+mean_co+sd_co'),
                    data=df, family = 'binomial')


M.step.t <- function(X, y, phi, probs, K) {
  
  soma1_1 <- 0
  soma1_2 <- 0
  soma2_1 <- 0
  soma2_2 <- 0
  for(i in 1:nrow(X)){
    soma1_1 <- soma1_1 + (X[i,1]-
                            phi$beta[[1]][1]*y[i])*probs[i,1]
    soma1_2 <- soma1_2 + (X[i,2]-
                            phi$beta[[1]][2]*y[i])*probs[i,1]
    soma2_1 <- soma2_1 + (X[i,1]-
                            phi$beta[[2]][1]*y[i])*probs[i,2]
    soma2_2 <- soma2_2 + (X[i,2]-
                            phi$beta[[2]][2]*y[i])*probs[i,2]
  }
  mus <- list(c(soma1_1/colSums(probs)[1],soma1_2/colSums(probs)[1]),
              c(soma2_1/colSums(probs)[2],soma2_2/colSums(probs)[2]))
  
  soma1_1 <- 0
  soma1_2 <- 0
  soma2_1 <- 0
  soma2_2 <- 0
  for(i in 1:nrow(X)){
    soma1_1 <- soma1_1 + (X[i,1]-phi$mu[[1]][1])*probs[i,1]
    soma1_2 <- soma1_2 + (X[i,2]-phi$mu[[1]][2])*probs[i,1]
    soma2_1 <- soma2_1 + (X[i,1]-phi$mu[[2]][1])*probs[i,2]
    soma2_2 <- soma2_2 + (X[i,2]-phi$mu[[2]][2])*probs[i,2]
  }
  betas <- list(c(soma1_1/sum(y*probs[,1]),soma1_2/sum(y*probs[,1])),
                c(soma2_1/sum(y*probs[,2]),soma2_2/sum(y*probs[,2])))
  
  x_beta <- list(cbind(phi$beta[[1]][1]*y, 
                       phi$beta[[1]][2]*y),
                 cbind(phi$beta[[2]][1]*y, 
                       phi$beta[[2]][2]*y))
  
  covs <- lapply(1:K, function(i) cov.wt((X-x_beta[[i]]), probs[,i]))
  phi$sig <- lapply(covs, "[[", "cov")
  
  phi$alpha <- colMeans(probs)
  
  phi$beta <- betas
  phi$mu <- mus
  
  return(phi)
}
