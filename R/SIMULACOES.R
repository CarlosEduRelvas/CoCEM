In_EM <- function(X, K, y){
  library(mclust)
  beta <- array()
  X_lag <- matrix(nr = nrow(X), nc = ncol(X))
  for(i in 1:ncol(X)){
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
  end <- ncol(X)*ncol(X)
  for(i in 1:K){
    variances[[i]] <- matrix(vars[ini:end], nrow=ncol(X))
    ini <- end+1
    end <- end+ncol(X)*ncol(X)
  }
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
  x_beta <- list()
  mus <- list()
  betas <- list()
  for(i in 1:K){
    x_beta[[i]] <- sweep(do.call(cbind, 
                                 rep(list(y),ncol(X))),
                         MARGIN=2, 
                         phi$beta[[i]],
                         `*`)
    
    mus[[i]] <- colSums((X-x_beta[[i]])*probs[,i])/colSums(probs)[i]
    
    betas[[i]] <- colSums(
      sweep(X, 2, phi$mu[[i]])*y*as.vector(matrix(
        rep(probs[,i],ncol(X)), nc=ncol(X), byrow = F)))/
      sum(y^(2)*probs[,i])
  }
  
  covs <- lapply(1:K, function(i) cov.wt((X-x_beta[[i]]), probs[,i]))
  phi$sig <- lapply(covs, "[[", "cov")
  
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
  phi <- In_EM(X,K = K, y=y)
  
  for(i in 1:max_iter) {
    oldphi <- phi
    probs <- E.step(X, y, phi, K)
    phi <- M.step(X, y, phi, probs, K)
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


generic_simulations <- function(seed, n=100, P=2, K=2){
  library(MASS)
  set.seed(seed)
  x = rnorm(n*K, mean=runif(1,0,10), sd = runif(1,0,5))
  
  set.seed(seed*100)
  random_betas <- c(runif(500,2,5), runif(500,-5,-2))
  random_means <- c(-10,-5,0,5,10)
  
  betas <- list()
  means <- list()
  sigma <- list()
  sssss <- list()
  for(i in 1:K){
    sssss[[i]] <- runif(1,1,10)
    sigma[[i]] <- diag(rep(sssss[[i]], P))
    betas[[i]] <- sample(random_betas, P, replace = T)
    means[[i]] <- sample(random_means, P, replace = T)
  }
  
  df = matrix(nrow=n*K, ncol=(P+1))
  cont <- 1
  set.seed(seed)
  for(j in 1:n){
    for(i in 1:K){
      df[cont,] = cbind(matrix(mvrnorm(1,
                                       means[[i]]+betas[[i]]*x[cont],
                                       sigma[[i]]), nr=1), i)
      cont<-cont+1
    }
  }
  df <- data.frame(df)
  output <- list(x, means, sssss, betas, df)
  names(output) <- c("x", "means", "s", "betas","df") 
  return(output)
}


perf_simulations <- function(n=200, K=2, P=2, KSIM=2){
  output <- matrix(nrow = n, ncol = 21)
  for(i in 1:n){
    dataset <- generic_simulations(seed=i,P = P, K=KSIM)
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
    
    distBeta.prop <- dist_parameters(proposal$beta, dataset$betas)
    distBeta.curr <- dist_parameters(rep(list(c(0,0)), 2), dataset$betas)
    distBeta.alll <- dist_parameters(rep(list(c(0,0)), 2), dataset$betas)
    distBeta.part <- dist_parameters(partial$beta, dataset$betas)
    
    distMu.prop <- dist_parameters(proposal$mu, dataset$means)
    distMu.curr <- dist_parameters(current$mu, dataset$means)
    distMu.alll <- dist_parameters(lapply(alldata$mu, 
                                          function(x) x[1:(length(x)-1)]), 
                                   dataset$means)
    distMu.part <- dist_parameters(partial$mu, dataset$means)
    
    output[i,] <- c(mean(y),
                    sd(y),
                    paste0(dataset$means, collapse = ","),
                    paste0(dataset$s, collapse = ","),
                    paste0(dataset$betas, collapse = ","),
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
                        'means',
                        'sigmas',
                        'betas',
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

### K = 2, P = 2
res <- perf_simulations()

jaccs <- res[,c(7,9,11,13)]
jaccs <- data.frame(jaccs)
jaccs <- apply(jaccs, c(1,2), as.numeric)

par(mfrow=c(2,2))
plot(jaccs[,2], jaccs[,1], xlab='CEM', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index')
abline(a=0,b=1)
plot(jaccs[,3], jaccs[,1], xlab='All', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index')
abline(a=0,b=1)
plot(jaccs[,4], jaccs[,1], xlab='2step', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index')
abline(a=0,b=1)



### K = 3, P = 2
res <- perf_simulations(K = 3,P = 2,KSIM = 3)

jaccs <- res[,c(7,9,11,13)]
jaccs <- data.frame(jaccs)
jaccs <- apply(jaccs, c(1,2), as.numeric)

par(mfrow=c(2,2))
plot(jaccs[,2], jaccs[,1], xlab='CEM', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index')
abline(a=0,b=1)
plot(jaccs[,3], jaccs[,1], xlab='All', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index')
abline(a=0,b=1)
plot(jaccs[,4], jaccs[,1], xlab='2step', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index')
abline(a=0,b=1)



### K = 5, P = 2
res <- perf_simulations(K = 5,n = 50,P = 2,KSIM = 5)

jaccs <- res[,c(7,9,11,13)]
jaccs <- data.frame(jaccs)
jaccs <- apply(jaccs, c(1,2), as.numeric)

par(mfrow=c(2,2))
plot(jaccs[,2], jaccs[,1], xlab='CEM', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index', xlim=c(0,1), ylim=c(0,1))
abline(a=0,b=1)
plot(jaccs[,3], jaccs[,1], xlab='All', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index', xlim=c(0,1), ylim=c(0,1))
abline(a=0,b=1)
plot(jaccs[,4], jaccs[,1], xlab='2step', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index', xlim=c(0,1), ylim=c(0,1))
abline(a=0,b=1)



### K = 4, P = 3
res <- perf_simulations(K = 4,n=50,P = 3,KSIM = 4)

jaccs <- res[,c(7,9,11,13)]
jaccs <- data.frame(jaccs)
jaccs <- apply(jaccs, c(1,2), as.numeric)

par(mfrow=c(2,2))
plot(jaccs[,2], jaccs[,1], xlab='CEM', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index', xlim=c(0,1), ylim=c(0,1))
abline(a=0,b=1)
plot(jaccs[,3], jaccs[,1], xlab='All', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index', xlim=c(0,1), ylim=c(0,1))
abline(a=0,b=1)
plot(jaccs[,4], jaccs[,1], xlab='2step', ylab = 'COCEM', pch=16,
     main='Adjusted Rand Index', xlim=c(0,1), ylim=c(0,1))
abline(a=0,b=1)





### K = 2, P = 2
### Initial Guess


dataset <- generic_simulations(seed=4325,P = 2, K=2)
X <- dataset$df[,-ncol(dataset$df)]
resp <- dataset$df[,ncol(dataset$df)]
y <- dataset$x

In_EM_random <- function(X, K, y){
  betas <- list(c(runif(1,-5,5), runif(1,-5,5)),
                c(runif(1,-5,5), runif(1,-5,5)))
  
  mus <- list(c(runif(1,-20,20), runif(1,-20,20)), 
              c(runif(1,-20,20), runif(1,-20,20)))
  
  variances <- list(diag(runif(2,0,30)), diag(runif(2,0,30)))
  
  alpha <- c(0.5,0.5)
  
  inicialization <- list(mu = mus, sig = variances, alpha = alpha, beta = betas)
  return(inicialization)
}

run.em.random <- function(X, y, K, max_iter=500) {
  phi <- In_EM_random(X,K = K, y=y)
  
  for(i in 1:max_iter) {
    oldphi <- phi
    probs <- E.step(X, y, phi, K)
    phi <- M.step(X, y, phi, probs, K)
    if((log.like(X, y, phi, K) - log.like(X, y, oldphi, K)) < 0.01)
      break
  }
  return(phi)
}

output <- matrix(nrow = 200, ncol = 3)
for(i in 1:200){
  cat(i, "\n")
  proposal <- run.em.random(X, y, 2)
  proposal.perf <- clusters_perf(resp, E.step(X, y, proposal, 2))
  output[i,] <- c(i,
                  proposal.perf[1],
                  proposal.perf[2])
}

par(mfrow=c(1,1))
plot(output[,2], pch=16,main='Adjusted Rand Index')

