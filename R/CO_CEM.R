### Simulation

library(MASS)
idade = sort(rep(seq(30,65,by=5), 100))
betas1 = c(0.2, 0.05)
betas2 = c(0.05, -0.05)
dfs = list()
cont <- 1
for(i in unique(idade)){
  dfs[[cont]] = cbind(mvrnorm(50,
                              c(23+betas1[1]*i, 30+betas1[2]*i),
                              matrix(c(1,0,0,1), nrow = 2)), 1, i)
  cont<-cont+1
  dfs[[cont]] = cbind(mvrnorm(50,
                              c(30+betas2[1]*i, 35+betas2[2]*i),
                              matrix(c(1,0,0,1), nrow = 2)), 2, i)
  cont<-cont+1
}
df = do.call(rbind, dfs)
plot(df[,1], df[,2], col=df[,3], pch=16, 
     xlab='X', ylab='Y', main='Clustering')

library(mixtools)
ellipse_bvn <- function(bvn, alpha, cor){
  Xbar <- apply(bvn,2,mean)
  S <- cov(bvn)
  ellipse(Xbar, S, alpha = alpha, col=cor)
}

for(i in unique(df[,4])){
  png(paste0('Idade_',i, ".png"), width = 720, height = 480)
  aux <- df[df[,4]==i,]
  plot(aux[,1], aux[,2], col=aux[,3], pch=16, 
       xlab='X', ylab='Y', main=paste0('Idade: ', i),
       xlim=c(25,40), ylim=c(29,36))
  ellipse_bvn(aux[1:50,c(1,2)], 0.2, "black")
  ellipse_bvn(aux[51:100,c(1,2)], 0.2, "red")
  dev.off()
}


library(rgl)
plot3d(df[,1], df[,2], df[,4], col=df[,3],
       xlab='X', ylab='Y',zlab='Idade')

library(mclust)
fit <- Mclust(df[,c(1,2)], G = 2)
fit$parameters
tab = table(fit$classification, df[,3])
tab
max(tab[1,1]+tab[2,2],
    tab[2,1]+tab[1,2])/nrow(df)


fit <- Mclust(df[,c(1,2,4)], G = 2)
tab = table(fit$classification, df[,3])
tab
max(tab[1,1]+tab[2,2],
    tab[2,1]+tab[1,2])/nrow(df)

fit <- In_EM(df[,c(1,2)], 2, df[,4])
probs <- E.step(df[,c(1,2)], df[,4], fit, 2)
class <- apply(probs, 1, which.max)
tab = table(class, df[,3])
tab
max(tab[1,1]+tab[2,2],
    tab[2,1]+tab[1,2])/nrow(df)


proposal <- run.em(df[,c(1,2)], df[,4], 2, max_iter=500)
probs <- E.step(df[,c(1,2)], df[,4], proposal, 2)
class <- apply(probs, 1, which.max)
tab = table(class, df[,3])
tab
max(tab[1,1]+tab[2,2],
    tab[2,1]+tab[1,2])/nrow(df)




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
                               phi$sig[[grupo]])
    }
  }
  return(prob/rowSums(prob))
}


## M-Step

M.step <- function(X, y, phi, probs, K) {
  mu <- list(
    colSums((X-unlist(phi['beta'][[1]][[1]])*matrix(rep(y, K), nrow=nrow(X), ncol=K))*probs[,1])/
      colSums(probs)[1], 
    colSums((X-unlist(phi['beta'][[1]][[2]])*matrix(rep(y, K), nrow=nrow(X), ncol=K))*probs[,2])/
      colSums(probs)[2]        
              )

  # beta <- list(
  #   colSums((X-matrix(phi['mu'][[1]][[1]], nrow=nrow(X), ncol=K, byrow = T))*probs[,1])/
  #     sum(y*probs[,1]), 
  #   colSums((X-matrix(phi['mu'][[1]][[2]], nrow=nrow(X), ncol=K, byrow = T))*probs[,2])/
  #     sum(y*probs[,2]))
  
  beta <- list(c(0,0), c(0,0))
  
#   SUM1 <- matrix(0,nrow=K, ncol=K)
#   SUM2 <- matrix(0,nrow=K, ncol=K)
#   covariate <- matrix(rep(y, K), nrow=nrow(X), ncol=K)
#   for(i in 1:nrow(X)){
#     N1 <- (X[i,]-unlist(phi['beta'][[1]][[1]])*covariate[i,]-
#              matrix(phi['mu'][[1]][[1]], nrow=1, ncol=K, byrow = T))
#     N2 <- (X[i,]-unlist(phi['beta'][[1]][[2]])*covariate[i,]-
#              matrix(phi['mu'][[1]][[2]], nrow=1, ncol=K, byrow = T))
#     N1 <- t(N1)%*%N1
#     N2 <- t(N2)%*%N2
#     SUM1 <- SUM1+probs[i,1]*N1
#     SUM2 <- SUM2+probs[i,2]*N2
#   }
#   variances <- list(SUM1/colSums(probs)[1], SUM2/colSums(probs)[2])
  
  covariate <- matrix(rep(y, K), nrow=nrow(X), ncol=K)
  covs <- lapply(1:K, function(i) cov.wt(X-phi['beta'][[1]][[i]]*covariate, probs[,i]))
#  mu <- lapply(covs, "[[", "center")
  variances <- lapply(covs, "[[", "cov")
  
  alpha <- colMeans(probs)
  #alpha <- c(0.5,0.5)
  
  return(list(mu = mu, sig = variances, alpha = alpha, beta = beta))
}


## Log - like

log.like <- function(X, y, phi, K) {
  prob <- matrix(nrow=nrow(X), ncol=K)
  for(grupo in 1:K){
    for(i in 1:nrow(X)){
      prob[i,grupo] <- phi$alpha[grupo]*dmvnorm(X[i,],
                               phi$mu[[grupo]]+phi$beta[[grupo]]*idade[i],
                               phi$sig[[grupo]])
    }
  }
  sum(log(rowSums(probs)))
}


## Run EM

run.em <- function(X, y, K, max.iter = 30) {
  phi <- In_EM(X,K)
  
  for(i in 1:100) {
    oldphi <- phi
    probs <- E.step(X, y, phi, K)
    phi <- M.step(X, y, phi, probs, K)
    #if((log.like(X, y, phi, K) - log.like(X, y, oldphi, K)) < 0.01)
    #  break
  }
  return(list(phi = phi, aic = 2*3*N - log.like(X, y, phi, N)))
}



library(MASS)
idade = sort(rep(seq(20,35,by=1), 10))
betas1 = c(0.2, 0.2)
betas2 = c(0.01, 0.01)
dfs = list()
cont <- 1
for(i in unique(idade)){
  dfs[[cont]] = cbind(mvrnorm(5,
                              c(20+betas1[1]*i, 30+betas1[2]*i),
                              matrix(c(1,0,0,1), nrow = 2)), 1)
  cont<-cont+1
  dfs[[cont]] = cbind(mvrnorm(5,
                              c(25+betas2[1]*i, 35+betas2[2]*i),
                              matrix(c(1,0,0,1), nrow = 2)), 2)
  cont<-cont+1
}
df = do.call(rbind, dfs)


X <- df[,-3]
y <- idade

learning_rate = 1
phi <- In_EM(X,K = K, y=y)
for(i in 1:100) {
  oldphi <- phi
  probs <- E.step(X, y, phi, K)
  probs <- t(apply(probs, 1, function(x) 0.5 + (x - 0.5) * (1-learning_rate) / (1+learning_rate)))
  probs <- probs/rowSums(probs)
  phi <- M.step(X, y, phi, probs, K)
  learning_rate = max(learning_rate - 0.01,0)
  #if((log.like(X, y, phi, K) - log.like(X, y, oldphi, K)) < 0.01)
}

phi.old <- In_EM(X,K,y)
for(i in 1:100){
  oldphi.old <- phi.old
  probs.old <- E.step.old(X, oldphi.old, K)
  phi.old <- M.step.old(X, probs.old, K)
}


pair <- comembership_table(df[,3], apply(probs.old, 1, which.max))
c((pair[[1]][1]+pair[[4]][1])/sum(unlist(pair)),
  adjustedRandIndex(apply(probs.old, 1, which.max), 
                    df[,ncol(df)]))

pair <- comembership_table(df[,3], apply(probs, 1, which.max))
c((pair[[1]][1]+pair[[4]][1])/sum(unlist(pair)),
  adjustedRandIndex(apply(probs, 1, which.max), 
                    df[,ncol(df)]))


phi_init <- In_EM(X,K = K, y=y)
probs_init <- E.step(X, y, phi_init, K)

pair <- comembership_table(df[,3], apply(probs_init, 1, which.max))
c((pair[[1]][1]+pair[[4]][1])/sum(unlist(pair)),
  adjustedRandIndex(apply(probs_init, 1, which.max), 
                    df[,ncol(df)]))



phi$alpha
phi.old$alpha








## EM algorithm implementation
## Joshua Moller-Mara
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

## log.likelihood
## Given points and parameters, how well does our model fit?
## This also gives us a terminating condition for the EM algorithm.
## We stop if we don't improve all that much.
## This is also used for AIC for model selection (choosing a value of k)
log.like <- function(X, phi, N) {
  probs <- 
    with(phi, do.call(cbind,
                      lapply(1:N, function(i)
                        alpha[i] * dmvnorm(X, mu[[i]], sig[[i]]))))
  sum(log(rowSums(probs)))
}

covs <- replicate(N, list(cov.wt(X[sample(nrow(X), 30),])))
mu <- lapply(covs, "[[", "center")
sig <- lapply(covs, "[[", "cov")
alpha <- rep(1/N, N)
phi <<- list(mu = mu, sig = sig, alpha = alpha, beta = list(c(0,0), c(0,0)))


run.em <- function(X, y, K, max.iter = 30) {
  for(i in 1:100) {
    oldphi <- phi
    probs <- E.step(X, phi, K)
    phi <- M.step(X, probs, K)
    #if((log.like(X, y, phi, K) - log.like(X, y, oldphi, K)) < 0.01)
    #  break
  }
  return(list(phi = phi, aic = 2*3*N - log.like(X, y, phi, N)))
}


teste<-function(learning_rate, prob){
  return(exp(learning_rate*prob)/(1+exp(learning_rate*prob)))
}

teste(0.000000001, 0.1)
