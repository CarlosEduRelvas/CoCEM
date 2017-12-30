library(MASS)
tempo = 1:30
time <- list()
dfs = list()
cont <- 1
for(i in unique(tempo)){
  dfs[[cont]] = cbind(mvrnorm(20,
                              c(20+0.4*i, 30+0.8*i),
                              matrix(c(1,0,0,1), nrow = 2)), 1)
  cont<-cont+1
  dfs[[cont]] = cbind(mvrnorm(20,
                              c(25+0.03*i^2-0.7*i, 35+0.08*i^2-1.9*i),
                              matrix(c(1,0,0,1), nrow = 2)), 2)
  cont<-cont+1
  time[[i]] <- rep(i, 40)
}
df = do.call(rbind, dfs)
time = unlist(time)
df <- cbind(df, time)

png(paste0('ALL_times', ".png"), width = 720, height = 480)
plot(df[,1], df[,2], col=df[,3], pch=16, 
     xlab='X', ylab='Y', main='Clustering',
     xlim=c(18,34), ylim=c(20,58))
dev.off()

library(mixtools)
ellipse_bvn <- function(bvn, alpha, cor){
  Xbar <- apply(bvn,2,mean)
  S <- matrix(c(1,0,0,1), nrow = 2)
  ellipse(Xbar, S, alpha = alpha, col=cor)
}

for(i in unique(df[,4])){
  png(paste0('Time_',i, ".png"), width = 720, height = 480)
  aux <- df[df[,4]==i,]
  plot(aux[,1], aux[,2], col=aux[,3], pch=16, 
       xlab='X', ylab='Y', main=paste0('Time: ', i),
       xlim=c(18,34), ylim=c(20,58))
  ellipse_bvn(aux[1:20,c(1,2)], 0.2, "black")
  ellipse_bvn(aux[21:40,c(1,2)], 0.2, "red")
  dev.off()
}




N = model.matrix(~as.factor(time)-1)
f1 = list(rep(0,30), rep(0,30))
f2 = list(rep(0,30), rep(0,30))


R = diag(rep(2/3,28))
Q = matrix(0, nr=30, nc=28)
for(i in 1:27){
  R[i+1,i] = 1/6
  R[i,i+1] = 1/6
  Q[i+1,i] = 1
  Q[i,i+1] = 1
}
J = Q%*%solve(R)%*%t(Q)
lambda = 10
I = diag(1, nr=30, nc=30)
PESO = solve(I + lambda*J)


In_EM_a <- function(X, K, Y){
  library(mclust)
  
  p <- ncol(Y)
  f1 = rep(list(c(0,0)),p)
  f2 = rep(list(c(0,0)),p)
  X_lag <- matrix(nr = nrow(X), nc = ncol(X))
  for(i in 1:K){
    lin_mod <- lm(X[,i]~Y- 1)  
    for(j in 1:p){
      f1[[j]][[i]] = lin_mod$coefficients[j]
      f2[[j]][[i]] = lin_mod$coefficients[j]
    }
    X_lag[,i] <- X[,i] - lin_mod$fitted
  }
  
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
  
  f <- list(f1, f2)
  
  inicialization <- list(mu = mus, sig = variances, alpha = alpha, beta = f)
  return(inicialization)
}


In_EM <- function(X, K, Y){
  library(mclust)
  
  p <- ncol(Y)
  f1 = rep(list(c(0,0)),p)
  f2 = rep(list(c(0,0)),p)
  
  fit <- Mclust(X, G=K)
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
  
  f <- list(f1, f2)
  
  inicialization <- list(mu = mus, sig = variances, alpha = alpha, beta = f)
  return(inicialization)
}

## E-Step : Calculating probabilities of each cluster

E.step <- function(X, Y, phi, K) {
  prob <- matrix(nrow=nrow(X), ncol=K)
  for(grupo in 1:K){
    for(i in 1:nrow(X)){
      media <- phi$mu[[grupo]]
      for(p in 1:ncol(Y)){
        media <- media + phi$beta[[grupo]][[p]]*Y[i,p]
      }
      prob[i,grupo] <- dmvnorm(X[i,],
                               media,
                               phi$sig[[grupo]])+10^(-10)
    }
  }
  return(prob/rowSums(prob))
}


## M-Step

xbeta_calculus <- function(X, Y, phi, K){
  x_beta <- list()
  for(i in 1:K){
    x_beta[[i]] <- matrix(0,nrow(X), ncol(X))
  }
  for(i in 1:K){
    for(j in 1:ncol(Y)){
      x_beta[[i]] <- x_beta[[i]] + sweep(
        do.call(cbind, rep(list(Y[,j]),ncol(X))),
        MARGIN=2, 
        phi$beta[[i]][[j]],
        `*`)
    }
  }
  return(x_beta)
}

M.step <- function(X, Y, phi, probs, K) {
  x_beta <- xbeta_calculus(X, Y, phi, K)
  
  for(i in 1:K){
    for(j in 1:ncol(Y)){
      phi$beta[[i]][[j]] <- colSums(
        (sweep(X, 2, phi$mu[[i]])-x_beta[[i]]+sweep(
          do.call(cbind, rep(list(Y[,j]),ncol(X))),
          MARGIN=2, 
          phi$beta[[i]][[j]],
          `*`))*Y[,j]*as.vector(matrix(
            rep(probs[,i],ncol(X)), nc=ncol(X), byrow = F)))/
        sum(Y[,j]^(2)*probs[,i])
      x_beta <- xbeta_calculus(X, Y, phi, K)
    }
    phi$mu[[i]] <- colSums((X-x_beta[[i]])*probs[,i])/colSums(probs)[i]
  }
  
  covs <- lapply(1:K, function(i) cov.wt((X-x_beta[[i]]), probs[,i]))
  phi$sig <- lapply(covs, "[[", "cov")
  
  phi$alpha <- colMeans(probs)
  
  return(phi)
}


log.like <- function(X, Y, phi, K) {
  prob <- matrix(nrow=nrow(X), ncol=K)
  for(grupo in 1:K){
    for(i in 1:nrow(X)){
      media <- phi$mu[[grupo]]
      for(p in 1:ncol(Y)){
        media <- media + phi$beta[[grupo]][[p]]*Y[i,p]
      }
      prob[i,grupo] <- phi$alpha[grupo]*dmvnorm(X[i,],
                                                media,
                                                phi$sig[[grupo]])
    }
  }
  sum(log(rowSums(prob)))
}


run.em <- function(X, y, K, max_iter=500) {
  library(splines)
  Y <- bs(y, degree=3, knots=seq(5,25,by=10), intercept=FALSE)
  phi <- In_EM(X,K = K, Y)
  
  for(i in 1:max_iter) {
    oldphi <- phi
    probs <- E.step(X, Y, phi, K)
    phi <- M.step(X, Y, phi, probs, K)
    cat(sum(unlist(phi$beta)), "\n")
    if((log.like(X, Y, phi, K) - log.like(X, Y, oldphi, K)) < 0.001)
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
  library(mclust)
  library(clusteval)
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


generic_simulations <- function(seed, n=30, replica = 5, P=2, K=2){
  set.seed(seed)
  random_betas <- c(runif(500,-1,0), runif(500,0,1))
  random_betas2 <- c(runif(500,-0.1,0), runif(500,0,0.1))
  random_betas3 <- c(runif(500,-0.01,0), runif(500,0,0.01))
  random_means <- c(-10,-5,0,5,10)
  
  betas <- list()
  betas2 <- list()
  betas3 <- list()
  means <- list()
  sigma <- list()
  sssss <- list()
  for(i in 1:K){
    sssss[[i]] <- runif(1,1,10)
    sigma[[i]] <- diag(rep(sssss[[i]], P))
    betas[[i]] <- sample(random_betas, P, replace = T)
    betas2[[i]] <- sample(random_betas2, P, replace = T)
    betas3[[i]] <- sample(random_betas3, P, replace = T)
    means[[i]] <- sample(random_means, P, replace = T)
  }
  
  
  tempo = 1:n
  time <- list()
  dfs = list()
  cont <- 1
  for(i in unique(tempo)){
    for(j in 1:K){
      dfs[[cont]] = cbind(mvrnorm(replica,
                                  means[[j]]+betas[[j]]*i+
                                    betas2[[j]]*i^2+betas3[[j]]*i^3,
                                  sigma[[j]]), j)
      cont<-cont+1
    }
    time[[i]] <- rep(i, K*replica)
  }
  df = do.call(rbind, dfs)
  time = unlist(time)
  
  output <- list(time, means, sssss, betas, betas2, betas3, df)
  names(output) <- c("time", "means", "s", "betas",
                     "betas2", "betas3", "df") 
  return(output)
}


perf_simulations <- function(size=100, replica=5, K=2, P=2, KSIM=2){
  output <- matrix(nrow = size, ncol = 17)
  for(i in 1:size){
    cat(i, "\n")
    dataset <- generic_simulations(seed=i,P = P, K=KSIM)
    X <- dataset$df[,-ncol(dataset$df)]
    resp <- dataset$df[,ncol(dataset$df)]
    time <- dataset$time
    Y <- bs(time, degree=3, 
            knots=seq(5,25,by=10), intercept=FALSE)
    
    proposal <- run.em(X, time, K)
    current <- run.em.old(X, K)
    alldata <- run.em.old(cbind(X,time), K)
    partial <- In_EM_a(X,K,Y)
    
    proposal.perf <- clusters_perf(resp, E.step(X, Y, proposal, K))
    current.perf <- clusters_perf(resp, E.step.old(X, current, K))
    all.perf <- clusters_perf(resp, E.step.old(cbind(X,time), alldata, K))
    partial.perf <- clusters_perf(resp, E.step(X, Y, partial, K))
    
    distMu.prop <- dist_parameters(proposal$mu, dataset$means)
    distMu.curr <- dist_parameters(current$mu, dataset$means)
    distMu.alll <- dist_parameters(lapply(alldata$mu, 
                                          function(x) x[1:(length(x)-1)]), 
                                   dataset$means)
    distMu.part <- dist_parameters(partial$mu, dataset$means)
    
    output[i,] <- c(mean(time),
                    sd(time),
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
                        'mu_dist_cocem',
                        'mu_dist_cem',
                        'mu_dist_all',
                        'mu_dist_mix')
  return(output)
}

### K = 2, P = 2
res <- perf_simulations(size = 200)

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
