library(MASS)
idade = sort(rep(seq(30,40,by=5), 100))
betas1 = c(0.05, 0.1)
betas2 = c(-0.1, -0.1)
sigmas = matrix(c(0.2, 0.01 ,0.1, 0.04), nr=2)
gammas = matrix(c(0.0001, 0.0001, 0.0001, 0.0001), nr=2)
dfs = list()
cont <- 1
for(i in unique(idade)){
  var <- exp(sigmas+i*gammas)
  dfs[[cont]] = cbind(mvrnorm(50,
                              c(20+betas1[1]*i, 25+betas1[2]*i),
                              var), 1, i)
  cont<-cont+1
  dfs[[cont]] = cbind(mvrnorm(50,
                              c(30+betas2[1]*i, 35+betas2[2]*i),
                              var), 2, i)
  cont<-cont+1
}
df = do.call(rbind, dfs)
plot(df[,1], df[,2], col=df[,3], pch=16, 
     xlab='X', ylab='Y', main='Clustering')

for(i in unique(idade)){
  plot(df[idade==i,1], df[idade==i,2], col=df[,3], pch=16, 
       xlab='X', ylab='Y', main='Clustering')
}


X <- df[,c(1,2)]
y <- idade
K=2

x_test <- c(matrix(c(0.02, 0.0001,0.0001,0.04), nr=2),
            matrix(c(0.05, 0.0001, 0.0001, 0.05), nr=2))
x_test <- x_test[c(1,3,5,7,2,4,6,8)]

define_function(x=x_test, list(K=2,
                               y=y,
                               X=X,
                               media = c(20,25), 
                               beta = c(0.05,0.1),
                               probs[,1]))

multiroot(define_function, start = x_test, 
          parms = list(K=2,
                       y=y,
                       X=X,
                       media = c(20,25), 
                       beta = c(0.05,0.1),
                       probs[,1]))$root

nleqslv(fn = define_function, x = x_test, 
        parametros =list(K=2,
                         y=y,
                         X=X,
                         media = c(20,25), 
                         beta = c(0.05,0.1),
                         probs[,1]), control=list(allowSingular=T))



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
    variances[[i]] <- matrix(c(log(vars[ini:end]+0.1), rep(0, K^2)), nrow=2*K, byrow = T)
    ini <- end+1
    end <- end+K*K
  }
  
  inicialization <- list(mu = mus, sig = variances, alpha = alpha, beta = betas)
  return(inicialization)
}


## E-Step : Calculating probabilities of each cluster


E.step <- function(X, y, phi, K) {
  prob <- matrix(nrow=nrow(X), ncol=K)
  for(grupo in 1:K){
    for(i in 1:nrow(X)){
      media <- phi$mu[[grupo]]+phi$beta[[grupo]]*y[i]
      vars <- create_E(ncol(X), y[i], phi$sig[[grupo]])
      prob[i,grupo] <- dmvnorm(X[i,],
                               media,
                               vars)+10^(-10)
    }
  }
  return(prob/rowSums(prob))
}


create_E <- function(P, yi, x){
  E <- matrix(nr=P, nc=P)
  for(i in 1:P){
    for(j in i:P){
      E[i,j] <- exp(x[i,j]+x[i+P,j]*yi)
      E[j,i] <- E[i,j]
    }
  }
  return(E)
}

define_function <- function(x, parametros){
  library(MASS)
  
  K = parametros[[1]]
  y = parametros[[2]] 
  X = parametros[[3]]
  media = parametros[[4]]
  beta = parametros[[5]]
  prob = parametros[[6]]
  
  P <- ncol(X)
  f <- matrix(0,nr=2*P, nc=P)
  
  mult2 = matrix(2, nr=P, nc=P)
  for(i in 1:K){
    mult2[i,i] <- 1
  }
  
  for(i in 1:length(y)){
    D = X[i,] - media - beta*y[i]
    C = D%*%t(D)
    E = create_E(ncol(X), y[i], matrix(x, nc=P, byrow = F))
    
    invE <- ginv(E)
    Eequal <- E*mult2
    I <- -invE%*%C%*%invE*Eequal
    II <- Eequal*t(invE)
    
    f[1:K,1:K] = f[1:K,1:K]-0.5*prob[i]*(I+II)
    
    f[(K+1):(2*K),1:K] = f[(K+1):(2*K),1:K]-0.5*prob[i]*y[i]*(I+II)
  }
  return(f)
}


estimate_VAR <- function(X, y, phi, probs, K){
  vars <- list()
  for(l in 1:K){
    media <- phi$mu[[l]]
    beta <- phi$beta[[l]]
    # vars[[l]] <- matrix(nleqslv(fn = define_function, x = phi$sig[[l]],
    #                             parametros = list(K, y, X, media, beta, probs[,l]),
    #                             control=list(allowSingular=T))$x,
    #                     nc = ncol(X), byrow = F)
    vars[[l]] <- matrix(multiroot(define_function, start = phi$sig[[l]],
                                  parms = list(K, y, X, media, beta, probs[,l]))$root,
                        nc = ncol(X), byrow = F)
  }
  return(vars)
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
  
  phi$beta <- betas
  phi$mu <- mus
  phi$sig <- estimate_VAR(X, y, phi, probs, K)
  
  phi$alpha <- colMeans(probs)
  
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
  library(mvtnorm)
  for(i in 1:max_iter) {
    oldphi <- phi
    probs <- E.step(X, y, phi, K)
    phi <- M.step(X, y, phi, probs, K)
    #if((log.like(X, y, phi, K) - log.like(X, y, oldphi, K)) < 0.01)
    #  break
  }
  return(phi)
}

