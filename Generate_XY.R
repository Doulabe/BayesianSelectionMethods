library(mvtnorm)

Generate.X <- function(n,r,m,f) {
  ## n : number of individuals
  ## r : pairwise correlation between SNPs
  ## m : number of SNPs 
  ## f : MAF

  un <- rep(1,m)
  V <- r*un%*%t(un)+(1-r)*diag(un) 
  Z1 <- rmvnorm(n,sigma=V)
  Z2 <- rmvnorm(n,sigma=V)
  U1 <- pnorm(Z1)
  U2 <- pnorm(Z2)
  G1 <- matrix(as.numeric(U1<f),ncol=m)
  G2 <- matrix(as.numeric(U2<f),ncol=m)
  G <- G1+G2
  X <- scale(G)
return(X)
}

Generate.Y <- function(n,X,beta,sigma){(cbind(rep(1),X)%*%beta+rnorm(n,mean=0,sd=sigma))[,1]}

Power_Fun <- function(beta,n,alpha,p=0){
z <- qnorm(1-alpha/2)
1-pnorm(z,mean=sqrt(n)*beta)+pnorm(-z,mean=sqrt(n)*beta)-p
}

