PF.smooth <- function(N, y, par.init, par.trans, par.em){
  # N: number of particles
  # y=(y_1,...,y_T): observations
  # param: parameters for initial/trasition/emission distributions (in "super Gaussian" case c(mean.x, sigma.x, sigma.y))
  Nobs <- length(y)
  xtilde <- r.init(N, par.init) 
  W1 <- d.em(y[1], xtilde, par.em) 
  W1 <- W1/sum(W1)
  xtilde <- sample(xtilde, N, prob=W1, replace=TRUE)
  Xtilde <- matrix(0, Nobs, N)
  Xtilde[1,] <- xtilde
  W <- matrix(0, Nobs, N)
  W[1,] <- W1
  for(t in 2:Nobs){
    Xtilde[t,] <- r.trans(N, Xtilde[t-1,], par.trans) 
    W[t,] <- d.em(y[t], Xtilde[t,], par.em)
    W[t,] <- W[t,]/sum(W[t,])
    index <- sample(1:N, N, prob=W[t,], replace=TRUE)
    Xtilde <- Xtilde[,index] 
  }
  return(list(ksi=Xtilde, omega=W))
}
# out <- PF.smooth(10, c(1,3,2,3,3,7,10,5,6,11), c(1,10,4))
# out

PF.filter <- function(N, y, par.init, par.trans, par.em){
  # N: number of particles
  # y=(y_1,...,y_Nobs): observations
  # param: parameters for initial/trasition/emission distributions (in "super Gaussian" case c(mean.x, sigma.x, sigma.y))
  Nobs <- length(y)
  xtilde <- r.init(N, par.init)
  W1 <- d.em(y[1], xtilde, par.em)
  W1 <- W1/sum(W1)
  Xtilde <- matrix(0, Nobs, N)
  Xtilde[1,] <- xtilde
  W <- matrix(0, Nobs, N)
  W[1,] <- W1
  for(t in 1:(Nobs-1)){
    index <- sample(1:N, N, replace=TRUE, prob=W[t,])
    Xtilde[t+1,] <- r.trans(N, Xtilde[t,index], par.trans)
    W[t+1,] <- d.em(y[t+1], Xtilde[t+1,], par.em)
    W[t+1,] <- W[t+1,]/sum(W[t+1,])
  }
  return(list(ksi=Xtilde[Nobs,], omega=W[Nobs,]))
}
# out1<-PF.filter(10,c(1,3,2,3,3,7,10,5,6,11),c(0,1),c(0,1),c(0,1))
# out1

PF.iter <- function(ksi, omega, yt, par.trans, par.em){
  # ksi=(ksi_t^1,...,ksi_t^N): states at time t of each particle
  # omega=(omega_t^1,...,omega_t^N): weights at time t for each particle
  # yt=y_t: observation at time t
  # param: parameters for initial/trasition/emission distributions (in "super Gaussian" case c(mean.x, sigma.x, sigma.y))
  N <- length(ksi)
  index <- sample(1:N, N, replace=TRUE, prob=omega)
  Xtilde <- r.trans(N, ksi, par.trans)
  W <- d.em(yt, Xtilde, par.em)
  W <- W/sum(W)
  return(list(ksi=Xtilde, omega=W))
}
# ksi=out1$ksi
# omega=out1$omega
# out2 <- PF.iter(10, ksi, omega, 1, c(0,1), c(0,1))
# out2

