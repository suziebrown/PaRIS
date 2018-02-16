r.transition <- function(N, mu, sigma){ # choose any transition distribution for x_t|x_{t-1}
  # N number of particles
  # other inputs are parameters of the distribution
  rnorm(N, mu, sigma)
}

d.emission <- function(yt, mu, sigma){ # choose emission distribution for y_t|x_t
  # yt=y_t: observation at time t
  # other inputs are parameters of the distribution
  dnorm(y, mu, sigma)
}

r.initial_distr <- function(N, mu, sigma){ # choose any initial distribution for x_1
  # N number of particles
  # other inputs are parameters of the distribution
  rnorm(N, mu, sigma)
}

PF.smooth <- function(N, y, param){
  # N: number of particles
  # y=(y_1,...,y_T): observations
  # param: parameters for initial/trasition/emission distributions (in "super Gaussian" case c(mean.x, sigma.x, sigma.y))
  Nobs <- length(y)
  xtilde <- r.initial_distr(N, param[1], param[2]) 
  W1 <- d.emission(y[1], xtilde, param[3]) 
  W1 <- W1/sum(W1)
  xtilde <- sample(xtilde, N, prob=W1, replace=TRUE)
  Xtilde <- matrix(0, Nobs, N)
  Xtilde[1,] <- xtilde
  W <- matrix(0, Nobs, N)
  W[1,] <- W1
  for(t in 2:(Nobs-1)){
    Xtilde[t,] <- r.transition(N, Xtilde[t-1,], param[2]) 
    W[t,] <- d.emission(y[t], Xtilde[t,], param[3])
    W[t,] <- W[t,]/sum(W[t,])
    index <- sample(1:N, N, prob=W[t,], replace=TRUE)
    Xtilde <- Xtilde[,index] 
  }
  Xtilde[Nobs,] <- r.transition(N, Xtilde[Nobs-1,], param[2]) 
  W[Nobs,] <- d.emission(y[Nobs], Xtilde[Nobs,], param[3])
  W[Nobs,] <- W[Nobs,]/sum(W[Nobs,])
  return(list(ksi=Xtilde, omega=W))
}
out <- PF.smooth(10, c(1,3,2,3,3,7,10,5,6,11), c(1,10,4))
out

PF.filter <- function(N, y, param){
  # N: number of particles
  # y=(y_1,...,y_Nobs): observations
  # param: parameters for initial/trasition/emission distributions (in "super Gaussian" case c(mean.x, sigma.x, sigma.y))
  Nobs <- length(y)
  xtilde <- r.initial_distr(N, param[1], param[2])
  W1 <- d.emission(y[1], xtilde, param[3])
  W1 <- W1/sum(W1)
  Xtilde <- matrix(0, Nobs, N)
  Xtilde[1,] <- xtilde
  W <- matrix(0, Nobs, N)
  W[1,] <- W1
  for(t in 1:(Nobs-1)){
    index <- sample(1:N, N, replace=TRUE, prob=W[t,])
    Xtilde[t+1,] <- r.transition(N, Xtilde[t,index], param[2])
    W[t+1,] <- d.emission(y[t+1], Xtilde[t+1,], param[3])
    W[t+1,] <- W[t+1,]/sum(W[t+1,])
  }
  return(list(ksi=Xtilde[Nobs,], omega=W[Nobs,]))
}
out1 <- PF.filter(10, c(1,3,2,3,3,7,10,5,6,11), c(1,10,4))
out1

PF.iter <- function(N, ksi, omega, yt, param){
  # N: number of particles
  # ksi=(ksi_t^1,...,ksi_t^N): states at time t of each particle
  # omega=(omega_t^1,...,omega_t^N): weights at time t for each particle
  # yt=y_t: observation at time t
  # param: parameters for initial/trasition/emission distributions (in "super Gaussian" case c(mean.x, sigma.x, sigma.y))
  index <- sample(1:N, N, replace=TRUE, prob=omega)
  Xtilde <- r.transition(N, ksi, param[2])
  W <- d.emission(yt, Xtilde, param[3])
  W <- W/sum(W)
  return(list(ksi=Xtilde, omega=W))
}
ksi=out1$ksi
omega=out1$omega
out2 <- PF.iter(10, ksi, omega, 1, c(1,10,4))
out2

