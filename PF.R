r.transition <- function(N, media, sigma){ # choose any f transition dist
  rnorm(N, mu, sigma)
}
#
d.emission <- function(y, media, sigma){ # choose g(y|x) emission dist
  dnorm(y, mu, sigma)
}

r.initial_distr <- function(N, media, sigma){ # choose any mu initial dist
  rnorm(N, mu, sigma)
}

PF.smooth <- function(N, y, param){
  # N: number of particles
  # y=(y_1,...,y_T): observations
  # param: parameters for initial/trasition/emission distributions
  # in my case c(mean.x, sigma.x, sigma.y)
  T <- length(y)
  xtilde <- r.initial_distr(N, param[1], param[2]) 
  W1 <- d.emission(y[1], xtilde, param[3]) 
  W1 <- W1/sum(W1)
  xtilde <- sample(xtilde, N, prob=W1, replace=TRUE)
  Xtilde <- matrix(0, T, N)
  Xtilde[1,] <- xtilde
  W <- matrix(0, T, N)
  W[1,] <- W1
  for(t in 2:(T-1)){
    Xtilde[t,] <- r.transition(N, Xtilde[t-1,], param[2]) 
    W[t,] <- d.emission(y[t], Xtilde[t,], param[3])
    W[t,] <- W[t,]/sum(W[t,])
    index <- sample(1:N, N, prob=W[t,], replace=TRUE)
    Xtilde <- Xtilde[,index] 
  }
  Xtilde[T,] <- r.transition(N, Xtilde[T-1,], param[2]) 
  W[T,] <- d.emission(y[T], Xtilde[T,], param[3])
  W[T,] <- W[T,]/sum(W[T,])
  return(list(ksi=Xtilde, omega=W))
}
out <- PF(10, c(1,3,2,3,3,7,10,5,6,11), c(1,10,4))
out




PF.filter <- function(N, y, param){
  # N: number of particles
  # y=(y_1,...,y_T): observations
  # param: parameters for initial/trasition/emission distributions
  # in my case c(mean.x, sigma.x, sigma.y)
  T <- length(y)
  xtilde <- r.initial_distr(N, param[1], param[2])
  W1 <- d.emission(y[1], xtilde, param[3])
  W1 <- W1/sum(W1)
  #xtilde <- sample(xtilde, N, prob=W1, replace=TRUE)
  Xtilde <- matrix(0, T, N)
  Xtilde[1,] <- xtilde
  W <- matrix(0, T, N)
  W[1,] <- W1
  for(t in 1:(T-1)){
    index <- sample(1:N, N, replace=TRUE, prob=W[t,])
    Xtilde[t+1,] <- r.transition(N, Xtilde[t,index], param[2])
    W[t+1,] <- d.emission(y[t+1], Xtilde[t+1,], param[3])
    W[t+1,] <- W[t+1,]/sum(W[t+1,])
  }
  return(list(ksi=Xtilde[T,], omega=W[T,]))
}
out1 <- PF.filter(10, c(1,3,2,3,3,7,10,5,6,11), c(1,10,4))
out1


PF.iter <- function(N, ksi, omega, y, param){
  index <- sample(1:N, N, replace=TRUE, prob=omega)
  Xtilde <- r.transition(N, ksi, param[2])
  W <- d.emission(y, Xtilde, param[3])
  W <- W/sum(W)
  return(list(ksi=Xtilde, omega=W))
}
ksi=out1$ksi[10,]
omega=out1$omega[10,]
out2 <- PF.iter(10, ksi, omega, 1, c(1,10,4))
