r.transition <- function(N, media, sigma){ # choose any f transition dist
  rnorm(N, media, sigma)
}
#
d.emission <- function(y, media, sigma){ # choose g(y|x) emission dist
  dnorm(y, media, sigma)
}

r.initial_distr <- function(N, media, sigma){ # choose any mu initial dist
  rnorm(N, media, sigma)
}

PF <- function(N, y, param){
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
  for(t in 2:T){
    Xtilde[t,] <- r.transition(N, Xtilde[t-1,], param[2]) 
    W[t,] <- d.emission(y[t], Xtilde[t,], param[3])
    W[t,] <- W[t,]/sum(W[t,])
    index <- sample(1:N, N, prob=W[t,], replace=TRUE)
    Xtilde <- Xtilde[,index] 
  }
  return(list(ksi=Xtilde[T,], omega=W[T,]))
}
PF(10, c(1,3,2,3,3,7,10,5,6,11), c(1,10,4))





# PF <- function(N, y, param){
#   # N: number of particles
#   # y=(y_1,...,y_T): observations
#   # param: parameters for initial/trasition/emission distributions
#   # in my case c(mean.x, sigma.x, sigma.y)
#   T <- length(y)
#   xtilde <- r.initial_distr(N, param[1], param[2]) 
#   W1 <- d.emission(y[1], xtilde, param[3]) 
#   W1 <- W1/sum(W1)
#   xtilde <- sample(xtilde, N, prob=W1, replace=TRUE)
#   Xtilde <- matrix(0, T, N)
#   Xtilde[1,] <- xtilde
#   W <- matrix(0, T, N)
#   W[1,] <- W1
#   for(t in 1:(T-1)){
#     index <- sample(1:N, N, replace=TRUE, prob=W[t,])
#     Xtilde[t+1,] <- r.transition(N, Xtilde[t,index], param[2])
#     W[t+1,] <- d.emission(y[t+1], Xtilde[t+1], param[3])
#     W[t+1,] <- W[t+1,]/sum(W[t+1,])
#   }
#   return(list(ksi=Xtilde[T,], omega=W[T,]))
# }


