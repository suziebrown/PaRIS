# generate some observations of Y:
StochVol <- function(nobs, phi, sigma.x, beta, x0=0){
  X <- numeric(nobs)
  Y <- numeric(nobs)
  X[1] <- rnorm(1,x0, sigma.x)
  for (t in 1:nobs){
    Y[t] <- beta*exp(X[t]/2)*rnorm(1)
    X[t+1] <- rnorm(1,phi*X[t], sigma.x)
  }
  Y
}

#StochVol(2000,0.975,0.16,0.63)

# forward transition density for X:
q <- function(x1,x2,sigma.x){
  dnorm(phi*x1, x2, sigma.x) ###???
}

# increment for "h1" in paper (section 4.2)
htilde2 <- function(x, x_new){
  x_new^2
}

