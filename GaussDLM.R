# generate some observations of Y:
GaussDLM <- function(nobs, a, b, sigma.x, sigma.y){
  X <- numeric(nobs)
  Y <- numeric(nobs)
  X[1] <- rnorm(1,0, sigma.x)
  for (t in 1:nobs){
    Y[t] <- rnorm(1,b*X[t], sigma.y)
    X[t+1] <- rnorm(1,a*X[t], sigma.x)
  }
  Y
}

# forward transition density for X:
q <- function(x1,x2){
  dnorm(a*x1, x2, sigma.x) #first and second arguments swapped to less obvious order for vectorisation
}