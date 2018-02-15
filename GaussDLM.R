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
