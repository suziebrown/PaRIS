# test of PaRIS

a <- 0.9
b <- 1
sigma.x <- 1
sigma.y <- 1
x0 <- 0

Y <- GaussDLM(10, a, b, sigma.x, sigma.y, x0)
paris(Y, param=c(x0, sigma.x, sigma.y), N=5, Ntilde=2)
