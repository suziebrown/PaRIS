# test of PaRIS

# a <- 0.9
# b <- 1
# sigma.x <- 1
# sigma.y <- 1
# x0 <- 0
# 
# Y <- GaussDLM(10, a, b, sigma.x, sigma.y, x0)

Y <- c(1,3,2,3,3,7,10,5,6,11)
paris(Y, N=5, Ntilde=2, c(0,1), c(0.452,1), c(0.476,1))
forwardFFBSm(Y, N=5, c(0,1), c(0.452,1), c(0.476,1))