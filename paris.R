## One iteration of PaRIS, going from (ksi(t), omega(t), tau(t)) to (ksi(t+1), omega(t+1), tau(t+1)):

paris.iter <- function(ksi, omega, tau, yt, Ntilde, par.trans, par.em){
  ## ksi, omega, tau passed in from previous iteration
  ## yt the most recent observation (scalar)
  ## param to be passed to PF.iter
  ## Ntilde number of indices to sample in Paris step (precision parameter)
  
  # check inputs okay:
  check.dim<-c(length(ksi),length(omega),length(tau))
  if(all(check.dim==check.dim[1])==FALSE) stop("The length of ksi, omega and tau must be the same!")
  # get number of particles:
  N <- length(ksi)
  tau_new <- numeric(N)
  
  # standard particle update:
  PFout <- PF.iter(ksi, omega, yt, par.trans, par.em)
  ksi_new <- PFout$ksi
  omega_new <- PFout$omega
  
  for (i in 1:N){
    # choose indices:
    J <- sample(1:N, Ntilde, replace=T, prob = omega*d.trans(ksi, ksi_new[i], par.trans)) 
    # update particle estimate of h:
    tau_new[i] <- (sum(tau[J]) + sum(htilde(ksi[J], ksi_new[i]))) / Ntilde 
  }
  
  # return:
  list(ksi=ksi_new, omega=omega_new, tau=tau_new)
}



## Wrapper for full PaRIS algorithm when implemented offline:

paris <- function(y, N, Ntilde, par.init, par.trans, par.em){
  ## y=(y_1, ..., y_Nobs) vector of observations from a HMM
  ## param to be passed to PF.iter
  ## N number of particles
  ## Ntilde number of indices to sample in Paris step (precision parameter)
  
  # number of time steps:
  Nobs <- length(y)
  ksi <- matrix(NA, Nobs, N)
  omega <- matrix(NA, Nobs, N)
  tau <- matrix(NA, Nobs, N)
  
  # intialisation:
  ksi[1,] <- r.init(N, par.init)
  omega[1,] <- d.em(y[1], ksi[1,], par.em)
  omega[1,] <- omega[1,]/sum(omega[1,]) 
  tau[1,] <- rep(0,N)
  
  # iterate the ol' Paris:
  for (i in 1:(Nobs-1)){
    parisout <- paris.iter(ksi[i,], omega[i,], tau[i,], y[i], Ntilde, par.trans, par.em)
    ksi[i+1,] <- parisout$ksi
    omega[i+1,] <- parisout$omega
    tau[i+1,] <- parisout$tau
  }
  
  # return:
  list(ksi=ksi, omega=omega, tau=tau)
}



