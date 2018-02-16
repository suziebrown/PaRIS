# skeleton code for increment of PaRIS algorithm
# Requires functions: PF, d.transition, htilde
# d.transition and htilde must be compatible with vectorised first argument.

paris.iter <- function(ksi, omega, tau, y, param, Ntilde){
  # check inputs okay:
  check.dim<-c(length(ksi),length(omega),length(tau))
  if(all(check.dim==check.dim[1])==FALSE) stop("The length of ksi, omega and tau must be the same!")
  # get number of particles:
  N <- length(ksi)
  tau_new <- numeric(N)
  
  # standard particle update:
  PFout <- PF.iter(N, ksi, omega, y, param)
  ksi_new <- PFout$ksi
  omega_new <- PFout$omega
  
  for (i in 1:N){
    # choose indices:
    J <- sample(1:N, Ntilde, replace=F, prob = omega*d.transition(ksi, ksi_new[i])) 
    # update particle estimate of h:
    tau_new[i] <- (sum(tau[J]) + sum(htilde(ksi[J], ksi_new[i]))) / Ntilde 
  }
  
  # return:
  list(ksi=ksi_new, omega=omega_new, tau=tau_new)
}


#~~~~~~~~~~~~
# wrapper for all iterations if using offline


paris <- function(y, param, N, Ntilde){
  # number of time steps:
  Nobs <- length(y)
  ksi <- matrix(NA, Nobs, N)
  omega <- matrix(NA, Nobs, N)
  tau <- matrix(NA, Nobs, N)
  
  # intialisation:
  ksi[1,] <- r.initial_distr(N, param[1], param[2])
  omega[1,] <- d.emission(y[1], ksi[1,], param[3])
  tau[1,] <- rep(0,N)
  
  # iterate the ol' Paris:
  for (i in 1:(Nobs-1)){
    parisout <- paris.iter(ksi[i,], omega[i,], tau[i,], y[i], param, Ntilde)
    ksi[i+1,] <- parisout$ksi
    omega[i+1,] <- parisout$omega
    tau[i+1,] <- parisout$tau
  }
  
  # return:
  list(ksi=ksi, omega=omega, tau=tau)
}



