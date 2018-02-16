# skeleton code for increment of PaRIS algorithm
# Requires functions: PF, q, htilde. q and htilde must be compatible with vectorised first argument.

paris.iter <- function(ksi, omega, tau, Ntilde){
  N <- length(ksi) #number of particles (could be passed in as argument)
  list(ksi_new, omega_new) <- PF.iter(N, ksi, omega) #invoke the vanilla SMC update
  for (i in 1:N){
    J <- sample(1:N, Ntilde, replace=F, prob = omega*q(ksi, ksi_new[i])) #choose indices
    tau_new[i] <- (sum(tau[J]) + sum(htilde(ksi[J], ksi_new[i]))) / Ntilde #update particle estimate of h
  }
  list(ksi=ksi_new, omega=omega_new, tau=tau_new)
}



# wrapper for all iterations if using offline
paris <- function(){
  
}