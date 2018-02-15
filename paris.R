paris <- function(ksi, omega, tau, Ntilde){
  # skeleton code for PaRIS algorithm
  # Requires functions: PF, q, htilde
  N <- length(ksi)
  c(ksi_new, omega_new) <- PF(ksi, omega) #invoke the vanilla SMC update
  for (i in 1:N){
    for (j in 1:Ntilde){
      J <- sample(1:N, Ntilde, replace=F, prob = omega*q(ksi, ksi_new))
    }
    tau_new <- Ntilde^(-1)*(sum(tau[J])+ sum(htilde(ksi[J], ksi_new)))
  }
}