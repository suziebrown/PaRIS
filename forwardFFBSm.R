###Forward-only implementation of Forward-filtering backward-smoothing 

forwardFFBSm.iter<-function(ksi,omega,tau){
  
  check.dim<-c(length(ksi),length(omega),length(tau))
  if(all(check.dim==check.dim[1])==FALSE) stop("The length of ksi, omega and tau must be the same!")
  
  N<-length(ksi)
  list(ksi_new, omega_new)<-PF(N, ksi, omega, y, param)
  
  transition<-numeric(N)
  tau_new<-numeric(N)
  
  for (i in 1:N){
    transition<-omega*d.transition(ksi,ksi_new[i])
    tau_new[i]<-sum(transition/sum(transition)*(tau+htilde(ksi,ksi_new[i])))
  }
  
  #smoothing.distrib<-sum(omega/sum(omega)*tau_new)
  #return("ParticleSample"=new.particle.sample,"Auxiliary"=tau_new,"SmoothingDistribution"=smoothing.distrib)
  return(list(ksi=ksi_new, omega=omega_new, tau=tau_new))
}

