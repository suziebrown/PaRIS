###Forward-only implementation of Forward-filtering backward-smoothing 

forwardFFBSm.iter<-function(ksi,omega,tau,yt,par.trans,par.em){
  
  check.dim<-c(length(ksi),length(omega),length(tau))
  if(all(check.dim==check.dim[1])==FALSE) stop("The length of ksi, omega and tau must be the same!")
  
  N<-length(ksi)
  out<-PF.iter(ksi, omega, yt, par.trans, par.em)
  ksi_new<-out$ksi
  omega_new<-out$omega
  # N: number of particles
  # yt: observation at time t
  # param: parameters for initial/trasition/emission distributions
  
  transition<-numeric(N)
  tau_new<-numeric(N)
  
  for (i in 1:N){
    transition<-omega*d.trans(ksi,ksi_new[i], par.trans)
    tau_new[i]<-sum(transition/sum(transition)*(tau+htilde(ksi,ksi_new[i])))
  }
  return(list(ksi=ksi_new, omega=omega_new, tau=tau_new))
}

# forwardFFBSm.old<-function(ksi,omega,tau,y,param){
#   
#   check.dim<-c(length(ksi),length(omega),length(tau))
#   if(all(check.dim==check.dim[1])==FALSE) stop("The length of ksi, omega and tau must be the same!")
#   
#   Nobs<-length(y)
#   smoothing.distrib<-matrix(NA,ncol=length(ksi),nrow=Nobs)
#   colnames(smoothing.distrib)<-paste("Particle",1:length(ksi))
#   rownames(smoothing.distrib)<-paste("Time",1:Nobs)
#   for (i in 1:Nobs){
#     list(ksi, omega, tau)<-forwardFFBSm.iter(ksi,omega,tau,y[i],param)
#     smoothing.distrib[i,]<-sum(omega/sum(omega)*tau)
#   }
#   return(smoothing.distrib)
# }



forwardFFBSm<-function(y, N, par.init, par.trans, par.em){
  #N: number of particles
  
  Nobs <- length(y)
  ksi <- matrix(NA, Nobs, N)
  omega <- matrix(NA, Nobs, N)
  tau <- matrix(NA, Nobs, N)
  smoothing.distrib <- matrix(NA, nrow=Nobs, ncol=N)
  colnames(smoothing.distrib)<-paste("Particle",1:N)
  rownames(smoothing.distrib)<-paste("Time",1:Nobs)
  
  # intialisation:
  ksi[1,] <- r.init(N, par.init)
  omega[1,] <- d.em(y[1], ksi[1,], par.em)
  omega[1,] <- omega[1,]/sum(omega[1,]) 
  tau[1,] <- rep(0,N)
  smoothing.distrib[1,] <- rep(0,N)   ###because tau[1,] is zero
  
  # iterate the forward FFBSm:
  for (i in 1:(Nobs-1)){
    FFBSm.out <- forwardFFBSm.iter(ksi[i,], omega[i,], tau[i,], y[i],par.trans,par.em)
    ksi[i+1,] <- FFBSm.out$ksi
    omega[i+1,] <- FFBSm.out$omega
    tau[i+1,] <- FFBSm.out$tau
    smoothing.distrib[i+1,]<-sum(omega[i+1,]/sum(omega[i+1,])*tau[i+1,])
  }
  
  # return:
  list(ksi=ksi, omega=omega, tau=tau, smoothing.distrib=smoothing.distrib)
}



