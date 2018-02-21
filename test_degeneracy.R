N.paris<-10
Ntilde<-2
#N.FFBSm<-30
par.init<-c(0,1)
par.trans<-c(0.7,0.2)
par.em<-c(1,1)
Nobs <- 20

y<-obs.sim(Nobs, par.init, par.trans, par.em)

htilde <- function(x, x_new){
  rep(x_new, length(x))
}

out <- paris(y, N.paris, Ntilde, par.init, par.trans, par.em)
tau <- out$tau
out2 <- PF.smooth(N.paris, y, par.init, par.trans, par.em)
ksi <- out2$ksi

par(mfrow=c(1,1))

xlim1 <- c(1,Nobs)

ylim2 <- range(ksi)
plot(NA,ylim=ylim2, xlim=xlim1, xlab="t", ylab="particles for X(0:t)")
for (i in 1:Nobs){
  points(rep(i, N.paris), ksi[i,], pch=1, col="purple")
}
for (i in 1:N.paris){
  lines(ksi[,i], col="purple")
}

#ylim1 <- range(tau)
ylim1 <- c(-4,4)
plot(NA, ylim=ylim1, xlim=xlim1, xlab="t", ylab="particles for h(X(0:t))")
for (i in 1:Nobs){
  points(rep(i, N.paris), tau[i,], pch=1, col="purple")
}
for (i in 1:N.paris){
  lines(tau[,i], col="purple")
}


#~~~~~~~~~~~

# par(mfrow=c(1,1))
# # intialisation:
# ksi <- r.init(N.paris, par.init)
# omega <- d.em(y[1], ksi, par.em)
# omega <- omega/sum(omega) 
# tau <- rep(0,N.paris)
# 
# plot(NA, ylim=ylim1, xlim=xlim1, xlab="t", ylab="particles for tau")
# for (i in 1:Nobs){
#   out <- paris.iter(ksi, omega, tau, y[i], Ntilde, par.trans, par.em)
#   tau <- out$tau
#   omega <- out$omega
#   ksi <- out$ksi
#   points(rep(i, N.paris), tau, pch=16, col="purple")
# }
# 
# 
# 
# 





