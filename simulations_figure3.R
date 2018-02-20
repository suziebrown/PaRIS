##~~~~~~~~~~~~~
## FIGURE 3

N.paris<-500
Ntilde<-2
N.FFBSm<-500
par.init<-c(0,1)
par.trans<-c(0.7,0.2)
par.em<-c(1,1)
Nobs <- 50
n.rep <- 1

y<-obs.sim(Nobs, par.init, par.trans, par.em)

htilde <- function(x, x_new){
  rep(x_new, length(x))
}

var.FFBSm <- matrix(NA,nrow=Nobs-1,ncol=n.rep)
var.paris <- matrix(NA,nrow=Nobs-1,ncol=n.rep)
time.FFBSm <- matrix(NA,nrow=Nobs-1,ncol=n.rep)
time.paris <- matrix(NA,nrow=Nobs-1,ncol=n.rep)

for (i in 2:Nobs){
  for (j in 1:n.rep){
    time.FFBSm[i-1,j] <- system.time(tau.FFBSm<-forwardFFBSm(y[1:i], N.FFBSm, par.init, par.trans, par.em)$tau)[3]
    var.FFBSm[i-1,j] <- var(tau.FFBSm[i,])
    time.paris[i-1,j] <- system.time(tau.paris<-paris(y[1:i], N.paris, Ntilde, par.init, par.trans, par.em)$tau)[3]
    var.paris[i-1,j] <- var(tau.paris[i,])
  }
}

eff.FFBSm <- var.FFBSm^(-1)*time.FFBSm^(-1)
Eff.FFBSm <- apply(eff.FFBSm, 1, mean)
eff.paris <- var.paris^(-1)*time.paris^(-1)
Eff.paris <- apply(eff.paris, 1, mean)

par(mfrow=c(1,1))
plot(Eff.paris, type='l', xlab='t', ylab='efficiency')
lines(Eff.FFBSm, lty=2)
legend("topright", c("PaRIS", "FFBSm"), lty=c(1,2))
