N.paris<-50
N.FFBSm <- 50
Ntilde<-2
#N.FFBSm<-30
par.init<-c(0,1)
par.trans<-c(0.7,0.2)
par.em<-c(1,1)

n.rep <- 20
Nobs<-100*(1:20)

y<-obs.sim(2000, par.init, par.trans, par.em)

htilde <- function(x, x_new){
  rep(x_new, length(x))
}

time.paris <- matrix(NA,nrow=length(Nobs),ncol=n.rep)
time.FFBSm <- matrix(NA,nrow=length(Nobs),ncol=n.rep)
rownames(time.paris) <- as.character(Nobs)
rownames(time.FFBSm) <- as.character(Nobs)

for (i in 1:length(Nobs)){
  for (j in 1:n.rep){
    time.paris[i,j] <- system.time(paris(y[1:(Nobs[i])], N.paris, Ntilde, par.init, par.trans, par.em))[3]
    time.FFBSm[i,j] <- system.time(forwardFFBSm(y[1:(Nobs[i])], N.FFBSm, par.init, par.trans, par.em))[3]
  }
}

T.paris <- apply(time.paris,1,mean)
T.FFBSm <- apply(time.FFBSm,1,mean)

par(mfrow=c(1,1))
plot(Nobs, T.paris, ylab="computation time", pch=16)
points(Nobs, T.FFBSm, col=2, pch=16)
legend("topleft", c("PaRIS", "FFBSm"), pch=16, col=c(1,2))