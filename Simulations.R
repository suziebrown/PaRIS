###Reproducing the results from Olsson&Westerborn, Section 4.1

N.paris<-150
Ntilde<-2
N.FFBSm<-50
par.init<-c(0,1)
par.trans<-c(0.7,0.2)
par.em<-c(1,1)

# increment for "h1" in paper
htilde <- function(x, x_new){
  x_new
}

###generate a y
Nobs=5000
y<-obs.sim(Nobs, par.init, par.trans, par.em)

###
h.FFBSM<-matrix(0,nrow=50,ncol=5)
colnames(h.FFBSM)<-c("1000","2000","3000","4000","5000")
h.paris<-matrix(0,nrow=50,ncol=5)
colnames(h.paris)<-c("1000","2000","3000","4000","5000")


for(j in 1:50){
  for (k in 1:5){
    tau.FFBSm<-forwardFFBSm(y[1:(1000*k)], N, par.init, par.trans, par.em)$tau
    h.FFBSM[j,k]<-mean(tau.FFBSm[dim(tau.FFBSm)[1],])
    tau.paris<-paris(y[1:(1000*k)], N, Ntilde, par.init, par.trans, par.em)$tau
    h.paris[j,k]<-mean(tau.paris[dim(tau.FFBSm)[1],])
  }
}

par(mfrow=c(1,2))
boxplot(t(t(h.FFBSM)/((1:5)*1000)),xlab="t",ylab="estimate of h/t", main="FFBSm",ylim=c(-0.015,0.015))
boxplot(t(t(h.paris)/((1:5)*1000)),xlab="t",ylab="estimate of h/t", main="PaRIS",ylim=c(-0.015,0.015))

#~~

X <- numeric(Nobs) 
X[1] <- r.init(1, par.init)
for (t in 1:Nobs){
  X[t+1] <- r.trans(1, X[t], par.trans)
}
ts.plot(X)