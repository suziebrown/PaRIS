### Reproducing the results from Olsson & Westerborn, Section 4.1

##~~~~~~~~~~~~~
## FIGURE 2

#pdf(width=12, height=16)
#par(mfrow=c(3,2))

N.paris<-150
Ntilde<-2
N.FFBSm<-50
par.init<-c(0,1)
par.trans<-c(0.7,0.2)
par.em<-c(1,1)

## increment for "h1" in paper (sort of)
htilde <- function(x, x_new){
  x_new
}

## generate a y
Nobs=5000
y<-obs.sim(Nobs, par.init, par.trans, par.em)

###
h.FFBSM<-matrix(0,nrow=50,ncol=5)
colnames(h.FFBSM)<-c("1000","2000","3000","4000","5000")
h.paris<-matrix(0,nrow=50,ncol=5)
colnames(h.paris)<-c("1000","2000","3000","4000","5000")


for(j in 1:50){
  for (k in 1:5){
    tau.FFBSm<-forwardFFBSm(y[1:(1000*k)], N.FFBSm, par.init, par.trans, par.em)$tau
    h.FFBSM[j,k]<-mean(tau.FFBSm[dim(tau.FFBSm)[1],])
    tau.paris<-paris(y[1:(1000*k)], N.paris, Ntilde, par.init, par.trans, par.em)$tau
    h.paris[j,k]<-mean(tau.paris[dim(tau.FFBSm)[1],])
  }
}

boxplot(t(t(h.FFBSM)/((1:5)*1000)),xlab="t",ylab="estimate of h1/t", main="FFBSm",ylim=c(-0.015,0.015))
boxplot(t(t(h.paris)/((1:5)*1000)),xlab="t",ylab="estimate of h1/t", main="PaRIS",ylim=c(-0.015,0.015))

##~~~
## increment for "h2" in paper (sort of)
htilde <- function(x, x_new){
  x_new^2
}

###
h.FFBSM<-matrix(0,nrow=50,ncol=5)
colnames(h.FFBSM)<-c("1000","2000","3000","4000","5000")
h.paris<-matrix(0,nrow=50,ncol=5)
colnames(h.paris)<-c("1000","2000","3000","4000","5000")


for(j in 1:50){
  for (k in 1:5){
    tau.FFBSm<-forwardFFBSm(y[1:(1000*k)], N.FFBSm, par.init, par.trans, par.em)$tau
    h.FFBSM[j,k]<-mean(tau.FFBSm[dim(tau.FFBSm)[1],])
    tau.paris<-paris(y[1:(1000*k)], N.paris, Ntilde, par.init, par.trans, par.em)$tau
    h.paris[j,k]<-mean(tau.paris[dim(tau.FFBSm)[1],])
  }
}

boxplot(t(t(h.FFBSM)/((1:5)*1000)),xlab="t",ylab="estimate of h2/t", main="FFBSm",ylim=c(-0.015,0.015))
boxplot(t(t(h.paris)/((1:5)*1000)),xlab="t",ylab="estimate of h2/t", main="PaRIS",ylim=c(-0.015,0.015))

##~~~
## increment for "h3" in paper (sort of)
htilde <- function(x, x_new){
  x*x_new
}

###
h.FFBSM<-matrix(0,nrow=50,ncol=5)
colnames(h.FFBSM)<-c("1000","2000","3000","4000","5000")
h.paris<-matrix(0,nrow=50,ncol=5)
colnames(h.paris)<-c("1000","2000","3000","4000","5000")


for(j in 1:50){
  for (k in 1:5){
    tau.FFBSm<-forwardFFBSm(y[1:(1000*k)], N.FFBSm, par.init, par.trans, par.em)$tau
    h.FFBSM[j,k]<-mean(tau.FFBSm[dim(tau.FFBSm)[1],])
    tau.paris<-paris(y[1:(1000*k)], N.paris, Ntilde, par.init, par.trans, par.em)$tau
    h.paris[j,k]<-mean(tau.paris[dim(tau.FFBSm)[1],])
  }
}

boxplot(t(t(h.FFBSM)/((1:5)*1000)),xlab="t",ylab="estimate of h3/t", main="FFBSm",ylim=c(-0.015,0.015))
boxplot(t(t(h.paris)/((1:5)*1000)),xlab="t",ylab="estimate of h3/t", main="PaRIS",ylim=c(-0.015,0.015))

#dev.off()




##~~~~~~~~~~~~~
## FIGURE 3


##~~~~~~~~~~~~~
## FIGURE 4


##~~~~~~~~~~~~~
## FIGURE 5


##~~~~~~~~~~~~~
## FIGURE 6

