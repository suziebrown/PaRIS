### Reproducing the results from Olsson & Westerborn, Section 4.1

##~~~~~~~~~~~~~
## FIGURE 2

N.paris<-25
Ntilde<-2
N.FFBSm<-50
par.init<-c(0,1)
par.trans<-c(0.7,0.2)
par.em<-c(1,1)

## generate a y
Nobs=5001
y<-obs.sim(Nobs, par.init, par.trans, par.em)



## increment for "h1" in paper (sort of)
htilde <- function(x, x_new){
  x_new
}

###
h1.FFSBm<-matrix(0,nrow=50,ncol=5)
colnames(h1.FFSBm)<-c("1000","2000","3000","4000","5000")
h1.paris<-matrix(0,nrow=50,ncol=5)
colnames(h1.paris)<-c("1000","2000","3000","4000","5000")


for(j in 1:50){
  for (k in 1:5){
    tau.FFBSm<-forwardFFBSm(y[1:(1000*k)], N.FFBSm, par.init, par.trans, par.em)$tau
    h1.FFSBm[j,k]<-mean(tau.FFBSm[dim(tau.FFBSm)[1],])
    tau.paris<-paris(y[1:(1000*k)], N.paris, Ntilde, par.init, par.trans, par.em)$tau
    h1.paris[j,k]<-mean(tau.paris[dim(tau.paris)[1],])
  }
}

##~~~
## increment for "h2" in paper (sort of)
htilde <- function(x, x_new){
  x_new^2
}

###
h2.FFSBm<-matrix(0,nrow=50,ncol=5)
colnames(h2.FFSBm)<-c("1000","2000","3000","4000","5000")
h2.paris<-matrix(0,nrow=50,ncol=5)
colnames(h2.paris)<-c("1000","2000","3000","4000","5000")


for(j in 1:50){
  for (k in 1:5){
    tau.FFBSm<-forwardFFBSm(y[1:(1000*k)], N.FFBSm, par.init, par.trans, par.em)$tau
    h2.FFSBm[j,k]<-mean(tau.FFBSm[dim(tau.FFBSm)[1],])
    tau.paris<-paris(y[1:(1000*k)], N.paris, Ntilde, par.init, par.trans, par.em)$tau
    h2.paris[j,k]<-mean(tau.paris[dim(tau.paris)[1],])
  }
}

system.time(for(j in 1:50){
  for (k in 1:5){
    tau.FFBSm<-forwardFFBSm(y[1:(1000*k)], N.FFBSm, par.init, par.trans, par.em)$tau
    h2.FFSBm[j,k]<-mean(tau.FFBSm[dim(tau.FFBSm)[1],])
  }
})

system.time(for(j in 1:50){
  for (k in 1:5){
    tau.paris<-paris(y[1:(1000*k)], N.paris, Ntilde, par.init, par.trans, par.em)$tau
    h2.paris[j,k]<-mean(tau.paris[dim(tau.paris)[1],])
  }
})

##~~~
## increment for "h3" in paper (sort of)
htilde <- function(x, x_new){
  x*x_new
}

h3.FFSBm<-matrix(0,nrow=50,ncol=5)
colnames(h3.FFSBm)<-c("1000","2000","3000","4000","5000")
h3.paris<-matrix(0,nrow=50,ncol=5)
colnames(h3.paris)<-c("1000","2000","3000","4000","5000")


for(j in 1:50){
  for (k in 1:5){
    tau.FFBSm<-forwardFFBSm(y[1:(1000*k)], N.FFBSm, par.init, par.trans, par.em)$tau
    h3.FFSBm[j,k]<-mean(tau.FFBSm[dim(tau.FFBSm)[1],])
    tau.paris<-paris(y[1:(1000*k)], N.paris, Ntilde, par.init, par.trans, par.em)$tau
    h3.paris[j,k]<-mean(tau.paris[dim(tau.paris)[1],])
  }
}

system.time(for(j in 1:50){
  for (k in 1:5){
    tau.paris<-paris(y[1:(1000*k)], N.paris, Ntilde, par.init, par.trans, par.em)$tau
    h3.paris[j,k]<-mean(tau.paris[dim(tau.paris)[1],])
  }
})

system.time(for(j in 1:50){
  for (k in 1:5){
    tau.FFBSm<-forwardFFBSm(y[1:(1000*k)], N.FFBSm, par.init, par.trans, par.em)$tau
    h3.FFSBm[j,k]<-mean(tau.FFBSm[dim(tau.FFBSm)[1],])
  }
})

################### Frappi's h3 with 5001 obs
h3.FFSBm.5001<-matrix(0,nrow=50,ncol=5)
colnames(h3.FFSBm.5001)<-c("1000","2000","3000","4000","5000")
h3.paris.5001<-matrix(0,nrow=50,ncol=5)
colnames(h3.paris.5001)<-c("1000","2000","3000","4000","5000")
for(j in 1:50){
  for (k in 1:5){
    tau.FFBSm.5001<-forwardFFBSm(y[1:((1000*k)+1)], N.FFBSm, par.init, par.trans, par.em)$tau
    h3.FFSBm.5001[j,k]<-mean(tau.FFBSm.5001[dim(tau.FFBSm.5001)[1],])
    tau.paris.5001<-paris(y[1:((1000*k)+1)], N.paris, Ntilde, par.init, par.trans, par.em)$tau
    h3.paris.5001[j,k]<-mean(tau.paris.5001[dim(tau.paris.5001)[1],])
  }
}
model_y.5001<- SSModel(y ~ -1 +SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),H = Ht)
distsmooth_y.5001 <- KFS(model_y.5001,smoothing="disturbance")  
alpha.h3.5001<-distsmooth_y.5001$att*c(0,distsmooth_y.5001$att[-length(distsmooth_y.5001$att)])
set.ylim3<-range(c(t(h3.FFSBm.5001)/((1:5)*1000+1)),c(t(h3.paris.5001)/((1:5)*1000+1)),cumsum(alpha.h3.5001)[(1:5)*1000]/((1:5)*1000+1))*c(0.9,1.1)
boxplot(t(t(h3.FFSBm.5001)/((1:5)*1000+1)),xlab="t",ylab=expression(paste("estimate of h3"[t],"/t")),cex.lab=1.2, ylim=set.ylim3)
points(cumsum(alpha.h3.5001)[(1:5)*1000]/((1:5)*1000+1),pch=8,col=2)
boxplot(t(t(h3.paris)/((1:5)*1000+1)),xlab="t",ylab=expression(paste("estimate of h3"[t],"/t")),cex.lab=1.2, ylim=set.ylim3)
points(cumsum(alpha.h3)[(1:5)*1000-1]/((1:5)*1000),pch=8,col=2)


###Producing the disturbance smoother (exact value of the smoothed sufficient statistic)
###refer to https://cran.r-project.org/web/packages/KFAS/vignettes/KFAS.pdf
library(KFAS)
Zt <- par.em[1]
#Ht <- matrix(NA)
Ht <- 1 
Tt <- par.trans[1]
Rt <- par.trans[2]
#Qt <- matrix(NA)
Qt <- 1
a1 <- 0
P1 <- 1
P1inf <- 0
model_y<- SSModel(y ~ -1 +SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),H = Ht)
#model_y <- fitSSM(model_y, c(log(var(y)), log(var(y))),
#                  method = "BFGS")$model
distsmooth_y <- KFS(model_y)#,smoothing="disturbance")  
points(cumsum(distsmooth_y$att)[(1:5)*1000]/((1:5)*1000),pch=8,col=3)

###Produce the plot

pdf(width=12, height=16)
par(mfrow=c(3,2))

set.ylim1<-range(c(t(h1.FFSBm)/((1:5)*1000)),c(t(h1.paris)/((1:5)*1000)))*c(1.1,1.1)
boxplot(t(t(h1.FFSBm)/((1:5)*1000)),xlab="t",ylab=expression(paste("estimate of h1"[t],"/t")),cex.lab=1.2, main="Forward-only FFBSm",ylim=set.ylim1)
points(cumsum(distsmooth_y$att)[(1:5)*1000]/((1:5)*1000),pch=8,col=2)
boxplot(t(t(h1.paris*2)/((1:5)*1000)),xlab="t",ylab=expression(paste("estimate of h1"[t],"/t")),cex.lab=1.2, main="PaRIS",ylim=set.ylim1)
points(cumsum(distsmooth_y$att)[(1:5)*1000]/((1:5)*1000),pch=8,col=2)

set.ylim2<-range(c(t(h2.FFSBm)/((1:5)*1000)),c(t(h2.paris*2)/((1:5)*1000)),cumsum(distsmooth_y$att^2)[(1:5)*1000]/((1:5)*1000))*c(0.9,1.1)
boxplot(t(t(h2.FFSBm)/((1:5)*1000)),xlab="t",ylab=expression(paste("estimate of h2"[t],"/t")),cex.lab=1.2, ylim=set.ylim2)
points(cumsum(distsmooth_y$att^2)[(1:5)*1000]/((1:5)*1000),pch=8,col=2)
boxplot(t(t(h2.paris*2)/((1:5)*1000)),xlab="t",ylab=expression(paste("estimate of h2"[t],"/t")),cex.lab=1.2, ylim=set.ylim2)
points(cumsum(distsmooth_y$att^2)[(1:5)*1000]/((1:5)*1000),pch=8,col=2)

alpha.h3<-distsmooth_y$alphahat*c(0,distsmooth_y$alphahat[-length(distsmooth_y$alphahat)])
set.ylim3<-range(c(t(h3.FFSBm)/((1:5)*1000)),c(t(h3.paris)/((1:5)*1000)),cumsum(alpha.h3)[(1:5)*1000-1]/((1:5)*1000))*c(0.9,1.1)
boxplot(t(t(h3.FFSBm)/((1:5)*1000)),xlab="t",ylab=expression(paste("estimate of h3"[t],"/t")),cex.lab=1.2, ylim=set.ylim3)
points(cumsum(alpha.h3)[(1:5)*1000-1]/((1:5)*1000),pch=8,col=2)
boxplot(t(t(h3.paris)/((1:5)*1000)),xlab="t",ylab=expression(paste("estimate of h3"[t],"/t")),cex.lab=1.2, ylim=set.ylim3)
points(cumsum(alpha.h3)[(1:5)*1000-1]/((1:5)*1000),pch=8,col=2)

dev.off()

##~~~~~~~~~~~~~
## FIGURE 3


##~~~~~~~~~~~~~
## FIGURE 4


##~~~~~~~~~~~~~
## FIGURE 5


##~~~~~~~~~~~~~
## FIGURE 6

