
N.paris<-25
Ntilde<-2
N.FFBSm<-50
par.init<-c(0,1)
par.trans<-c(0.7,0.2)
par.em<-c(1,1)

## generate a y
Nobs=500
y<-obs.sim(Nobs, par.init, par.trans, par.em)



## increment for "h1" in paper (sort of)
htilde <- function(x, x_new){
  x_new^2
}

###
h1.FFBSm<-matrix(0,nrow=50,ncol=5)
colnames(h1.FFBSm)<-c("100","200","300","400","500")
h1.paris<-matrix(0,nrow=50,ncol=5)
colnames(h1.paris)<-c("100","200","300","400","500")


for(j in 1:50){
  for (k in 1:5){
    tau.FFBSm<-forwardFFBSm(y[1:(100*k)], N.FFBSm, par.init, par.trans, par.em)$tau
    h1.FFBSm[j,k]<-mean(tau.FFBSm[dim(tau.FFBSm)[1],])
    tau.paris<-paris(y[1:(100*k)], N.paris, Ntilde, par.init, par.trans, par.em)$tau
    h1.paris[j,k]<-mean(tau.paris[dim(tau.paris)[1],])
  }
}

library(KFAS)
Zt <- par.em[1]
Ht <- par.em[2]^2
Tt <- par.trans[1]
Rt <- 1
Qt <- par.trans[2]^2
a1 <- par.init[1]
P1 <- par.init[2]^2
P1inf <- 0
model_y<- SSModel(y ~ -1 +SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),H = Ht)
model_y <- fitSSM(model_y, c(log(var(y)), log(var(y))),method = "BFGS")$model
distsmooth_y <- KFS(model_y, smoothing="disturbance")

par(mfrow=c(1,2))

#set.ylim1<-range(c(t(h1.FFBSm)/((1:5)*100)),c(t(h1.paris)/((1:5)*100)))*c(1.1,1.1)
set.ylim1 <- c(0, 0.05)
boxplot(t(t(h1.FFBSm)/((1:5)*100)),xlab="t",ylab=expression(paste("estimate of h1"[t],"/t")),cex.lab=1.2, main="Forward-only FFBSm",ylim=set.ylim1)
points(cumsum(distsmooth_y$att^2)[(1:5)*100]/((1:5)*100),pch=8,col=4)
boxplot(t(t(h1.paris)/((1:5)*100)),xlab="t",ylab=expression(paste("estimate of h1"[t],"/t")),cex.lab=1.2, main="PaRIS",ylim=set.ylim1)
points(cumsum(distsmooth_y$att^2)[(1:5)*100]/((1:5)*100),pch=8,col=4)
