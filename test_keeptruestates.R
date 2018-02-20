# test retaining the true states for comparison

obs.sim <- function(Nobs, par.init, par.trans, par.em){ 
  X <- numeric(Nobs) 
  Y <- numeric(Nobs)
  X[1] <- r.init(1, par.init)
  for (t in 1:Nobs){
    Y[t] <- r.em(1, X[t], par.em)
    X[t+1] <- r.trans(1, X[t], par.trans)
  }
  list(X,Y)
}

r.trans <- function(N, xt, par.trans){ # transition distribution for x_t|x_{t-1}
  rnorm(N, par.trans[1]*xt, par.trans[2])
}
d.trans <- function(x1, x2, par.trans){ # forward transition density for X
  dnorm(x1, par.trans[1]*x2, par.trans[2])
}
d.em <- function(yt, xt, par.em){ # emission density for y_t|x_t
  dnorm(yt, par.em[1]*xt, par.em[2])
}
r.em <- function(N, xt, par.em){ # emission distribution for y_t|x_t
  rnorm(N, par.em[1]*xt, par.em[2])
}
r.init <- function(N, par.init){ # initial distribution for x_1
  rnorm(N, par.init[1], par.init[2])
}


N.paris<-25
Ntilde<-2
N.FFBSm<-50
par.init<-c(0,1)
par.trans<-c(0.7,0.2)
par.em<-c(1,1)

Nobs=500
XY<-obs.sim(Nobs, par.init, par.trans, par.em)
y <- XY[[2]]

htilde <- function(x, x_new){
  x_new
}

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
Ht <- matrix(NA)
Tt <- par.trans[1]
Rt <- par.trans[2]
Qt <- matrix(NA)
a1 <- 0
P1 <- 1
P1inf <- 0
model_y<- SSModel(y ~ -1 +SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),H = Ht)
model_y <- fitSSM(model_y, c(log(var(y)), log(var(y))),method = "BFGS")$model
distsmooth_y <- KFS(model_y)


par(mfrow=c(1,2))

set.ylim1<-range(c(t(h1.FFBSm))/((1:5)*100)),c(t(h1.paris)/((1:5)*100)))*c(2,2)
#set.ylim1 <- c(-0.03, 0.07)
boxplot(t(t(h1.FFBSm)/((1:5)*100)),xlab="t",ylab=expression(paste("estimate of h1"[t],"/t")),cex.lab=1.2, main="Forward-only FFBSm",ylim=set.ylim1)
points(cumsum(XY[[1]])[(1:5)*100]/((1:5)*100),pch=8,col=2)
points(cumsum(distsmooth_y$alphahat)[(1:5)*100]/((1:5)*100),pch=16,col=4)

boxplot(t(t(h1.paris)/((1:5)*100)),xlab="t",ylab=expression(paste("estimate of h1"[t],"/t")),cex.lab=1.2, main="PaRIS",ylim=set.ylim1)
points(cumsum(XY[[1]])[(1:5)*100]/((1:5)*100),pch=8,col=2)
points(cumsum(distsmooth_y$alphahat)[(1:5)*100]/((1:5)*100),pch=16,col=4)
