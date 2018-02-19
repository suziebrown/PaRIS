## generic function to simulate observations Y

obs.sim <- function(Nobs, par.init, par.trans, par.em){ 
  X <- numeric(Nobs) 
  Y <- numeric(Nobs)
  X[1] <- r.init(1, par.init)
  for (t in 1:Nobs){
    Y[t] <- r.em(1, X[t], par.em)
    X[t+1] <- r.trans(1, X[t], par.trans)
  }
  Y
}


##~~~~~~~~~~
## functions for stupid model

r.trans <- function(N, xt, par.trans){ # transition distribution for x_t|x_{t-1}
  xt
}
d.trans <- function(x1, x2, par.trans){ # forward transition density for X
  ifelse(x1==x2, 1, 0)
}
d.em <- function(yt, xt, par.em){ # emission density for y_t|x_t
  dnorm(yt, par.em[1]*xt, par.em[2])
}
r.em <- function(N, xt, par.em){ # emission distribution for y_t|x_t
  rnorm(N, par.em[1]*xt, par.em[2])
}
r.init <- function(N, par.init){ # initial distribution for x_1
  1
}


#~~~~~~~~~~~~~
## test

N.paris<-25
Ntilde<-2
N.FFBSm<-50
par.init<-NA
par.trans<-NA
par.em<-c(1,1)
Nobs<-500

htilde <- function(x, x_new){
  x_new
}

y<-obs.sim(Nobs, par.init, par.trans, par.em)

h.FFBSM<-matrix(0,nrow=50,ncol=5)
colnames(h.FFBSM)<-c("1000","2000","3000","4000","5000")
h.paris<-matrix(0,nrow=50,ncol=5)
colnames(h.paris)<-c("1000","2000","3000","4000","5000")


for(j in 1:50){
  for (k in 1:5){
    tau.FFBSm<-forwardFFBSm(y[1:(100*k)], N.FFBSm, par.init, par.trans, par.em)$tau
    h.FFBSM[j,k]<-mean(tau.FFBSm[dim(tau.FFBSm)[1],])
    tau.paris<-paris(y[1:(100*k)], N.paris, Ntilde, par.init, par.trans, par.em)$tau
    h.paris[j,k]<-mean(tau.paris[dim(tau.paris)[1],])
  }
}

boxplot(t(t(h.FFBSM)/((1:5)*100)),xlab="t",ylab="estimate of h1/t", main="FFBSm")
boxplot(t(t(h.paris)/((1:5)*100)),xlab="t",ylab="estimate of h1/t", main="PaRIS")



