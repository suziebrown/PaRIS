##~~~~~~~~~~~~~
## FIGURE 4

N.paris<-50
#N.FFBSm<-50
par.init<-c(0,1)
par.trans<-c(0.7,0.2)
par.em<-c(1,1)
Nobs <- 500
n.rep <- 100

y<-obs.sim(Nobs, par.init, par.trans, par.em)[[2]]

htilde <- function(x, x_new){
  x_new
}

var.paris1<-matrix(0,ncol=50,nrow=n.rep)
colnames(var.paris1)<-10*1:50

Ntilde <- 1


for (j in 1:n.rep){
  tau.paris<-paris(y, N.paris, Ntilde, par.init, par.trans, par.em)$tau
  for (i in 1:50){
    var.paris1[j,i] <- var(tau.paris[i*10,])
  }
}


Var1 <- apply(var.paris1,2,mean)

par(mfrow=c(1,1))
plot(Var1)

## in the paper they used the same particles for each value of Ntilde, and only repeated the backward sampling. Our implememntation doesn't easily allow this