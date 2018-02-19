###Test about computational time

Nvec<-c(10,20,50,100,150)
Ntilde<-2
N.FFBSm<-50
par.init<-c(0,1)
par.trans<-c(0.7,0.2)
par.em<-c(1,1)

htilde <- function(x, x_new){
  x_new
}

Nobs=500
y<-obs.sim(Nobs, par.init, par.trans, par.em)

###
time.FFBSM<-numeric(length(Nvec))
time.paris<-numeric(length(Nvec))

for(i in 1:length(Nvec)){
    time.FFBSM[i]<-system.time(forwardFFBSm(y, Nvec[i], par.init, par.trans, par.em))[3]   ###we want the elapsed time
}

for(i in 1:length(Nvec)){
  time.paris[i]<-system.time(paris(y, Nvec[i], Ntilde, par.init, par.trans, par.em))[3]   ###we want the elapsed time
}

rbind(time.FFBSM,time.paris)
