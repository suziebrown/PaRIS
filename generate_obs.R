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
## functions for Gaussian DLM

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


##~~~~~~~~~~
## functions for Stochastic Volatiliy

r.trans <- function(N, xt, par.trans){ # transition distribution for x_t|x_{t-1}
  rnorm(N, par.trans[1]*xt, par.trans[2])
}
d.trans <- function(x1, x2, par.trans){ # forward transition density for X
  dnorm(x1, par.trans[1]*x2, par.trans[2])
}
d.em <- function(yt, xt, par.em){ # emission density for y_t|x_t
  dnorm(yt, 0, par.em[1]^2 * exp(xt))
}
r.em <- function(N, xt, par.em){ # emission distribution for y_t|x_t
  par.em[1]*exp(xt/2)*rnorm(N)
}
r.init <- function(N, par.init){ # initial distribution for x_1
  rnorm(N, par.init[1], par.init[2])
}


