source("~/Desktop/Research/scalemixture/scalemix_utils.R")
source("~/Desktop/Research/scalemixture/scalemix_likelihoods.R")
source("~/Desktop/Research/scalemixture/scalemix_priors.R")
source("~/Desktop/Research/scalemixture/generic_samplers.R")
source("~/Desktop/Research/scalemixture/scalemix_sampler_02.R")

library(fields)   # For rdist

# ------------ 1. Simulation settings -------------
n.s <- 200        # Number of sites
n.t <- 20         # Number of time points
tau <-9      # Nugget SD
delta <- 0.7      # For R
lambda <- 0.5     # Powered exponential range
gamma <-  1       # Powered exponential smoothness
rho <- 0.1

# Threshold for fitting
thresh <- 11


# -------------- 2. Generate fake data -----------------
S     <- cbind(seq(0, 1, length=n.s), seq(0, 1, length=n.s))
# Cor   <- corr.fn(rdist(S), lambda = lambda, gamma = gamma)
Cor   <- corr.fn(rdist(S), rho)
C.Cor <- chol(Cor)

set.seed(3333)
u<-rep(NA,n.t)
R<-rep(NA,n.t)
for(t in 1:n.t){
  u[t] <-runif(1,0,1)
  R[t]<-pow(1/(1-u[t]),delta/(1-delta))
}


X <- matrix(NA, n.s, n.t)
X.s <- matrix(NA, n.s, n.t)
for(t in 1:n.t) {
  Z.t<-crossprod(C.Cor, rnorm(n.s))
  Z.to.W.s<-1/(1-pnorm(Z.t))
  X.s[ ,t] <- R[t]*Z.to.W.s
  X[ ,t] <- X.s[ ,t] + sqrt(tau)*rnorm(n.s)
  
}


prob.below <- 0.8
theta.gpd <- c(thresh, 1, 0)

thresh.X <- qmixture.me.interp(prob.below, tau_sqd = tau, delta = delta) 
sum(X < thresh.X) / length(X)
cen <- X < thresh.X  



## ------------ 3. Marginal transformation -----------------
library(evd)

Y <- X
Y[cen] <- NA
Y[!cen] <- scalemix.me.2.gpd(x = X[!cen], tau_sqd = tau, delta = delta, theta.gpd = theta.gpd, prob.below=prob.below)




## --------------- 4. Running Metropolis -------------------
initial.values <- list(delta = delta, rho=rho, tau=tau, theta.gpd=theta.gpd, prob.below=prob.below, X.s=X.s, R=R)
n.updates <- 4000
thin <- 10
echo.interval <- 50
true.params <- list(delta = delta, rho=rho, tau=tau, theta.gpd=theta.gpd, prob.below=prob.below, X.s=X.s, R=R)

scalemix.sampler.02(Y=Y, S=S, cen=cen, thresh=thresh,
                                initial.values=initial.values,
                                n.updates=n.updates, thin=thin,
                                experiment.name="Huser-wadsworth",
                                echo.interval=echo.interval,
                                sigma.m=NULL, prop.Sigma=NULL, 
                                true.params=true.params, n.pl=1, lower.prob.lim=0.5)
