Rcpp::sourceCpp('~/Desktop/Research/scalemixture/integration.cpp')
Rcpp::sourceCpp('~/Desktop/Research/scalemixture/likelihood.cpp')

## ------------- 1. Marginal CDF ---------------
## Plot the CDF function
delta <- 0.8
tau <- 1   # Always variance = std^2
x_vals <- seq(1.01, 20, length = 10000)
system.time(cdf_vals<-pmixture_me(x_vals, tau, delta))  #0.363
plot(x_vals, cdf_vals, type="l")
grid()
legend("topleft",pch=20,legend=expression(delta==0.8))

## Compare with asymptotic estimate
x_vals <- seq(1.01, 50, length = 20000)
system.time(cdf_vals<-pmixture_me(x_vals, tau, delta))  #0.277
equiv<-function(x)   (delta/(2*delta-1))*x^{(delta-1)/delta} - (1-delta)*x^{-1}/(2*delta-1)
cdf_est<-equiv(x_vals)
plot(x_vals, 1-cdf_vals, type="l", col="grey",lwd=3, ylab="1-F(x)", 
     xlab=expression(x), main=expression(delta==0.8))
lines(x_vals, cdf_est,col="red")
grid(col="grey50")
legend('topright', lwd=c(3,1), col=c("grey","red"),legend = c("1-F(x)","Asymptotic estimate"))





## ------------- 2. Marginal PDF ---------------
## Compute the PDF function
delta <- 0.8
tau <- 1   # Always variance = std^2
x_vals <- seq(1.01, 50, length = 20000)
system.time(cdf_vals<-pmixture_me(x_vals, tau, delta))  #0.363
system.time(pdf_vals<-dmixture_me(x_vals, tau, delta))  #0.183

## Compare with the numerical results from CDF values
library(pspline)
deriv<-predict(sm.spline(x_vals, cdf_vals), x_vals, 1)


## Compare with the asymptotic results 
equiv.pdf<-function(x) {(1-delta)/(2*delta-1)*(x^{-1/delta}-x^{-2})}
pdf_est<-equiv.pdf(x_vals)


## Plot all the results together
layout(matrix(c(1,1,1,1,1,1,2,2), nrow = 4, ncol = 2, byrow = TRUE))
par(mai=c(0.1,0.62,0.82,0.32))
plot(x_vals, pdf_vals, type="l",col="grey30",lwd=3,ylim=c(min(pdf_vals),max(pdf_est)),
     xlab="",ylab="Marginal density", main=expression(delta==0.8), xaxt="n")
grid(col="grey40")
lines(x_vals,deriv,col="red")
lines(x_vals,pdf_est, col="blue")
legend('topright', lwd=c(3,1,1), col=c("grey30","red","blue"),legend = c("f(x)","Numerical derivative","Asymptotic estimate"))

par(mai=c(0.6,0.62,0.1,0.32))
plot(x_vals,pdf_est-pdf_vals,type='n',ylab="Diffrence",xlab=expression(x))
grid(ny=0,col="grey40")
abline(h=0)
lines(x_vals,pdf_est-pdf_vals,lwd=2,col="darkred")

layout(matrix(c(1,1,1,1), nrow = 2, ncol = 2, byrow = TRUE))
par(mai=c(1.02,0.82,0.82,0.42))



## ------------- 3. Find xrange ---------------
delta <- 0.3
tau <- 1   # Always variance = std^2
x_vals <- seq(1.01, 40, length = 20000)
system.time(cdf_vals<-pmixture_me(x_vals, tau, delta))  #0.363
plot(x_vals, cdf_vals, type="l")
grid()
legend("bottomright",pch=20,legend=expression(delta==0.3))

left<-which(cdf_vals<0.61&cdf_vals>0.59)[1]
right<-which(cdf_vals<0.82&cdf_vals>0.81)[1]
lines(c(x_vals[left],x_vals[left]),c(0,cdf_vals[left]),lty=2)
lines(c(x_vals[right],x_vals[right]),c(0,cdf_vals[right]),lty=2)

fd<-find_xrange_pmixture_me(cdf_vals[left], cdf_vals[right], c(20,30), tau, delta)
lines(c(fd[1],fd[1]),c(0,pmixture_me(fd[1], tau, delta)),col="red")
lines(c(fd[2],fd[2]),c(0,pmixture_me(fd[2], tau, delta)), col="red")




## ------------- 4. Quantile function for marginal x ---------------
delta <- 0.3
tau <- 1   # Always variance = std^2
x_vals <- seq(1.01, 40, length = 20000)
system.time(cdf_vals<-pmixture_me(x_vals, tau, delta))  #0.363
plot(x_vals, cdf_vals, type="l")
legend("bottomright",pch=20,legend=expression(delta==0.3))

p.vals<-seq(0.51,0.91,by=0.1)
q.vals<-qmixture.me.interp(p=p.vals, tau_sqd = tau, delta=delta)
abline(h=p.vals,lty=2,col="grey")
abline(v=q.vals,lty=3,col="blue")




## ------------- 5. Marginal transformation ---------------
Y<-scalemix.me.2.gpd(c(20,30),1,0.3,c(11,1,0))
gpd.me.2.scalemix(Y, tau_sqd=1, delta=0.3, theta.gpd=c(11,1,0))




## ------------- 6. Fix issues for large quantile --------------
pmixture_me_uni(4, tau, delta)
pmixture_me_old(4, tau, delta)
pmixture_me(4, tau, delta)



## ------------- 7. Use interpolation for PDF --------------
delta <- 0.3
tau <- 1   # Always variance = std^2
x_vals <- seq(1.01, 50, length = 20000)
system.time(pdf_vals<-dmixture_me(x_vals, tau, delta))  #0.089
system.time(pdf_interp<-dmixture.me(x_vals, tau, delta))  #0.002


## Plot all the results together
plot(x_vals, pdf_vals, type="l",col="grey30",lwd=3,
     xlab="",ylab="Marginal density", main=expression(delta==0.3), xaxt="n")
grid(col="grey40")
lines(x_vals,pdf_interp,col="red")
legend('topright', lwd=c(3,1), col=c("grey30","red"),legend = c("f(x)","Numerical derivative"))



## ------------- 8. Use RcppArmadillo for likelihood / multivariate update --------------
eig2logdet_c(1:4)
eig2logdet(1:4)

library(fields)
S     <- cbind(seq(0, 1, length=200), rep(1, 200))
Cor   <- corr.fn(rdist(S), lambda = 0.5, gamma = 1)
eig.Sigma <- eigen(Cor, symmetric=TRUE)
V <- eig.Sigma$vectors
d <- eig.Sigma$values

eig2inv.quadform.vector(V, 1/d, 1:200)
eig2inv_quadform_vector(V, 1/d, 1:200)

eig2inv.times.vector(V, 1/d, 1:200)
eig2inv_times_vector(V, 1/d, 1:200)

M<-matrix(rnorm(600),ncol=3)
dmvn.eig(M, V, 1/d)
dmvn_eig(M, V, 1/d)




n.s<-200; n.t<-2; tau<-2; delta<-0.6
X <- matrix(NA, n.s, n.t)
X.s <- matrix(NA, n.s, n.t)
C.Cor <- chol(Cor)
u<-rep(NA,n.t)
R<-rep(NA,n.t)
for(t in 1:n.t){
  u[t] <-runif(1,0,1)
  R[t]<-pow(1/(1-u[t]),delta/(1-delta))
}
for(t in 1:n.t) {
  Z.t<-crossprod(C.Cor, rnorm(n.s))
  Z.to.W.s<-1/(1-pnorm(Z.t))
  X.s[ ,t] <- R[t]*Z.to.W.s
  X[ ,t] <- X.s[ ,t] + sqrt(tau)*rnorm(n.s)
}

X_s_likelihood_conditional(X.s[,1], R[1], V, d)
X.s.likelihood.conditional(X.s[,1], R[1], V, d)

X_s_likelihood_conditional_on_X(X.s[,1], X[,1], R[1], V, d, 0.04)
X.s.likelihood.conditional.on.X(X.s[,1], X[,1], R[1], V, d, 0.04)

system.time(X_s<-var_at_a_time_update_X_s(X[,1], R[1], V, d, 0.04,n_chain = 300)) #3.352  6.562  9.819 
system.time(X_s<-var.at.a.time.update.X.s(X[,1], R[1], V, d, 0.04,n.chain = 300)) #4.431  8.627  12.811
