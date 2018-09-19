Rcpp::sourceCpp('~/Desktop/Research/Rcpp/integration.cpp')

## ------------- 1. Marginal CDF ---------------
## Plot the CDF function
delta <- 0.3
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
     xlab=expression(x), main=expression(delta==0.3))
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

