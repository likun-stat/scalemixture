Rcpp::sourceCpp('~/Desktop/Research/Rcpp/integration.cpp')


################################################################################
################################################################################
## Approximates the marginal quantile function by taking values of the
## marginal CDF of X and doing linear interpolation.  If no values of the CDF
## are supplied, it computes n.x of them, for x in (lower, upper).
##
## tau_sqd ........................... The measurement variance standard
##                                     deviation
## delta ............................. Parameters that get passed to dhuser.wadsworth
## mu ................................ The mean of the normal that gets mixed
##                                     This should probably always be zero
##                                     if the model is used as a copula.
## cdf.vals .......................... If the CDF has already been calculated,
##                                     then they can be passed in and re-used.
## x.vals ............................ Design points for the numerical
##                                     approximation.  Probably best to keep
##                                     as NULL and let the function figure out
##                                     good design points.
## n.x ............................... The number of design points to evaluate
##                                     the CDF.
## lower ............................. The best guess as the smallest design
##                                     point that is necessary.  If this is
##                                     too large, it will get fixed in the 
##                                     function.
## upper ............................. The best guess as the largest design
##                                     point that is necessary.  If this is
##                                     too small, it will get fixed in the 
##                                     function.
##
qmixture.me.interp <- function(p, tau_sqd, delta, mu=0, cdf.vals = NULL, x.vals = NULL,
                               n.x=200, lower=5, upper=20) {
  
  if (is.null(x.vals)) {
    x.range <- find_xrange_pmixture_me(min(p), max(p), c(lower, upper), 
                                       tau_sqd, delta, relerr = 1e-10)
    # x.vals <- seq(x.range[1], x.range[2], length=n.x)
    x.vals <- exp(seq(log(x.range[1]), log(x.range[2]), length=n.x))
    cdf.vals <- pmixture_me(x.vals, tau_sqd, delta)
  } else {
    if (is.null(cdf.vals)) {
      cdf.vals <- pmixture_me(x.vals, tau_sqd, delta)
    }
  }
  q.vals <- spline(x=cdf.vals, y=x.vals, xout=p)$y
  return(q.vals)
  
}
##
################################################################################
################################################################################






################################################################################
################################################################################
## The thing that gets integrated dr to result in the marginal CDF of X
##
## 
##
## tau_sqd ........................... The measurement variance standard
##                                     deviation
## delta ............................. Parameters that get passed to dhuser.wadsworth
##
mix.distn.integrand <- function(x, xval, tau_sqd, delta) {
  prod <- rep(0, length(x))

  half_result = (delta/(2*delta-1))*(xval-x)^(-(1-delta)/delta)-((1-delta)/(2*delta-1))*(xval-x)^(-1)
  prod <- dnorm(x, 0.0, sqrt(tau_sqd)) * half_result
  return(prod)
}
##
################################################################################
################################################################################






################################################################################
################################################################################
## Integrates mix.distn.integrand dr to result in the marginal CDF of X
##
## tau_sqd ........................... The measurement variance standard
##                                     deviation
## delta ............................. Parameters that get passed to dhuser.wadsworth
##
pmixture.uni <- function(xval, tau_sqd, delta) {
  integ <- integrate(mix.distn.integrand, -Inf, xval-1,  xval=xval, tau_sqd = tau_sqd, delta = delta, rel.tol = 1e-10)
  return(pnorm(xval-1, 0.0, sqrt(tau_sqd))-integ$value)
}
pmixture <- Vectorize(pmixture.uni, "xval")

##
################################################################################
################################################################################






################################################################################
################################################################################
## Transforms observations from a Gaussian scale mixture to a GPD
##
## x ................................. A vector of observations from a HOT
##                                     scale mixture
## tau_sqd ........................... The measurement variance standard
##                                     deviation
## delta ............................. Parameters that get passed to dhuser.wadsworth
## theta.gpd ......................... A vector, (thresh, scale, shape)
## 
##
scalemix.me.2.gpd <- function(x, tau_sqd, delta, theta.gpd) {
  require(evd)
  
  thresh <- theta.gpd[1]
  scale <- theta.gpd[2]
  shape <- theta.gpd[3]
  
  unifs <- pmixture_me(x, tau_sqd, delta)
  gpds <- qgpd(unifs, loc=thresh, scale=scale, shape=shape)
  
  return(gpds)
}
##
################################################################################
################################################################################







################################################################################
################################################################################
## Transforms observations from a GPD to a Gaussian scale mixture
##
## y ................................. A vector of observations from a GPD
## tau_sqd ........................... The measurement variance standard
##                                     deviation
## delta ............................. Parameters that get passed to dhuser.wadsworth
## theta.gaussian .................... A parameter, mu
## 
##
gpd.me.2.scalemix <- function(y, tau_sqd, delta, theta.gpd) {
  require(evd)
  
  thresh <- theta.gpd[1]
  scale <- theta.gpd[2]
  shape <- theta.gpd[3]
  
  unifs <- pgpd(y, loc=thresh, scale=scale, shape=shape)
  scalemixes <- qmixture.me.interp(unifs, tau_sqd = tau_sqd, delta = delta, n.x=500)
  
  return(scalemixes)
}
##
################################################################################
################################################################################
