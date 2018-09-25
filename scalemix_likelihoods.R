################################################################################
## The log likelihood of the data, where the data comes from a HOT scale mixture
## of Gaussians, transformed to GPD
##
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X ................................. a matrix of Y, transformed to HOT 
##                                     scale mixture of Gaussian
## X.s ............................... X, but without the measurement error
## cen 
## prob.below
## theta.gpd
## delta
## tau_sqd
##
marg.transform.data.mixture.me.likelihood <- function(Y, X, X.s, cen, prob.below,
                                                      theta.gpd, delta,
                                                      tau_sqd, thresh.X=NULL) {
  if (!is.matrix(Y)) Y <- matrix(Y, ncol=1)
  if (!is.matrix(X)) X <- matrix(X, ncol=1)
  if (!is.matrix(X.s)) X.s <- matrix(X.s, ncol=1)
  
  ll <- matrix(0, nrow(Y), ncol(Y))
  
  loc <- theta.gpd[1]
  scale <- theta.gpd[2]
  shape <- theta.gpd[3]
  
  if (is.null(thresh.X)) thresh.X <- qmixture.me.interp(prob.below, tau_sqd = tau_sqd, delta = delta)
  
  if (sum(cen) > 0)  {
    ll[cen] <- pnorm(thresh.X, mean=X.s[cen], sd=sqrt(tau_sqd), log.p=TRUE)
  }
  if (sum(!cen) > 0) {
    ll[!cen] <- dnorm(X[!cen], mean=X.s[!cen], sd=sqrt(tau_sqd), log=TRUE) +
      dgpd(Y[!cen], loc=loc, scale=scale, shape=shape, log=TRUE)  -
      dmixture.me(X[!cen], tau_sqd = tau_sqd, delta = delta, log=TRUE) - 0.5*log(tau_sqd)
  }
  
  return(sum(ll)) 
}

#                                                                              #
################################################################################





################################################################################
## The log likelihood of X.s, when it is conditioned on scaling factor R and Σ(λ,γ)
##
## X.s ............................... X, but without the measurement error
## R 
## V ................................. eigen vectors of covariance Σ(λ,γ)
## d ................................. a vector of eigenvalues
##

X.s.likelihood.conditional<-function(X.s, R, V, d){
  X.s.to.Z <- qnorm(1-R/X.s)
  loglik <- -0.5*eig2inv.quadform.vector(V, 1/d, X.s.to.Z)+0.5*sum(X.s.to.Z^2)-2*sum(log(X.s))
  return(loglik)
}


X.s.likelihood.conditional.on.X<-function(X.s, X, R, V, d, tau_sqd){
  if(any(X.s<R)) return(-Inf) else{
    X.s.to.Z <- qnorm(1-R/X.s)
    loglik <- -0.5*eig2inv.quadform.vector(V, 1/d, X.s.to.Z)+0.5*sum(X.s.to.Z^2)-2*sum(log(X.s))
    loglik <- loglik - 0.5*sum((X.s-X)^2)/tau_sqd
    return(loglik)
  }
}

var.at.a.time.update.X.s <- function(X, R, V, d, tau_sqd, v.q=0.5, n.chain=100){
  X <- as.vector(X)
  n.s <- length(X)
  X.s <- X
  X.s[which(X.s<R)] <- R + abs(rnorm(length(which(X.s<R)),sd=0.5)) # X.s has to be larger than R
  accept <- rep(0, n.s)
  
  for(i in 1:n.chain){
    for(iter in 1:n.s){
      X.s.update<-X.s
      X.s.update[iter]<-X.s[iter]+rnorm(1,0,sd=v.q)
      log.num <- X.s.likelihood.conditional.on.X(X.s.update, X=X, R=R, V=V, d=d, tau_sqd=tau_sqd)
      log.denom <- X.s.likelihood.conditional.on.X(X.s, X=X, R=R, V=V, d=d, tau_sqd=tau_sqd)
      
      r <- exp(log.num - log.denom)
      
      if(runif(1) < r){    
        X.s <- X.s.update
        accept[iter]<-accept[iter]+1
      }
    }
  }
  
  return(list(X.s=X.s, accept=accept))
}

#                                                                              #
################################################################################





################################################################################
## For the generic Metropolis sampler
## Samples from the scaled Gaussian process (update the smooth process).
## The mixing distribution comes from from the Huser-wadsworth scale mixing distribution.
## The PROPOSAL will be the conjugate update from the model that treats the
## X process (i.e. X.s, but with measurement error) as the response, ignoring
## the marginal transformation.  Then the Metropolis part either rejects or
## accepts the entire draw.
##
## data............................... a n.t vector of scaling factors
## params............................. a 2-vector of parameters (theta in dhuser.thibaud):
##                                     theta[1] = gamma (from H-O-T (2017))
##                                     theta[2] = beta (from H-O-T(2017))
## Y ................................. a (n.s x n.t) matrix of data that are
##                                     marginally GPD, and conditionally
##                                     independent given X(s)
## X.s ............................... the latent Gaussian process, without the
##                                     measurement error
## cen 
## prob.below
## theta.gpd
## theta.mix
## theta.gaussian
## V ................................. eigen vectors of covariance Σ(λ,γ)
## d ................................. a vector of eigenvalues
##

X.s.update.mixture.me <- function(R, Y, X, X.s, cen, 
                                  prob.below, theta.gpd, delta,
                                  tau_sqd, V, d, v.q=0.5, n.chain=100,
                                  thresh.X=NULL) {
  
  n.s <- nrow(Y)
  n.t <- ncol(Y)
  
  accepted <- rep(FALSE, n.t)  
  
  
  if (is.null(thresh.X)) thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau_sqd, delta = delta)
  
  for(t in 1:n.t){
    # Proposal
    # Treat the X process as the response, ignoring the marginal transformation
    # i.e. prop.X.s ~ X.s | X, R, other.params
    
    metropolis <- var.at.a.time.update.X.s(X=X[,t], R=R[t], V=V, d=d, tau_sqd=tau_sqd, v.q = v.q, n.chain = n.chain)
    prop.X.s <- metropolis$X.s
    
    
    # M-H ratio    
    log.rat <- 
      X.s.likelihood.conditional.on.X(X.s[,t], X=X[,t], R=R[t], V=V, d=d, tau_sqd=tau_sqd)+ # Prop density of current value
      marg.transform.data.mixture.me.likelihood(Y[ ,t], X[ ,t], prop.X.s,               # Likelihood of proposal
                                                cen[ ,t], prob.below,
                                                theta.gpd, delta,
                                                tau_sqd, thresh.X=thresh.X) +
      X.s.likelihood.conditional(prop.X.s, R[t], V, d) -                                # Likelihood of proposal
      X.s.likelihood.conditional.on.X(prop.X.s, X=X[,t], R=R[t], V=V, d=d, tau_sqd=tau_sqd)-# Prop density of proposal
      marg.transform.data.mixture.me.likelihood(Y[ ,t], X[ ,t], X.s[ ,t],               # Likelihood of current value
                                                cen[ ,t], prob.below,
                                                theta.gpd, delta,
                                                tau_sqd, thresh.X=thresh.X) -
      X.s.likelihood.conditional(X.s[ ,t], R[t], V, d)                                  # Likelihood of current value
    
    
    
    
    if (runif(1) < exp(log.rat)) {
      X.s[ ,t] <- prop.X.s
      accepted[t] <- TRUE
    }
  }
  
  return(list(X.s=X.s, accepted=accepted))
}

#                                                                              #
################################################################################





################################################################################
## Updates the censored Xs, just by drawing from the log likelihood of the
## data, censored at the threshold (on the Y scale), where the data comes
## from a Huser-wadsworth scale mixture
## of transformed Gaussians, transformed to GPD
##
## X.s ............................... X, but without the measurement error
## cen 
## prob.below
## theta.mix
## theta.gaussian
##
update.censored.obs.mixture.me <- function(X.s, cen, tau_sqd, thresh.X) {
  
  ll <- matrix(0, nrow(X.s), ncol(X.s))
  sd=sqrt(tau_sqd)
  
  # Draws a collection of univariate truncated normals
  B <- pnorm(thresh.X, X.s[cen], sd=sd)
  U <- B * runif(sum(cen))
  X.update <- qnorm(U, mean=X.s[cen], sd=sd)
  
  return(X.update) 
}

#                                                                              #
################################################################################



