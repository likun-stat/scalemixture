library(fields)
library(scales)

###################################################################################
## Main sampler

## Y ........................................... censored observations on GPD scale
## S ............................................................. location indices
## cen ........................................................... indicator matrix
## thresh .................................................... threshold on Y scale
## initial.values .................. a list: delta, rho, tau, theta.gpd, prob.below
##                                             X.s, R
## n.updates .................................................... number of updates
## thin ............................................. number of runs in each update
## experiment.name
## echo.interval ......................... echo process every echo.interval updates
## sigma.m
## prop.Sigma
## true.params ..................... a list: delta, rho, tau, theta.gpd, prob.below
##                                              X.s, R
## lower.prob.lim ......................................... for updating prob.below

scalemix.sampler.02 <- function(Y, S, cen, thresh,
                                initial.values,
                                n.updates, thin=10,
                                experiment.name="Huser-wadsworth",
                                echo.interval=50,
                                sigma.m=NULL, prop.Sigma=NULL, 
                                true.params=NULL, n.pl=1, lower.prob.lim=0.5) {
  
  library(doParallel)
  library(foreach)
  save.bit <- TRUE
  
  # Constants to control how many Metropolis steps get executed at each
  # iteration of the main loop
  n.metr.updates.delta <- 2
  n.metr.updates.rho <- 2
  n.metr.updates.theta.gpd <- 4
  n.metr.updates.prob.below <- 4
  n.metr.updates.tau <- 2
  n.metr.updates.R <- 2
  
  # Constants to control adaptation of the Metropolis sampler
  c.0 <- 10
  c.1 <- 0.8
  k <- 3  # the iteration offset
  metr.opt.1d <- 0.41
  metr.opt.2d <- 0.35
  
  # Hyper parameters for the prior of the mixing distribution parameters and 
  # the correlation parameters
  hyper.params.delta <- 1
  hyper.params.rho <- 2
  hyper.params.theta.gpd <- 0
  hyper.params.tau <- c(0.1,0.1)
  # hyper.params.prob.below <- c(0, 10)
  
  # A small number
  eps <- 1e-06
  
  # Bookkeeping
  n.s <- nrow(Y)
  n.t <- ncol(Y)
  h <- rdist(S)
  diag(h)  <- 0
  registerDoParallel(cores=n.t)
  
  # Load initial values
  delta <- initial.values$delta
  log.rho <- log(initial.values$rho)
  tau <- initial.values$tau
  scale <- initial.values$theta.gpd[2]
  shape <- initial.values$theta.gpd[3]
  prob.below <- initial.values$prob.below
  # logit.prob.below <- logit(prob.below, c(lower.prob.lim, 1))
  X.s <- initial.values$X.s
  R <- initial.values$R
  
  X.s.accept <- rep(1, n.t)
  
  theta.gpd <- c(thresh, scale, shape)
  
  
  #**********   Eigendecomposition of the correlation matrix   ****************#
  Sigma <- corr.fn(h, log.rho)
  eig.Sigma <- eigen(Sigma, symmetric=TRUE)
  V <- eig.Sigma$vectors
  d <- eig.Sigma$values
  
  # Initialize X and X.s
  thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau, delta = delta)
  X <- Y
  X[!cen] <- gpd.2.scalemix.me(Y[!cen], tau_sqd = tau, delta = delta, theta.gpd = theta.gpd, prob.below=prob.below)
  # X[cen] <- update.censored.obs.mixture.me(X.s, cen = cen, tau_sqd = tau, thresh.X = thresh.X)
  
  # Initialize trace objects
  accepted <- rep(0, n.updates)
  
  X.s.trace <- array(NA, dim=c(n.updates, n.s, n.t))
  X.trace <- array(NA, dim=c(n.updates, n.s, n.t))
  X.s.accept.trace <- matrix(NA, n.updates, n.t)
  log.rho.trace <- rep(NA, n.updates)
  tau.trace <- rep(NA, n.updates)
  delta.trace <- rep(NA, n.updates)
  theta.gpd.trace <- matrix(NA, n.updates, 2)
  # prob.below.trace <- rep(NA, n.updates)
  R.trace <- matrix(NA, n.updates, n.t)
  
  sd.ratio.trace <- rep(NA, n.updates)
  
  # Column names for trace objects
  colnames(theta.gpd.trace) <- c("scale", "shape")
  
  # Fill in trace objects with initial values
  X.s.trace[1, , ] <- X.s
  X.trace[1, , ] <- X
  X.s.accept.trace[1, ] <- X.s.accept
  log.rho.trace[1] <- log.rho
  tau.trace[1] <- tau
  delta.trace[1] <- delta
  theta.gpd.trace[1, ] <- theta.gpd[c(2,3)]
  # prob.below.trace[1] <- prob.below
  R.trace[1, ] <- R
  if(is.null(prop.Sigma$theta.gpd))  prop.Sigma$theta.gpd<-diag(2)
  sd.ratio.trace[1] <- prop.Sigma$theta.gpd[2, 2]
  
  # For tuning Metropolis updates of theta
  if (is.null(sigma.m$delta)) sigma.m$delta <- 1
  if (is.null(sigma.m$rho)) sigma.m$rho <- 1
  if (is.null(sigma.m$theta.gpd)) sigma.m$theta.gpd <- 1
  # if (is.null(sigma.m$prob.below)) sigma.m$prob.below <- 1
  if (is.null(sigma.m$tau)) sigma.m$tau <- 1
  if (is.null(sigma.m$R)) sigma.m$R <- rep(1, n.t)
  r.hat.delta <- NA
  r.hat.rho <- NA
  r.hat.theta.gpd <- NA
  r.hat.prob.below <- NA
  r.hat.tau <- NA
  # r.hat.R <- rep(NA, n.t)
  
  sd.ratio <- 7
  
  
  for (i in 2:n.updates) {
    
    ################################################################
    ## Update Metropolis adaptation parameters
    ################################################################
    gamma1 <- c.0 / (i + k)^(c.1)
    gamma2 <- 1 / (i + k)^(c.1)
    
    for (j in 1:thin) {
    
    ################################################################
    ## Update X -- branching: use X.s
    ################################################################
    X[cen] <- update.censored.obs.mixture.me(X.s = X.s, cen = cen, tau_sqd = tau, thresh.X = thresh.X)
    X[!cen] <- gpd.2.scalemix.me(Y[!cen], tau_sqd = tau, delta = delta, theta.gpd = theta.gpd, prob.below=prob.below)

    ################################################################
    ## Update X.star  *Parallel
    ################################################################
    X.s.res <- X.s.update.mixture.me.par(R, Y, X, X.s, cen, 
                                         prob.below, theta.gpd, delta,
                                         tau, V, d, v.q=2, n.chain = 100, thresh.X = thresh.X)
    X.s <- X.s.res$X.s
    X.s.accept <- X.s.res$accepted
    
    
    
    ################################################################
    ## Update rho
    ################################################################
    
    
    metr.out.rho <- static.metr(z = R, starting.theta = log.rho,
                                likelihood.fn = rho.update.mixture.me.likelihood, prior.fn = log.rho.prior,
                                hyper.params = hyper.params.rho, n.updates = n.metr.updates.rho, prop.Sigma = 1, sigma.m=sigma.m$rho, verbose=FALSE,
                                X.s = X.s, R = R, S = S)
    r.hat.rho <- metr.out.rho$acc.prob
    log.rho <- metr.out.rho$trace[n.metr.updates.rho]
    sigma.m$rho <- exp(log(sigma.m$rho) + gamma1*(r.hat.rho - metr.opt.1d))
    
        
    ## Re-create covariance matrix and eigenvectors/eigenvalues
    Sigma   <- corr.fn(h, log.rho)
    eig.Sigma <- eigen(Sigma, symmetric=TRUE)
    V <- eig.Sigma$vectors
    d <- eig.Sigma$values

  
      
    ################################################################
    ## Update R  *Parallel
    ################################################################
     Metr.R<-foreach(t = 1:n.t, .combine = "rbind") %dopar% {
       
       metr.out.R <- static.metr(z = R, starting.theta = R[t],
                                   likelihood.fn = Rt.update.mixture.me.likelihood, prior.fn = huser.wadsworth.prior,
                                   hyper.params = delta, n.updates = n.metr.updates.R, prop.Sigma = 1, sigma.m=sigma.m$R[t], verbose=FALSE,
                                   X.s = X.s[,t], delta = delta, V = V, d = d)
       c(metr.out.R$trace[n.metr.updates.R],
       exp(log(sigma.m$R[t]) + gamma1*(metr.out.R$acc.prob - metr.opt.1d)))
     }
    
    R <- Metr.R[, 1]
    sigma.m$R <- Metr.R[, 2]
    
    
    
    ################################################################
    ## Update tau
    ################################################################
    metr.out.tau <-  static.metr(z = R, starting.theta = tau, 
                                 likelihood.fn = tau.update.mixture.me.likelihood, prior.fn = tau.sqd.prior,
                                 hyper.params = hyper.params.tau, n.updates = n.metr.updates.tau, prop.Sigma = 1, sigma.m=sigma.m$tau, verbose=FALSE, 
                                 Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, delta = delta, theta.gpd = theta.gpd)
    tau <- metr.out.tau$trace[n.metr.updates.tau]
    r.hat.tau <- metr.out.tau$acc.prob
    sigma.m$tau <- exp(log(sigma.m$tau) + gamma1*(r.hat.tau - metr.opt.1d))
    
    
    
    ################################################################
    ## Update delta
    ################################################################
    metr.out.delta <- static.metr(z = R, starting.theta = delta, 
                                     likelihood.fn = delta.update.mixture.me.likelihood, prior.fn = interval.unif,
                                     hyper.params = hyper.params.delta, n.updates = n.metr.updates.delta, prop.Sigma = 1, sigma.m=sigma.m$delta, verbose=FALSE, 
                                     Y = Y, X.s = X.s, cen = cen, prob.below = prob.below, theta.gpd = theta.gpd, tau_sqd = tau)
    delta <- metr.out.delta$trace[n.metr.updates.delta]
    r.hat.delta <- metr.out.delta$acc.prob
    sigma.m$delta <- exp(log(sigma.m$delta) + gamma1*(r.hat.delta - metr.opt.1d))
    

    
    ################################################################
    ## Update prob.below
    ################################################################
#    metr.out.prob.below <- static.metr(R, logit.prob.below,
#                                       logit.prob.below.update.mixture.me.likelihood,
#                                       normal.scalar,
#                                       hyper.params=hyper.params.prob.below,
#                                       n.updates=n.metr.updates.prob.below,
#                                       prop.Sigma=1,
#                                       sigma.m=sigma.m$prob.below,
#                                       Y=Y, X.s=X.s, cen=cen,
#                                       theta.mix=theta.mix, theta.gaussian=theta.gaussian,
#                                       theta.gpd=theta.gpd, lower.prob.lim=lower.prob.lim)
#    logit.prob.below <- metr.out.prob.below$trace[n.metr.updates.prob.below]
#    prob.below <- ilogit(logit.prob.below, c(lower.prob.lim, 1))
#    r.hat.prob.below <- metr.out.prob.below$acc.prob
#    sigma.m$prob.below <- exp(log(sigma.m$prob.below) +
#                    gamma1*(r.hat.prob.below - metr.opt.1d))
    
    # Re-calculate the threshold on the X-scale    
    thresh.X <- qmixture.me.interp(p = prob.below, tau_sqd = tau, delta = delta)

    
    
    ################################################################
    ## Update theta.gpd
    ################################################################
    metr.out.theta.gpd <- static.metr(z = R, starting.theta = theta.gpd[2:3],
                                        likelihood.fn = theta.gpd.update.mixture.me.likelihood, prior.fn = gpd.prior,
                                        hyper.params = hyper.params.theta.gpd, n.updates = n.metr.updates.theta.gpd, prop.Sigma = prop.Sigma$theta.gpd, sigma.m=sigma.m$theta.gpd, verbose=FALSE, 
                                        Y =Y, X.s = X.s, cen = cen, prob.below = prob.below, delta = delta, tau_sqd = tau, loc = thresh, thresh.X=thresh.X)
    
    theta.gpd[c(2,3)] <- metr.out.theta.gpd$trace[n.metr.updates.theta.gpd, ]
    r.hat.theta.gpd <- metr.out.theta.gpd$acc.prob
    sigma.m$theta.gpd <- exp(log(sigma.m$theta.gpd) + gamma1*(r.hat.theta.gpd - metr.opt.2d))
    }
    
    
    
    # ---------------- Attempt to adapt proposal covariance ------------------
    if (r.hat.theta.gpd > 0) {
      sd.ratio.hat <- sd(metr.out.theta.gpd$trace[ ,1]) / sd(metr.out.theta.gpd$trace[ ,2])
    } else {
      sd.ratio.hat <- 1
    }
    sd.ratio <- exp(log(sd.ratio) + gamma1*(log(sd.ratio.hat) - log(sd.ratio)))
    prop.Sigma$theta.gpd <-  matrix(c(1, prop.Sigma$gpd.corr/sd.ratio, prop.Sigma$gpd.corr/sd.ratio, 1/sd.ratio^2), 2, 2)

   
    
    # ------------------------- Fill in trace objects ---------------------------
    X.s.trace[i, , ] <- X.s
    X.trace[i, , ] <- X
    X.s.accept.trace[i, ] <- X.s.accept
    log.rho.trace[i] <- log.rho
    tau.trace[i] <- tau
    delta.trace[i] <- delta
    theta.gpd.trace[i, ] <- theta.gpd[c(2,3)]
    # prob.below.trace[i] <- prob.below
    R.trace[i, ] <- R
    sd.ratio.trace[i] <- sd.ratio


    
    # ----------------------------- Echo progress --------------------------------
    cat("Done with", i, "updates,\n")
    if ((i %% echo.interval) == 0) {
      cat(" Acceptance rate for X^* is", mean(X.s.accept.trace[(i-echo.interval):i, ]), "\n")
      pdf(file=sprintf("%s_progress.pdf", experiment.name))
      par(mfrow=c(3,2))
      plot(log.rho.trace, type="l", ylab=expression(rho))
      if (!is.null(true.params)) abline(h=log(true.params$rho), lty=2, col=2, lwd=3)
      plot(tau.trace, type="l", ylab=expression(tau^2))
      if (!is.null(true.params)) abline(h=true.params$tau, lty=2, col=2, lwd=3)
      plot(delta.trace, type="l", ylab=expression(delta))
      if (!is.null(true.params)) abline(h=true.params$delta, lty=2, col=2, lwd=3)
      # plot(prob.below.trace, type="l", ylab=expression(pi))
      # if (!is.null(true.params)) abline(h=true.params$prob.below, lty=2, col=2, lwd=3)
      plot(theta.gpd.trace[ ,1], type="l", ylab="scale")
      if (!is.null(true.params)) abline(h=true.params$theta.gpd[2], lty=2, col=2, lwd=3)
      plot(theta.gpd.trace[ ,2], type="l", ylab="shape")
      if (!is.null(true.params)) abline(h=true.params$theta.gpd[3], lty=2, col=2, lwd=3)
      par(mfrow=c(3, 3))
      for (j in 1:min(18, n.t)) {
        plot(R.trace[ , j], type="l", main=paste("Replication", j), 
             ylab=paste0("R[", j, "]"),  ylim=range(R.trace[ , 1:min(18, n.t)], na.rm=TRUE))
        if (!is.null(true.params)) abline(h=true.params$R[j], lty=2, lwd=3, col="gray80")
      }
      # For each year, plot the location with an exceedance that happens to be first in the matrix
      for (j in 1:min(18, n.t)) { 
        plot.loc <- which(!cen[ ,j])[1]
        if (!is.na(plot.loc)) {
          plot(X.s.trace[ , plot.loc, j], type="l",
               ylab=paste0("X^*[", plot.loc, ",", j, "]"))
          if (!is.null(true.params)) abline(h=true.params$X.s[plot.loc, j], lty=2, lwd=3, col="gray80")
        } else {
          plot(1, 1, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
        }
      }
      # For each year, plot the location with a non-exceedance that happens to be first in the matrix
      for (j in 1:min(18, n.t)) { 
        plot.loc <- which(cen[ ,j])[1]
        if (!is.na(plot.loc)) {
          plot(X.s.trace[ , plot.loc, j], type="l",
               ylim=c(min(X.s.trace[ , plot.loc, j], na.rm=TRUE), thresh.X),
               ylab=paste0("X^*[", plot.loc, ",", j, "]"))
          abline(h=thresh.X, lty=1, lwd=3, col="gray80")
        } else {
          plot(1, 1, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
        }
      }
      dev.off()
      
      # lines(curr.draw, col=alpha("gray80", 0.05))
      # # points(X.imp[ ,j], col=2, pch=20, cex=0.5)
      # Sys.sleep(0.1)
      state <- list(Y=Y, cen=cen, S=S, thresh=thresh,
                    experiment.name="Huser-wadsworth",
                    i=i, sigma.m=sigma.m, prop.Sigma=prop.Sigma, sd.ratio.trace=sd.ratio.trace,
                    X=X, X.s=X.s, theta.gpd=theta.gpd, prob.below=prob.below,
                    delta=delta, R=R, tau=tau, log.rho=log.rho)
      out.obj <- list(Y=Y, cen=cen, S=S, thresh=thresh,
                      i=i,
                      X.s.trace=X.s.trace,
                      X.trace=X.trace,
                      X.s.accept.trace=X.s.accept.trace,
                      log.rho.trace=log.rho.trace,
                      tau.trace=tau.trace,
                      delta.trace=delta.trace,
                      theta.gpd.trace=theta.gpd.trace,
                      # prob.below.trace=prob.below.trace,
                      R.trace=R.trace,
                      sd.ratio.trace=sd.ratio.trace)
      save(state, out.obj, file=sprintf("%s_progress_%1.1s.RData", experiment.name, save.bit))
      save.bit <- !save.bit
    }
  }
  
  return(out.obj)
  
  
}



  
  
  
