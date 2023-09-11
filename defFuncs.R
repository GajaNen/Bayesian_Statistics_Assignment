
################################################################################
#### estimation, MH, convergence, interpretation of estimates and intervals ####
################################################################################


GibbsRegression <- function(iv_all, iter, priors, d, burn, mu_MH=0, sigma_MH=0.1, 
                            rw, CI_level=0.05, chains, n){ 
  # iv_all: list of vectors of initial values, d: data, priors: vector of prior parameters,
  # burn: number of burn-in, iter: n of iterations, mu_MH & sigma_MH for MH step,
  # rw: whether random walk should be used for MH step, CI_level: credible intervals,
  # n: number of parameters
  
  N <- nrow(d)
  val <- array(rep(NA, iter*(n+1)*chains), dim=c(iter, n+1, chains)) # storage for samples
  X <- as.matrix(cbind(rep(1, N), d[,1:(n-2)]))
  y <- d[,(n-1)]
  
  # separate chains
  for (c in 1:chains){
    
    iv <- iv_all[[c]]
    reg_coef <- as.matrix(iv[1:(n-1)]) # current intercept and regression coefficients
    tau <- iv[n] # current tau
    
    # sample from the conditional posterior iter times 
    for (i in 1:iter){
      
      # tau, posterior shape and rate (a and b)
      a <- N/2 + priors[n,1]
      b <- sum((y - (X %*% reg_coef))**2)/2 + priors[n,2]
      tau <- val[i,n,c] <- rgamma(1, shape=a, rate=b)
      
      # intercept: posterior mu and tau
      tau01 <- N * tau + priors[1,2]
      mu01 <- (sum(y - (X[,-1] %*% reg_coef[-1]))*tau + priors[1,1]*priors[1,2])/tau01
      reg_coef[1] <- val[i,1,c] <- rnorm(1, mu01, sqrt(1/tau01))
      
      # first regression coefficient - MH step
      tau11 <- sum(X[,2]**2) * tau + priors[2,2]
      mu11 <- (sum(X[,2] * (y - (X[,-2] %*% reg_coef[-2]))) * tau + priors[2,1]*priors[2,2])/tau11
      
      # if random walk use the latest value of the coefficients as mean, otherwise
      # use user-specified mean
      mu_MH <- rw*reg_coef[2] + (1-rw)*mu_MH
      sampled_MH <- MetropolisHastings(reg_coef[2], mu=mu_MH, sigma=sigma_MH, mu_post = mu11, sd_post=1/sqrt(tau11), random=rw)
      reg_coef[2] <- val[i,2,c] <- sampled_MH[1] # value sampled in this step
      val[i,n+1,c] <- sampled_MH[2] # accepted or not
      
      # for other coefficients compute posterior mean and tau estimated in this step as before
      for (j in 3:(n-1)){
        tauj1 <- sum(X[,j]**2) * tau + priors[j,2]
        muj1 <- (sum(X[,j]*(y - (X[,-j] %*% reg_coef[-j])))*tau + priors[j,1]*priors[j,2])/tauj1
        reg_coef[j] <- val[i,j,c] <- rnorm(1, muj1, sqrt(1/tauj1))
        
      }
      
    }
    
  }
  val[,n,] <- 1/sqrt(as.vector(val[,n,])) # tau to SD for all chains
  # remove burn in for each chain
  val <- val[c((burn+1):iter),,]
  # means, SD, CI, AR per chain (over rows, dim=1):
  means_c <- apply(val, c(2,3), mean) # column is chain, row is parameter
  means_chain <- means_c[1:n,] # remove acceptance column
  medians_chain <- apply(val, c(2,3), median)[1:n,] # remove accept
  AR_chain <- means_c[n+1,] # keep only accept
  SD_chain <- apply(val, c(2,3), sd)[1:n,] # remove accept
  MCerror_chain <- SD_chain/sqrt(iter-burn)
  CI_chain <- apply(val,c(2,3), function(x){quantile(x, probs=c(CI_level/2,1-CI_level/2))})[1:2,1:n,1:chains]
  # collapse chains by column (parameter) and convert to dataframe
  out <- as.data.frame(apply(val,2,function(y){y}))
  means <- apply(out[,1:n], 2, mean)
  medians <- apply(out[,1:n], 2, median)
  SD <- apply(out[,1:n], 2, sd)
  MCerror <- SD/sqrt((iter-burn)*chains) # niter is now higher
  CI <- out[,1:n] %>% map(quantile, probs=c(CI_level/2,1-CI_level/2))
  AR <- mean(out[,n+1]) # acceptance rate
  
  return(list(bychain=val,allsamples=out, M=means, Md=medians, SD=SD, MCerror=MCerror, CI=CI, Acceptance=AR, Mchain=means_chain, MdChain=medians_chain, SDchain=SD_chain, CIchain=CI_chain, MCerrorChain=MCerror_chain, AcceptanceChain=AR_chain))
  
}

###--------------------------------------------------------------------------###

MetropolisHastings <- function(previous, mu=NULL, sigma=NULL, mu_post, sd_post, random){ # by default it's a Metropolis sampler
  
  # sample from the proposal
  candidate <- rnorm(1, mu, sigma)
  # sample from uniform
  u <- runif(1)
  # acceptance ratio if random walk
  if (random){
    r <- dnorm(candidate, mu_post, sd_post)/dnorm(previous, mu_post, sd_post)
  }
  else { # AR if not random walk
    r <- dnorm(candidate, mu_post, sd_post)/dnorm(previous, mu_post, sd_post)*(dnorm(previous, mu, sigma)/dnorm(candidate, mu, sigma))
  }
  # if denominator is 0, do not accept, else accept with p = min(1, r)
  if (is.nan(r) | is.na(r)) {
    r <- 0} 
  accept <- (u <= r)
  current <- if (accept) candidate else previous
  return(c(current, accept))
  
}

###--------------------------------------------------------------------------###

AssessConvergence <- function(samplesChain, samplesALL, chainMean, nchains, maxlag, discard=NULL, nparam, MHindex){
  nsamples <- nrow(samplesChain[,,1])
  samplesChain <- samplesChain[,1:nparam,]
  # G-R statistic
  grandMean <- rowMeans(chainMean)
  Bvar <- sum((chainMean-grandMean)**2)*(nsamples/(nchains-1)) # between-chain variance per parameter
  Wvar <- matrix(NA, nparam, nchains) # storage for within chain variance per parameter
  AR <- matrix(NA, maxlag, nchains) # storage for autocorrelations per chain per parameter
  arlag1 <- matrix(NA, nrow=nchains) # storage for autocorrelations at lag1 per chain
  for (i in 1:nchains){
    Wvar[,i] <- colSums((samplesChain[,,i] - chainMean[,i])**2)*1/(nsamples-1) #  within chain variance per parameter
    ar <- myAutoCorrelation(nsamples, maxlag, samplesChain[,MHindex,i]) 
    AR[,i] <- ar[[1]] # autocorrelation per chain up to maxlag
    arlag1[i] <- ar[[2]] # autocorrelations at lag 1 per chain
  }
  for (j in 1:nparam){
    # plots by chain
    print(plot(x$bychain[,j,1], type="l", xlab="iteration", ylab="sampled value"))
    lines(samplesChain[,j,2], col="red")
  }
  W <- rowMeans(Wvar) # average within chain variance per parameter
  RG <- (((nsamples-1)/nsamples)*W + Bvar/nsamples)/W # Rubin-Gelman as ratio of weighted sum of within and between chain variance and within variance
  return(list(R=RG, AR=arlag1))
}

###--------------------------------------------------------------------------###

myAutoCorrelation <- function(nsamples, maxlag=30, samples){
  autocorr <- matrix(NA, nrow=maxlag)
  # compute autocorrelations for MH step coefficient
  for (i in 1:maxlag){
    autocorr[i] <- cor(samples[1:(nsamples-i)], samples[(i+1):nsamples])
  }
  return(list(autocorr, autocorr[1]))
}


################################################################################
########################## posterior predictive check ##########################
################################################################################

checkAssumption <- function(samples, X, y, maxlag, nparam){
  
  n.iter <- nrow(samples)
  N <- nrow(X)
  X <- as.matrix(cbind(rep(1, N), X))
  simulated <- matrix(NA, nrow = N, ncol = n.iter)
  res_sim <- matrix(NA, nrow = N, ncol = n.iter)
  res_obs <- matrix(NA, nrow = N, ncol = n.iter)
  t_sim <- matrix(NA, nrow = n.iter)
  t_obs <- matrix(NA, nrow = n.iter)
  for (j in 1:n.iter){
    parameters <- samples[j,]
    # predicted outcomes for each sample of parameters
    y_hat <- X %*% as.matrix(unlist(parameters[1:4]))
    # simulate N observations for each sample of parameters
    simulated[,j] <- rnorm(N, X %*% as.matrix(unlist(parameters[1:4])), sd = parameters[[5]])
    # calculate residuals for observed values
    res_obs[,j] <- y - y_hat
    # calculate residuals for simulated values in this step
    res_sim[,j] <- simulated[,j] - y_hat
    # get a distribution of residual autocorrelations for observed
    AR_o <- myAutoCorrelation(N, maxlag=maxlag, res_obs[,j])
    AR_obs <- AR_o[[1]]
    # a distro of residual AR for simulated
    AR_s <- myAutoCorrelation(N, maxlag=maxlag, res_sim[,j])
    AR_sim <- AR_s[[1]]
    # calculate discrepancy measure for observed data
    t_obs[j] <- mean(abs(AR_obs) > 0.1)
    # calculate discrepancy measure for simulated data
    t_sim[j] <- mean(abs(AR_sim) > 0.1)
  }
  # calculate two-sided ppp
  ppp <- mean(t_sim > t_obs)
  return(ppp)
  
}


################################################################################
############################### model selection ################################
################################################################################

# DIC

my_dic <- function(data, samples, means, nparam){
  N <- nrow(data)
  n.iter <- nrow(samples)
  n.pred <- nparam - 2
  y <- data[,3]
  X <- as.matrix(cbind(rep(1, N), data[,1:n.pred]))
  # deviance given posterior means
  D_mean <- my_likelihood(means, X, y, n.pred = n.pred)
  mean_D <- matrix(NA, n.iter)
  for (i in 1:n.iter){
    # for each sample of parameters calculate deviance
    mean_D[i] <- my_likelihood(samples[i,], X, y, n.pred = n.pred)
  }
  # average deviance across the posterior samples
  mean_D <- mean(mean_D)
  # effective number of parameters (penalty)
  pD <- mean_D - D_mean
  DIC <- pD + mean_D
  return(DIC)
}

###--------------------------------------------------------------------------###

my_likelihood <- function(parameters, X, y, n.pred){
  # compute loglikelihood of data given a vector of parameters summed over i
  L <- dnorm(y, mean = X %*% as.matrix(unlist(parameters[1:(n.pred+1)])), sd = parameters[[n.pred+2]])
  return(sum(-2*log(L)))
}

###--------------------------------------------------------------------------###