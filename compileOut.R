require(datasets)
require(tidyverse)
require(purr)
require(foreign)
require(haven)
require(logspline)
require(MASS)
require(pander)
require(car)

source(defFuncs.R)

################################################################################
################################## DATA ########################################
################################################################################

set.seed(1813544)
dat <- as.data.frame(read_spss("0_BA_Daten_Version3.1 final set.sav"))
# select relevant variables
dat_clean <- dat[,c(4:7)]
# center predictors
dat_clean[,c(1:3)] <- scale(dat_clean[,1:3], scale=F)
# fit a lm
m <- lm(TMD_all~WHO_all+SWL_all+FFA_all, data = dat_clean)
koef <- coef(m)
pander(psych::describe(dat_clean)[,c(3,4,5,8,9,11,12)], 
       caption = "Descriptive statistics of the data.")


################################################################################
#### estimation, MH, convergence, interpretation of estimates and intervals ####
################################################################################


# priors
prior_informative <- matrix(c(0, 0, 0, -0.39, 151, 0.001, 0.001, 0.001, 100, 0.151), 
                            ncol = 2)
prior_uninformative <- matrix(c(0, 0, 0, 0, 0.01, 0.001, 0.001, 0.001, 0.001, 0.01), 
                              ncol = 2)

# models with convergence, autocorrelations and table
x <- GibbsRegression(iv_all=list(c(0,0,0,0,1), c(10, -0.3, -0.4, -0.5, 0.5)), 
                     iter=10000, 
                     priors=prior_uninformative, 
                     d=dat_clean, 
                     burn=2000, 
                     sigma_MH = 0.1, 
                     chains=2, 
                     n=5, 
                     rw=1)
convx <- AssessConvergence(x$bychain, x$allsamples, x$Mchain, 2, 10, nparam=5, MHindex=2)

x2 <- GibbsRegression(iv_all=list(c(0,0,0,0,1), c(10, -0.3, -0.4, -0.5, 0.5)), 
                      iter=10000, 
                      priors=prior_informative, 
                      d=dat_clean, 
                      burn=2000, 
                      sigma_MH = 0.1, 
                      chains=2, 
                      n=5, 
                      rw=1)

convx2 <- AssessConvergence(x2$bychain, 
                            x2$allsamples, 
                            x2$Mchain, 
                            2, 
                            10, 
                            nparam=5, 
                            MHindex=2)

out2 <- as.data.frame(cbind(x2$M, t(matrix(unlist(x2$CI), ncol=5, nrow=2)), convx2$R, x2$MCerror))
colnames(out2) <- c("Mean", "95% CI lwrbd", "95% CI uppbd", "R", "MC error")
rownames(out2) <- c("intercept", "well-being", "satisfaction with life", "mindfulness", "sigma")

myAutoCorrelation(8000,samples = x2$bychain[,2,1])
myAutoCorrelation(8000,samples = x2$bychain[,2,2])
x2$AcceptanceChain

x4 <- GibbsRegression(iv_all=list(c(0,0,0,0,1), c(10, -0.3, -0.4, -0.5, 0.5)), iter=10000, priors=prior_informative, d=dat_clean, burn=2000, mu_MH=koef[2], sigma_MH = 0.11821, chains=2, n=5, rw=1)
convx4 <- AssessConvergence(x4$bychain, x4$allsamples, x4$Mchain, 2, 10, nparam=5, MHindex=2)

myAutoCorrelation(8000, samples=x4$bychain[,2,1])
myAutoCorrelation(8000, samples=x4$bychain[,2,2])
x4$AcceptanceChain

x3 <- GibbsRegression(iv_all=list(c(0,0,0,0,1), c(10, -0.3, -0.4, -0.5, 0.5)), iter=10000, priors=prior_informative, d=dat_clean, burn=2000, mu_MH=0, sigma_MH = 1, chains=2, n=5, rw=0)
convx3 <- AssessConvergence(x3$bychain, x3$allsamples, x3$Mchain, 2, 10, nparam=5, MHindex=2)

myAutoCorrelation(8000, samples=x3$bychain[,2,1])
myAutoCorrelation(8000, samples=x3$bychain[,2,2])
x3$AcceptanceChain

pander(out2, caption = "Estimates, CI and convergence for the model with informative prior.")


################################################################################
########################## posterior predictive check ##########################
################################################################################

checkAssumption(x$bychain[,1:5,1], dat_clean[,1:3], dat_clean[,4], maxlag = 40)

checkAssumption(x$bychain[,1:5,1], dat_clean[,1:3], dat_clean[,4], maxlag = 10)

durbinWatsonTest(m)

################################################################################
############################### model selection ################################
################################################################################

# DIC

DIC_uninf <- my_dic(dat_clean, x$bychain[,,1], x$M, 5)
DIC_inf <- my_dic(dat_clean, x2$bychain[,,1], x2$M, 5)

h1m <- x$allsamples[which(x$allsamples[,2,] < 0 & x$allsamples[,4,1] < 0),]
meansh1 <- apply(h1m, 2, mean)
h2m <- x$allsamples[which(x$allsamples[,2,] > 0 & x$allsamples[,4,1] < 0),]
meansh2 <- apply(h2m, 2, mean)

DIC_h1 <- my_dic(dat_clean, h1m, meansh1, 5)
DIC_h2 <- my_dic(dat_clean, h2m, meansh2, 5)

dic <- matrix(c("Hu informative prior", "Hu uninformative", "H1", "H2", 
                DIC_inf, DIC_uninf, DIC_h1, DIC_h2), ncol=2)
pander(dic)

# BAYES FACTOR 

# BAYES FACTOR 

samples <- x$allsamples
b_wb <- koef[2] # mean for posterior
b_tmd <- koef[4]
BF <- matrix(NA, nrow = 10, ncol=3)
for (i in 1:10){
  J <- i
  b <- J/nrow(dat_clean) # fraction of data to be used
  X <- as.matrix(dat_clean[,c(1,3)]) # selected predictors involved in the hypothesis
  sigma2 <- 70.7281 # residual variance of the model
  # covariance matrix for regression coefficients
  cov_m <- sigma2*solve(t(X)%*%X)
  # sample from prior of unconstrained
  samples_pr <- mvrnorm(10000, c(0,0), cov_m/b)
  # sample from posterior of unconstrained
  samples_post <- mvrnorm(10000, c(b_wb, b_tmd), cov_m)
  # complexity for H1: proportion of prior in agreement with H1
  ci <- mean((samples_pr[,1] < 0) & (samples_pr[,2] < 0))
  # fit for H1: proportion of posterior in agreement with H1
  fi <- mean((samples_post[,1] < 0) & (samples_post[,2] < 0))
  BF[i,1] <- fi/ci
  # complexity for H2: proportion of prior in agreement with H2
  ci2 <- mean((samples_pr[,1] > 0) & (samples_pr[,2] < 0))
  # fit for H2: proportion of posterior in agreement with H2
  fi2 <- mean((samples_post[,1] > 0) & (samples_post[,2] < 0))
  BF[i,2] <- fi2/ci2
  BF[i,3] <- fi/ci/(fi2/ci2)
}

colnames(BF) <- c("BF1u", "BF2u", "BF12")
rownames(BF) <- c(1:10)

pander(BF)