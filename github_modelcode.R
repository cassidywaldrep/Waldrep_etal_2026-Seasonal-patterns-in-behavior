
# Code: Time–activity budgets of eastern mallards: Seasonal patterns and differences among migratory strategies
# Author: Cassidy Waldrep

# Load packages and functions

library(dplyr) 
library(forcats)
library(nimble)
library(nimbleHMC)
library(coda)
library(parallel)

# read in your data. Data associated with this manuscript is summarized per season to protect the raw data. 
# this model uses daily values.
# data <- read.csv()

#-------------------------------------------------------------------------#
###### BAYESIAN MODEL IN NIMBLE  ###################
#-------------------------------------------------------------------------#

# indexing ---------------------------------------------------------------------
# total number of data points
nobs_bhvr <- nrow(data)

# number of sub-stages (seasons)
nss <- length(unique(data$substage_index))

# number of populations (resident, migratory, partial migrant)
npop <- length(unique(data$pop_index))

# number of individuals in study
nind_bhvr <- length(unique(data$bandnum_index))

# number of behaviors
nbhvr <- ncol(dplyr::select(data, starts_with("num_fixes")))

# behaviors (for reference: order matches order in beta_bhvr)
colnames(dplyr::select(data, starts_with("num_fixes")))


# specify GLMMs in BUGS language -----------------------------------------------
model_code <- nimbleCode({
  ## priors --------------------------------------------------------------------
  
  ### behaviors ----------------------------------------------------------------
  for(b in 1:(nbhvr - 1)){ # behaviors
    # overall intercept (i.e., mean probability of each behavior for reference level [logit-scale])
    alpha_bhvr[b] ~ dnorm(0, sd = 1)
    
    # sub-stage and population-level intercept
    for(s in 1:nss){ # sub-stages
      for(k in 1:(npop - 1)){ # populations
        beta_bhvr[k, s, b] ~ dnorm(0, sd = 1)
      } # k
    } # s
    
    for(s in 1:(nss - 1)){ # sub stages
      beta_bhvr[npop, s, b] ~ dnorm(0, sd = 1)
    } # j
    
    # reference level = resident, stationary, winter 
    beta_bhvr[npop, 8, b] <- 0  
    
    # random effect for individual
    for(i in 1:nind_bhvr){ # individuals
      eps_ind_bhvr[i, b] ~ dnorm(0, sd = sd_ind_bhvr[b])
    } # m
    sd_ind_bhvr[b] ~ dunif(0, 1)
  }

  ## likelihoods ---------------------------------------------------------------
  
  ### behaviors --------------------------------------------------------------
  for(m in 1:nobs_bhvr){ # data points
    for(b in 1:(nbhvr - 1)){ # all behaviors except one (stationary)
      # covariate effects on probability of each behavior
      pi_bhvr[m, b] <- exp(alpha_bhvr[b] + 
                             beta_bhvr[pop_bhvr[m], ss_bhvr[m], b] +
                             eps_ind_bhvr[ind_bhvr[m], b])
      
      # taking propensity scores and making them be between 0 and 1 (probability)
      p_bhvr[m, b] <- pi_bhvr[m, b] / (1 + sum(pi_bhvr[m, 1:(nbhvr - 1)]))
    } # b
    
    # calculate probability of reference category (take all probabilities and substract from 1)
    p_bhvr[m, nbhvr] <- 1 - sum(p_bhvr[m, 1:(nbhvr - 1)])
    
    # number of ACC fixes classified as each behavior out of the total
    # p_bhvr is a vector of ACC fixes probabilities
    num_bhvr_fixes[m, 1:nbhvr] ~ dmulti(p_bhvr[m, 1:nbhvr], total_fixes[m]) 
  } # m
})

# set up NIMBLE model ----------------------------------------------------------
# constants
nimble_constants <- list(
  nobs_bhvr = nobs_bhvr,
  nss = nss,
  npop = npop,
  nind_bhvr = nind_bhvr,
  nbhvr = nbhvr,
  ss_bhvr = data$substage_index,
  pop_bhvr = data$pop_index,
  ind_bhvr = data$bandnum_index
)

# data
nimble_data <- list(
  total_fixes = data$total_num_fixes,
  num_bhvr_fixes = dplyr::select(data, starts_with("num_fixes"))
)

# function to set beta initial values with reference level = 0
beta_bhvr_inits <- function(npop, nss, nbhvr){
  beta <- array(
    rlogis(npop * nss * nbhvr, 0, 1),
    dim = c(npop, nss, nbhvr)
  )
  beta[npop, 8,] <- 0
  beta
}

# initial values
inits_function <- function() {
  list(
    alpha_bhvr = rnorm(nbhvr - 1, 0, 1),
    beta_bhvr = beta_bhvr_inits(npop = npop, nss = nss, nbhvr = nbhvr - 1),
    eps_ind_bhvr = matrix(rnorm(nind_bhvr * (nbhvr - 1), 0, 1), nrow = nind_bhvr, ncol = (nbhvr - 1)),
    sd_ind_bhvr = runif((nbhvr - 1), 0.1, 1)
  )
}

inits <- inits_function()

# run NIMBLE model in parallel
# number of chains
nc <- 3

# set up cluster
cl <- makeCluster(nc)

clusterExport(cl, c("model_code", "nimble_constants", "nimble_data"))

for (j in seq_along(cl)) {
  nimble_inits <- inits_function() 
  clusterExport(cl[j], "nimble_inits")
}

# run MCMC in parallel with HMC sampling
start_time <- Sys.time()
out <- clusterEvalQ(cl, {
  
  .libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")
  
  library(nimble)
  library(nimbleHMC)
  library(coda)
  library(parallel)
  
  # build model
  model <- nimbleModel(code = model_code, 
                       constants = nimble_constants,  
                       dat =  nimble_data, 
                       inits = nimble_inits, 
                       buildDerivs = T)
  
  # configure MCMC
  mcmc_Conf  <- configureHMC(model)
  
  # build MCMC
  modelMCMC  <- buildMCMC(mcmc_Conf)
  
  # compile model and MCMC
  Cmodel     <- compileNimble(model)
  CmodelMCMC <- compileNimble(modelMCMC, project = model)
  
  # run model
  CmodelMCMC$run(niter = 11000,
                 nburnin = 2200,
                 thin = 10)
  
  # extract samples
  samples <- as.matrix(CmodelMCMC$mvSamples)
  
  return(samples)
})
end_time <- Sys.time()
(run_time <- end_time - start_time)  # model run time
stopCluster(cl)

# save output
save(out, file = "results.RData")

