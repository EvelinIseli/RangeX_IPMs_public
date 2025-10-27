# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# DATA ANALYSIS SCRIPT: BOOTSTRAPPING LAMBDA ON CLUSTER ----

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Author              : Evelin Iseli
### Data used           : IPM_bootstrappedVR_para_cluster.csv
### Date last modified  : 27.10.2025
### Purpose             : Use the bootstrapped vital rate variables to calculate multiple lambdas per species * site * treatment combination to use them later for 
###                       getting confidence intervals around lambda. Additionally, go through all produced kernels to calculate necessary starting population size.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


rm(list = ls())

## PACKAGES etc. ---------------------------------------------------------------

library(dplyr); library(tidylog); library(janitor); library(tidyverse); library(data.table) # data manipulation
library(parallel) # parallel computing
library(Matrix) # calculating lambda

## FUNCTIONS -------------------------------------------------------------------

# function to get lambda via an iterative method (faster for 1000x1000 or bigger matrices)
lambda.iter = function(k, tol = 1e-8) {
  # k: IPM full kernel
  # tol: tolerance, i.e. track convergence while iteratively getting lambda
  
  if(sum(k == 0) == length(k)) { return(list(lambda = NA, w = NA)) }
  else {
    qmax = 10*tol
    lam = 1
    x = rep(1,nrow(k))   
    k = Matrix(k)
    while(qmax > tol) {
      x1 = k%*%x
      qmax = sum(abs(x1-lam*x))  
      lam = sum(x1)
      # lambda is 0 when population is 0
      if(lam == 0) { x = 0; qmax = 0 }
      else x = x1/lam # one individual in total
    } 
    # lambda: converged population growth rate
    # w: stablized size distribution
    return(list(lambda = lam,w = x/sum(x)))
  }
} 

# make the IPM kernel
mk.kernel <- function(L, U, par, n, ceiling_growth = FALSE, ceiling_seeds = FALSE) {
  # n: number of meshpoints
  # par: vital rates
  # L, U: lower and upper limit of size
  # fun: full kernel 
  # mesh points 
  h <- (U - L)/n
  meshpts <- L + ((1:n) - 1/2) * h
  
  # no ceiling
  if(ceiling_growth == TRUE & ceiling_seeds == TRUE) {
    # ceiling
    G <- h * (outer(meshpts, meshpts, G.t0t_ceiling, params = par))
    R <- h * (outer(meshpts, meshpts, R.t0t_ceiling, params = par))
  } else if (ceiling_growth == TRUE & ceiling_seeds == FALSE){
    G <- h * (outer(meshpts, meshpts, G.t0t_ceiling, params = par))
    R <- h * (outer(meshpts, meshpts, R.t0t, params = par))
  } else if (ceiling_growth == FALSE & ceiling_seeds == TRUE) {
    G <- h * (outer(meshpts, meshpts, G.t0t, params = par))
    R <- h * (outer(meshpts, meshpts, R.t0t_ceiling, params = par))
  } else {
    G <- h * (outer(meshpts, meshpts, G.t0t, params = par))
    R <- h * (outer(meshpts, meshpts, R.t0t, params = par))
  }
  K <- G+R
  return(list(meshpts = meshpts, G=G, R=R, K = K))
}

# 1. probability of surviving (logistic)
survival.t0 <- function(z, params) {
  if (is.na(params$surv.slope)) { # check whether the mdoel is intercept only
    # if yes, define the linear predictor using only the intercept
    u <- params$surv.int
  } else {
    # if not, define the linear predictor using both intercept and slope
    u <- params$surv.int + params$surv.slope * z
  }
  
  # inverse logit transformation
  survival.p <- exp(u)/(1 + exp(u))
  survival.p[u>700] <- 1 # exp(710) gives Inf values
  # return the probability of survival
  return(survival.p)
}

# 2.1 growth function (gaussian)
growth.t0t <- function(z1, z, params) {
  # mean size at t based on size at t0
  if (is.na(params$growth.slope)) {
    mu <- params$growth.int  # use only  intercept
  } else {
    mu <- params$growth.int + params$growth.slope * z  # use both intercept and slope
  }
  
  # sd around mean
  sig <- params$growth.sigma
  
  ## debugging
  #print(paste("z:", z, "mu:", mu, "sig:", sig))
  
  # probability density of new biomass at time t for current biomass at t0
  pd <- dnorm(z1, mean = mu, sd = sig)
  # return probability density
  return(pd)
}


# 2.2 growth function (ceiling)
growth.t0t_ceiling <- function(z1, z, params) {
  # mean size at t based on size at t0
  if (is.na(params$growth.slope)) {
    mu <- params$growth.int  # use only  intercept
  } else {
    mu <- params$growth.int + params$growth.slope * z  # use both intercept and slope
  }
  
  # ceiling
  mu[mu > params$U] = params$U
  print(sum(mu > params$U)) # only mention if ceiling was used
  
  # sd around mean
  sig <- params$growth.sigma
  
  ## debugging
  #print(paste("z:", z, "mu:", mu, "sig:", sig))
  
  # probability density of new biomass at t for current biomass at t0
  pd <- dnorm(z1, mean = mu, sd = sig)
  # return probability density
  return(pd)
}

# 2.3 growth function (gaussian)
# to predict size
growth.t0 <- function(z, params) {
  # mean size at t based on size at t0
  if (is.na(params$growth.slope)) {
    biomass_t1_pred <- params$growth.int  # use only  intercept
  } else {
    biomass_t1_pred <- params$growth.int + params$growth.slope * z  # use both intercept and slope
  }
  return(biomass_t1_pred)
}


# 3. flowering probability (logistic)
flowering.t0 <- function(z, params) {
  if (is.na(params$flow.slope)) { # check whether the model is intercept only
    # if yes, define the linear predictor using only the intercept
    u <- params$flow.int
  } else {
    # if not, define the linear predictor using both intercept and slope
    u <- params$flow.int + params$flow.slope * z
  }
  
  # probability of flowering, inverse logit
  flowering.p <- exp(u)/(1+exp(u))
  flowering.p[u>700] <- 1 # exp(710) gives Inf values
  # return probability of flowering
  return(flowering.p) 
}

# 4.1 seed production (poisson)
seeds.t0 <- function(z, params) {
  # mean seed production at t based on size at t0
  if (is.na(params$seed.slope)) {
    u <- params$seed.int  # use only  intercept
  } else {
    u <- params$seed.int + params$seed.slope * z  # use both intercept and slope
  }
  
  # convert linear predictor into mean seed production
  mu <- exp(u) 
  # generate mean seed count based on Poisson distribution (only integers and no negatives)
  seeds <- pmax(0, rpois(length(mu), lambda = mu))  
  # return seed production
  return(seeds)
}

# 4.2 seed production (ceiling)
seeds.t0_ceiling <- function(z, params) {
  # mean seed production at t based on size at t0
  if (is.na(params$seed.slope)) {
    u <- params$seed.int  # use only  intercept
  } else {
    u <- params$seed.int + params$seed.slope * z  # use both intercept and slope
  }
  
  # convert linear predictor into mean seed production
  mu <- exp(u) 
  # generate mean seed count based on Poisson distribution (only integers and no negatives)
  seeds <- pmax(0, rpois(length(mu), lambda = mu))  
  
  # add a ceiling
  seeds[seeds > params$max_fecundity] <- params$max_fecundity
  # return seed production
  return(seeds)
}

# 4.3 mean seed production (poisson)
seeds.t0_mean <- function(z, params) {
  # mean seed production at t based on size at t0
  if (is.na(params$seed.slope)) {
    u <- params$seed.int  # use only  intercept
  } else {
    u <- params$seed.int + params$seed.slope * z  # use both intercept and slope
  }
  
  # convert linear predictor into mean seed production
  mu <- exp(u) 
  
  return(mu)
}

# 5. probability of seed germination (logistic)
germination <- function(params) {
  # define linear predictor (intercept)
  u <- params$germ.int
  
  # inverse logit transformation
  germination.p <- exp(u) / (1 + exp(u))
  
  # return probability
  return(germination.p)
}

# 6. probability of seedling establishment (logistic)
establishment <- function(params) {
  # define linear predictor (intercept)
  u <- params$est.int
  
  # inverse logit transformation
  establishment.p <- exp(u) / (1 + exp(u))
  
  # return probability
  return(establishment.p)
}


# 7. seedling size (constant)
seedling.t <- function(z1, params) {
  # probability density of recruits
  rpd <- dnorm(z1, mean = params$r_mean, sd = params$r_sd)
  # return probability density
  return(rpd)
}


# growth kernel: survival-growth
G.t0t <- function(z1, z, params) {
  # combine survival and growth
  surv_growth <- survival.t0(z, params) * growth.t0t(z1, z, params)
  # return
  return(surv_growth)
}

# growth kernel: survival-growth, but with ceiling
G.t0t_ceiling <- function(z1, z, params) {
  # combine survival and growth
  surv_growth <- survival.t0(z, params) * growth.t0t_ceiling(z1, z, params)
  # return
  return(surv_growth)
}

# reproduction kernel: reproduction
R.t0t <- function(z1, z, params) {
  # calculate fecundity kernel
  fecundity <- flowering.t0(z, params) * 
    seeds.t0_mean(z, params) * 
    germination(params) *
    establishment(params) *
    seedling.t(z1, params)
  # return
  return(fecundity)
}

# reproduction kernel: reproduction, but with seed production ceiling
R.t0t_ceiling <- function(z1, z, params) {
  # calculate fecundity kernel
  fecundity <- flowering.t0(z, params) * 
    seeds.t0_ceiling(z, params) * 
    germination(params) *
    establishment(params) *
    seedling.t(z1, params)
  # return
  return(fecundity)
}


# function to project an IPM kernel forward (output: matrix (if n > 1) with the number of individuals in different size classes over time)
project.ipm <- function(n0, k, nstep){
  # k: value of mk.kernel, with meshpoints and kernel
  # nstep: number of steps to project
  # n0: number of individuals in different size classes
  
  # project only one step
  if(nstep==1) { 
    nt <- k %*% n0
    re <- nt
  }
  # project n steps
  else if(nstep >=2) {
    ntmat <- matrix(0, nrow=length(n0), ncol = nstep)
    ntmat[,1] <- n0; # each column is nt in one year
    for(t in 2:nstep) { 
      ntmat[,t] <- k %*% ntmat[,t-1] 
    }
    re <- ntmat
  }
  return(re)
}

# function to calculate the number of individuals in the smallest size class required to reach a target pop. size after nstep years
required_sc1_simple <- function(target_pop_size, kernel, nstep, sc1_index = 1) { # target_pop_size: wanted population size 
  # kernel: IPM kernel K 
  # nstep: number of steps/ years of wanted projection 
  # sc1_index: which position in state vector corresponds to smallest size class 
  
  m <- nrow(kernel)  # infer number of meshpoints 
  n0 <- c(1, rep(0, m-1))  # 1 individual in smallest class 
  
  proj <- project.ipm(n0 = n0, k = kernel, nstep = nstep) 
  gT <- sum(proj[, nstep]) # per-1-sc1 yield at year nstep 
  tibble(gT = gT, required_sc1 = target_pop_size / gT) 
}

# function to calculate the required number of seeds for a target pop. size after nstep years
required_seeds_from_params <- function(par, n_mesh, target_pop_size = 5000, nstep = 10, sc1_index = 1) {
  # par: list of vital rate parameters
  # target_pop_size: wanted population size
  # nstep: number of steps/ years of wanted projection
  # sc1_index: which position in state vector corresponds to smallest size class
  # n_mesh: number of meshpoints
  
  K_list <- mk.kernel(L = par$L, U = par$U, par = par, n = n_mesh)
  out_sc1 <- required_sc1_simple(target_pop_size = target_pop_size, kernel = K_list$K,
                                 nstep = nstep, sc1_index = sc1_index)
  
  p_seed_to_sc1 <- plogis(par$germ.int) * plogis(par$est.int) # seed->sc1 probability implied by the same params (logit intercepts)
  
  tibble(
    gT            = out_sc1$gT,
    required_sc1  = out_sc1$required_sc1,
    required_seeds = out_sc1$required_sc1 / p_seed_to_sc1,
    p_seed_to_sc1 = p_seed_to_sc1
  )
}


# function to bootstrap and get a confidence interval around the lambdas as well as find the starting population size with starting population  n = 1
lambda.startpop.boot <- function(a, params_list, start_pop = TRUE) {
  # a = index of the params_list
  
  da <- params_list[[a]]  # get parameters for the specific combination
  
  # generate the kernel with an appropriate number of mesh points (n)
  n_meshpoints <- 500 
  ipm <- mk.kernel(L = as.numeric(da$L), U = as.numeric(da$U), par = da, n = n_meshpoints, ceiling_growth = TRUE, ceiling_seeds = FALSE)
  
  # check for NA values in the kernel (shouldn't have any)
  if (sum(is.na(ipm$K)) > 0) { 
    return(list(lambda = NA, n0 = NA))  # return list with NA if there's an issue with the kernel
  } else {
    # calculate lambda using the lambda.iter function
    cal.a <- lambda.iter(ipm$K)
    lam.a <- cal.a$lambda
    
    # calculate starting population size only for lambda > 1
    if (lam.a > 1) {
      out_sc1 <- required_seeds_from_params(
        par = da,
        n_mesh = n_meshpoints,
        target_pop_size = 5000,
        nstep = 10,
        sc1_index = 1
      )
    } else {
      out_sc1 <- NA  # no need to calculate for lambda <= 1
    }
    
    # return lambda and required starting population (without kernel)
    if (lam.a > 1) {
      return(list(lambda = lam.a, required_sc1 = out_sc1$required_sc1, required_seeds = out_sc1$required_seeds, gT = out_sc1$gT, p_seed_to_sc1 = out_sc1$p_seed_to_sc1))  # consistent list return
    } else {
      return(list(lambda = lam.a, required_sc1 = NA, required_seeds = NA, gT = NA, sp_seed_to_sc1 = NA))  # always return a list with lambda and n0
    }
  }
}


## LOAD DATA -------------------------------------------------------------------

bootpara_wide <- read_csv("/cluster/home/eviseli/ipms/d20251027/IPM_bootstrappedVR_para_cluster.csv")
#bootpara_wide <- read_csv("/Users/eviseli/Documents/GitHub/RangeX_IPMs/cluster/d20251027/IPM_bootstrappedVR_para_cluster.csv")
#bootpara_wide <- read_csv("/Users/eviseli/Documents/GitHub/RangeX_IPMs/cluster/d20250205/IPM_bootstrappedVR_para_cluster.csv")


## ANALYSIS -------------------------------------------------------------------

### LAMBDAS

# generate a list containing different parameters for making full kernels

# initialize an empty list to store results
params_list_bootpara <- list()

# loop through each row of the reshaped data
for (i in 1:nrow(bootpara_wide)) { # nrow(bootpara_wide)
  list_name <- paste(bootpara_wide$species[i], bootpara_wide$site[i], bootpara_wide$treat_combi[i], "bootstrap", bootpara_wide$bootstrap[i], sep = "_")
  
  params_list_bootpara[[list_name]] <- list(
    r_mean = bootpara_wide$r_mean[i],
    r_sd = bootpara_wide$r_sd[i],
    surv.int = bootpara_wide$intercept_survival[i],
    surv.slope = bootpara_wide$slope_survival[i],
    growth.int = bootpara_wide$intercept_size_t1_log[i],
    growth.slope = bootpara_wide$slope_size_t1_log[i],
    growth.sigma = bootpara_wide$sigma[i],
    flow.int = bootpara_wide$intercept_flower_status[i],
    flow.slope = bootpara_wide$slope_flower_status[i],
    seed.int = bootpara_wide$intercept_number_seeds[i],
    seed.slope = bootpara_wide$slope_number_seeds[i],
    germ.int = bootpara_wide$intercept_germ_rate[i],
    est.int = bootpara_wide$intercept_establishment[i],
    L = bootpara_wide$L[i],
    U = bootpara_wide$U[i]
  )
}


# number of cores
n.cores <- 24 # for laptop: detectCores() - 1 
#n.cores <- 4

# calculate lambda and starting population for each entry in params_list_boot in parallel
lambda_results_wstartpop <- mclapply(seq_along(params_list_bootpara), 
                                     function(a) {
                                       tryCatch({
                                         result <- lambda.startpop.boot(a, params_list = params_list_bootpara, start_pop = TRUE)
                                         
                                         # debugging: check the class of the result
                                         if (!is.list(result)) {
                                           message(paste("Iteration", a, "returned a non-list result:", class(result)))
                                         }
                                         
                                         return(result)
                                       }, error = function(e) {
                                         message(paste("Error in iteration", a, ":", e$message))
                                         return(list(lambda = NA_real_,
                                                     gT_sc1 = NA_real_, required_sc1 = NA_real_, p_seed_to_sc1 = NA_real_, required_seeds_sc1 = NA_real_))  # return a consistent list with NAs in case of error
                                       })
                                     }, 
                                     mc.cores = n.cores)
names(lambda_results_wstartpop) <- names(params_list_bootpara)

# prepare a data frame out of the lambda_results list
lambda_results_df <- bind_rows(lambda_results_wstartpop, .id = "parameter_name")


## SAVE OUTPUT -----------------------------------------------------------------


# save lambda_results_df
write.csv(lambda_results_df, "/cluster/home/eviseli/ipms/d20251027/results/IPM_bootstrappedLambda_para_cluster.csv")




