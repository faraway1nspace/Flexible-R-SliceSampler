############################
# Demonstration of using the Flexible R Slice Sampler for a random-effects model
# and comparison with JAGS
# 
# A Poisson regression model with: 
# - individual random-effects for the intercepts
# - individual random-effects for slopes
# - hafl-Student-T priors on sigmas

source('../src/flexible_slice_sampler.R')
library(rjags)
library(coda)

# student t distribution
dstudent_t <- function (x, df, mu = 0, sigma = 1, log=FALSE) {
    if (log) {
        return(dt((x - mu)/sigma, df = df, log = TRUE) - log(sigma))
    } else {
        return(dt((x - mu)/sigma, df = df)/sigma)
    }}

###############
# GENERATE DATA: Poisson regression with n=40 individuals, and T=4 observations y each.
# this is like taking a repeated count 4 times for each i individual
generate_fake_poisson_data <- function(seed=42,
                                       n=40, # grouping term (individuals)
                                       T=4, # reps per indi
                                       beta0=log(1.5), #intercept
                                       beta1=0.4, # slope
                                       sigma0=0.4, # root-variance in random-intercepts
                                       sigma1=0.05 # root-variance in random-slopes
){
    # generates fake data for a poissoin regression with grouping
    # by n=40 individuals and T=4 measurements
    set.seed(seed)
    # true betas - fixed effects (intercept and slope)
    beta_true <- c(beta0, beta1)
    # betas in matrix form
    BETA <- matrix(beta_true, nrow=2,ncol=n)
    # model matrix
    X <<- matrix(c(rep(1,T), 0:(T-1)), ncol=2)
    # random effects: random intercepts per n=40 individuals
    epsilion0 <- rnorm(n, 0, sd=sigma0)
    # random effects: random slope per n=40 individuals
    epsilion1 <- rnorm(n, 0, sd=sigma1)
    # random effects in matrix format
    EPSILON <- matrix(c(epsilion0, epsilion1), 2, ncol=length(epsilion0),byrow=TRUE)
    # poisson lambda
    LOGLAMBDA_true <- X%*%BETA +  X%*%EPSILON
    # Y: true counts per n x T
    y_true <- rpois(
        n=prod(dim(LOGLAMBDA_true)),
        lambda = exp(LOGLAMBDA_true)
    )
    return(list(
        Y=matrix(y_true, T, n),
        X=X,
        n=n, # grouping term (individuals)
        T=T, # reps per indi
        beta0=beta0, #intercept
        beta1=beta1, # slope
        sigma0=sigma0, # root-variance in random-intercepts
        sigma1=sigma1 # root-variance in random-slo        
    ))
}

# generate the dat
simulate_data <- generate_fake_poisson_data(
    seed=42, n=40, T=4, beta0=log(1.5), beta1=0.4, sigma0=0.4, sigma1=0.05
)
# assign to workspace
for(varnm in names(simulate_data)){ assign(varnm, simulate_data[[varnm]]) }

print(Y)
print(X)


##################################
# PART I: Analysis with Slice Sampler
# ... later we'll compare with Jags

# NAMES OF VARIABLES
names_fixedeffects = c('beta0','beta1')
names_randomeffects = c(paste0('eps0_',1:n), paste0('eps1_',1:n))
names_sigmas <- c('sigma0', 'sigma1')

# INITIAL ESTIMATES
x.init.beta <- setNames( c(log(mean(Y)),0), nm=names_fixedeffects)
x.init.eps <- setNames( c(rep(0,n), rep(0,n)), nm=names_randomeffects)
x.init.sigmas <- setNames( c(0.05,0.05), nm=names_sigmas)

x.init <- c(x.init.beta, x.init.eps, x.init.sigmas)

# DATA
data_likelihood <- list(
    Y = Y, # count series
    X=X, # model matrix
    n = n, # number of groups
    T = T, # number of observations per groups
    names_fixedeffects = names_fixedeffects,
    names_randomeffects = names_randomeffects,
    names_sigmas = names_sigmas
)

#######################################
# Posterior Functions: 
# - `log_posterior_REregression_dnorm_prior` - normal prior on fixed effects
# - `log_posterior_dnormlike_halfStudentt_hyperprior` half-student T on the sigmas that govern the random-effects
# - `log_posterior_REregression_dnorm_REprior` - the random-effects themselves; notice that he priors

# log posterior with log-normal density on lambda
log_posterior_REregression_dnorm_prior <- function(x_target,
                                           x_all,
                                           data_likelihood,
                                           prior_parameters
                                    ){
    # data for likelihood: y values
    Y <- data_likelihood$Y
    X <- data_likelihood$X
    n <- data_likelihood$n
    names_fixedeffects <- data_likelihood$names_fixedeffects
    names_randomeffects <- data_likelihood$names_randomeffects

    # get fixed effects
    betas_vector <- x_all[names_fixedeffects]
    epsilons_vector <- x_all[names_randomeffects]
    
    # matrices for fixed and random effects
    BETAS <- matrix(betas_vector, nrow=2, ncol=n)
    EPS <- matrix(epsilons_vector, 2, ncol=n,byrow=TRUE)

    # expected values of count-distribution (combine fixed and random effects)
    LAMBDA <- exp(X%*%BETAS +  X%*%EPS)

    # loglikelihood
    loglike <- sum(dpois(
        x=as.numeric(Y), # vectorize observations
        lambda=as.numeric(LAMBDA), # vectorize expectations
        log=TRUE
    ))
    
    # normal prior on betas
    log_prior <- dnorm(
        x=x_target,
        mean=prior_parameters[['prior_mean']],
        sd=prior_parameters[['prior_sd']], # note the sqrt(1/tau) for JAGS-compatibility,
        log=TRUE
    )
    return(loglike + log_prior)}

#
log_posterior_REregression_dnorm_REprior <- function(x_target,
                                            x_all,
                                            data_likelihood,
                                            prior_parameters
                                            ){
    # data for likelihood: y values
    Y <- data_likelihood$Y
    X <- data_likelihood$X
    n <- data_likelihood$n
    names_fixedeffects <- data_likelihood$names_fixedeffects
    names_randomeffects <- data_likelihood$names_randomeffects
    names_sigmas <- data_likelihood$names_sigmas
    
    # get fixed effects
    betas_vector <- x_all[names_fixedeffects]
    # get random effects
    epsilons_vector <- x_all[names_randomeffects]
    # get names of sigmas (hyperprior)
    sigmas <- x_all[names_sigmas]
    
    # matrices for fixed and random effects
    BETAS <- matrix(betas_vector, nrow=2, ncol=n)
    EPS <- matrix(epsilons_vector, 2, ncol=n,byrow=TRUE)

    # expected values of count-distribution (combine fixed and random effects)
    LAMBDA <- exp(X%*%BETAS +  X%*%EPS)
    
    # loglikelihood
    loglike <- sum(dpois(
        x=as.numeric(Y), # vectorize observations
        lambda=as.numeric(LAMBDA), # vectorize expectations
        log=TRUE
    ))

    # normal prior on betas
    log_prior <- dnorm(
        x=x_target,
        mean=0,
        sd=sigmas[prior_parameters$which_sigma_idx],
        log=TRUE
    )
    return(loglike + log_prior)}


# a Half-Student-T hyperprior for the sigmas that govern the random effects
# 
log_posterior_dnormlike_halfStudentt_hyperprior <- function(x_target,
                                            x_all,
                                            data_likelihood,
                                            prior_parameters
                                            ){
    if(x_target<0){ return(-999999999999999999999999)}

    # names of randomeffects
    names_randomeffects <- data_likelihood$names_randomeffects
    names_sigmas <- data_likelihood$names_sigmas

    # get randomeffects
    epsilons_vector <- x_all[names_randomeffects]
    # get hyperpriors
    sigmas <- x_all[names_sigmas]
    
    # reshape ranomd effects into matrix
    EPS <- matrix(epsilons_vector, 2, ncol=n,byrow=TRUE)
    
    # density of random-intercepts
    logdensity_0 <- sum(dnorm(EPS[1,], mean=0, sd=sigmas[0], log=TRUE))
    # density of random-slopes
    logdensity_1 <- sum(dnorm(EPS[2,], mean=0, sd=sigmas[1], log=TRUE))
    
    # half-t prior on sigmas
    logprior <- dstudent_t(
        x=x_target,
        mu=0,
        df=prior_parameters[['prior_df']], # 
        sigma=prior_parameters[['prior_sd']],
        log=TRUE
    )
    return(logdensity_0 + logdensity_1 + logprior)}


# EMPTY STORAGE CONTAINTER FOR POSTERIOR FUNCTIONS
# ... each list has a named-entry for every variable
list_of_log_posteriors <- setNames(vector(mode='list',length=length(x.init)),nm=names(x.init))
x.lowb <- x.uppb <- setNames(vector(mode='numeric', length=length(x.init)),nm=names(x.init)) # storage containter for upper/lower bounds
list_of_priors <- setNames(vector(mode='list',length=length(x.init)),nm=names(x.init)) # storage containter for prior parameters
w.slice <- setNames(vector(mode='numeric', length=length(x.init)),nm=names(x.init)) # slice parameter

# FIXED EFFECTS - load posteriors, priors, boundaries, and w.slice
list_of_log_posteriors[names_fixedeffects] <- list(log_posterior_REregression_dnorm_prior)
x.lowb[names_fixedeffects] <- -10 # reasonable lower bounds
x.uppb[names_fixedeffects] <- 10 # reasonable upper bounds
w.slice[names_fixedeffects] <- 0.1 # slice-step size initial guess
list_of_priors[names_fixedeffects] <- list(list(
        'prior_mean'=0, # normal prior mean, on fixed effects beta0 and beta1
        'prior_sd'=1, # normal prior sd, on fixed effects beta0 and beta1
    ))

# HYPERPRIORS (sigma)- load posteriors, priors, boundaries, and w.slice
# sigma0 and sigma1 control the variance of random-effects
# - sigma0: root-variance of random-effect intercepts
# - sigma1: root-variance of random-effect slopes
list_of_log_posteriors[names_sigmas] <- list(log_posterior_dnormlike_halfStudentt_hyperprior)
x.lowb[names_sigmas] <- 0 # lower bound >0
x.uppb[names_sigmas] <- 5 # upper bound on sigmas
w.slice[names_sigmas] <- c(0.1,0.03) # step size initial guess
list_of_priors[names_sigmas] <- list(
    'sigma0'=list('prior_df'=3,'prior_sd'=0.05), # half-student t hyperprior
    'sigma1'=list('prior_df'=3,'prior_sd'=0.001) # # half-student t hyperprior
)

# RANDOM EFFECTS - load posteriors, priors, boundaries, and w.slice
# - by the definition of 'random effect', the prior parameters on the random
# effects are themselves variables included as targest of slice-sampling
# therefore, the "prior parameters" are just an index to grab the correct variable (sigma0 or sigma1) from the vector of variables being sampled
list_of_log_posteriors[names_randomeffects] <- list(log_posterior_REregression_dnorm_REprior)
x.lowb[names_randomeffects] <- -10 # lower bound on random-effects
x.uppb[names_randomeffects] <- 10  # upper bound on random-effects
w.slice[names_randomeffects] <- 0.1 # step-size initial guess
list_of_priors[names_randomeffects] <- list(c(
    rep(list(which_sigma_idx=1),n), # prior is just an index to grab sigma from x
    rep(list(which_sigma_idx=2),n) # prior is just an index to grab sigma from x
))

# TEST: fixed effects the functions: ensure no erros
for(x.nm in names(x.init.beta)){
    posterior_func <- list_of_log_posteriors[[x.nm]]
    print(posterior_func(
        x_target = x.init[[x.nm]],
        x_all = x.init,
        data_likelihood = data_likelihood,
        prior_parameters = list_of_priors[[x.nm]]
    ))
}

# test the random 
for(x.nm in names(x.init.eps)){
    posterior_func <- list_of_log_posteriors[[x.nm]]
    print(posterior_func(
        x_target = x.init[[x.nm]],
        x_all = x.init,
        data_likelihood = data_likelihood,
        prior_parameters = list_of_priors[[x.nm]]
    ))
}

# test the sigmas
for(x.nm in names(x.init.sigmas)){
    posterior_func <- list_of_log_posteriors[[x.nm]]
    print(posterior_func(
        x_target = x.init[[x.nm]],
        x_all = x.init,
        data_likelihood = data_likelihood,
        prior_parameters = list_of_priors[[x.nm]]
    ))
}

# RUN THE SLICE SAMPLE
out <- slice.sample(
    x.init, # initial estimates of variables
    list_of_log_posteriors=list_of_log_posteriors, # list of posterior-densities
    data_likelihood=data_likelihood, # y data and model-matrices
    prior_parameters_list=list_of_priors, # hyperparameters for priors
    nslice=5000, # number of samples to collect
    thin = 5, # thinning the MCMC
    x.lowb=x.lowb, # sensible lower bounds to all variables
    x.uppb=x.uppb, # sensible upper bounds to all variables
    w=w.slice, # W slice step
    m=12, # steps to shift W
    pass.counter=2000, # rejection sampling fallback n-steps
    w_auto_adjust=TRUE, # whether to auto-adjust w
    w_auto_adjust_factor=0.5, # auto-adjustment factor
    print_interval=50  # print interval                       
)
# the slice sampler mixes well: 
# ... watch the w values: they should coverge and not be crazy-high

# get samples from Slice Samples
ss.samp <- out$samples

# Better mixing than JAGS!
plot(mcmc(ss.samp),ask=TRUE) # ... beautiful mixing

# get posterior means and SE:
ss_estimates <- apply(ss.samp[,c('beta0','beta1','sigma0','sigma1')], 2, function(x){ c('ss-mean'=mean(x),'ss-se'=sd(x))})
print(ss_estimates)
#       beta0      beta1     sigma0      sigma1
#ss-mean 0.3112876 0.42600474 0.16347697 0.001090720
#ss-se   0.1008498 0.05173868 0.03150973 0.001232339

# compare to true values
print(c(beta0, beta1)) # 0.4054651 0.4000000



####################################
# PART II: JAGS POISSON RANDOM-EFFECTS COMPARISON
####################################

jags_model_syntax <- "model{
   # fixed effects (intercept and trend)
   beta[1] ~ dnorm(0,1) # intercept
   beta[2] ~ dnorm(0,1) # trend
   # hyperprior: half-normal prior on random-effects dispersion
   #sigma[1] ~ dnorm(0,pow(0.03,-2)) T(0,)
   #sigma[2] ~ dnorm(0,pow(0.03,-2)) T(0,)
   sigma[1] ~ dt(0,pow(prior_sigma0,-2), prior_df) T(0,)
   sigma[2] ~ dt(0,pow(prior_sigma1,-2), prior_df) T(0,)
   # likelihood
   for(i in 1:n){
      epsilon0[i] ~ dnorm(0,pow(sigma[1],-2))
      epsilon1[i] ~ dnorm(0,pow(sigma[2],-2))
      for(j in 1:T){
          Y[j,i] ~ dpois(exp(lambda[j,i]))
          lambda[j,i] <- beta[1] + (beta[2] + epsilon1[i])*(T-1) + epsilon0[i]
      }
   }
}"
cat(jags_model_syntax, file ='/tmp/jags_script.jags',append=FALSE)

# compile jags model
jags <- jags.model(
    '/tmp/jags_script.jags',
    data = list(
        "Y" = Y,
        "n"=n,
        "T"=T,
        "prior_df"=3,
        "prior_sigma0"=0.05,
        "prior_sigma1"=0.001
        ),
    n.chains = 1, n.adapt = 40000, 
)

# jags mixes poorly for random-effects models
plot(mcmc.samples,ask=TRUE)

# run jags
mcmc.samples <- coda.samples(jags,c('beta','sigma'),2000000,thin=100)

# get jags estimate
jags_estimates <- apply(mcmc.samples[[1]], 2, function(x){ c('jags-mean'=mean(x),'jags-se'=sd(x))})
# compare jags estimates to true values
print(jags_estimates) # 
#            beta[1]   beta[2]   sigma[1]    sigma[2]
#jags-mean 0.1046751 0.3185370 0.41466436 0.001118042 
#jags-se   0.9649312 0.3222891 0.07132055 0.001363403 ** EXTREME HIGH VARIANCE
print(c(beta0, beta1)) # 0.4054651 0.4000000



######################################
# CONCLUSIONS: 
# - the Flexible R Slice Sampler has much better MCMC mixing that the JAGS random effects.
# - the Flexible R Slice Sampler more accurate estimates than JAGS for random-effects, also with lower variance.
# - JAGS works well for lots of data and fixed-effects, but not for random-effects with difficult distributions (like sparse count)
