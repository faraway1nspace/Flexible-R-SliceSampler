# Data Imputation - demo of interweaving other processes between Slice-Sampling

# One of the benefits of a pure-R MCMC-sampler is that we can easily interwave
# other native-R functions inbetween calls to sampling from the posterior.

# We'll do data-imputation of missing data to demonstrate the interwoven functions:
# we'll do i) one round of slice-sampling the posteriors, then ii) one round of
# random data-imputation. Theoretically, this is called "Monte Carlo" integration.
# Importantly, we can tailor our favourite data-imputation methods, or a mixture of many.

# SCENARIO: a common problem with getting field-data in ecology (like counts of
# animals) is uncertainty. We've encountered situations where some of the recorded
# counts are not integers by "min, best-guess, max" values.

# We'll do a poisson regression, but some of our count-data will be a mix of
# integers and min/best-guess/max.

set.seed(42)
library(coda)
source("../src/flexible_slice_sampler.R", chdir=TRUE) # load sampler
source("../src/example_log-posteriors.R") # import the log-posteior

###############
# SET-UP X-DATA

n<-20
# true betas
beta_names = c('x1','x2','x3')
betas_true <- c(
    '(Intercept)'=log(10), 'x1'=-1, 'x2'=0,'x3'=1
)
# dataframe: x1,x2,x3
x_data <- as.data.frame(setNames(lapply(1:3,function(x){ runif(n,-2,2)}),beta_names))
# model matrix
mm <- model.matrix(~x1+x2+x3, data=x_data)
# true y variable
y_true <- rpois(n, lambda=exp(mm%*%betas_true))

#######################
# SIMULATE MISSING DATA

# the imputed data will NOT be missing-at random. Like in ecology, there will be 
# more missingness with higher-counts

# sample missing data ~ counts (higher-counts more likely missing)
idx_missing <- sort(sample(1:n, size=4,prob=y_true))
print(y_true[idx_missing])

# n x 4 matrix including min/best/max values
y_minbestmax <- matrix(NA, n,4,dimnames=list(1:n, c('count','min','best','max')))
y_minbestmax[,'count'] <- y_true
y_minbestmax[idx_missing,'count']<- NA # missing data

# add min max best values for y
for(idx in idx_missing){
    # random min best max values centred around the y-true
    y_minbestmax[idx,c('min','best','max')] <- sort(c(
        max(0, sample(round(seq(y_true[idx]*0.80,y_true[idx]*0.95)),size=1)),
        max(0, sample(round(seq(y_true[idx]*0.90,y_true[idx]*1.1)),size=1)),
        max(0, sample(round(seq(y_true[idx]*1.1,y_true[idx]*1.3)),size=1))
    ))
}

print('done simulating the data')

####################
# IMPUTATION-FUNCTION

# truncated rpoisson: sample poisson centred on lambda, between [min,max]
rpois_truncated <- function(lambda, min, max){
    # rpois, truncated
    x <- seq(min,max,1)
    x_density <- dpois(x, lambda=lambda)
    return(sample(x, size=1, prob = x_density/sum(x_density)))}

# mixture imputation
impute <- function(min, best, max, p=0.5){
    # complex imputation function: mixture of 'best' and a rpois(truncated)
    best_or_rpois <- runif(1) # toggle between
    # just return the best-guess
    if(best_or_rpois < p){ return(best) }
    # truncated poission imputation: sampled on 'best'
    return( rpois_truncated(lambda=best, min=min, max=max) )}


y_imputed <- y_minbestmax
for(idx in idx_missing){
    y_imputed[idx,'count'] <- impute(
        min=y_minbestmax[idx,'min'],
        best=y_minbestmax[idx,'best'],
        max=y_minbestmax[idx, 'max'],
        0.3
    )
}
print(y_imputed)

################################
# log posterior: poisson regression

log_posterior_poisson_regression_betas <- function(x_target,
                                                   x_all,
                                                   data_likelihood,
                                                   prior_parameters
                                                   ){
    # poisson regression with normal-density on betas
    # expects the prior_parmaters to have named-entries 'beta_mean' and 'beta_sigma'
    
    # variables
    beta_names <- colnames(data_likelihood$mm) # model matrix colunm-names
    betas <- x_all[beta_names]
    
    # log likelihood
    mu <- exp(data_likelihood$mm%*%betas)
    log_like <- sum(dpois(data_likelihood$y, lambda=mu, log=TRUE))
    
    # prior
    log_prior <- dnorm(
        x_target,
        mean=prior_parameters[['beta_mean']],
        sd=prior_parameters[['beta_sigma']],
        log=TRUE
    )
    return(log_like+log_prior)}


list_log_posteriors <- list(
    '(Intercept)'=log_posterior_poisson_regression_betas,
    'x1'=log_posterior_poisson_regression_betas,    
    'x2'=log_posterior_poisson_regression_betas,
    'x3'=log_posterior_poisson_regression_betas
)
    
##################################
# PRIORS:
list_prior_parameters <- list()
# prior on intercept
list_prior_parameters[["(Intercept)"]] <- list('beta_mean'=log(10), 'beta_sigma'=10)
# prior on regression coefficients
for(beta_name in beta_names){
    list_prior_parameters[[beta_name]] <- list('beta_mean'=0, 'beta_sigma'=1)
}

##################
# INITIAL VALUES FOR SAMPLER
x.init <- setNames(rep(0,ncol(mm)), nm=colnames(mm))
x.init['(Intercept)'] <- log(mean(y_minbestmax[,1],na.rm=TRUE))

###################
# PARAMETERS FOR SLICE SAMPLER

n_mcmc <- 4000
w <- setNames(c(0.5,0.5,0.5,0.5), nm=names(x.init))
x.uppb <- c(log(100),10,10,10)
x.lowb <- c(log(0.001),-10,-10,-10)

########################
# MCMC: WITH DYNAMIC DATA-IMPUTATION INTERWOVEN IN MCMC
# Two steps: i) impute new data; ii) sample posteriors

mcmc_samples <- matrix(NA, n_mcmc, length(x.init), dimnames=list(1:n_mcmc, names(x.init)))

# MCMC
x_star <- x.init
for(j in 1:n_mcmc){

    # step 1: data imputation : sample min/best/max
    y_imputed <- y_minbestmax
    # loop through missing values to imput
    for(idx in idx_missing){

        # impute y using min/best/max function
        y_imputed[idx,'count'] <- impute(
            min=y_minbestmax[idx,'min'],
            best=y_minbestmax[idx,'best'],
            max=y_minbestmax[idx, 'max'],
            0.3)
    }

    # (dynamic) data with imputed values
    data_likelihood <- list(y = y_imputed[,'count'],mm = mm)

    # step2: slice-sample
    slice_samps <- slice.sample(
        x_star, # initial estimates of variables
        list_log_posteriors, # list of log-posterior densities per variable
        data_likelihood, # y data and model-matrices
        list_prior_parameters, # hyperparameters for priors
        nslice=4, # number of slices
        thin=1, # thin
        x.lowb=x.lowb, # lower safety bounds on variables
        x.uppb=x.uppb, # upper safety bounds on variables
        w=w, # W hyperparameter governing Slice Sampler (see Details)
        m=12, # number of steps of W (see Details)
        do_checks=FALSE
        )
    
    # update x_star with last sample from slice-sampler
    x_star <- tail(slice_samps$samples, 1)[1,]

    # simple take the smoothed mean of w
    w_star <- slice_samps$w
    w <- w*0.8 + 0.2*w_star # moving average of w
    # print the w to monitor its converge
    if(j%%20==0){ print(as.numeric(w))}
    
    # store MCMC sample
    mcmc_samples[j,]<- x_star
}

# inspect converge
plot(mcmc(mcmc_samples))

# get posterior means and compare to true
print(colMeans(mcmc_samples))
print(betas_true)

print(apply(mcmc_samples,2,sd))
# 95%CI
print(apply(mcmc_samples,2,function(x){quantile(x,c(0.025,0.975))}))


###########
# COMPARE TO STATIC-IMPUTATION 
# how do the values compare to just imputing the "best" values only


y_imputed <- y_minbestmax
for(idx in idx_missing){ y_imputed[idx,'count'] <- y_minbestmax[idx,'best']}
data_likelihood <- list(y = y_imputed[,'count'],mm = mm)

slice_samps <- slice.sample(
        x.init, # initial estimates of variables
        list_log_posteriors, # list of log-posterior densities per variable
        data_likelihood, # y data and model-matrices
        list_prior_parameters, # hyperparameters for priors
        nslice=4000, # number of slices
        x.lowb, # lower safety bounds on variables
        x.uppb, # upper safety bounds on variables
        w=w, # W hyperparameter governing Slice Sampler (see Details)
        m=12 # number of steps of W (see Details)
)

# posterior means of static-imputation
print(colMeans(slice_samps$samples)) # static
# compare to dynamic imputation
print(colMeans(mcmc_samples)) # dynamic

# posterior variance/sd
print(apply(slice_samps$samples,2,sd)) # static
print(apply(mcmc_samples,2,sd)) # dynamic



