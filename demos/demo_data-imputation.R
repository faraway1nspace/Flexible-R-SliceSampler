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

source("../src/flexible_slice_sampler.R")
source("../src/example_log-posteriors.R"

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
                                                   prior_parameters,
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
        mean=prior_parameters['beta_mean'],
        sd=prior_parameters['beta_sigma'],
        log=TRUE
    )
    return(log_like+log_prior)}

##################################
# PRIORS:
list_prior_parameters <- list()

# prior on intercept
list_prior_parameters["(Intercept)"] <- list('beta_mean'=0, 'beta_sigma'=10)
# prior on regression coefficients
