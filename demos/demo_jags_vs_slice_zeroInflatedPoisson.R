# A simple comparison of JAGS syntax vs the Flexible R Slice Sampler

source("../src/flexible_slice_sampler.R")

library(coda)
library(rjags)

# zero-inflated count data (we'll model with a Zero-inflated poisson
y <- c(1,0,0,0,10,0,3,0,0,0, 0,0,0,30)

#######
# JAGS

# jags syntax for zero inflated count data
jags_model_syntax <- "model{
   # prior on (log)lambda
   loglambda ~ dnorm(prior_lambda_mean, prior_lambda_tau)
   lambda <- exp(loglambda)
   # prior on psi
   psi ~ dbeta(prior_psi_a,prior_psi_b)
  
   # likelihood
   for(i in 1:N){
      
      # zero-inflation process:
      zip[i] ~ dbern(psi)
      
      # toggle-able lambda based on zero-inflation
      lambda1[i] <- lambda*zip[i] + 0.0001
      
      # poisson likelkhood
      y[i] ~ dpois(lambda1[i])

   }
}
"
cat(jags_model_syntax, file ='/tmp/jags_script.jags',append=FALSE)

# compile jags model
jags <- jags.model(
    '/tmp/jags_script.jags',
    data = list("y" = y, "N"=length(y),
                "prior_lambda_mean"=2.9,
                "prior_lambda_tau"=1/10,
                "prior_psi_a"=1, "prior_psi_b"=1
                ),
    n.chains = 1, n.adapt = 100
)
# run jags
mcmc.samples <- coda.samples(jags,c('psi','lambda'),4000)

# get jags estimate
jags_estimates <- apply(mcmc.samples[[1]], 2, function(x){ c('jags-mean'=mean(x),'jags-se'=sd(x))})


##############
# Flexible-R-Slice-Sampler

#### DATA
data_likelihood <- list(y = y)

#### PRIORS
prior_parameters_list <- list(
    # lambda poisson process
    'lambda' = c("prior_lambda_mean"=2.9, "prior_lambda_tau"=1/10),
    # psi zero-inflation process
    'psi' = c( "prior_psi_a"=1, "prior_psi_b"=1)
)

# Instead of a JAGS model{}, we specify our posteriors in an R function
# .. with arguments x_target, x_all, data_likelihood, and prior_parameters
list_of_log_posteriors <-list('lambda'=NA, 'psi'=NA)

# log posterior with log-normal density on lambda
log_posterior_ziplambda <- function(x_target,
                                    x_all,
                                    data_likelihood,
                                    prior_parameters
                                    ){
    # data for likelihood: y values
    y <- data_likelihood$y

    # (naive) loglikelihood for zip poisson
    likelihoods_poisson_unconditional <- dpois(y,lambda=x_all['lambda'],log=FALSE)
    likelihoods_zip <- (1-x_all['psi'])*(y==0) + x_all['psi']*likelihoods_poisson_unconditional
    loglike <- sum(log(likelihoods_zip))

    # normal prior on log-lambda
    log_prior <- dnorm(
        log(x_target),
        mean=prior_parameters[['prior_lambda_mean']],
        sd=sqrt(1/prior_parameters[['prior_lambda_tau']]), # note the sqrt(1/tau) for JAGS-compatibility,
        log=TRUE
    )
    return(loglike + log_prior)}


log_posterior_ziplambda(2.3, c('lambda'=2.3, 'psi'=0.5), data_likelihood, prior_parameters_list[['lambda']])


# log-posterior with Beta prior on psi
log_posterior_zippsi <- function(x_target,
                                 x_all,
                                 data_likelihood,
                                 prior_parameters
                                 ){
    # data for likelihood: y values
    y <- data_likelihood$y

    # (naive) loglikelihood
    likelihoods_poisson_unconditional <- dpois(y,lambda=x_all['lambda'],log=FALSE)
    likelihoods_zip <- (1-x_all['psi'])*(y==0) + x_all['psi']*likelihoods_poisson_unconditional
    loglike <- sum(log(likelihoods_zip))
    # logprior density on psi with Beta prior
    log_prior <- dnorm(
        x_target,
        shape1=prior_parameters[['prior_beta_a']],
        shape2=prior_parameters[['prior_beta_b']],        
        log=TRUE
    )
    return(likelihood_zip + log_prior)}
    
