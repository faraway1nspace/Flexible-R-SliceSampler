
log_posterior_poisson_regression_betas <- function(x_target,
                                                   x_all,
                                                   data_likelihood,
                                                   prior_parameters
                                                   ){
    # poisson regression with normal-density on betas
    # expects the prior_parmaters to have named-entries 'beta_mean' and 'beta_sigma'
    # the model-matrix (aka X-data) should be an element of `data_likelihood`
    
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


           
