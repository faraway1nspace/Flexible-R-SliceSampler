# demonstration of how to do a multivariate regression via Slice Sampling
#
# In particular, how to specify the log-densities (log-ikelihood + log-prior)
# To induce 'sparsity' we'll use a student-T prior
# We'll compare it to a Ridge Regression

# import Slice Sampler
source("../src/flexible_slice_sampler.R")

library(glmnet) # to compare against Ridge Regression

# load sparse example
data(QuickStartExample)
x <- QuickStartExample$x;
y <- QuickStartExample$y

# least squares (for comparison)
ls = glm(y~., data=cbind(y=y,as.data.frame(x)))
# fit a ridge regression 
ridge = glmnet(x, y, alpha=0, lambda=5) # ridge regression

#########
# slice sampler set-up: unnormalized posterior

log_posterior_regression_with_student_t <- function(x_target,
                                                    x_all,
                                                    data_likelihood,
                                                    prior_parameters
                                                    ){
    # default density: a regression with student-t priors


    # data for likelihood: model matrix
    mm <- data_likelihood$mm
    # data for likelihood: y values
    y <- data_likelihood$y

    # variable-names
    beta_names <- names(mm)
    betas <- x_all[beta_names]
    sigma <- x_all['sigma']
    
    # log likelihood
    log_like <- dnorm(y,mean=mm%*%betas,sd=sigma, log=TRUE)
    # prior
    log_prior <- log_studentt(
        x_target,
        df=prior_parameters[['df']],
        mu=0,
        sigma=prior_parameters[['tau']])
    return(log_like+log_prior)}

    
