# demonstration of how to do a multivariate regression via Slice Sampling
#
# In particular, how to specify the log-densities (log-ikelihood + log-prior)
# To induce 'sparsity' we'll use a student-T prior
# We'll compare it to a Ridge Regression (l2-regularization) and the Lasso (l1-regularization)

# import Slice Sampler
source("../src/flexible_slice_sampler.R")

library(glmnet) # to compare against Ridge Regression
library(coda) # for MCMC visuals
#### IMPORT DEMO DATA

# load sparse example
data(QuickStartExample)
x <- QuickStartExample$x;
y <- QuickStartExample$y
X_df = data=cbind(y=y,as.data.frame(x)) # X as data.frame

# least squares (for comparison to Bayesian)
maxlike_model = glm(y~., data=X_df)
# fit a ridge regression l2-regularization (for comparison to Bayesin)
ridge = glmnet(x, y, alpha=0, lambda=2) # ridge regression
# fit a lasso regression l1-regularization (for comparison to Bayesin)
lasso = glmnet(x, y, alpha=1, lambda=0.1) # lasso regression

# in this demo, we'll see how Student-T priors induce a variable-slection like the Lasso


############################################################
### Set-Up Slice Sampler

#### DATA 

# data for likelihood (response variable and X-matrix aka model-matrix)
X_model_matrix = model.matrix(y~., data=X_df)

data_likelihood <- list(
    y = y, # y data
    mm = X_model_matrix # X-data, aka model-matrix
)

#### PRIORS
prior_parameters_list <- list()

# prior on intercent (uniform prior with min and max)
prior_parameters_list[["(Intercept)"]]<-list('min'=-10,'max'=10)

# priors on beta regression parameters (with student-t parameters df and sigma)
# these encourage sparsity i.e. a Bayesian type of regularization
for(beta_i in colnames(data_likelihood$mm)){
    if(beta_i!='(Intercept)'){
        # long-tailed student-t distribution: peaked at zero
        prior_parameters_list[[beta_i]] <- list('tau'=0.02, 'df'=1)
    }
}

# prior on root-Variance sigma (uniform prior with min=0 and max
prior_parameters_list[["sigma"]]<-list('min'=0.0000001,'max'=10)

### INITIAL-VALUES 
# like all MCMC routines, we need relatively decent initial values to start the MCMC chains

n_par <-length(prior_parameters_list)
x.init <- rep(0, n_par) # initialize to zero
names(x.init)<-names(prior_parameters_list)

# use MLE estimates as initial estimates
for(beta_mle in names(coef(maxlike_model))){ x.init[[beta_mle]] <- coef(maxlike_model)[beta_mle] }
x.init['sigma']<-sigma(maxlike_model)

#### POSTERIOR DENSITIES (unnormalized)
# the great thing about an R-based Bayesian slice-sampler is that all posteriors can be 
#   specified as pure-R density functions, like dnorm(..., log=TRUE) dt(,...), etc.
# Also, the posteriors (likelihood * dprior) can be unnormalized, meaning we don't need a 
#   tractable normalizing constant, freeing us to specify models we way we want to
# 
# The slice sampler requires a log-postior for each variable to sample.
# Each posterior is a function with the same form
# `x_target` : float, candidate value of the target variable
# `x_all`: numeric vector, all other variables necessary for likelihood calculation
# `data_likelihood`: list, with (x,y)-data needed for likelihood calculation
# `prior_parameters`: list, with entries of parameters specifying the prior distribution on x_target

# list of posteriors (one function for each variable)
list_of_log_posteriors <- list()

# posterior for regression with student-t priors on regression betas
log_posterior_regression_with_student_t_prior <- function(x_target,
                                                    x_all,
                                                    data_likelihood,
                                                    prior_parameters
                                                    ){
    # default density: a regression with student-t priors
    mm <- data_likelihood$mm
    # data for likelihood: y values
    y <- data_likelihood$y

    # variable-names
    beta_names <- colnames(mm)
    betas <- x_all[beta_names]
    sigma <- x_all['sigma']
    
    # log likelihood
    log_likes <- dnorm(y,mean=mm%*%betas,sd=sigma, log=TRUE)
    log_like <- sum(log_likes)
    # prior on x_target: a student t distribution
    log_prior <- log_studentt(
        x_target, df=prior_parameters[['df']], mu=0,
        sigma=prior_parameters[['tau']]
    )
    return(log_like+log_prior)}

# posterior for regression with uniform prior (for Intercept)
log_posterior_regression_with_uniform_prior <- function(x_target,
                                                    x_all,
                                                    data_likelihood,
                                                    prior_parameters
                                                    ){
    # expects prior parameters called 'min' and 'max' for dunif
    # variables
    beta_names <- colnames(data_likelihood$mm)
    betas <- x_all[beta_names]
    sigma <- x_all['sigma']
    # log likelihood
    log_like <- sum(dnorm(
        data_likelihood$y,
        mean=data_likelihood$mm%*%betas,
        sd=sigma,
        log=TRUE
    ))
    # prior
    log_prior <- dunif(
        x_target,
        min=prior_parameters[['min']],
        max=prior_parameters[['max']],
        log=TRUE
    )
    return(log_like+log_prior)}
 
# insert the log-posteriors into list, one per parameters
# WARNING: must be in the same order as `x_init`

# log-posterior of intercept
list_of_log_posteriors[['(Intercept)']] <- log_posterior_regression_with_uniform_prior
# log-posteriors of beta-regressoin parameters
for(beta_param in paste0("V",1:20)){
    list_of_log_posteriors[[beta_param]] <- log_posterior_regression_with_student_t_prior
}
# log-posterior of sigma 
list_of_log_posteriors[['sigma']] <- log_posterior_regression_with_uniform_prior

####################
# OTHER PARAMETERS FOR SLICE (upper bounds, lower bounds)
x.uppb <- rep(30,n_par)
x.lowb <- rep(-30,n_par); x.lowb['sigma']<-0.0000001

# W is a key hyperparameter governing the efficiency of the Slice Sampler see details
w.slice <- rep(0.40, n_par)



##########
# SAMPLE

samples_mcmc <- slice.sample(
    x.init, # initial estimates of variables
    list_of_log_posteriors, # list of log-posterior densities per variable
    data_likelihood, # y data and model-matrices
    prior_parameters_list, # hyperparameters for priors
    nslice=4000, # number of slices
    x.lowb, # lower safety bounds on variables
    x.uppb, # upper safety bounds on variables
    w=w.slice, # W hyperparameter governing Slice Sampler (see Details)
    m=10, # number of steps of W (see Details)
    )

# inspect the MCMC time-series (no-rejection-sampling necessary)
plot(mcmc(samples_mcmc$samples))

# inspect the W parameters
# 

# get estimates of the posteriors
beta_hats_bayesian <- apply(
    samples_mcmc$samples,2,function(x){
        c(
            'mean'=mean(x),
          'median'=median(x),
          'se'=sd(x),
          'lcl95'=quantile(x, 0.025,names=FALSE),
          'ucl95'=quantile(x, 0.025,names=FALSE))
    })

# PLOT: compare Ridge Regression to Bayesian Student-T priors
beta_names <- paste0('V',1:20)
plot(x=coef(ridge)[beta_names,], y=beta_hats_bayesian['mean',beta_names],
     main='Student-T priors vs Ridge Regression',
     xlab='Ridge estimates', ylab='Bayesian Estimates'
)
abline(0,1)

# PLOT: compare Lasso Regression to Bayesian Student-T priors
plot(x=coef(lasso)[beta_names,], y=beta_hats_bayesian['mean',beta_names],
     main='Student-T priors vs Lasso Regression',
     xlab='Lasso estimates', ylab='Bayesian Estimates'
)
abline(0,1)

# PLOT: compare to MLEs Bayesian Student-T priors
plot(x=coef(maxlike_model)[beta_names], y=beta_hats_bayesian['mean',beta_names],
     main='Student-T priors vs Ridge Regression',
     xlab='MLE', ylab='Bayesian Estimates'
)
abline(0,1)


## CONCLUSION: 
