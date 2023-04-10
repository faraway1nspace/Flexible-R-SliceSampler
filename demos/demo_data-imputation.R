# Data Imputation - demo of interweaving other processes between Slice-Sampling

# One of the benefits of a pure-R MCMC-sampler is that we can easily interwave
# other native-R functions inbetween calls to sampling from the posterior.

# We'll do data-imputation of missing data to demonstrate the interwoven functions:
# we'll do i) one round of slice-sampling the posteriors, then ii) one round of
# random data-imputation. Theoretically, this is called "Monte Carlo" integration.

# SCENARIO: a common problem with getting field-data in ecology (like counts of
# animals) is uncertainty. We've encountered situations where some of the recorded
# counts are not integers by "min, best-guess, max" values.

# We'll do a poisson regression, but some of our count-data will be a mix of
# integers and min/best-guess/max.

set.seed(42)

source("../src/flexible_slice_sampler.R")


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
y <- rpois(n, lambda=exp(mm%*%betas_true))

#######################
# SIMULATE MISSING DATA










