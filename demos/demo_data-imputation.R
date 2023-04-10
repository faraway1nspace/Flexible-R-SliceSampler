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

# the imputed data will NOT be missing-at random. Like in ecology, there will be 
# more missingness with higher-counts

# sample missing data ~ counts (higher-counts more likely missing)
idx_missing <- sort(sample(1:n, size=4,prob=sqrt(y)))

# n x 4 matrix including min/best/max values
y_minbestmax <- matrix(NA, n,4,dimnames=list(1:n, c('count','min','best','max')))
y_minbestmax[,'count'] <- y
y_minbestmax[idx_missing,'count']<- NA # missing data

# add min max best values for y
for(idx_miss in idx_missing){

}







