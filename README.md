# Flexible-R-SliceSampler

Just another Slice Sampler for Bayesian Sampling -- but highly efficient and flexible in pure R (no more JAGS/BUGS/STAN). 

Why would you want a pure-R Gibbs Sampler (and specifically a "Slice Sampler") when JAGS/BUGs exists? We used in Rankin et al 2020 for a complex Bayesian analysis of Dugongs that wouldn't work in JAGS. 

## Motivation: Why use Flexible-R-SliceSampler (vs. JAGS or BUGS)
You should use this flexibl R-based Slice Sampler for Bayesian analysis if:
- **No 'rejection' sampling** - don't waste compute-resources throwing out samples a la Metropolis-Hasting.
- **Matrix Operations** - these are cumbersome in JAGS but can be easy and optimized in R  
- **Novel Distributions** - if you need a distribution that is not supported in JAGS/BUGS/STAN  
    - whatever probability distribution you code in R, you can use in this Slice Sampler
- **Intractible Priors** - you have an intractible Prior distribution (and likewise intractible posterior) that can't be normalized  
    - for example, one of our priors was a _mixture distribution_ based on the outputs from  *another* independent model. You can pipe-in any external info you want, so long as you can code-it in R and can sample from it with a `dfunction`-type of log-density.
- **Penalities on Derived Quantitis** - you want to calculate derived quantities, and place priors on them
    - you want to penalize some quantity that is a *derivative* of your process-parameters
	- e.g., in ecology's mark-recapture, we often have prior information on derivatives of our statistical models (like we know the population-size isn't 0 nor is it 10000) but these quantities are derivatives, often *not* direct parameters in our likelihood.
	- regulation can also be thought-of as a penalty on derviative quantities (the $\ell_1$-norm of parameters). See the demos.
- **Missing Data Imputation** - we can do subsampling or random-inputation, per MCMC iteration.  
    - changing the data is illegal in JAGS, but here we can flexibly subsample or impute the data for each MCMC iteration.  
	- because everything is R-based, you can run the slice sampler for a few iterations, then run a subsampling/imputation step outside the function, and then continue with the slice-sampler. This R-based Gibbs-routine opens-up lots of possibilities.
- **Interweave Complex Processes** -you want to interweave complex R-processes in between vanilla MCMC  
    - e.g. let's say a part of your Gibbs-sampling depends on a Hidden-Markov-Model process that is unavailable in JAGS/BUGS, like sampling of HMM-states. You can sample from these states outside of the slice sampler using conventional HMM R-code, then run the Slice Sampler conditional on those states as "data", and repeat.
- **BUGS/JAGS is annoying** - why code in BUGS/JAGS if you can do it better & faster in R?

### Why you Shouldn't use Flexible-R-SliceSampler (vs. JAGS or BUGS)
- You don't understand likelihoods and can't compute them yourself
- You have a simple model that is easily run in JAGS 

## Syntax Comparison to JAGS

Let's run a Zero-Inflated poisson model, with poisson variable `lambda` and zero-inflation variable `psi`. Let's say the data is simply: `y <- c(1,0,0,0,10,0,3,0,0,0,0,0,0,30)`. In JAGS, the model syntax would be:

```R
jags_model_syntax <- "model{
   # prior on psi
   psi ~ dbeta(prior_psi_a,prior_psi_b)
   
   # prior on (log)lambda
   tau <- 1/(prior_lambda_sigma * prior_lambda_sigma)
   loglambda ~ dnorm(prior_lambda_mean, tau)
   lambda <- exp(loglambda)

   # likelihood
   for(i in 1:N){
      
      # zero-inflation process:
      zip[i] ~ dbern(psi)
      
      # toggle-able lambda based on zero-inflation
      lambda1[i] <- lambda*zip[i] + 0.000001 
      # .. NOTICE THE CONSTANT 0.000001 that is necessary to stabilize the ZIPpoisson!
      
      # poisson likelkhood
      y[i] ~ dpois(lambda1[i])
   }
}
"
```

Flexible-R-SliceSampler, 




## Files

- `src/flexible_slice_sampler.R` - source code with flexible slice sampler
- `demo/demo_multivariate_regression_with_student-t_priors.R` - script to demonstrate sparse Bayesian regression using Student-T priors -- in comparison with the Lasso (l1-regularization)
- `demo/demo_jags_vs_slice_zeroInflatedPoisson.R` - script to compare JAGS, for a Zero-Inflated Poisson model


## Citation

If you use this slice-sampler, please cite the following study:

```
Rankin RW, Marsh H. 2020. Technical Appendices: 8-11 in Marsh H, Collins K. Grech A., Miller R and Rankin RW. 2020. *An assessment of the distribution and abundance of dugongs and in-water, large marine turtles along the Queensland coast from Cape York to Hinchinbrook Island.* A report to the Great Barrier Reef Marine Park Authority, April 2020.
```
