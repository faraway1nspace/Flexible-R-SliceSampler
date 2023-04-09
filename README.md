# Flexible-R-SliceSampler

Just another Slice Sampler for Bayesian Sampling -- but highly efficient and flexible in pure R (i.e., no more weird JAGS/BUGS/STAN synatx ). 

Why would you want a pure-R Gibbs Sampler (and specifically a "Slice Sampler") when JAGS/BUGs exists? We used in Rankin et al 2020 for a complex Bayesian analysis of Dugongs which could be run in JAGS... here was missing data, complex mixture priors, and a hierarchical process that required a custom Gibbs sampler. 

The Slice Sampler (based on Neal 2003) worked so well, I decided to open-source it for others

Be sure to check out our examples in `demos/`. And be sure to read above the key-parameters below.


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

### Reasons NOT to use Flexible-R-SliceSampler (vs. JAGS/BUGS)
- You have a simple model with simple data that is easily run in JAGS/BUGS
- You're unfamiliar with likelihoods and have difficulty making R-likelihood functions

## Syntax Comparison to JAGS

Let's run a Zero-Inflated poisson model, with poisson variable `lambda` and zero-inflation variable `psi`. Let's say the data is simply: `y <- c(1,0,0,0,10,0,3,0,0,0,0,0,0,30)`. 

In JAGS, the model syntax would be:

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
      
      # poisson likelihood
      y[i] ~ dpois(lambda1[i])
   }
}
"
```

### Syntax in Flexible-R-SliceSampler

We will now re-write the above JAGS ZIPpoisson model as log-posterior **R-functions**. 

In particular, for each variable (`lambda` and `psi`), we must write a log-posterior function that computes the likelihood and the prior (logged). The posteriors can be **unnormalized**, making it relatively simple to compute (basically, we just the sum of the log-likelihood and log-priors).

Here are the log-posterior functions for `lambda` and `psi`, which correspond to the above JAGS/BUGS model. Notice that we make use of native R density functions like `dnorm`, `dpois`, etc.:

```R
# log posterior with log-normal density on lambda
log_posterior_ziplambda <- function(x_target,
                                    x_all,
                                    data_likelihood,
                                    prior_parameters
                                    ){
    # data for likelihood: y values
    y <- data_likelihood$y

    # (naive) loglikelihood for zip poisson (works, but primariy for readability
    likelihoods_poisson_unconditional <- dpois(y,lambda=x_all['lambda'],log=FALSE)
    likelihoods_zip <- (1-x_all['psi'])*(y==0) + x_all['psi']*likelihoods_poisson_unconditional
    loglike <- sum(log(likelihoods_zip))

    # normal prior on log-lambda
    log_prior <- dnorm(
        log(x_target),
        mean=prior_parameters[['prior_lambda_mean']],
        sd=prior_parameters[['prior_lambda_sigma']], # note the sqrt(1/tau) for JAGS-compatibility,
        log=TRUE
    )
    return(loglike + log_prior)}
	
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
    log_prior <- dbeta(
        x_target,
        shape1=prior_parameters[['prior_psi_a']],
        shape2=prior_parameters[['prior_psi_b']],        
        log=TRUE
    )
    return(loglike + log_prior)}
	

# collect log-posterior functions in a named-list
list_of_log_posteriors = list(
    "lambda" = log_posterior_ziplambda,
	"psi" = log_posterior_zippsi
)	
```
Each study will have it's own likelihoood and priors, so you'll need to custom-write these functions anew each time. This may seem cumbersome, but it is *flexible*.

Each log-posterior density function must have the same arguments:
- `x_target`: the candidate value of the variable for that posterior function
- `x_all`: a named numeric variable with _all_ variables in a vector, needed to compute the likelihood.
- `data_likelihood`: a named-list with all the data necessary to compute the likellihood (like the  JAGS argument `data`)
- `prior parameters`: a named-list with the prior-parameters for each variable to sample


The log-posterior functions are collected in a named-list `list_of_log_posteriors`, with an entry for each variable to sample (lambda and psi). The list is passed as an argument `slice.sample` function. 

Read below for more about the other arguments of `slice.sample`. Or, have a look at the following demo files.



## Files

- `src/flexible_slice_sampler.R` - source code with flexible slice sampler
- `demo/demo_multivariate_regression_with_student-t_priors.R` - script to demonstrate sparse Bayesian regression using Student-T priors -- in comparison with the Lasso (l1-regularization)
- `demo/demo_jags_vs_slice_zeroInflatedPoisson.R` - script to compare JAGS, for a Zero-Inflated Poisson model


## Slice Sampling Overview.

Neal (2003) describes a rejection-free method of MCMC sampling for complex multivariate joint-posteriors. _Slice-Sampling also has relatively few "tunable" hyperparameters, unlike some Metropolis-Hastings-like methods that fail to converge if the hyperparameters are poorly set.

**The most important hyperparameter is the STEP-SIZE parametre `W`** -- it is used to laterally "step" arond the univariate-parameter interval of `x` to find good regions of the posterior-density to sample `x`. You'll need to read Neal 2003 to get a better understanding.

In Flexible-R-SliceSampler, we have a good method of automatically adjusting the step-size `W` parameters. However, **the `W` parameters must be carefully watched** to ensure they converge on acceptable values, and the sampler doesn't hang, choke, or mix badly. 

Usually, if the `W` step-sizes don't converge, there is something _else_ wrong with the model.

### Slice Sampling in Brief

Given a point $x$ and posterior-density function $f(x)$, we:
1. sample $z$ uniformly along vertical from 0 to f(x): $z\sim\mathcal{U}(0,f(x))$
2. find the left-most point ($L$) at which $f(L) = z$
3. find the right-most point ($R$) at which $f(R) == z$
4. sample a new $x_{\text{new}}$ uniformly between $L$ and $R$ : $x_{\text{new}}\sim\mathcal{U}(L,R)$
5. If $f(x_{\text{new}}) > z$ then accept $x = x_{\text{new}}$, and move-on to next variable
6. If not then "shrink" the boundaries L and R:  
    - If $x_{\text{new}} < x$, then $L \leftarrow x_{\text{new}}$ 
	- Elif $x_{\text{new}} > x$, then $R \leftarrow x_{\text{new}}$ 
    - Repeat steps 4 to 6 until $f(x_{\text{new}}) > z$
	
	
Notice there is some hand-waiving for the steps 2 & 3 ("find the left-most point ($L$) and right-most point $L$ where $f=z$. This is the key step called by Neal the "stepping-out procedure. Basically, we move $x$ outward by a step-size $W$ until $f(x-m*W) < z$. The following picture from Wikipedia does a great job describing the stepping-out procedure, and the "shrinkage" procedure of step 6.

**The critical point is that the hyperparameter $W$ must be carefully set**. This is the key thing the user must initially set.

```note
We have found that the best step-size W is the long-run average of R-L over the entire density. Therefore, we can monitor R and L, and slowly adjust W to the long-run average of R-L
```



## Citation

If you use this slice-sampler, please cite the following study:

```
Rankin RW, Marsh H. 2020. Technical Appendices: 8-11 in Marsh H, Collins K. Grech A., Miller R and Rankin RW. 2020. *An assessment of the distribution and abundance of dugongs and in-water, large marine turtles along the Queensland coast from Cape York to Hinchinbrook Island.* A report to the Great Barrier Reef Marine Park Authority, April 2020.
```

The original Slice Samplier paper by Neal 2003 should be cited as: 

```
Neal, R. M. 2003. Slice sampling. The Annals of Statistics 31:705–767.
```
