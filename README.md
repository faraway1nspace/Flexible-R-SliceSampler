# Flexible-R-SliceSampler

Just another Slice Sampler for Bayesian Sampling -- but highly efficient and flexible in pure R (no more JAGS/BUGS/STAN). 

Why would you want a pure-R Gibbs Sampler (and specifically a "Slice Sampler") when JAGS/BUGs exists? We used in Rankin et al 2020 for a complex Bayesian analysis of Dugongs that wouldn't work in JAGS. 

You should use this flexibl R-based Slice Sampler for Bayesian analysis if:
- you need a distributions not in JAGS/BUGS/STAN  
- you have an intractible Prior distribution (and likewise intractible posterior)  
    - for example, our prior was a mixture distribution based on *another* model
- you want interweave complex R-processes in between vanilla MCMC  
    - e.g. let's say a part of your Gibbs-sampling depends on a Hidden-Markov-Model process that is unavailable in JAGS/BUGS
- you want to calculate derived quantities, and place priors on them
    - e.g., lets say you want to penalize some quantity that is a *derivative* of your process-parameters
- **Missing data imputation** - like want do to subsampling or random-inputation per MCMC iteration  
    - changing the data is illegal in JAGS, but here we can flexibly subsample or impute the data for each MCMC iteration -- 

## Citation

If you use this slice-sampler, please cite the following study:

```
Rankin RW, Marsh H. 2020. Technical Appendices: 8-11 in Marsh H, Collins K. Grech A., Miller R and Rankin RW. 2020. *An assessment of the distribution and abundance of dugongs and in-water, large marine turtles along the Queensland coast from Cape York to Hinchinbrook Island.* A report to the Great Barrier Reef Marine Park Authority, April 2020.
```
