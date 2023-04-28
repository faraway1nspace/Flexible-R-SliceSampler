# Demo Files for `slice.sample`

I've provided some interesting examples how to specify log-posterior density functions, and use the slice sampler.

- `demos/demo_jags_vs_slice_zeroInflatedPoisson.R` - Here, we benchmark the slice-sampler to JAGS, using a Zero-Inflated Poisson model.
- `demos/demo_multivariate_regression_with_student-t_priors.R` - A multiple-regression example. In particular, I show how the Bayesian slice-sampled estimates compare to MLE, as well as $\ell_1$ and $\ell_2$-regularized estimates (aka Lasso and Ridge-regression estimates). A fascinating result is that Bayesian regression with long-tailed Student-t priors results in a very similar estimates as the Lasso ($\ell_1$-regularization). I.e., the Student-t priors induce a type of automatic variable-selection and shrinkage almost identical to lasso, but all in a Bayesian framework!
- `demos/demo_random-effects.R` - Random Effects aka Hierarchical Bayes model, compared to JAGS. Here, we have individual random effects for a Poisson regression, including random intercepts and random slopes. Our simple slice sampler mixes much better than JAGS

More to come.

