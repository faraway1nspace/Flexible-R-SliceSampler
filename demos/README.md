# Demo Files for `slice.sample`

I've provided some interesting examples how to specify log-posterior density functions, and use the slice sampler.

- `demo/demo_jags_vs_slice_zeroInflatedPoisson.R` - Here, we benchmark the slice-sampler to JAGS, using a Zero-Inflated Poisson model.
- `demo/demo_multivariate_regression_with_student-t_priors.R` - A multiple-regression example. In particular, we should how the Bayesian slice-sampled estimates compare to MLE, as well as $\ell_1$ and $\ell_2$ regression. A fascinating result is that Bayesian regression with long-tailed Student-t priors results in a very similar estimates as the Lasso ($\ell_1$-regularization). I.e., with Student-t priors, we can induce a type of automatic variable-selection and shrinkage almost identical to lasso, but all in a Bayesian framework!

More to come.

