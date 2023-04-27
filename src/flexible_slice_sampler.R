
# default functions
slice.sample <- function(x.init, # initial estimates of variables
                         list_of_log_posteriors, #log.densities.list, per variable in x
                         data_likelihood, # y data and model-matrices
                         prior_parameters_list, # hyperparameters for priors, per variable in x
                         nslice=2000,
                         thin=1,
                         x.lowb=rep(-10^6, length(x.init)),
                         x.uppb=rep(10^6, length(x.init)),
                         w=rep(3.0,length(x.init)),
                         m=10, # steps to shift W
                         pass.counter=1000,
                         w_auto_adjust=TRUE, # whether to auto-adjust w
                         w_auto_adjust_factor=0.8, # auto-adjustment factor
                         print_interval=100  # print interval                       
                         ){
    #' Gibbs Slice sampler, native to R
    #' x.init: numeric vector of initial estimates
    #' data: y-data and model-matrix for ingestion by likelihood function
    #' list_of_log_posteriors: a list of (log) density functions for each parameter in x, typically the log-prior + log-likelihood. See details
    #' prior_parameters_list: list of hyperparameters 
    #' nlice: int, number of slices to sample / number of iterations
    #' x.lowb: numeric float, lower-bounds of x.init
    #' x.uppb: numeric float, upper-nounds of x.init
    #' w: window parameters: the most important hyperparameter to guide behaviour: see Details
    #' m: see Slice Sampling Paper
    #' pass.counter: last-ditch rejection sampling if things go badly
    #' print_interval: how often to print w-slices estimates
    #' w_auto_adjust: bool, whether to auto-adjust the w-slice intervals
    #' w_auto_adjust_factor: float, smoothing-factor to mix new w estimates with past estimates
    
    # number of parameters to sample
    npar <- length(x.init)

    # store the samples
    samples_mcmc <- matrix(0,nslice,npar)
    colnames(samples_mcmc) <- names(list_of_log_posteriors)

    # initialize `x` with initial values of `x.init`
    x <-x.init # x is our current estimate of the variables

    # monitor the slice size: this is import to tune/auto-adjust the values of `w`
    w_monitor_slice_size <- numeric(npar) # monitor the slice-size

    # loop through n-slices/samples
    for(i in 1:nslice){

        for(j in 1:thin){

        # loop through each parameter p: we sample univariately, conditional on other parameters
            for(p in 1:npar){ # 

                # target density for this parameter p
                log_density_for_p = list_of_log_posteriors[[p]]
            
                # get (posterio) density at x=x
                log_f_at_x <- log_density_for_p(
                    x_target=x[[p]],                
                    x_all=x,
                    data_likelihood,
                    prior_parameters_list[[p]]
                )
            
                # sample random slice height < f(x) y.r
                random_slice_height <- log_f_at_x - rexp(1,1)

                # sample window size (adjustment factor to grow L)
                uv <- runif(2,0,1)

                # R and L are boundaries within-which we sample a new value for x
                # ... we initiaize R and L to the most recent values of ALL parameters
                L <- R <- x # notice we condition on the other parameters in

                # move p's L to the left by w[p]*uv[1] (i.e., move outside density)
                L[p] <- pmax( x[p]-w[p]*uv[1], x.lowb[p]);
                # log-density at L[p]
                f_at_L <- log_density_for_p(
                    x_target=L[[p]],
                    x_all=L,
                    data_likelihood,
                    prior_parameters_list[[p]]
                )

                # move p's R to the right by the window size w[p]
                R[p] <- pmin(L[[p]]+w[[p]], x.uppb[[p]])
                J <- floor(m*uv[2]) # J steps to move L
                K <- (m-1)-J # K steps to move R

                # GOAL: adjust L[p] and R[p] to be at the boundards of the log-density
                # step 1: incrementally move L[p] outside of density until f_at_L < slice_height
                while(
                    (random_slice_height < f_at_L) & (J>0) & (L[p]>x.lowb[p])
                ){
                    # decrease L[p] by w[p] 
                    L[p] <- pmax(L[p]-w[p], x.lowb[p])
                    # continue for J steps
                    J <- J-1
                    # new log-density at L
                    f_at_L <- log_density_for_p(
                        x_target=L[p],
                        x=L,
                        data_likelihood,
                        prior_parameters_list[[p]]
                    )
                } # end when f_at_L < slice_height or J==0

                # step 2: incrementally move R[p] outside of density until f_at_R < slice_height
                f_at_R <- log_density_for_p(
                    x_target=R[[p]],
                    x_all=R,
                    data_likelihood,
                    prior_parameters_list[[p]]
                )

                # while loop to adjust R[p] by W[p]
                while(
                    (random_slice_height < f_at_R) & (K>0) & (R[p]<x.uppb[p])
                ){
                    # increment R[p] by w[p] (i.e., push rightward until outside of density
                    R[p] <- pmin(R[p]+w[p], x.uppb[p])
                    # decrement K
                    K <- K-1
                    # new log-density at R
                    f_at_R <- log_density_for_p(
                        x=R, x_target=R[p], data_likelihood, prior_parameters_list[[p]]
                    )
                } # end when f_at_R < slice_height or K==0

                # Step 3: sample between L and R
                xstar <- x 
                # insert new sample at p: sample xstar[p] from Uniform(L[p], R[p])
                xstar[p] <- runif(1,L[p],R[p])

                # log-denisty at the new x-star
                f_at_xstar <- log_density_for_p(
                    x_target=xstar[[p]],
                    x_all=xstar,
                    data_likelihood,
                    prior_parameters_list[[p]]
                )
                # we accept x if log_density_at_x is > random_slice_height            
                accept_f_at_xstar <- (f_at_xstar >= random_slice_height)
                
                # if f is not > random_slice_height, then we adjust L and R until it is
                pass.counter.i <- pass.counter
                pass <-  accept_f_at_xstar & pass.counter.i>0
                # if not pass, we micro-adjust L and R
                while(!pass){

                    # new L is position of X[p]
                    L[p] <- xstar[p]*(xstar[p]<x[p])+ L[p]*(xstar[p]>=x[p])
                    # new R is position of X[p]
                    R[p] <- xstar[p]*(xstar[p]>x[p])+ R[p]*(xstar[p]<=x[p])
                    
                    # sample new x[p] from between L and R
                    xstar[p] <- runif(1,L[p],R[p])

                    # log-density at start
                    f_at_xstar <- log_density_for_p(
                        x_target=xstar[p],
                        x_all=xstar,
                        data_likelihood,
                        prior_parameters_list[[p]]
                    )
                    # we accept x if log_density_at_x is > random_slice_height            
                    accept_f_at_xstar <- (f_at_xstar >= random_slice_height)
                                        # decrement pass.counter.i
                    pass <- accept_f_at_xstar & pass.counter.i>0
                    pass.counter.i <- pass.counter.i-1
                    if(pass.counter.i<0){stop(sprintf('failed to get xstar; something wrong with %d',p))}
                } # done

                # update xstar to x
                x <- xstar

                # monitor our slice of R[p]-L[p] (we'll use it to dynamically adjust w)
                w_monitor_slice_size[p] <- R[p]-L[p]
                
            } # thin
            
        } # p
        # done loopting through all parameters p

        # accept x into our storage container
        samples_mcmc[i,] <- x

        # update values of W
        if(w_auto_adjust){
            # moving average that stabilizes with increasing i
            q<-w_auto_adjust_factor
            w <- q^(8.0/i)*w + (1-q^(8.0/i))*w_monitor_slice_size
        }
        
        # print monitor the w
        if (i%%print_interval==0){
            # here, we print the W-slice widths. These should stabilize
            print(sprintf("Iter:%d/%d; W=%s", i,nslice, paste(round(w,3),collapse=',')))
        }
    }
    return(
        list(
            samples=samples_mcmc,
            w=w_monitor_slice_size
        ))
}


# DEMO log density function
logdensity_univariatemean_norm <- function(
                                           x_target, # target variable to sample (here mean)
                                           x_all, # all other variables necessary for likelihood fucntion
                                           data_likelihood, # y data for the likelihood function
                                           prior_par # prior parmeters
                                           ){
    # likelihood and prior for a univariate 
    loglike <- sum(dnorm(
        data_likelihood$y,
        mean = x_target,
        sd= x_all[2],
        log=TRUE
    ))

    # prior on mean
    logdensity_prior <- dnorm(
        x_target,prior_par$mu,
        prior_par$disperion,
        log=TRUE
    )
    return (loglike + logdensity_prior)}


# student t distribution
log_studentt<- function(x, df=1, mu=0, sigma=1){
    return( log(1)-log(sigma) + dt((x - mu)/sigma, df, log=TRUE) )
}
