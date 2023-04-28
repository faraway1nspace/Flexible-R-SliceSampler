#########
# Check data and arguments prior to running `slice.sample` core algorithm

check_args_and_data <- function(x.init, # initial estimates of variables
                         list_of_log_posteriors, #log.densities.list, per variable in x
                         data_likelihood, # y data and model-matrices
                         prior_parameters_list, # hyperparameters for priors, per variable in x
                         nslice=2000, # number of samples to collect
                         thin=1, # thinning (accept only every `thin` sample, discard rest 
                         x.lowb=rep(-10^6, length(x.init)), # sensible lower bounds on variables
                         x.uppb=rep(10^6, length(x.init)),  # sensible upper bounds on variables
                         w=rep(3.0,length(x.init)),
                         m=10, # steps to shift W
                         pass.counter=1000,
                         w_auto_adjust=TRUE, # whether to auto-adjust w
                         w_auto_adjust_factor=0.8, # auto-adjustment factor
                         print_interval=100  # print interval                       
                         ){
    #' A Native R Slice sampler based on "Neal RM 2003 "Slice Sampling". doi:10.1214/aos/1056562461
    #' x.init: numeric vector of initial estimates
    #' data: y-data and model-matrix for ingestion by likelihood function
    #' list_of_log_posteriorss: a list of (log) density functions for each parameter in x, typically the log-prior + log-likelihood. See details
    #' prior_parameters_list: list of hyperparameters 
    #' nlice: int, number of slices to sample / number of iterations
    #' thin: int, number of steps to run sampler and discard samples, improves mixing
    #' x.lowb: numeric float, lower-bounds of x.init
    #' x.uppb: numeric float, upper-nounds of x.init
    #' w: slice step-size parameter: the most important hyperparameter to guide behaviour: see Details
    #' m: default number of steps of step-size W. see Slice Sampling Paper by Neal 2003
    #' pass.counter: last-ditch rejection sampling if things go badly
    #' print_interval: how often to print w-slices estimates
    #' w_auto_adjust: bool, whether to auto-adjust the w-slice intervals
    #' w_auto_adjust_factor: float, smoothing-factor to mix new w estimates with past estimates
    
    # number of parameters to sample
    npar <- length(x.init)

                                        # CHECK 1: classes of arguments
    check <- try(stopifnot(class(list_of_log_posteriors)=='list'))
    if(class(check) == "try-error"){
        print(sprintf('`list_of_log_posteriors` should be a list of functions, not %s',class(list_of_log_posteriors)))
    }
    
    for(i in range(length(list_of_log_posteriors))){
        check <- try(stopifnot(class(list_of_log_posteriors[[i]])=='function'))
        if(class(check) == "try-error"){
            print(sprintf('`list_of_log_posteriors[[%d]]` should be a function', i))
        }
    }
    
    check <- try(stopifnot(class(x.init)=='numeric'))
    if(class(check) == "try-error"){
        print('`x.init` should be a numeric vector of initial guesses of the values of the variables you are trying to sample')
    }
    
    check <- try(stopifnot(class(data_likelihood)=='list'))
    if(class(check) == "try-error"){ print('`data_likelihood` should be list') }

    check <- try(stopifnot(class(prior_parameters_list)=='list'))
    if(class(check) == "try-error"){ print('`prior_parameters_list` should be list of prior parametres') }
    
    check <- try(stopifnot(as.integer(nslice)==nslice))
    if(class(check) == "try-error"){ print('`nslice` should be an integer, number of samples to take of x') }


    check <- try(stopifnot(as.integer(thin)==thin))
    if(class(check) == "try-error"){ print('`thin` should be an integer, number of MCMC steps in between accepting a sample') }

    check <- try(stopifnot(class(x.lowb)=='numeric'))
    if(class(check) == "try-error"){ print('`x.lowb` should be a numeric vector, the sensible lower bounds of values for each variable in x') }    

    check <- try(stopifnot(class(x.uppb)=='numeric'))
    if(class(check) == "try-error"){ print('`x.uppb` should be a numeric vector, the sensible upper bounds of values for each variable in x') }    

    check <- try(stopifnot(class(w)=='numeric'))
    if(class(check) == "try-error"){ print('`w` should be a numeric vector, the stepwise per parametres') }

    check <- try(stopifnot(as.integer(m)==m))
    if(class(check) == "try-error"){ print('`m` should be an integer, the number of steps to move a stepsize of `w` to find boundaries of density') }
    
}
# TODO check that arguments of each log-posterior function are x_target, x_all, except    
    
                             
                             

