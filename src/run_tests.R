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
                         pass.counter=1000,  #
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

    ## CHECK 1: classes of arguments
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
        stop('`x.init` should be a numeric vector of initial guesses of the values of the variables you are trying to sample')
    }
    
    check <- try(stopifnot(class(data_likelihood)=='list'))
    if(class(check) == "try-error"){ stop('`data_likelihood` should be list') }

    check <- try(stopifnot(class(prior_parameters_list)=='list'))
    if(class(check) == "try-error"){ stop('`prior_parameters_list` should be list of prior parametres') }
    
    check <- try(stopifnot(as.integer(nslice)==nslice))
    if(class(check) == "try-error"){ stop('`nslice` should be an integer, number of samples to take of x') }

    check <- try(stopifnot(as.integer(thin)==thin))
    if(class(check) == "try-error"){ stop('`thin` should be an integer, number of MCMC steps in between accepting a sample') }

    check <- try(stopifnot(class(x.lowb)=='numeric'))
    if(class(check) == "try-error"){ stop('`x.lowb` should be a numeric vector, the sensible lower bounds of values for each variable in x') }    

    check <- try(stopifnot(class(x.uppb)=='numeric'))
    if(class(check) == "try-error"){ stop('`x.uppb` should be a numeric vector, the sensible upper bounds of values for each variable in x') }    

    check <- try(stopifnot(class(w)=='numeric'))
    if(class(check) == "try-error"){ stop('`w` should be a numeric vector, the stepwise per parametres') }

    check <- try(stopifnot(as.integer(m)==m))
    if(class(check) == "try-error"){ stop('`m` should be an integer, the number of steps to move a stepsize of `w` to find boundaries of density') }

    check <- try(stopifnot((m<=20 & m>=8)))
    if(class(check) == "try-error"){ warning(sprintf('Good values for `m` tend to be in the range of 10-12, please double check your value for `m`. Proceeding with value %d',m)) }
    
    check <- try(stopifnot(w_auto_adjust_factor>0 & w_auto_adjust_factor<1.0))
    if(class(check) == "try-error"){ stop('`w_auto_adjust_factor` should be >0 and <1') }

    ## CHECK 2: LENGTHS of variables
    check <- try(stopifnot(length(x.init)==length(list_of_log_posteriors)))
    if(class(check) == "try-error"){ stop('`x.init` and `list_of_log_posteriors` should have the same length (an entry for each variable to estimate)')}

    check <- try(stopifnot(length(x.init)==length(prior_parameters_list)))
    if(class(check) == "try-error"){ stop('`x.init` and `prior_parameters_list` should have the same length (an entry for each variable to estimate)')}

    check <- try(stopifnot(length(x.init)==length(x.lowb)))
    if(class(check) == "try-error"){ stop('`x.init` and `w.lowb` should have the same length (an entry for each variable to estimate)')}

    check <- try(stopifnot(length(x.init)==length(x.uppb)))
    if(class(check) == "try-error"){ stop('`x.init` and `w.uppb` should have the same length (an entry for each variable to estimate)')}    

    check <- try(stopifnot(length(x.init)==length(w)))
    if(class(check) == "try-error"){ stop('`x.init` and `w` should have the same length (an entry for each variable to estimate)')}

    ## check for NAs
    check <- try(stopifnot(all( all(all(!is.na(x.init)), all(!is.na(x.uppb)), all(!is.na(x.lowb)), all(!is.na(w))))))
    if(class(check) == "try-error"){ stop('NAs in x.init or x.lowb or x.uppb. Please set with real values') }

    ## check min and max bounds on x
    check <- try(stopifnot( all(x.lowb < x.uppb) ))
    if(class(check) == "try-error"){ stop('`x.lowb` should be much lower than `x.uppb`')}

    check <- try(stopifnot( all(x.init <= x.uppb) ))
    if(class(check) == "try-error"){ stop('`x.init` should be less than or equal to `x.uppb`')}

    check <- try(stopifnot( all(x.lowb <= x.init) ))
    if(class(check) == "try-error"){ stop('`x.lowb` should be less than or equal to `x.init`')}

    ## check slice sizes
    check <- try(stopifnot( all(w < (x.uppb - x.lowb))))
    if(class(check) == "try-error"){ stop('`w` should be much lower than `x.uppb-x.lowb`')}
    
    ## check slice sizes
    check <- try(stopifnot( all(w < (x.uppb - x.lowb)/2)))
    if(class(check) == "try-error"){
        is_high_values_of_w <- which(w >= (x.uppb - x.lowb)/2)
        nm_haigh_values_of_w <- names(is_high_values_of_w)
        if(!is.null(nm_haigh_values_of_w)){
            offending_variables <- paste0(nm_haigh_values_of_w,collapse=',')
            warning(sprintf('Values for `w` seem TOO LARGE for these variables:%s. Please reduce either `w` or check x.lowb and x.upp have a large enough spread.',offending_variables))
        } else {
            offending_variables <- paste0(is_high_values_of_w,collapse=',')
            warning(sprintf('Values for `w` seem TOO LARGE for these elements:%s. Please reduce either `w` or check x.lowb and x.upp have a large enough spread.',offending_variables))
        }
    }

    ## Check the argument names of the log-posteriors
    ## Ensure each log-postior has argument names: x_target,x_all,data_likelihood,prior_parameters
    names_of_expected_arguments <- c('x_target','x_all','data_likelihood','prior_parameters')
    which_misnamed <- numeric(0)
    for(i in 1:length(x.init)){
        names_of_arguments <- names(formals(list_of_log_posteriors[[i]]))        
        check <- try(stopifnot( all(names_of_expected_arguments%in%names_of_arguments) ))
        if(class(check) == "try-error"){ which_misnamed<-c(which_misnamed, i)}
    }
    if(length(which_misnamed)>0){
        stop(sprintf("Unexpected argument names for elements %s in `list_of_log_posteriors`: expected %s", paste(which_misnamed,collapse=','), paste(names_of_expected_arguments,collapse=', ')))
    }

    ## Check each posterior function at x.init
    points_to_check <- x.init
    which_failed <- c()
    for(i in 1:length(points_to_check)){
        check <- try( out<-list_of_log_posteriors[[i]]( x_target=points_to_check[[i]], x_all=points_to_check, data_likelihood=data_likelihood, prior_parameters=prior_parameters_list[[i]] ))
        if( (class(check) == "try-error") | is.na(out) | is.null(out) ){ which_failed <- c(which_failed, i); print(sprintf('list_of_log_posteriors[[%d]] fails to compute at x.init[[%d]]: please check',i,i)) }
    }
    if (length(which_failed)){
        stop(sprintf("Cannot compute `list_of_log_posteriors` at values in `x.init`. Check elements: %s", paste(which_failed,collapse=',')))
    }

    ## Check each posterior function at x.lowb
    points_to_check <- x.lowb
    which_failed <- c()
    for(i in 1:length(points_to_check)){
        check <- try( out<-list_of_log_posteriors[[i]]( x_target=points_to_check[[i]], x_all=points_to_check, data_likelihood=data_likelihood, prior_parameters=prior_parameters_list[[i]] ))
        if( (class(check) == "try-error") | is.na(out) | is.null(out) ){ which_failed <- c(which_failed, i); print(sprintf('list_of_log_posteriors[[%d]] fails to compute at x.lowb[[%d]]: please check',i,i)) }
    }
    if (length(which_failed)){
        stop(sprintf("Cannot compute `list_of_log_posteriors` at values in `x.lowb`. Check elements: %s", paste(which_failed,collapse=',')))
    }

    ## Check each posterior function at x.uppb
    points_to_check <- x.uppb
    which_failed <- c()
    for(i in 1:length(points_to_check)){
        check <- try( out<-list_of_log_posteriors[[i]]( x_target=points_to_check[[i]], x_all=points_to_check, data_likelihood=data_likelihood, prior_parameters=prior_parameters_list[[i]] ))
        if( (class(check) == "try-error") | is.na(out) | is.null(out) ){ which_failed <- c(which_failed, i); print(sprintf('list_of_log_posteriors[[%d]] fails to compute at lower-bounds x.lowb[[%d]]: please check',i,i)) }
    }
    if (length(which_failed)){
        stop(sprintf("Cannot compute `list_of_log_posteriors` at values in upper-bounds `x.uppb`. Check elements: %s", paste(which_failed,collapse=',')))
    }
    
    ## SUCCESS exit
    print('SUCCESS: Passed checks on inputs and priors. Be sure to set `do_checks=FALSE` for future runs on the same data and priors.')
    
}
# TODO check that arguments of each log-posterior function are x_target, x_all, except    
    
                             
                             

