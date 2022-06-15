# make sure that gauseR is loaded
library('gauseR')

name_match <- function(fixed = NULL,
                       Nsp = NULL,
                       spname = NULL) {
  #' @title name_match
  #' 
  #' @description Used in the gause_wrapper_fixed function to match provided
  #' fixed model parameters to model coefficients
  #' 
  #' @param fixed A named vector with named elements, indicating model
  #'  parameters to be fixed. Intrinsic growth rates are denoted as r1, r2, etc. 
  #'  and species interaction coefficients are denoted as a11, a12, a22, etc.
  #' @param Nsp An integer indicating the number of species present. Please 
  #' don't use more than 10 species!!!
  #' @param spname A character vector indicating the species names
  
  # create empty coefficient matrix to match structure in the lm
  coef_matrix <- matrix(nrow = Nsp + 1, ncol = Nsp)
  rownames(coef_matrix) <- c('(Intercept)', paste0('lag_abund', spname))
  colnames(coef_matrix) <- spname
  
  # get the names and split into r and a and 1 and 2 etc.
  fixed_names <- names(fixed) # get provided names
  params <- substring(fixed_names, 1, 1) # get letters i.e. r or a
  ids <- substring(fixed_names, 2) # get numbers i.e. 1, 2, 11, 12, 21, 22, etc.
  
  # intrinsic growth rate parameter locations
  loc_r <- (1:(Nsp*2))[params %in% 'r'] # which fixed parameters are r?
  r_sp <- as.numeric(ids[loc_r]) # get their location
  coef_matrix['(Intercept)',r_sp] <- fixed[loc_r] # put values into matrix
  
  # intrinsic growth rate parameter locations
  loc_a <- (1:(Nsp*2))[params %in% 'a'] # which fixed parameters are a?
  a_sp1 <- as.numeric(substring(ids[loc_a], 1, 1)) # first species id 
  a_sp2 <- as.numeric(substring(ids[loc_a], 2, 2)) # second species id
  # need to do it this way so matrix doesn't fill in whole rows with duplicates
  for (i in 1:length(loc_a)) {
    coef_matrix[a_sp1[i] + 1, a_sp2[i]] <- fixed[loc_a[i]] # values into matrix
  }
  
  return(coef_matrix)
} 

# function to optimise in the log-space ripped straight from gauseR and modified a bit
lv_log_optim <- function (pars, opt_data, parm_signs, odefun = lv_interaction_log) {
  #' @title lv_log_optim
  #' 
  #' @description A slightly modified copy of the lv_optim from gauseR to 
  #' optimise the data in the log-space rather than normal space
  #' 
  #' @param pars A vector of parameter values in log space to be optimized. 
  #' Must include a logged starting abundance for each species, followed by the
  #' logged absolute values of the growth rates, followed by the logged 
  #' absolute value of the elements of the interaction matrix
  #' @param opt_data Abundance data for optimization. Must include one column 
  #' labeled 'time' with time steps, and a column for each species abundance.
  #' @param parm_signs A vector that provides the desired sign of each parameter 
  #' (i.e. -1 or 1). If value is zero, then the term is held at zero (but 
  #' should be left out of the pars vector).
  #' @param odefun The function to use to simulate the ODE - defaults to 
  #' lv_interaction_log
  #' 
  #' @return Squared error between model fits for given parameter values 
  #' and observations
  
  nsp <- ncol(opt_data) - 1
  logn0 <- pars[1:nsp]
  logparms <- pars[-c(1:nsp)]
  n0 <- exp(logn0)
  times <- opt_data$time
  par_opt_use <- numeric(length(parm_signs))
  par_opt_use[parm_signs != 0] <- exp(logparms) * parm_signs[parm_signs != 
                                                               0]
  out <- deSolve::ode(y = log(n0), times = times, func = odefun, 
                      parms = par_opt_use)
  predictions <- (out[, -1, drop = FALSE])
  observations <- log(opt_data[, -1, drop = FALSE])
  for (i in 1:ncol(observations)) {
    mo <- mean(observations[, i], na.rm = T)
    predictions[, i] <- predictions[, i]/mo
    observations[, i] <- observations[, i]/mo
  }
  return(sum((predictions - observations)^2, na.rm = T))
}

# function definition
gause_wrapper_fixed <- function(time, 
                                species, 
                                fixed = NULL, 
                                logfit = FALSE, 
                                ...) {
  #' @title gause_wrapper_fixed
  #' 
  #' @description Like the gause_wrapper function in the gauseR package, this 
  #' function is used to automatically obtain lotka-volterra model parameters 
  #' with certain settings and with the ability to keep certain parameters fixed
  #' as well as optimise according to the log-transformed data. 
  #' 
  #' @param time Vector of time steps corresponding to observations in species 
  #' data.frame. This can be numeric, Date, or POSIXct.
  #' @param species A data.frame with one column per species to be fitted. 
  #' Note - column names cannot include white spaces or non-standard 
  #' special characters.
  #' @param fixed A named vector with named elements, indicating model
  #'  parameters to be fixed. Intrinsic growth rates are denoted as r1, r2, etc. 
  #'  and species interaction coefficients are denoted as a11, a12, a22, etc.
  #' @param logfit A logical indicating whether the fitting should be optimised 
  #' using the normal-space time-series data or the log-space time-series data.
  #' 
  #' @return A list with simulated time series (out), paramter estimates 
  #' (parameter_intervals), optimizer output (optout), and raw data used for 
  #' fitting (rawdata).
  
  # check validity of provided fixed parameters
  if (!is.null(fixed)) { # if anything provided for fixed
    if (!is.vector(fixed)) { # check that it's a vector
      stop('provide a named vector for fixed parameters')
    } else if (is.null(names(fixed))) { # check that they have names
      stop('provide a named vector for fixed parameters')
    } 
  }
  
  # check validity of provided species data
  if (is.null(dim(species))) { # if only one species (i.e. it is a vector)
    if (is.null(names(species))) { # if it doesn't have a name
      nm <- 'N1' # give it a name
    } else { # otherwise
      nm <- names(species)[1] # use its given name
    }
    species <- data.frame(species) # turn into data frame
    colnames(species) <- nm # using species names decided earlier
  }
  
  # some basic info, transformed into more useful forms
  # how many species present?
  Nsp <- ncol(species)
  
  # get species names
  spname <- colnames(species)
  
  # create list of lagged time-series data
  laglst <- apply(species, 2, get_lag, time = time)
  # get values at 
  
  # create list + matrix of per-capita growth rates
  percap_lst <- lapply(laglst, function(i) {percap_growth(i$x, i$laggedx, i$dt)})
  percap_mat <- as.matrix(as.data.frame(percap_lst))
  
  # get lagged abundances from laglst
  lag_abund <- as.data.frame(sapply(laglst, `[`, 'laggedx'), col.names = spname)
  lag_abund <- as.matrix(lag_abund)
  
  # create variable name map to detect input parameters and keep them fixed
  if (!is.null(fixed)) {
    name_map <- c('Intercept', paste0('lag_abund', spname))
    name_mat <- matrix(NA, nrow = Nsp + 1, ncol = Nsp, dimnames = list(4:6,1:2))
    name_mat <- name_match(fixed, Nsp, spname)
  }
  
  # # create multivariate linear regression as initial guess for parameters
  # mod <- lm(percap_mat ~ 1+lag_abund) 
  # mod_coef <- c(coef(mod))
  # names(mod_coef) <- rownames(mod)
  
  # I think the solution is to manually make everything on the LHS and RHS using
  # paste and then we can specify the offsets individually
  # this necessitates several linear models for each species
  
 modList <- lapply(1:Nsp, function(i) {
    
    # create data frame to evaluate lm in later
    modDF <- list2DF(percap_lst)
    
    # get number of observations
    Nobs <- nrow(percap_mat) 
    
    # write out LHS
    
    list2env(percap_lst)#, envir = environment())
    LHS <- paste(spname[i], '~ ', collapse = '')
    
    # get offset data
    whichOffset <- !is.na(name_mat[,i])
    varNames <- row.names(name_mat)
    
    # get lag data
    lagList <- lag_abund %>% as.data.frame %>% as.list
    names(lagList) <- paste0('lag_abund', names(lagList))
    lagNames <- names(lagList)
    modDF <- cbind(modDF, list2DF(lagList))#, environment())
    
    # create formula matrix incl RHS elements
    formulaMat <- matrix(NA, nrow = length(whichOffset), ncol = 5)
    formulaMat[,1] <- ifelse(whichOffset, 'offset(', '')
    formulaMat[,2] <- c(paste('rep(1, ', Nobs, ')', collapse = ''), lagNames)
    formulaMat[,3] <- ifelse(whichOffset, ' * ', '')
    formulaMat[,4] <- ifelse(whichOffset, name_mat[whichOffset,i], '')
    formulaMat[,5] <- ifelse(whichOffset, ')', '')
    
    # collapse into single line formula
    modFormula <- paste(
      paste(LHS), # LHS
      paste(apply(formulaMat, 1 , paste0, collapse = ''), collapse=' + ')) # RHS
    
    # do lm
    mod <- lm(modFormula, modDF)
    
    # return
    return(mod)
  })
  
  # now we need to optimise to fit the growth data
  names(modList) <- spname

  # TODO CONTINUE THIS!!!
  
}
