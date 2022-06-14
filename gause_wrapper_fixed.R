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

# Demo data
fixed <- c(
  'a22' = gause_wrapper(
    subset(gause_1934_science_f02_03, Treatment == 'Pa')$Day,                     
    subset(gause_1934_science_f02_03, Treatment == 'Pa')$Volume_Species2, 
    doplot = FALSE)$parameter_intervals['a11','mu'],
  'r2' = gause_wrapper(
    subset(gause_1934_science_f02_03, Treatment == 'Pa')$Day,                     
    subset(gause_1934_science_f02_03, Treatment == 'Pa')$Volume_Species2, 
    doplot = FALSE)$parameter_intervals['r1','mu'],
  'a11' = gause_wrapper(
    subset(gause_1934_science_f02_03, Treatment == 'Pc')$Day,                     
    subset(gause_1934_science_f02_03, Treatment == 'Pc')$Volume_Species1, 
    doplot = FALSE)$parameter_intervals['a11','mu'],
  'r1' = gause_wrapper(
    subset(gause_1934_science_f02_03, Treatment == 'Pc')$Day,                     
    subset(gause_1934_science_f02_03, Treatment == 'Pc')$Volume_Species1, 
    doplot = FALSE)$parameter_intervals['r1','mu'])
  
time <- gause_1934_science_f02_03$Day
species <- cbind(
  'caudatum' = subset(gause_1934_science_f02_03, Treatment == 'Mixture')$Volume_Species1,
  'aurelia' = subset(gause_1934_science_f02_03, Treatment == 'Mixture')$Volume_Species2)
logfit = FALSE

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
  
  
  name_mat[]
  
  lv
  
  
}

function (time, species, N_starting = NULL, r_starting = NULL, 
          A_starting = NULL, doplot = TRUE, keeptimes = FALSE, parm_signs = NULL, 
          doopt = TRUE, ...) 
{
  
    lagged_abund <- sapply(laglst, function(dat) dat$laggedx)
    colnames(lagged_abund) <- colnames(species)
    modlst <- NULL
    for (i in 1:Nsp) {
      moddat <- data.frame(dNNdt = dNNdtlst[[i]], lagged_abund)
      modlst[[i]] <- eval(parse(text = paste("lm(dNNdt~", 
                                             paste(colnames(lagged_abund), collapse = "+"), 
                                             ", data=moddat)", sep = "")))
    }
    r_start <- unname(sapply(modlst, function(dat) {
      stats::coef(dat)["(Intercept)"]
    }))
    A_start <- unname(sapply(modlst, function(dat) {
      stats::coef(dat)[-1]
    }))
    if (!is.null(r_starting)) {
      r_start[!is.na(r_starting) & r_starting == 0] <- 0
    }
    if (!is.null(A_starting)) {
      A_start[!is.na(A_starting) & A_starting == 0] <- 0
    }
  }
  else {
    r_start <- r_starting
    A_start <- A_starting
  }
  parms <- c(r_start, A_start)
  if (is.null(N_starting)) {
    initialN <- species[which.min(time), ]
    for (i in 1:Nsp) {
      if (is.na(initialN[i]) || initialN[i] == 0) {
        initialN[i] <- mean(species[, i][species[, i] != 
                                           0], na.rm = T) * 0.01
      }
    }
  }
  else {
    initialN <- N_starting
  }
  if (doopt) {
    opt_data <- data.frame(time = time, species)
    if (is.null(parm_signs)) {
      parm_signs <- sign(parms)
    }
    pars <- c(log(initialN), log(abs(parms[parms != 0])))
    optout <- stats::optim(par = pars, fn = lv_optim, hessian = TRUE, 
                           opt_data = opt_data, parm_signs = parm_signs, ...)
    parms <- numeric(length(parm_signs))
    parms[parm_signs != 0] <- exp(optout$par[-c(1:length(initialN))]) * 
      parm_signs[parm_signs != 0]
    initialN <- exp(optout$par[1:Nsp])
  }
  else {
    optout <- NA
    parms <- c(r_start, A_start)
    initialN <- unlist(species[1, ])
  }
  if (keeptimes) {
    timessim <- time
  }
  else {
    timessim <- seq(min(time), max(time), length = 100)
  }
  out <- deSolve::ode(y = initialN, times = timessim, func = lv_interaction, 
                      parms = parms, ...)
  if (doplot) {
    graphics::matplot(out[, 1], out[, -1], type = "l", col = 1:Nsp, 
                      lty = 1:Nsp, xlab = "time", ylab = "N", lwd = 2, 
                      ylim = range(c(out[, -1], species), na.rm = T))
    graphics::matpoints(time, species, col = 1:Nsp, pch = 1:Nsp)
    graphics::legend("topleft", colnames(species), col = 1:Nsp, 
                     lwd = 2, lty = 1:Nsp, pch = 1:Nsp, bty = "n")
  }
  if (doopt) {
    fisher_info <- unname(solve(-optout$hessian))
    optout$par_sd <- sqrt(abs(diag(fisher_info)))
    parm_signs_sp <- c(rep(1, ncol(opt_data) - 1), parm_signs)
    parameter_intervals <- data.frame(lower_sd = numeric(length(parm_signs_sp)), 
                                      mu = numeric(length(parm_signs_sp)), upper_sd = numeric(length(parm_signs_sp)))
    parameter_intervals[parm_signs_sp != 0, ] <- data.frame(lower_sd = exp(optout$par - 
                                                                             optout$par_sd) * parm_signs_sp[parm_signs_sp != 
                                                                                                              0], mu = exp(optout$par) * parm_signs_sp[parm_signs_sp != 
                                                                                                                                                         0], upper_sd = exp(optout$par + optout$par_sd) * 
                                                              parm_signs_sp[parm_signs_sp != 0])
    tmp1 <- parameter_intervals$lower_sd
    tmp2 <- parameter_intervals$upper_sd
    ps <- which(parameter_intervals$lower_sd > parameter_intervals$upper_sd)
    parameter_intervals$lower_sd[ps] <- tmp2[ps]
    parameter_intervals$upper_sd[ps] <- tmp1[ps]
  }
  else {
    parameter_intervals <- matrix(ncol = 1, data = c(initialN, 
                                                     parms))
  }
  row.names(parameter_intervals) <- c(paste(colnames(species), 
                                            "0", sep = ""), paste("r", 1:Nsp, sep = ""), paste("a", 
                                                                                               rep(1:Nsp, each = Nsp), rep(1:Nsp, Nsp), sep = ""))
  return(list(out = out, parameter_intervals = parameter_intervals, 
              optout = optout, rawdata = data.frame(time, species)))
}
