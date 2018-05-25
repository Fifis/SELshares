# Universal function for all estimation methods (OLS/GMM/SEL)
# Data generating process for a linear model
generate.data <- function(N = 150, # Number of observations
                          params = c(1,1,1,1), # Vector: intercept, slope coefficient for linear model, sd(log(X)), sd(U)
                          boundary = 1.4, # Threshold value that separates the strata
                          Pl = c(0.9, 0.3), # Vector: probabilities of retaining a value from the first and second strata respectively
                          strat.var = "X", # "X" or "Y" for exogenous or endogenous sampling respectively
                          heteroskedastic = FALSE, # FALSE for homoskedastic design and TRUE for heteroskedastic error term
                          seed = 1 # Seed for Monte Carlo
) {
  set.seed(seed)
  beta <- params[1:2]
  sdXs <- params[3]
  sdUs <- params[4]
  
  # Initialise the model, variables drawn from the target population
  Xs <- sort(rlnorm(N, sd = sdXs))
  Us <- rnorm(N, sd = sdUs)
  if (heteroskedastic) Us <- Us*sqrt(.1+.2*Xs+.3*Xs^2)
  Ys <- cbind(1, Xs) %*% beta + Us
  
  # Variables on which selection is made
  S <- if (strat.var == "X") Xs else if (strat.var == "Y") Ys else stop("strat.var should be either 'X' or 'Y'")
  D <- ifelse(S < boundary, 0, 1) # Label of the strata, 0 (lower) or 1 (upper) 
  T <- runif(N) 		 # T is the variable to decide whether the observation is retained or not 
  include <- ifelse(D==0, T < Pl[1], T < Pl[2])
  weight.i <- ifelse(D==0, 1/Pl[1], 1/Pl[2])
  Zs <- data.frame(Ys, Xs, D, T, include, weight.i)
  
  # Uncomment below to check for correct proportions
  # aggregate(Zs$include, by=list(Zs$D), FUN=mean)
  
  # Retained sample
  Z <- Zs[include, c(1:3, ncol(Zs))]
  # Since the distribution is heavy-tailed log-normal, we trim off a couple of points to get rid of outliers and improve all estimators
  # Z is sorted on X values, so the worst values are at the end of the dataframe; floor(log10(N)) seems a very sparing value
  # For N<100, it will trim 1 point, for 100<N<1000, 2 points etc.
  Z <- Z[1:(nrow(Z)-floor(log10(nrow(Z)))),]
  colnames(Z) <- c("Y", "X", "D", "weight")
  
  return(Z)
}

# Function for case 1: no shares
# Wrapper for Owen's cemplik for estimation of coefficients without shares
smoothEmplik <- function(z, sel.weights){
  sel <- apply(sel.weights, MARGIN = 1, function(w) suppressWarnings(cemplik(z, ct=w)))
  return(sel)
}

# Function for case 2: with shares
# Wrapper for Owen's cemplik2 for estimation of coefficients and shares
smoothEmplik2 <- function(z, sel.weights, shift) {
  sel <- apply(sel.weights, MARGIN = 1, function(w) suppressWarnings(cemplik2(z, ct=w, shift=shift)))
  return(sel)
}

# Function for case 1: no shares
# A function that takes parameters and returns the likelihood
linearSmoothEmplik <- function(beta, Z, sel.weights, trim = NULL) {
  if (is.null(trim)) trim <- rep(1, nrow(Z))
  rho <- (Z$Y - cbind(1, Z$X) %*% beta)*Z$weight # = 0; conditional moment restriction
  empliklist <- smoothEmplik(rho, sel.weights)  # returns a list; One for each conditioning vector point
  logsemplik <- trim %*% unlist(lapply(empliklist, '[[', 1))
  return(logsemplik)
}

# Function for case 1: no shares
# Modified function in the inner loop: takes a parameter (intercept) and a hypothesised slope and returns the likelihood
# Here, beta = beta0 and slope = beta1
linearSmoothEmplikTest <- function(beta, slope, Z, sel.weights, trim = NULL) {
  theta <- c(beta, slope)
  if (is.null(trim)) trim <- rep(1, nrow(Z))
  rho <- (Z$Y - cbind(1, Z$X) %*% theta)*Z$weight # = 0; conditional moment restriction
  empliklist <- smoothEmplik(rho, sel.weights)  # returns a list; One for each conditioning vector point
  logsemplik <- trim %*% unlist(lapply(empliklist, '[[', 1))
  return(logsemplik)
}

# Function for case 2: with shares
# Original function in the inner loop: takes parameters, mu and returns the likelihood
# Here, theta = c(beta0, beta1, Q)
linearMuSmoothEmplik <- function(mu, theta, Z, sel.weights, trim = NULL) {
  beta <- theta[1:2]
  Q <- theta[3]
  if (is.null(trim)) trim <- rep(1, nrow(Z))
  rho1 <- (Z$Y - cbind(1, Z$X) %*% beta) * Z$weight	# = 0; conditional moment restriction
  rho2 <- ((1 - Z$D) - Q) * Z$weight
  empliklist <- smoothEmplik2(rho1, sel.weights, shift = mu*rho2) # returns a list; one for each conditioning element
  logsemplik <- trim %*% unlist(lapply(empliklist, '[[', 1))
  return(logsemplik)
}

# Function for case 2: with shares
# Modified function in the inner loop: takes parameters, mu and the hypothesised slope and returns the likelihood
# Here, theta = c(beta0, Q) and slope = beta1
linearMuSmoothEmplikTest <- function(mu, theta, slope, Z, sel.weights, trim = NULL) {
  beta <- c(theta[1], slope)
  Q <- theta[2]
  if(is.null(trim)) trim <- rep(1, nrow(Z))
  rho1 <- (Z$Y - cbind(1, Z$X) %*% beta) * Z$weight
  rho2 <- ((1 - Z$D) - Q) * Z$weight
  empliklist <- smoothEmplik2(rho1, sel.weights, shift = mu*rho2)
  logsemplik <- trim %*% unlist(lapply(empliklist, '[[', 1))
  return(logsemplik)
}

# Function for case 2: with shares
# Modified function in the inner loop: takes mu and the hypothesised coefficients and returns the likelihood
# It is called 'Wald' for comparability with OLS/GMM Wald test results for two coefficients
# Here, theta = Q and coefs = c(beta0, beta1)
linearMuSmoothEmplikWald <- function(mu, theta, coefs, Z, sel.weights, trim = NULL) {
  beta <- c(coefs[1], coefs[2])
  Q <- theta
  if(is.null(trim)) trim <- rep(1, nrow(Z))
  rho1 <- (Z$Y - cbind(1, Z$X) %*% beta) * Z$weight
  rho2 <- ((1 - Z$D) - Q) * Z$weight
  empliklist <- smoothEmplik2(rho1, sel.weights, shift = mu*rho2)
  logsemplik <- trim %*% unlist(lapply(empliklist, '[[', 1))
  return(logsemplik)
}

# Function for case 2: with shares
linearSharesSmoothEmplik <- function(theta, Z, sel.weights, trim = NULL) {
  muopt <- optim(0, linearMuSmoothEmplik, theta=theta, Z=Z, sel.weights=sel.weights, method="Brent", lower=-10, upper=10)$value
  return(muopt)
}

# Function for case 2: with shares
linearSharesSmoothEmplikTest <- function(theta, slope, Z, sel.weights, trim = NULL) {
  muopt <- optim(0, linearMuSmoothEmplikTest, theta=theta, slope=slope, Z=Z, sel.weights=sel.weights, method="Brent", lower=-10, upper=10)$value
  return(muopt)
}

# Function for case 2: with shares
linearSharesSmoothEmplikWald <- function(theta, coefs, Z, sel.weights, trim = NULL) {
  muopt <- optim(0, linearMuSmoothEmplikWald, theta=theta, coefs=coefs, Z=Z, sel.weights=sel.weights, method="Brent", lower=-10, upper=10)$value
  return(muopt)
}

# Function for both cases: find good initial values for SEL using GMM as a reference for order of magnitude
# We do not trust GMM because it is (1) more sensitive to outliers and (2) less efficient
# In order to speed up convergence and reduce maximisation time, we check a grid of points not far away from GMM initval and pick the best one
# This function is completely agnostic about the true values for unrestricted estimation, covers a wide grid of possible points
# Receives the dataset from the external function and does the computaion
maximise.SEL.noshares <- function(Z,
                                  params, # Used only if size and power are being tested
                                  sel.weights,
                                  powersize = TRUE,
                                  wald = TRUE,
                                  verbose = FALSE,
                                  optmethod = "Nelder-Mead"
) {
  traceval <- if (verbose) 1 else 0
  tic0 <- Sys.time()
  if (verbose) print("Searching for a good initial value...")
  start.values <- lm(Y~X, data=Z, weights=Z$weight)$coeff
  start.value.multipliers <- c(-1, -0.5, 0, 0.5, 0.75, 1, 1.25, 1.5, 2)
  start.values.grid.b0 <- sort(start.value.multipliers*start.values[1]) # A reasonable range for initial guesses
  start.values.grid.b1 <- sort(start.value.multipliers*start.values[2])
  start.values.grid.b0 <- c(min(start.values.grid.b0)-1, start.values.grid.b0, max(start.values.grid.b0)+1) # If any of the ends is 0, multiplication will not do anything; addition needed
  start.values.grid.b1 <- c(min(start.values.grid.b1)-1, start.values.grid.b1, max(start.values.grid.b1)+1)
  # The user must know how many initial values have been checked!
  if (verbose) pb <- txtProgressBar(style=3, width=50)
  test.start.grid <- matrix(NA, nrow=length(start.values.grid.b0), ncol=length(start.values.grid.b1))
  progress <- 0 # We are not using Vectorize and outer for the sake of keeping track of progress
  for (i in 1:length(start.values.grid.b0)) {
    for (j in 1:length(start.values.grid.b1)) {
      test.start.grid[i,j] <- linearSmoothEmplik(c(start.values.grid.b0[i], start.values.grid.b1[j]), Z = Z, sel.weights = sel.weights)
      progress <- progress+1
      if (verbose) setTxtProgressBar(pb, progress/prod(dim(test.start.grid)))
    }
  }
  if (verbose) close(pb)
  #
  best.index <- which(test.start.grid==max(test.start.grid), arr.ind = TRUE)
  start.values <- c(beta0=start.values.grid.b0[best.index[1]], beta1=start.values.grid.b1[best.index[2]])
  tic.initval.ur <- Sys.time()
  diff.initval <- as.numeric(difftime(tic.initval.ur, tic0, units = "mins"))
  if (verbose) print(paste0("Searched [", paste(round(range(start.values.grid.b0), 2), collapse=","), "]x[", paste(round(range(start.values.grid.b1), 2), collapse=","),
                            "] in ", round(diff.initval, 1), " mins; SEL(", paste(round(start.values[1:2], 2), collapse=","), ")=", round(max(test.start.grid), 2), "."))
  
  if (verbose) print("Conducting optimisation of unrestricted SEL without shares...")
  optim.unrestricted <- optim(start.values, fn=function(vpar) {linearSmoothEmplik(vpar, Z, sel.weights = sel.weights)}, control=list(fnscale=-1, trace=traceval, REPORT=1), method=optmethod)
  tic.ur <- Sys.time()
  diff.ur <- as.numeric(difftime(tic.ur, tic.initval.ur, units = "mins"))
  thetahat <- optim.unrestricted$par
  SELur <- optim.unrestricted$value
  if (verbose) print(paste0("Unconstrained optimisation finished in ", round(diff.ur, 1), " mins; SEL(", paste(round(thetahat, 2), collapse=","), ")=", round(SELur, 2), "."))
  
  if (powersize) {
    tic.r0 <- Sys.time()
    # Measuring the power and size of the test, i.e. H01: slope = 0 and H02: slope = true value
    start.values0 <- lm(Y~1, data=Z, weights = Z$weight)$coeff # Restricted model with slope = 0
    start.values1 <- lm(I(Y-params[2]*X)~1, data=Z, weights = Z$weight)$coeff # Restricted model with slope = true slope
    range0 <- if (start.values0<0) c(min(2*start.values0, start.values0-3), max(0.5*start.values0, start.values0+3)) else c(min(0.5*start.values0, start.values0-3), max(2*start.values0, start.values0+3))
    range1 <- if (start.values1<0) c(min(2*start.values1, start.values1-2), max(0.5*start.values1, start.values1+2)) else c(min(0.5*start.values1, start.values1-2), max(2*start.values1, start.values1+2))
    
    if(verbose) print(paste0("Performing optimisation of restricted SEL without shares under false H0: slope=0..."))
    optim.restricted0 <- optim(start.values0, fn=function(vpar) {linearSmoothEmplikTest(vpar, slope=0, Z, sel.weights = sel.weights)}, control=list(fnscale=-1, trace=traceval, REPORT=1), method="Brent", lower=range0[1], upper=range0[2])
    tic.r0 <- Sys.time()
    diff.r0 <- as.numeric(difftime(tic.r0, tic.ur, units = "mins"))
    thetahat0 <- optim.restricted0$par
    SELr0 <- optim.restricted0$value
    print(paste0("Constrained optimisation 1 finished in ", round(diff.r0, 1), " mins; SEL(", round(thetahat0, 2), ",0)=", round(SELr0, 2), "."))
    
    if (verbose) print(paste0("Performing optimisation of restricted SEL without shares under true H0: slope=", params[2], "..."))
    optim.restricted1 <- optim(start.values1, fn=function(vpar) {linearSmoothEmplikTest(vpar, slope=params[2], Z, sel.weights = sel.weights)}, control=list(fnscale=-1, trace=traceval, REPORT=1), method="Brent", lower=range1[1], upper=range1[2])
    tic.r1 <- Sys.time()
    diff.r1 <- as.numeric(difftime(tic.r1, tic.r0, units = "mins"))
    thetahat1 <- optim.restricted1$par
    SELr1 <- optim.restricted1$value
    print(paste0("Constrained optimisation 2 finished in ", round(diff.r1, 1), " mins; SEL(", round(thetahat1, 2), ",", params[2],")=", round(SELr1, 2), "."))
  }
  
  results <- list(noshares.ur = c(thetahat, SELur=SELur))
  if (powersize) {
    results$noshares.h0zero <- c(beta0=thetahat0[1], beta1=0, SELr0=SELr0)
    results$noshares.h1true <- c(beta0=thetahat1[1], beta1=params[2], SELr1=SELr1)
  }
  if (wald) {
    SELwald <- as.numeric(linearSmoothEmplik(params[1:2], Z, sel.weights))
    results$noshares.h1wald <- c(beta0=params[1], beta1=params[2], SELwald=SELwald)
  }
  results$minutes <- c(noshares.ur.initval=diff.initval, noshares.ur.optim=diff.ur)
  if (powersize) results$minutes <- c(results$minutes, noshares.h0=diff.r0, noshares.h1=diff.r1)
  return(results)
}

# Function for case 1: no shares
stratSampleLinearSEL <- function(N, # See generate.data for arguments of the same name
                                  params, 
                                  boundary, 
                                  Pl, 
                                  strat.var, 
                                  heteroskedastic,
                                  seed = 1,
                                  trim = NULL, # Trimming function for trimmed SEL
                                  band = -1, # Multiple of Silverman's rot
                                  powersize = FALSE, # Whether SEL under the true and false H0 should be returned;
                                  wald = FALSE, # Estimate the SEL at the true parameter values for comparison with OLS and GMM Wald test
                                  # Leave FALSE to speed up calculations, just for estimation purposes
                                  verbose = FALSE, # Benchmarks time, appends the timing to the output list, and passes trace=1 to optim
                                  optmethod = "Nelder-Mead", # Method that is passed to optim()
                                  gridtransform = NULL # Function that is applied to sel.weights.grid
) { 
  tic0 <- Sys.time()
  if(seed%%10 == 0) cat(seed, " ")
  if(seed%%200 == 0) cat(seed, "\n")
  traceval <- if (verbose) 1 else 0
  
  Z <- generate.data(N=N, params=params, boundary=boundary, Pl=Pl, strat.var=strat.var, heteroskedastic=heteroskedastic, seed=seed)
  
  sel.weights.grid <- Z$X
  if (!is.null(gridtransform)) {
    ftransform <- get(gridtransform)
    sel.weights.grid <- ftransform(sel.weights.grid)
  }
  if (band < 0) {
    bn <- - band * 1.06*sd(sel.weights.grid)*length(sel.weights.grid)^{-1/5} # Bandwidth, Silverman's rule of thumb
  } else bn <- band
  sel.weights <- outer(sel.weights.grid, sel.weights.grid, function(x, y) dnorm((x-y)/bn) )
  sel.weight.den <- rowSums(sel.weights) # sum over 'j' to obtain marginal on 'i'
  sel.weights <- sel.weights / sel.weight.den # w_{ij} ~ p(j|i)
  
  results <- maximise.SEL.noshares(Z=Z, params=params, sel.weights=sel.weights, powersize=powersize, wald=wald, verbose=verbose, optmethod=optmethod)
  
  results$call <- c(N, params, boundary, Pl, strat.var, heteroskedastic, seed, band, bn, powersize, wald, verbose, optmethod)
  names(results$call) <- c("N", "beta0", "beta1", "sdX", "sdU", "boundary", "P1", "P2", "strat.var", "heteroskedastic", "seed", "band", "bn", "powersize", "wald", "verbose", "optmethod")
  
  return(results)
}

# Function for case 2: with shares
stratSampleLinearQSEL <- function(N, # See generate.data for arguments of the same name
                                  params, 
                                  boundary, 
                                  Pl, 
                                  strat.var, 
                                  heteroskedastic,
                                  seed = 1,
                                  trim = NULL, # Trimming function for trimmed SEL
                                  band = -1, # Multiple of Silverman's rot
                                  powersize = FALSE, # Whether SEL under the true and false H0 should be returned;
                                  # Leave FALSE to speed up calculations, just for estimation purposes
                                  wald = FALSE, # Estimate the SEL at the true parameter values for comparison with OLS and GMM Wald test
                                  verbose = FALSE, # Benchmarks time, appends the timing to the output list, and passes trace=1 to optim
                                  optmethod = "Nelder-Mead", # Method that is passed to optim()
                                  gridtransform = NULL # Function that is applied to sel.weights.grid
) { 
  tic0 <- Sys.time()
  if(seed%%10 == 0) cat(seed, " ")
  if(seed%%200 == 0) cat(seed, "\n")
  traceval <- if(verbose) 1 else 0
  
  Z <- generate.data(N=N, params=params, boundary=boundary, Pl=Pl, strat.var=strat.var, heteroskedastic=heteroskedastic, seed=seed)
  n <- nrow(Z) # observations in the sample
  M <- n - sum(Z$D) # observations on the sample from 1st stratum
  
  sel.weights.grid <- Z$X
  if (!is.null(gridtransform)) {
    ftransform <- get(gridtransform)
    sel.weights.grid <- ftransform(sel.weights.grid)
  }
  if (band < 0) {
    bn <- - band * 1.06*sd(sel.weights.grid)*length(sel.weights.grid)^{-1/5} # Bandwidth, Silverman's rule of thumb
  } else bn <- band
  sel.weights <- outer(sel.weights.grid, sel.weights.grid, function(x, y) dnorm((x-y)/bn) )
  sel.weight.den <- rowSums(sel.weights) # sum over 'j' to obtain marginal on 'i'
  sel.weights <- sel.weights / sel.weight.den # w_{ij} ~ p(j|i)
  
  Q.mm <- (M/Pl[1]) / (M/Pl[1] + (n-M)/Pl[2])
  
  initval <- maximise.SEL.noshares(Z=Z, params=params, sel.weights=sel.weights, powersize=powersize, wald=wald, verbose=verbose, optmethod=optmethod)
  
  start.values <- c(initval$noshares.ur[1:2], Q=Q.mm)
  # start.values <- c(params[1:2], Q=if(heteroskedastic)Q.het else Q.hom) # Uncomment if you want to start the search from true values (infeasible in practice; for testing only; Q.het and Q.hom must be in memory!)
  if (verbose) print("Conducting optimisation of unrestricted SEL with shares...")
  tic.initval.ur <- Sys.time()
  optim.unrestricted <- optim(start.values, fn=function(vpar) {linearSharesSmoothEmplik(vpar, Z, sel.weights = sel.weights)}, control=list(fnscale=-1, trace=traceval, REPORT=1), method=optmethod)
  tic.ur <- Sys.time()
  diff.ur <- as.numeric(difftime(tic.ur, tic.initval.ur, units = "mins"))
  thetahat <- optim.unrestricted$par
  SELur <- optim.unrestricted$value
  if (verbose) print(paste0("Unconstrained optimisation finished in ", round(diff.ur, 1), " mins; SEL(", paste(round(thetahat, 2), collapse=","), ")=", round(SELur, 2), "."))
  
  if (powersize) {
    # Measuring the power and size of the test, i.e. H01: slope = 0 and H02: slope = true value
    # As usually, first, we use the maxima of SEL without shares and then add shares and maximise the full function
    start.values0 <- c(initval$noshares.h0zero[1], Q=Q.mm)
    start.values1 <- c(initval$noshares.h1true[1], Q=Q.mm)
    # start.values0 <- c(params[1]+exp(params[3]^2/2), Q=if(heteroskedastic) Q.het else Q.hom) # Uncomment if you want to start the search from true values (infeasible in practice; for testing only; Q.het and Q.hom must be in memory!)
    # start.values1 <- c(params[1], Q=if(heteroskedastic) Q.het else Q.hom)
    tic.initval.r <- Sys.time()
    if (verbose) print(paste0("Performing optimisation of restricted SEL with shares under false H0: slope=0..."))
    optim.restricted0 <- optim(start.values0, fn=function(vpar) {linearSharesSmoothEmplikTest(vpar, slope=0, Z, sel.weights = sel.weights)}, control=list(fnscale=-1, trace=traceval, REPORT=1), method=optmethod)
    tic.r0 <- Sys.time()
    diff.r0 <- as.numeric(difftime(tic.r0, tic.initval.r, units = "mins"))
    thetahat0 <- optim.restricted0$par
    SELr0 <- optim.restricted0$value
    if (verbose) {
      print(paste0("Constrained optimisation 1 finished in ", round(diff.r0, 1), " mins; SEL(", round(thetahat0[1], 2), ",0,", round(thetahat0[2], 2), ")=", round(SELr0, 2), "."))
      print(paste0("Performing optimisation of restricted SEL with shares under true H0: slope=", params[2], "..."))
    }
    optim.restricted1 <- optim(start.values1, fn=function(vpar) {linearSharesSmoothEmplikTest(vpar, slope=params[2], Z, sel.weights = sel.weights)}, control=list(fnscale=-1, trace=traceval, REPORT=1), method=optmethod)
    tic.r1 <- Sys.time()
    diff.r1 <- as.numeric(difftime(tic.r1, tic.r0, units = "mins"))
    thetahat1 <- optim.restricted1$par
    SELr1 <- optim.restricted1$value
    if(verbose) print(paste0("Constrained optimisation 2 finished in ", round(diff.r1, 1), " mins; SEL(", round(thetahat1[1], 2), ",", params[2], ",", round(thetahat1[2], 2), ")=", round(SELr1, 2), "."))
  }
  
  if (wald) {
    start.valuesWald <- Q.mm # Does not really matter for Brent optimiser
    # start.valuesWald <- Q.het
    if (verbose) print(paste0("Performing optimisation of restricted SEL under true H0: coefs=(", paste(params[1:2], collapse = ","), ")."))
    ticWald0 <- Sys.time()
    optim.restrictedWald <- optim(start.valuesWald, fn=function(vpar) {linearSharesSmoothEmplikWald(vpar, coefs=params[1:2], Z, sel.weights = sel.weights)},
                                  control=list(fnscale=-1, trace=traceval, REPORT=1), method="Brent", lower=max(Q.mm/3, 0.02), upper=min(Q.mm*4, 0.98))
    ticWald1 <- Sys.time()
    diff.Wald <- as.numeric(difftime(ticWald1, ticWald0, units = "mins"))
    thetahatWald <- optim.restrictedWald$par
    SELWald <- optim.restrictedWald$value
    if(verbose) print(paste0("Constrained optimisation (Wald) finished in ", round(diff.Wald, 1), " mins; SEL(", params[1], ",", params[2], ",", round(thetahatWald[1], 2), ")=", round(SELWald, 2), "."))
  }
  
  results <- list(unrestricted = c(thetahat, SELur=SELur), Q.mm = Q.mm)
  if (powersize) {
    results$h0zero <- c(thetahat0[1], beta1=0, thetahat0[2], SELr0=SELr0)
    results$h1true <- c(thetahat1[1], beta1=params[2], thetahat1[2], SELr1=SELr1)
  }
  if (wald) results$h1Wald <- c(beta0=params[1], beta1=params[2], Q=thetahatWald, SELWald=SELWald)
  
  # Processing timestamps
  results$minutes <- c(shares.ur.initval=as.numeric(initval$minutes["noshares.ur.initval"]+initval$minutes["noshares.ur.optim"]), shares.ur.optim=diff.ur)
  if (powersize) results$minutes <- c(results$minutes, shares.r.initval=as.numeric(initval$minutes["noshares.h0"]+initval$minutes["noshares.h1"]), shares.h0=diff.r0, shares.h1=diff.r1)
  if (wald) results$minutes <- c(results$minutes, shares.hWald=diff.Wald)
  
  results$call <- c(N, params, boundary, Pl, strat.var, heteroskedastic, seed, band, bn, powersize, wald, verbose, optmethod)
  names(results$call) <- c("N", "beta0", "beta1", "sdX", "sdU", "boundary", "P1", "P2", "strat.var", "heteroskedastic",
                           "seed", "band", "bn", "powersize", "wald", "verbose", "optmethod")
  
  return(results)
}

StratSampOLS <- function(N, # See generate.data for arguments of the same name
                         params,
                         boundary, 
                         Pl,
                         strat.var = "X", 
                         heteroskedastic = FALSE,
                         seed=1,
                         method = "OLS", # "OLS" or "GMM"
                         powersize = TRUE # Whether t-stats under the true and false H0 should be returned;
                         # Set to FALSE to speed up calculations
                         # Also does Wald test for beta0=true value, beta1=true value
) {
  Z <- generate.data(N=N, params=params, boundary=boundary, Pl=Pl, strat.var=strat.var, heteroskedastic=heteroskedastic, seed=seed)
  the.weights <- if (method=="GMM") Z$weight else if (method=="OLS") rep(1, nrow(Z)) else stop("method must be either 'GMM' or 'OLS'")
  stra.OLS <- lm(Y ~ X, data = Z, weights = the.weights)
  
  # White's HC0 heteroskedasticity-robust VCOV computed by hand (6x times faster than sandwich::vcovHC)
  # Xtilde <- model.matrix(stra.OLS)
  # Xbread <- Xtilde * sqrt(stra.OLS$weights)
  # Xmeat <- Xtilde * (stra.OLS$weights * stra.OLS$residuals)
  # bread <- solve(crossprod(Xbread)) # E(XX'/b)
  # meat <- crossprod(Xmeat) # E(XX'U^2/b^2)
  # thetavar2 <- bread %*% meat %*% bread
  
  output <- list() # We need not carry over those huge lm-class objects; we only need coefficients, shares and, if hypotheses are tested, t-stats
  output$coefficients <- stra.OLS$coefficients
  
  if (powersize) {
    if (!("sandwich" %in% (.packages()))) stop("Please load library(sandwich) before invoking StratSampOLS with powersize=TRUE") 
    thetavarHC0 <- vcovHC(stra.OLS, type="HC0")
    thetavarHC3 <- vcovHC(stra.OLS, type="HC3")
    
    # Test if the hand-computed result is the same as the one provided by the package
    # all.equal(thetavarHC0, thetavar2)
    
    t0HC0 <- (stra.OLS$coefficients["X"] - 0)/sqrt(thetavarHC0[2,2]) # H0: theta[2] = 0
    t1HC0 <- (stra.OLS$coefficients["X"] - params[2])/sqrt(thetavarHC0[2,2]) # H0: theta[2] = true value
    t0HC3 <- (stra.OLS$coefficients["X"] - 0)/sqrt(thetavarHC3[2,2])
    t1HC3 <- (stra.OLS$coefficients["X"] - params[2])/sqrt(thetavarHC3[2,2])
    output$zstatHC0 <- as.numeric(c(t0HC0, t1HC0))
    output$zstatHC3 <- as.numeric(c(t0HC3, t1HC3))
    
    # Extra regressions incorporating the information from the hypotheses
    stra.OLS0 <- lm(Y~1, data=Z, weights = the.weights)
    stra.OLS1 <- lm(I(Y-params[2]*X)~1, data=Z, weights = the.weights)
    output$h0zero <- c(stra.OLS0$coefficients, 0)
    output$h1true <- c(stra.OLS1$coefficients, 1)
    
    # Wald test
    WaldHC0 <- as.numeric(t(stra.OLS$coefficients - params[1:2]) %*% solve(thetavarHC0) %*% (stra.OLS$coefficients - params[1:2]))
    WaldHC3 <- as.numeric(t(stra.OLS$coefficients - params[1:2]) %*% solve(thetavarHC3) %*% (stra.OLS$coefficients - params[1:2]))
    output$Wald <- c(WaldHC0 = WaldHC0, WaldHC3 = WaldHC3)
  }
  
  # Calculating the empirical shares
  if (strat.var == "X") {
    sum.1p1 <- sum(Z$X < boundary)/Pl[1]
    sum.1p2 <- sum(Z$X >= boundary)/Pl[2]
  } else if (strat.var == "Y") {
    sum.1p1 <- sum(Z$Y < boundary)/Pl[1]
    sum.1p2 <- sum(Z$Y >= boundary)/Pl[2]
  }
  output$shares <- (sum.1p1)/(sum.1p1 + sum.1p2)
  
  return(output)
}
