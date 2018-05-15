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
  Xs <- rlnorm(N, sd = sdXs)
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
  # library(plyr)
  # prop.check <- ddply(Zs, .(D), function(df) prop=mean(df$include))
  # print(prop.check)
  
  # Retained sample
  Z <- Zs[include, c(1:3, ncol(Zs))]
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
linearSharesSmoothEmplik <- function(theta, Z, sel.weights, trim = NULL) {
  muopt <- optim(0, linearMuSmoothEmplik, theta=theta, Z=Z, sel.weights=sel.weights, method="Brent", lower=-10, upper=10)$value
  return(muopt)
}

# Function for case 2: with shares
linearSharesSmoothEmplikTest <- function(theta, slope, Z, sel.weights, trim = NULL) {
  muopt <- optim(0, linearMuSmoothEmplikTest, theta=theta, slope=slope, Z=Z, sel.weights=sel.weights, method="Brent", lower=-10, upper=10)$value
  return(muopt)
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
  traceval <- if(verbose) 1 else 0
  
  Z <- generate.data(N=N, params=params, boundary=boundary, Pl=Pl, strat.var=strat.var, heteroskedastic=heteroskedastic, seed=seed)
  
  sel.weights.grid <- Z$X
  if (!is.null(gridtransform)) sel.weights.grid <- gridtransform(sel.weights.grid)
  if (band < 0) {
    bn <- - band * 1.06*sd(sel.weights.grid)*length(sel.weights.grid)^{-1/5} # Bandwidth, Silverman's rule of thumb
  } else bn <- band
  sel.weights <- outer(sel.weights.grid, sel.weights.grid, function(x, y) dnorm((x-y)/bn) )
  sel.weight.den <- rowSums(sel.weights) # sum over 'j' to obtain marginal on 'i'
  sel.weights <- sel.weights / sel.weight.den # w_{ij} ~ p(j|i)
  
  start.values <- lm(Y~X, data=Z, weights=Z$weight)$coeff
  # start.values <- params[1:2] # Uncomment if you want to start the search from true values (infeasible in practice; for testing only)
  
  if (verbose) print("Conducting optimisation of unrestricted SEL, details below")
  optim.unrestricted <- optim(start.values, fn=function(vpar) {linearSmoothEmplik(vpar, Z, sel.weights = sel.weights)}, control=list(fnscale=-1, trace=traceval, REPORT=1), method=optmethod )
  diff1 <- difftime(tic1 <- Sys.time(), tic0, units = "mins")
  if (verbose) print(paste0("Unconstrained optimisation finished in ", round(diff1, 1), " mins"))
  thetahat <- optim.unrestricted$par
  SELur <- optim.unrestricted$value
  
  if (powersize) {
    # Measuring the power and size of the test, i.e. H01: slope = 0 and H02: slope = true value
    start.values0 <- lm(Y~1, data=Z, weights = Z$weight)$coeff # Restricted model with slope = 0
    start.values1 <- lm(I(Y-params[2]*X)~1, data=Z, weights = Z$weight)$coeff # Restricted model with slope = true slope
    if(verbose) print(paste0("Performing optimisation of restricted SEL under false H0: slope=0"))
    optim.restricted0 <- optim(start.values0, fn=function(vpar) {linearSmoothEmplikTest(vpar, slope=0, Z, sel.weights = sel.weights)}, control=list(fnscale=-1, trace=traceval, REPORT=1), method="Brent", lower=-2, upper=15 )
    diff2 <- difftime(tic2 <- Sys.time(), tic1, units = "mins")
    if(verbose) {
      print(paste0("Constrained optimisation 1 finished in ", round(diff2, 1), " mins"))
      print(paste0("Performing optimisation of restricted SEL under true H0: slope=", params[2]))
    }
    optim.restricted1 <- optim(start.values1, fn=function(vpar) {linearSmoothEmplikTest(vpar, slope=params[2], Z, sel.weights = sel.weights)}, control=list(fnscale=-1, trace=traceval, REPORT=1), method="Brent", lower=-1, upper=5 )
    diff3 <- difftime(tic3 <- Sys.time(), tic2, units = "mins")
    if(verbose) print(paste0("Constrained optimisation 2 finished in ", round(diff3, 1), " mins"))
    thetahat0 <- optim.restricted0$par
    SELr0 <- optim.restricted0$value
    thetahat1 <- optim.restricted1$par
    SELr1 <- optim.restricted1$value
  }
  
  results <- list(unrestricted = c(thetahat, SELur=SELur))
  if (powersize) {
    results$h0zero <- c(thetahat0[1], 0, SELr0=SELr0)
    results$h1true <- c(thetahat1[1], 1, SELr1=SELr1)
  }
  if (wald) {
    SELwald <- as.numeric(linearSmoothEmplik(params[1:2], Z, sel.weights))
    results$SELwald <- SELwald
  }
  results$minutes <- if (powersize) as.numeric(c(diff1, diff2, diff3)) else as.numeric(diff1)
  
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
  if (!is.null(gridtransform)) sel.weights.grid <- gridtransform(sel.weights.grid)
  if (band < 0) {
    bn <- - band * 1.06*sd(sel.weights.grid)*length(sel.weights.grid)^{-1/5} # Bandwidth, Silverman's rule of thumb
  } else bn <- band
  sel.weights <- outer(sel.weights.grid, sel.weights.grid, function(x, y) dnorm((x-y)/bn) )
  sel.weight.den <- rowSums(sel.weights) # sum over 'j' to obtain marginal on 'i'
  sel.weights <- sel.weights / sel.weight.den # w_{ij} ~ p(j|i)
  
  Q.mm <- (M/Pl[1]) / (M/Pl[1] + (n-M)/Pl[2])
  start.values <- c(lm(Y~X, data=Z, weights=Z$weight)$coeff, Q=Q.mm)
  # start.values <- c(params[1:2], Q=Q.mm) # Uncomment if you want to start the search from true values (infeasible in practice; for testing only)
  
  if (verbose) print("Conducting optimisation of unrestricted SEL, details below")
  optim.unrestricted <- optim(start.values, fn=function(vpar) {linearSharesSmoothEmplik(vpar, Z, sel.weights = sel.weights)}, control=list(fnscale=-1, trace=traceval, REPORT=1), method=optmethod )
  diff1 <- difftime(tic1 <- Sys.time(), tic0, units = "mins")
  if (verbose) print(paste0("Unconstrained optimisation finished in ", round(diff1, 1), " mins"))
  thetahat <- optim.unrestricted$par
  SELur <- optim.unrestricted$value
  
  if (powersize) {
    # Measuring the power and size of the test, i.e. H01: slope = 0 and H02: slope = true value
    start.values0 <- c(lm(Y~1, data=Z, weights = Z$weight)$coeff, Q=Q.mm) # Restricted model with slope = 0
    start.values1 <- c(lm(I(Y-params[2]*X)~1, data=Z, weights = Z$weight)$coeff, Q=Q.mm) # Restricted model with slope = true slope
    if(verbose) print(paste0("Performing optimisation of restricted SEL under false H0: slope=0"))
    optim.restricted0 <- optim(start.values0, fn=function(vpar) {linearSharesSmoothEmplikTest(vpar, slope=0, Z, sel.weights = sel.weights)}, control=list(fnscale=-1, trace=traceval, REPORT=1), method=optmethod )
    diff2 <- difftime(tic2 <- Sys.time(), tic1, units = "mins")
    if(verbose) {
      print(paste0("Constrained optimisation 1 finished in ", round(diff2, 1), " mins"))
      print(paste0("Performing optimisation of restricted SEL under true H0: slope=", params[2]))
    }
    optim.restricted1 <- optim(start.values1, fn=function(vpar) {linearSharesSmoothEmplikTest(vpar, slope=params[2], Z, sel.weights = sel.weights)}, control=list(fnscale=-1, trace=traceval, REPORT=1), method=optmethod )
    diff3 <- difftime(tic3 <- Sys.time(), tic2, units = "mins")
    if(verbose) print(paste0("Constrained optimisation 2 finished in ", round(diff3, 1), " mins"))
    thetahat0 <- optim.restricted0$par
    SELr0 <- optim.restricted0$value
    thetahat1 <- optim.restricted1$par
    SELr1 <- optim.restricted1$value
  }
  
  results <- list(unrestricted = c(thetahat, SELur=SELur), Q.mm = Q.mm)
  if (powersize) {
    results$h0zero <- c(thetahat0[1], 0, thetahat0[2], SELr0=SELr0)
    results$h1true <- c(thetahat1[1], 1, thetahat1[2], SELr1=SELr1)
  }
  results$minutes <- if (powersize) as.numeric(c(diff1, diff2, diff3)) else as.numeric(diff1)
  
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
