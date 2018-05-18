rm(list=ls())
source("EL/scelcount.R")
source("EL/scelcount2.R")
source("functions.R")

# Most important switch: do you want to estimate coefficients and shares, or coefficients only?
shares <- TRUE

# Set to FALSE if you are planning to run this code on one multi-core machine; TRUE if you are working on a cluster
multimachine <- FALSE

##########################
#  MAIN
##########################

cat(as.character(start.time <- Sys.time() ), "\n")

# Create variables to pass to stratSampleLinearQSEL and call the function
# Unique for all simulations
beta0 <- 1 # Linear model intercept
beta1 <- 1 # Linear model slope
sdXs <- 1 # Regressor standard deviation (for log-normal distribution)
sdUs <- 1 # Error variance
MC <- 1000 # Number of Monte-Carlo simulations
boundary <- 1.4 # Boundary between two strata
Pl <- c(0.9, 0.3) # Probability of being included in the sample, per stratum
params <- c(beta0, beta1, sdXs, sdUs)
powersize <- TRUE # TRUE because empirical size and power need to be computed
wald <- TRUE
verbose <- TRUE
N <- 50 # 50, 150 or 500 was used in the paper
# Naïve bandwidths used in the paper for untransformed X: should be 0.3 for N=50, 0.4 for N=150 and 0.8 for N=500
# band <- if (N<100) 0.3 else if (N<200) 0.4 else 0.8
# Naïve bandwidths used in the paper for log-transformed X: should be 0.5 for N=50, 0.45 for N=150 and 0.35 for N=500
band <- if (N<100) 0.5 else if (N<200) 0.45 else 0.35
heteroskedastic <- TRUE
strat.var <- "Y"
optmethod <- "Nelder-Mead" # Passed to optim()
gridtransform <- "log" # "log" drastically improves the performance of the SEL estimator under log-normally distributed X; otherwise use 'NULL'

filename <- paste0("SELQmu-N", N, strat.var, ifelse(heteroskedastic, "-het", "-hom"), "-MC", MC, "-bw", band, gridtransform, ".RData")
cat("Output will be in ", filename, "\n")

if (shares) {
  SELofseed <- function(sid) stratSampleLinearQSEL(N=N, params=params, boundary=boundary, Pl=Pl, strat.var=strat.var, heteroskedastic=heteroskedastic, band=band, powersize=powersize, wald=wald, verbose=verbose, optmethod=optmethod, gridtransform=gridtransform, seed=sid) 
} else {
  SELofseed <- function(sid) stratSampleLinearSEL(N=N,  params=params, boundary=boundary, Pl=Pl, strat.var=strat.var, heteroskedastic=heteroskedastic, band=band, powersize=powersize, wald=wald, verbose=verbose, optmethod=optmethod, gridtransform=gridtransform, seed=sid)
}

if (multimachine) { # This part applies only for snow::parLapply implementation
  library(Rmpi)
  library(snow)
  cl <- makeMPIcluster(length(readLines(Sys.getenv("OAR_NODE_FILE"))))
  clusterExport(cl, c("generate.data", "SELofseed", "cemplik", "cemplik2",
                      "smoothEmplik", "linearSmoothEmplik", "linearSmoothEmplikTest",
                      "smoothEmplik2", "linearMuSmoothEmplik", "linearMuSmoothEmplikTest", "linearSharesSmoothEmplik", "linearSharesSmoothEmplikTest",
                      "linearMuSmoothEmplikWald", "linearSharesSmoothEmplikWald",
                      "stratSampleLinearSEL",  "stratSampleLinearQSEL",
                      "params", "N", "strat.var", "heteroskedastic", "boundary", "Pl", "band",
                      "powersize", "wald", "verbose", "optmethod", "gridtransform"))
  res <- parLapply(cl, 1:MC, SELofseed)
  stopCluster(cl)
} else { # This part applies only for parallel::mclapply implementation
  library(parallel)
  num.workers <- if (.Platform$OS.type=="windows") 1 else detectCores()
  res <- mclapply(1:MC, SELofseed, mc.cores=num.workers)
}

save(res, file=filename)
cat("It took", round(difftime(Sys.time(), start.time, units="hours"), 2), "hours to do", MC, "replications\n")

if (multimachine) mpi.quit()
  