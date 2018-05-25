rm(list=ls());
source("functions.R")
print(Sys.time())

library(parallel)
library(sandwich)
# library(plyr) # Uncomment only if you want to check for correct proportions in StratSampOLS

# Create variables to pass to StratSampOLS and call the function
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
num.workers <- if (.Platform$OS.type=="windows") 1 else detectCores()
# Varying across simulations
skedasis <- c(TRUE, FALSE) # Estimate homo- and heteroskedastic models
strat.vars <- c("Y", "X") # Estimate models with endo- and exogenous stratification
sample.sizes <- c(50, 150, 500) # Original sample size BEFORE stratification
methods <- c("OLS", "GMM")

# Calculate true shares by approximation
Nbig <- 999999 # With this value, the block below should run in under 1 minute
nrep <- 20
Q.hom <- mclapply(1:nrep, function(sid) {Zhom <- generate.data(N=Nbig, params=params, boundary=boundary, Pl=c(1,1), strat.var="Y", heteroskedastic=FALSE, seed=sid); return(sum(Zhom$Y<boundary)/nrow(Zhom))}, mc.cores = num.workers)
Q.het <- mclapply(1:nrep, function(sid) {Zhet <- generate.data(N=Nbig, params=params, boundary=boundary, Pl=c(1,1), strat.var="Y", heteroskedastic=TRUE,  seed=sid); return(sum(Zhet$Y<boundary)/nrow(Zhet))}, mc.cores = num.workers)
Q.hom <- mean(unlist(Q.hom))
Q.het <- mean(unlist(Q.het))
Q.X <- plnorm(boundary)
print("Estimated theoretical aggregate shares, starting simulations...")

##########################
# MAIN LOOPS
##########################

# These loops will run through all possible model designs w.r.t. skedasticity, stratification variable,
# OLS/GMM method and sample sizes, outputting the statistics into a table

# Uncomment this to try specific designs
# skedasis <- c(TRUE)
# strat.vars <- c("X")
# sample.sizes <- c(50, 150, 500)
# methods <- c("OLS", "GMM")

results <- data.frame(array(NA, dim=c(length(skedasis)*length(strat.vars)*length(sample.sizes)*length(methods), 4+6+3+6)))
names(results) <- c("Heteroskedastic", "Stratification", "N", "Method",
                    "beta0.Bias", "beta0.SE", "beta0.RMSE", "beta1.Bias", "beta1.SE", "beta1.RMSE", "Q.Bias", "Q.SE", "Q.RMSE",
                    "size1.HC0", "size1.HC3", "size2.HC0", "size2.HC3", "power1.HC0", "power1.HC3")
thisrow <- 1

for (heteroskedastic in skedasis) {
  for (strat.var in strat.vars) {
    for (N in sample.sizes) {
      for (method in methods) {
        start.time <- Sys.time() 
        results$Heteroskedastic[thisrow] <- heteroskedastic
        results$Stratification[thisrow] <- strat.var
        results$N[thisrow] <- N
        results$Method[thisrow] <- method
        
        res	<- mclapply(1:MC, function(sid) StratSampOLS(N = N, params = params, boundary = boundary, Pl = Pl, strat.var = strat.var, method = method, heteroskedastic = heteroskedastic, powersize = powersize, seed = sid), mc.cores = num.workers)
        coeff <- data.frame(do.call("rbind", (lapply(res, '[[', "coefficients"))) ) # Extract only the estimated coefficients
        shares <- do.call("c", (lapply(res, '[[', "shares"))) # Extract only estimated shares
        colnames(coeff) <- c("beta0", "beta1")
        
        results[thisrow, c("beta0.Bias", "beta1.Bias")] <- colMeans(sweep(coeff, 2, params[1:2]))
        results[thisrow, c("beta0.SE", "beta1.SE")] <- sapply(coeff, sd)
        results[thisrow, c("beta0.RMSE", "beta1.RMSE")] <- sqrt(colMeans(sweep(coeff, 2, params[1:2])^2))
        
        # Computing share statistics for the GMM case
        if (method=="GMM") {
          Q.true <- if (strat.var=="X") Q.X else (if(heteroskedastic) Q.het else Q.hom)
          results$Q.Bias[thisrow] <- mean(shares) - Q.true
          results$Q.SE[thisrow] <- sd(shares)
          results$Q.RMSE[thisrow] <- sqrt(mean((Q.true - shares)^2))
        }
        
        if (powersize) {
          zstatsHC0 <- data.frame(do.call("rbind", (lapply(res, '[[', "zstatHC0"))) )
          zstatsHC3 <- data.frame(do.call("rbind", (lapply(res, '[[', "zstatHC3"))) )
          Wald <- data.frame(do.call("rbind", (lapply(res, '[[', "Wald"))) )
          colnames(zstatsHC0) <- c("t0", "t1")
          colnames(zstatsHC3) <- c("t0", "t1")
          results[thisrow, c("size1.HC0", "size1.HC3", "size2.HC0", "size2.HC3", "power1.HC0", "power1.HC3")] <- 
          c(sum(abs(zstatsHC0$t1) >= qnorm(0.975))/MC, sum(abs(zstatsHC3$t1) >= qnorm(0.975))/MC,
            sum(Wald$WaldHC0 >= qchisq(0.95, 2))/MC, sum(Wald$WaldHC3 >= qchisq(0.95, 2))/MC,
            sum(abs(zstatsHC0$t0) >= qnorm(0.975))/MC, sum(abs(zstatsHC3$t0) >= qnorm(0.975))/MC)
        }
        
        # Uncomment to obtain density plots
        # den.const <- density(coeff$beta0, bw="ucv")
        # den.X <- density(coeff$beta1, bw="ucv")
        # pdf(paste0("coef-", (if(method=="GMM") "gmm" else "ols"), (if (heteroskedastic) "-het" else "-hom") , "-strat-", strat.var, "-", N, ".pdf"), 7, 5)
        # plot(den.X, lwd=2, xlim=range(den.X$x, den.const$x), ylim=range(0, den.X$y, den.const$y),
        # main=paste0("MC estimates with stratification by ", strat.var, ", N=", N), xlab=paste0("N = ", MC))
        # lines(den.const, lwd=2, lty=2)
        # legend("topright", c("Slope", "Intercept"), lwd=2, lty=c(1, 2))
        # dev.off()
        
        elapsed <- difftime(Sys.time(), start.time, units="secs")
        msg <- paste0((if (heteroskedastic) "Heterosk." else "Homosk."), ", strat. by ", strat.var, ", N=", N, ", ", method, ", MC=", MC, " done in")
        cat(msg, round(elapsed, 1), "seconds\n")
        
        thisrow <- thisrow+1
      }
    }
  }
}

write.csv(results, "results.csv", row.names = FALSE)
