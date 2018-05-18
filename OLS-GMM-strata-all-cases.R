rm(list=ls());
source("functions.R")
sink(file = "runlogs.txt", append = FALSE) # Flush the output file
print(Sys.time())
sink()

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
Nbig <- 1e7 # With this value, the block below should run in under 5 seconds
set.seed(20170808)
Xs <- rlnorm(Nbig, sd = sdXs)
Us <- rnorm(Nbig, sd = sdUs)
Xbeta <- cbind(1, Xs) %*% c(beta0, beta1)
Ys <- Xbeta + Us # Homoskedastic data
Q.hom <- sum(Ys<boundary)/Nbig
Us <- Us*sqrt(.1+.2*Xs+.3*Xs^2)
Ys <- Xbeta + Us # Heteroskedastic data
Q.het <- sum(Ys<boundary)/Nbig
rm(Xs, Xbeta, Us, Ys)
print("Estimated theoretical aggregate shares, starting simulations...")

##########################
# MAIN LOOPS
##########################

# These loops will run through all possible model designs w.r.t. skedasticity, stratification variable,
# OLS/GMM method and sample sizes, outputting the annotated statistics into runlogs.txt in the current working directory

# Uncomment this to try specific designs
# skedasis <- c(TRUE)
# strat.vars <- c("X")
# sample.sizes <- c(50, 150, 500)
# methods <- c("OLS", "GMM")

for (heteroskedastic in skedasis) {
  for (strat.var in strat.vars) {
    for (N in sample.sizes) {
      for (method in methods) {
        
        start.time <- Sys.time() 
        stime	<- mclapply(1:MC, function(sid) StratSampOLS(N = N, params = params, boundary = boundary, Pl = Pl, strat.var = strat.var, method = method, heteroskedastic = heteroskedastic, powersize = powersize, seed = sid), mc.cores = num.workers)
        coeff <- data.frame(do.call("rbind", (lapply(stime, '[[', "coefficients"))) ) # Extract only the estimated coefficients
        shares <- do.call("c", (lapply(stime, '[[', "shares"))) # Extract only estimated shares
        colnames(coeff)[1] <- "Intercept"
        
        # Computing share statistics for the GMM case
        if (method=="GMM") {
          Q.true <- if (strat.var=="X") plnorm(boundary) else (if(heteroskedastic) Q.het else Q.hom)
          share.bias <- mean(shares) - Q.true
          share.msd <- sd(shares)
          share.rmse <- sqrt(mean((Q.true - shares)^2))
          share.table <- c(share.bias, share.msd, share.rmse)
          names(share.table) <- c("Share bias", "Share SE", "Share RMSE")
        }
        
        # Computing coefficient statistics
        bias <- colMeans(sweep(coeff, 2, params[1:2]))
        msd <- sapply(coeff, sd)
        rmse <- sqrt(colMeans(sweep(coeff, 2, params[1:2])^2))
        
        # The calculated RMSE might experience precision issue: see with the fragment below
        # Uncomment these in order to check precision
        # rmse2 = sqrt(bias^2 + msd^2)
        # print((rmse2 - rmse)/rmse2) # Relative difference
        
        err.table <- t(rbind(bias, msd, rmse))
        rownames(err.table) <- c("Intercept", "X")
        colnames(err.table) <- c("Bias", "SE", "RMSE")
        err.table2 <- matrix(c(err.table[1, ], err.table[2, ]), nrow=1)
        colnames(err.table2) <- c("Intercept Bias", "Intercept SE", "Intercept RMSE", "Slope Bias", "Slope SE", "Slope RMSE")
        
        if (powersize) {
          zstatsHC0 <- data.frame(do.call("rbind", (lapply(stime, '[[', "zstatHC0"))) )
          zstatsHC3 <- data.frame(do.call("rbind", (lapply(stime, '[[', "zstatHC3"))) )
          Wald <- data.frame(do.call("rbind", (lapply(stime, '[[', "Wald"))) )
          colnames(zstatsHC0) <- c("t0", "t1")
          colnames(zstatsHC3) <- c("t0", "t1")
          h0HC0 <- sum(abs(zstatsHC0$t0) >= qnorm(0.975))/MC
          h1HC0 <- sum(abs(zstatsHC0$t1) >= qnorm(0.975))/MC
          h0HC3 <- sum(abs(zstatsHC3$t0) >= qnorm(0.975))/MC
          h1HC3 <- sum(abs(zstatsHC3$t1) >= qnorm(0.975))/MC
          WaldHC0 <- sum(Wald$WaldHC0 >= qchisq(0.95, 2))/MC
          WaldHC3 <- sum(Wald$WaldHC3 >= qchisq(0.95, 2))/MC
        }
        
        # Uncomment to obtain density plots
        # den.const <- density(coeff$Intercept, bw="ucv")
        # den.X <- density(coeff$X, bw="ucv")
        # pdf(paste0("coef-", (if(method=="GMM") "gmm" else "ols"), (if (heterosk) "-het" else "-hom") , "-strat-", strat.var, "-", N, ".pdf"), 7, 5)
        # plot(den.X, lwd=2, xlim=range(den.X$x, den.const$x), ylim=range(0, den.X$y, den.const$y),
        # main=paste0("MC estimates with stratification by ", strat.var, ", N=", N), xlab=paste0("N = ", MC))
        # lines(den.const, lwd=2, lty=2)
        # legend("topright", c("Slope", "Intercept"), lwd=2, lty=c(1, 2))
        # dev.off()
        
        msg <- paste0((if (heteroskedastic) "Heterosk." else "Homosk."), (if(method=="GMM") " GMM" else " OLS"), ", strat. by ", strat.var, ", N=", N, ", MC=", MC)
        elapsed <- difftime(Sys.time(), start.time, units="secs")
        sink(file = "runlogs.txt", append = TRUE)
        cat(msg, "\n")
        print(round(elapsed, 2))
        print(err.table2)
        if (method=="GMM") print(share.table)
        if (powersize) {
          cat(paste0("Empirical test size (nominal size 5%): the true H0: beta1=", params[2], " was rejected in ", round(h1HC0*100, 1), "% cases (HC0) and ", round(h1HC3*100, 1), "% cases (HC3)\n"))
          cat(paste0("Empirical test size (nominal size 5%): the true H0: beta1=beta2=true values was rejected in ", round(WaldHC0*100, 1), "% cases (HC0) and ", round(WaldHC3*100, 1), "% cases (HC3)\n"))
          cat(paste0("Empirical test power (nominal size 5%): the false H0: beta1=0 was rejected in ", round(h0HC0*100, 1), "% cases (HC0) and ", round(h0HC3*100, 1), "% cases (HC3)\n"))
        }
        cat("\n_______________________________________\n")
        sink()
        
        print(paste(msg, "in", round(elapsed, 1), "s"))
        
      }
    }
  }
}

# Check runlogs.txt in your current working directory to see the numbers!
