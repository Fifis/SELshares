# Contains 'cemplik2', a modification of the original 'cemplik' code by A. Owen

# This modification allows to pass an extra argument to the cemplik function
# Feature needed for the code to be used in a smoothed empirical likelihood (SEL) setting
# when imposing both conditional and unconditional moment restrictions

# Code modified by Antonio Cosma in September 2017

# Modification of scel.R to handle data with lots of ties.
# The data vectors are rows of z and the count variable
# indicates their multiplicity.
#
# Technical debt: 
#  a) this is copy/paste/edit of the underlying algorithm to allow weights.
#  b) there are lots of internal functions to avoid cluttering the namespace.
#     Those might induce a mild cost from reassignment at each call.
#
# The underlying algorithm is the one for
# self-concordant empirical likelihood for a vector mean,
# as described in:
#
# @article{owen:2013,
#  title={Self-concordance for empirical likelihood},
#  author={Owen, A. B.},
#  journal={Canadian Journal of Statistics},
#  volume={41},
#  number={3},
#  pages={387--397},
#  year={2013}
# }

# This file has two callable functions:
#   cemplik    does one EL calculation,
#   ctracelr   calls emplik on a trajectory from mu0 to mu1 in N+1 steps

# cemplik returns:
#  logelr     log empirical likelihood ratio
#  lam        Lagrange multiplier (vector of length d)
#  wts        observation weights/probabilities (vector of length n)
#  converged  TRUE if algorithm converged
#             FALSE usually means that mu is not in the convex hull of the data
#             you then get a very small likelihood (instead of zero)
#  iter       number of iterations taken
#  ndec       Newton decrement (see Boyd & Vandenberghe) 
#  gradnorm   norm of gradient of log empirical likelihood


# There are also two internal functions:
#   mllog     gives a family of alternatives to -log() and derivative thereof
#             in order to attain self-concordance.
#   svdlm     does least squares regression via the SVD

# These lengthen the emplik function, but making them internal keeps them
# out of the users' name spaces

# Art B. Owen, Feb 2014

cemplik2 = function( z,  # matrix with one data vector per row, a column vector is ok when d=1
                   ct,  # vector of counts
                shift,  # computes 1 + shift + z%*%lam instead of 1 + z%*%lam
                   mu,  # hypothesized mean, default (0 ... 0) in R^d
                  lam,  # starting lambda, default (0 ... 0)
                  eps,  # lower cutoff for -log( ), default 1/nrow(z)
                    M,  # upper cutoff for -log( ), default Inf
         thresh=1e-30,  # convergence threshold for log likelihood (default is aggressive)
          itermax=100,  # upper bound on number of Newton steps (seems ample)
        verbose=FALSE)  # controls printed output
{ 
# empirical likelihood test for whether
# mean of (rows of) z is mu

if( min(ct)<0 )
  stop("Negative weights not allowed.")
# NB  negative weights might be useful but they can destroy convexity or even boundedness.
# They also make the Newton step fail to be of least squares type.

wttol = 0.001 
# Warn if any nonzero weights are smaller than wttol
if( any(0<ct & ct<wttol) ){
  warning(paste("Positive counts below",wttol,"have been replaced by zero."))
  ct[ct<wttol]=0
}
nonz = (ct>0) # These are observations with nonzero weight.  Used in Newton step below.

if( sum(ct)<=0 )stop("Total weight must be positive.")
# NB  this still allows some zero counts.

# Internal function mllog, modified -log( ) with derivatives
  
mllog = function( x, eps, M, der=0 ){
# minus log and its first der derivatives, on  eps < x < M
# 4th order Taylor approx to left of eps and right of M
# der = 0 or 1 or 2
# 4th order is lowest that gives self concordance

if( missing(M) )
  M = Inf
if( eps>M )
  stop("Thresholds out of order")

lo = x < eps
hi = x > M
md = (!lo) & (!hi)

# Coefficients for 4th order Taylor approx below eps
coefs      = rep(0,5)
coefs[1]   = -log(eps)
coefs[2:5] = (-eps)^-(1:4)/(1:4)

# Coefficients for 4th order Taylor approx above M
Coefs      = rep(0,5)
Coefs[1]   = -log(M)
Coefs[2:5] = (-M)^-(1:4)/(1:4)

# degree 4 polynomial approx to log
h = function(y,cvals){ # cvals are coefs at eps, Coefs at M
  # sum c[t+1] y^t
  tee = 1:4
  ans = y*0
  ans = ans + cvals[1]
  for( j in tee )
    ans = ans + y^j*cvals[j+1]
  ans
}

# first derivative of h at y, from approx at pt
hp = function(y,pt){
  tee = 0:3
  ans = y*0
  for( j in tee )
    ans = ans + (-y/pt)^j
  ans = ans * (-pt)^-1
  ans
}

# second derivative of h at y, from approx at pt
hpp = function(y,pt){
  tee = 0:2
  ans = y*0
  for( j in tee )
    ans = ans + (j+1) * (-y/pt)^j
  ans = ans *(-pt)^-2 
  ans
}

# function value
f      = x*0
f[lo]  = h( x[lo]-eps, coefs )
f[hi]  = h( x[hi]-M,   Coefs )
f[md]  = -log(x[md])

if( der<1 )return(cbind(f))

# first derivative
fp     = x*0
fp[lo] = hp( x[lo]-eps, eps )
fp[hi] = hp( x[hi]-M, M )
fp[md] = -1/x[md]

if( der<2 )return(cbind(f,fp))

# second derivative
fpp     = x*0
fpp[lo] = hpp( x[lo]-eps, eps )
fpp[hi] = hpp( x[hi]-M, M )
fpp[md] = 1/x[md]^2

return( cbind(f,fp,fpp) )
# End of mllog()
} 

# Internal function to do a linear model via SVD
# Empirical likelihood's Newton steps are of
# least squares type.

svdlm = function(X,y){
# Linear model regression coefficient via SVD

  # Tolerances for generalized inverse via SVD
  RELTOL = 1e-9     
  ABSTOL = 1e-100

  # Get Xplus = generalized inverse of X
  # If svd algorithm failures are encountered
  # it sometimes helps to try svd(t(X)) and
  # translate back. First check to ensure that
  # X does not contain NaN or Inf or -Inf.
  svdX     = svd(X)
  d        = svdX$d
  lo       = d < (RELTOL * max(d) + ABSTOL)
  dinv     = 1/d
  dinv[lo] = 0
  Xplus    = svdX$v %*% diag(dinv,nrow=length(dinv)) %*% t(svdX$u)
             # taking care with diag when dinv is 1x1
             # to avoid getting the identity matrix of
             # size floor(dinv)

  # Return X^+ y
  Xplus %*% matrix(y,ncol=1)
}
# end of svdlm

  
# Backtracking line search parameters [Tweak only with extreme caution.] 
# See Boyd and Vandenberghe, pp464-466.
ALPHA = 0.3  # seems better than 0.01 on some 2d test data (sometimes fewer iters)
BETA  = 0.8
# We need  0 < ALPHA < 1/2   and 0 < BETA < 1

# Backtrack threshold: you can miss by this much.
BACKEPS = 0
# Consider replacing 0 by 1e-10 if backtracking seems to be
# failing due to round off.

if( is.vector(z) )
  z = matrix(z,ncol=1)

n = nrow(z)
d = ncol(z)

if( missing(mu)  )
  mu  = rep(0,d)
z = t( t(z)-mu ) # subtract mu from each z[i,]

if( missing(eps) )eps = 1/n
if( missing(M)   )M   = Inf

#
# Use lam = 0 or initial lam, whichever is best
#

init0 = mllog( rep(1,n) + shift, eps=eps, M=M, der=2 ) # i.e. fn and derivatives at lam = 0
# Weights only appear OUTSIDE of mllog
for( j in 1:ncol(init0))
  init0[,j] = init0[,j]*ct

if( missing(lam) ){
  init = init0
  lam  = rep(0,d)
}else{
  init = mllog( 1+z%*%lam + shift, eps=eps, M=M, der=2 )
  for( j in 1:ncol(init0))
    init[,j] = init[,j]*ct
  if( sum(init0[,1]) < sum(init[,1]) ){
    lam = rep(0,d)
    init = init0
  }
}

# Initial f, g
fold = sum(init[,1])
gold = apply( z * init[,2],2,sum )

converged = FALSE
iter      = 0
oldvals   = init
while( !converged ){
  iter = iter + 1

  # Get Newton Step
  rootllpp = sqrt(oldvals[,3])  # sqrt 2nd deriv of -llog lik
  zt = z
  for( j in 1:d )
    zt[,j] = zt[,j] * rootllpp
  yt   = oldvals[,2] / rootllpp
  step = -svdlm(zt[nonz,],yt[nonz])  #  SVD is more reliable than step = -lm( yt~zt-1 )$coef
                                     #  also we kick out the ones with ct=0 to avoid unnecessary NaNs
  backtrack = FALSE 
  s = 1   # usually called t, but R uses t for transpose
  while( !backtrack ){
    newvals = mllog( 1+z%*%(lam+s*step)+shift,eps=eps,M=M,der=2 )
    for( j in 1:ncol(newvals))
      newvals[,j] = newvals[,j]*ct
    fnew    = sum(newvals[,1])
    targ    = fold + ALPHA * s * sum( gold*step ) + BACKEPS # (BACKEPS for roundoff, should not be needed)
    if(  fnew <= targ ){
      # backtracking has converged
      backtrack = TRUE
      oldvals   = newvals
      fold      = fnew
      gold      = apply( z * oldvals[,2],2,sum )
      # take the step
      lam       = lam + s*step
    }else{
      s = s * BETA
    }
  }

  # Newton decrement and gradient norm
  ndec     = sqrt( sum( (step*gold)^2 ) )
  gradnorm = sqrt( sum(gold^2))

  if(verbose)print(c(fold,gradnorm,ndec,lam))

  converged = ( ndec^2 <= thresh)
  if( iter > itermax )break
}

wts    = (ct/sum(ct))/(1+z%*%lam) # NB numerator was 1/n, now has the weights.
logelr = sum( ct*mllog(1+z%*%lam+shift,eps=eps,M=M,der=0) )

list(logelr=logelr,lam=lam, wts=wts,
   converged=converged,iter=iter,ndec=ndec,gradnorm=gradnorm)
}



