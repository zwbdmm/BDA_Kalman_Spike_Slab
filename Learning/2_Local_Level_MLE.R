##################################################################
# 2_Local_Level_MLE.R
#
# The file 1_Local_Level used fixed values for disturbances,
#  sigsq.ep and sigsq.eta.
#
# Now I will estimate them via MLE
#
#################################################################
#==================================
# Regular MLE
#  Uses the "vanilla" log likelihood
# Not how book estimates the parameters in example (see below)
#==================================

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to make variance matrix omega (2.4),(2.5)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getOmega <- function(sigsq.ep, sigsq.eta, P1, n){
    # Note: When it comes to a symmetric matrix, I don't think
    #       apply methods are that much faster than a for-method
    # 
    # Init OMEGA by filling its diagonal
    omega <- diag(P1 + sigsq.ep + (1:n)*sigsq.eta  ,  n)
    # Fill every element
    for(ii in 1:n){
        jj <- 1
        while(jj < ii){
            omega[ii,jj] <- P1 + (jj-1) * sigsq.eta
            omega[jj,ii] <- omega[ii,jj] # symmetry
            jj <- jj + 1
        }
    }
    return(omega)
}

#~~~~~~~~~~~~~~~~
# LOG LIKELIHOOD
# This function's PROPORTIONAL TO the log likelihood
#  for the local level model
# We want to optim wrt sigsq.ep and sigsq.eta, which I enter as x
#~~~~~~~~~~~~~~~~
localLevelLogLikely <- function(x, Y, a1, P1){
    if(length(x)!=2){
        stop("Incorrect dim of x")
    }
    sigsq.ep <- x[1]
    sigsq.eta <- x[2]
    print(sigsq.eta)
    #
    omega <- getOmega(sigsq.ep = sigsq.ep , sigsq.eta = sigsq.eta, P1 = P1, n = length(Y))
    # pre-compute
    #    (as.numeric gets rid of the attributes)
    logdet <- as.numeric(determinant(omega, logarithm=TRUE)$modulus)
    inv <- chol2inv(chol(omega))
    # formula
    #  > the colsums thing is a quadratic form
    return(
        -0.5 * logdet - 0.5 * colSums( (Y-a1) * (inv %*% (Y-a1)) )
    )
}

# maximization using fnscale = -1
#
# Making P1 very large is like a diffuse prior (lots of uncertainty)
optim(par = c(15000, 1500), fn=localLevelLogLikely, control=list(fnscale = -1), Y=Nile[-1], a1=0, P1=1e7)


# =================================================
# Concentrated diffuse log likelihood
#   (1) Let P -> infinity  [making variance of first component huge...kind of like improper prior]
#   (2) Let sigsq.eta = q * sigsq.epsilon
#
# Condition (1) gives it the name "diffuse"
# Condition (2) gives it the name concentrated
#
# The purpose is that it reduces the size of the parameter space (to optimize over)
#  to one-dimensional, in just q
#
# Pages 36 and 37 of the text
# 
# =================================================

concDiffLogLikely <- function(x, Y){
    q <- x
    n <- length(Y)
    A <- numeric(n + 1)
    F <- numeric(n)
    P <- numeric(n + 1)
    v <- numeric(n)
    P[2] <- q + 1
    A[2] <- Y[1]
    for(tt in 2:n){
        v[tt] <- Y[tt] - A[tt]
        F[tt] <- P[tt] + 1
        A[tt+1] <- A[tt] + (P[tt]/F[tt]) * v[tt]
        P[tt + 1] <- P[tt] *((F[tt]-P[tt])/F[tt]) + q
    }
    sigsq.ep.hat <- (1/(n-1)) * sum( (v[-1]^2) / (F[-1])  )
    print(sigsq.ep.hat)

    return(
       -0.5*(n-1)*log(sigsq.ep.hat)-0.5*sum(log(F[-1]))
    )

    #
    # Compute full likelihood for comparison w/ table 2.1
    #
    #return(
    #   -0.5*n*log(2*pi) - 0.5*(n-1) -0.5*(n-1)*log(sigsq.ep.hat)-0.5*sum(log(F[-1]))
    #)
}

optim(par = 1, method = "Brent", lower = 0, upper = 20, fn=concDiffLogLikely, control=list(fnscale = -1), Y=as.numeric(Nile))
# Gives same siqsq.ep and q as the text.
