######################################################################################
# 3_Smoothing.R
#
#  Replication of Illustration 2.4.3 on pages 22 and 23
# 
# > Whereas filtering estimates a = E(alpha | Y_{t-1}),  
#  smoothing estimates         alpha_hat = E(alpha | Y_n)
# 
# > That is, we go back and refine our estimates of alpha given all the data.
# > It uses a weighted average (r) of all the future estimation errors (v),
#    and averages typically smooth out estimates, right?
# > In any case, a plot of alpha_hat will be smoother than a plot of a
######################################################################################

# =========================================================================
# Filter Stage (This is the function from 1_Local_Level...)
# =========================================================================
filter <- function(Y, sigsq.ep, sigsq.eta, a1, P1){
    n <- length(Y)
    # Init lists of the params we want
    A.filter <- numeric(n)
    P.filter <- numeric(n)
    A.t <- numeric(n + 1)
    P.t <- numeric(n + 1)
    #
    A.t[1] <- a1
    P.t[1] <- P1
    # 
    # Also save the prediction error, Vt , and the variance
    vt <- numeric(n)
    Ft <- numeric(n)

    for(tt in 1:n){
        # Prediction error, Pred variance, Kalman gain
        vt[tt] <- Y[tt] - A.t[tt]
        Ft[tt] <- P.t[tt] + sigsq.ep
        Kt <- P.t[tt]/Ft[tt]
        # Update variance matrices
        P.filter[tt] <- P.t[tt] * (1 - Kt)
        P.t[tt + 1] <- P.filter[tt] + sigsq.eta
        # Update mean matrices
        A.t[tt + 1] <- A.t[tt] + Kt*vt[tt]
        A.filter[tt] <-  A.t[tt + 1] ## These are the same for the local level model
    }

    return(list(
        A.t = A.t[1:n], A.filter = A.filter, P.t = P.t[1:n], P.filter = P.filter, vt = vt, Ft = Ft
    ))
}
# =========================================================================
# Filter Stage (This is the function from 1_Local_Level...)
# =========================================================================
smooth <- function(vt, Ft, at, Pt){
    n <- length(vt)
    # init targets
    alpha.hat <- numeric(n)
    V <- numeric(n)
    # init needed vars
    r <- numeric(n) # weighted sum of disturbances v
    N <- numeric(n) # weighted sum of inverse variances
    r[n] <- 0
    N[n] <- 0
    for(tt in n:2){ # need to handle tt = 1 separately, because r[0] isn't elegantly handled in R
        # Update smoothed mean
        Lt <- 1 - Pt[tt]/Ft[tt]
        r[tt - 1] <- vt[tt] / Ft[tt] + Lt * r[tt]
        alpha.hat[tt] <- at[tt] + Pt[tt]*r[tt - 1] 
        # Update smoothed variance
        N[tt - 1] <- 1/Ft[tt] + (Lt^2)*N[tt]
        V[tt] <- Pt[tt] * (1 - Pt[tt]*N[tt - 1])
    }
    # tt = 1
    # > alpha hat
    L1 <- 1 - Pt[1]/Ft[1]
    r0 <- vt[1] / Ft[1] + L1 * r[1]
    alpha.hat[1] <- at[1] + Pt[1]*r0 
    # > variance
    N0 <- 1/Ft[1] + (Lt^2)*N[1]
    V[1] <- Pt[1] * (1 - Pt[1]*N0)
    

    return(list(alpha.hat=alpha.hat, V = V))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Run It
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res.filt <- filter(Y = as.numeric(Nile), sigsq.ep = 15099, sigsq.eta = 1469.1, a1 = 0, P1 = 1e7)
res.smooth <- smooth(vt=res.filt$vt, Ft = res.filt$Ft, at = res.filt$A.t, Pt = res.filt$P.t)


plot(Nile, type = "p", pch = 20, main = "Nile, data and smoothed model", ylim = c(450, 1300))
lines(start(Nile)[1]:end(Nile)[1], res.smooth$alpha.hat)
lines(start(Nile)[1]:end(Nile)[1], res.smooth$alpha.hat + qnorm(0.05)*sqrt(res.smooth$V), col = "red")
lines(start(Nile)[1]:end(Nile)[1], res.smooth$alpha.hat + qnorm(0.95)*sqrt(res.smooth$V), col = "red")

# compare to filter only
dev.new() # linux version of "windows()" or "quartz()"
plot(Nile, type = "p", pch = 20, main = "Nile, data and filter")
lines(start(Nile)[1], end(Nile)[1], res.filt$A.t[-1])