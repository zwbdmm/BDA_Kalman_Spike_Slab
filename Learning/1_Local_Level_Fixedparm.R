####################################################
##
## 1_Local_Level_Fixedparm.R
##  
## Durbin and Koopman, page 16.
## Local level model (no trend, random walk states)
## Fixed variance parameters -- no MLE
##
## Nile data.  Described from source:
# ``` The data appearing in
# ``` Table 1 are measurements of the annual volume of discharge from the Nile River at Aswan
# ``` for the years 1871 to 1970. The measurements are meteorologically significant as evidence
# ``` of a possible abrupt change in rainfall regime near the turn of the last cent
##
## so we can think of the "state" alpha as the amount of water in the river
## and the observation y as the amount discharging from some point
##
##
####################################################

library(stats)
plot(Nile, type = "p", pch = 20, main = "Raw rainfall data, as used in Durbin and Koopman")
##
## Point of Kalman: 
##   (1) Find distribution of next state given all current observations
##   (2) Update distribution of next state when a new observation comes in.
##
## We know everything is normally distributed...just need mean and variance
## The means are called ```a``` and the variances are called ```P```
##
##
# =========================================================================
# Function to return a.filter (a_{t|t}), a.t (a_{t+1}), and P.filter, P.t
# =========================================================================
getKalmanParams <- function(Y, sigsq.ep, sigsq.eta, a1, P1){
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

#
# Run it
#
ret <- getKalmanParams(Y = as.numeric(Nile), sigsq.ep = 15099, sigsq.eta = 1469.1, a1 = 0, P1 = 1e7)

# Recreate the plot of P_t
plot(ret$P.t,  type = "l", ylim = c(5000, 17500), ylab = "P_t", main = "Filtered State Variance by year\n Stabilizes due to stationarity")
# Last Ft is not populated, because we only go through n+1
plot(ret$Ft,  type = "l", ylim = c(20000, 32500), ylab = "F_t", main = "Prediction Variance by year\n Stabilizes due to stationarity")
#
# Compare data to model
#
plot(Nile, type = "p", pch = 20, main = "Nile, data and local level model")
lines(1872:1970, ret$A.t[-1])
