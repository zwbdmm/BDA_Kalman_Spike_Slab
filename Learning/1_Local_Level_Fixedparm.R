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
getKalmanParams(Y, sigsq.ep, sisgsq.et, a1, P1){
    n <- length(Y)
    # Init lists of the params we want
    A.filter <- numeric(n)
    P.filter <- numeric(n)
    A.t <- numeric(n)
    P.t <- numeric(n)
    #
    A.t[1] <- a1
    P.t[1] <- P1
    # 
    # Also save the prediction error, Vt
    vt <- numeric(n)

    for(t in 1:n){
        vt <- Y[t] - A.t[t]
    }
}



