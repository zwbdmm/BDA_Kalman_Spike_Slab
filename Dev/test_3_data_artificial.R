rm(list=ls(all=TRUE))
library('tidyverse')
library('ggplot2') 
library('bsts')     # load the bsts package
# source: http://www.unofficialgoogledatascience.com/2017/07/fitting-bayesian-structural-time-series.html

x1<- c(rnorm(50,70,5),rnorm(50,40,5))
x2<- c(rnorm(20,70,5),rnorm(80,40,5))
x3<- c(rnorm(80,70,5),rnorm(20,40,5))
L_Tr <- c(seq(1:20),seq(19,10),seq(11,30),seq(29,20),seq(19,58))
plot.ts(L_Tr)
ii <- seq(1:100)
SSS <- 100*sin(2*pi*ii/7) # period is 7
plot.ts(SSS)
y<- 0.5*x1 + x2 + 0.8*x3 + rnorm(100,0,4)+SSS+L_Tr
D <- data.frame(x1=x1,x2=x2,x3=x3,y=y)

D$x4 <- rnorm(100,40,1)
D$x5 <- rnorm(100,20,1)
D$x6 <- rnorm(100,10,1)
D$x7 <- rnorm(100,5,1)
D$x8 <- rnorm(100,2,1)

#a<-initial.claims$iclaimsNSA
#View(as.data.frame(a))


ss <- AddLocalLinearTrend(list(), D$y)
ss <- AddSeasonal(ss, D$y, nseasons = 7)
model1 <- bsts(D$y,
               state.specification = ss,
               niter = 1000)

#names(model1)

plot(model1)
plot(model1, "components")

pred1 <- predict(model1, horizon = 12)
plot(pred1, plot.original = 156)

#### Addding regressors:
#### Addding regressors:
#### Addding regressors:
#### Addding regressors:
model2 <- bsts(y ~ .,
               state.specification = ss,
               niter = 1000,
               data = D)

model3 <- bsts(y ~ .,
               state.specification = ss,
               niter = 1000,
               data = D,
               expected.model.size = 3)  # Passed to SpikeSlabPrior.

plot(model2, "comp")
plot(model2, "coef")
plot(model3, "coef")
