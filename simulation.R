# Stan simulation - Non-local unit-information prior linear mixed models
## Load library
library(Matrix)
library(rstan)
library(bayesplot)

## Generate simulation data (input)
source("h_simdata.R")
datalist = h_simdata(M=96, lmd=5)

X = datalist$X
Z = datalist$Z
N = datalist$N
M = datalist$M
n = datalist$n

## Simulation parameters
sgm = runif(1, 0.1, 2)
s_1 = runif(1, 0.1, 2)
s_2 = runif(1, 0.1, 2)
be = as.matrix(c(runif(1, 40, 55), runif(dim(X)[2]-1, 3, 12)), ncol=1)

## Generate simulation data (output)
y = as.matrix(X%*%be + Z%*%c(rnorm(M, 0, s_1), rnorm(N, 0, s_2)) + rnorm(N*2, 0, sgm))

## Pre-compute statistics
source("h_stat.R")
statlist = h_stat(X, y, n)

## Stan
standata = list(
  p=dim(X)[2], N=N, M=M, n=n,
  m=rep(0, dim(X)[2]), q=2, dlt_q=0,
  X=X, y=c(y),
  XX=statlist$XX, yy=c(statlist$yy), Xy=c(statlist$Xy),
  XX_p=statlist$XX_p, yy_p=c(statlist$yy_p), Xy_p=c(statlist$Xy_p),
  XX_g=statlist$XX_g, yy_g=c(statlist$yy_g), Xy_g=statlist$Xy_g
)

fit = stan(file="nl-uip-lmm.stan", data=standata, iter=2000, warmup=1000, chains=4)

## Plot
print(fit)
mcmc_trace(fit)