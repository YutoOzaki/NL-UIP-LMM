# Stan simulation - Non-local unit-information prior linear mixed models
## Load library
library(rstan)
library(bayesplot)
library(rpart)

## Load data
datalist = readRDS("data.rds")

## Run analysis
stanfile = c("nl-uip-lmm.stan", "uip-lmm.stan")
stanfitlist = vector(mode="list", length=2)
lnZ = vector(mode="numeric", length=2)

q = 2

for(i in 1:2) {
  if(i == 1) {
    X = datalist$X
  } else {
    X = datalist$X[, -q]
  }
  n = datalist$n
  y = as.matrix(datalist$y, ncol=1)
  
  ## Pre-compute statistics
  source("h_stat.R")
  statlist = h_stat(X, y, n)
  
  ## Stan
  standata = list(
    p=dim(X)[2], N=sum(n), M=length(n), n=n,
    m=rep(0, dim(X)[2]),
    q=q, dlt_q=0,
    X=X, y=c(y),
    XX=statlist$XX, yy=c(statlist$yy), Xy=c(statlist$Xy),
    XX_p=statlist$XX_p, yy_p=c(statlist$yy_p), Xy_p=c(statlist$Xy_p),
    XX_g=statlist$XX_g, yy_g=c(statlist$yy_g), Xy_g=statlist$Xy_g
  )
  
  fit = stan(file=stanfile[i], data=standata, iter=2000, warmup=1000, chains=4)
  stanfitlist[[i]] = fit
  
  ## Log marginal likelihood
  if(i == 1) {
    df_beta = data.frame(
      beta_1=c(extract(fit, permuted=FALSE, pars="beta[1]")),
      beta_2=c(extract(fit, permuted=FALSE, pars="beta[2]")),
      beta_3=c(extract(fit, permuted=FALSE, pars="beta[3]")),
      beta_4=c(extract(fit, permuted=FALSE, pars="beta[4]"))
    )
  } else {
    df_beta = data.frame(
      beta_1=c(extract(fit, permuted=FALSE, pars="beta[1]")),
      beta_3=c(extract(fit, permuted=FALSE, pars="beta[2]")),
      beta_4=c(extract(fit, permuted=FALSE, pars="beta[3]"))
    )
  }
  
  u = cbind(
    data.frame(
      psi=-c(extract(fit, permuted=FALSE, pars="lp__")),
      u=c(extract(fit, permuted=FALSE, pars="u")),
      sgm=c(extract(fit, permuted=FALSE, pars="sgm")),
      s_1=c(extract(fit, permuted=FALSE, pars="s_1")),
      s_2=c(extract(fit, permuted=FALSE, pars="s_2"))
    ),
    df_beta
  )
  
  fit.tree = rpart(psi ~ ., data=u, method="anova", cp=0.008)
  idx_A = unique(fit.tree$where)
  K = length(idx_A)
  lnZ_k = matrix(0, nrow=K, ncol=1)
  for(k in 1:K) {
    u_k = u[fit.tree$where == idx_A[k], ]
    b_k = sapply(u_k, max)
    a_k = sapply(u_k, min)
    idx_psi = order(u_k$psi)
    S = sum(exp(u_k$psi - max(u_k$psi)))
    W = cumsum(exp(u_k$psi[idx_psi] - max(u_k$psi)))
    idx_med = min(which(W/S >= 0.5))
    c_k = u_k$psi[idx_psi[idx_med]]
    lnZ_k[k] = -c_k + sum(log(b_k - a_k))
  }
  lnZ[i] = max(lnZ_k) + log(sum(exp(lnZ_k - max(lnZ_k))))
}

lnBF = lnZ[1] - lnZ[2]

## Print
print(stanfitlist[[1]])
print(stanfitlist[[2]])
print(data.frame(lnZ_1=lnZ[1], lnZ_2=lnZ[2], lnBF=lnBF, BF=exp(lnBF)))