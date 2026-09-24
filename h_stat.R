h_stat <- function (X, y, n) {
  N = sum(n)
  M = length(n)
  
  XX = t(X)%*%X  
  yy = t(y)%*%y
  Xy = t(X)%*%y
  
  X_p = X[seq(from=1, to=2*N-1, by=2), ] + X[seq(from=2, to=2*N, by=2), ]
  y_p = y[seq(from=1, to=2*N-1, by=2), ] + y[seq(from=2, to=2*N, by=2), ]
  XX_p = t(X_p)%*%X_p
  yy_p = t(y_p)%*%y_p
  Xy_p = t(X_p)%*%y_p
  
  idx_g = c(1, sapply(2:M, function(i){2*sum(n[1:(i - 1)]) + 1}))
  X_g = t(sapply(1:M, function(i){colSums(X[idx_g[i]:(idx_g[i]+2*n[i]-1),])}))
  y_g = sapply(1:M, function(i){sum(y[idx_g[i]:(idx_g[i]+2*n[i]-1), ])})
  XX_g = array(0, dim=c(M, dim(X)[2], dim(X)[2]))
  for(i in 1:M) XX_g[i,,] = t(X_g[i, , drop=FALSE])%*%X_g[i, , drop=FALSE]
  yy_g = y_g^2
  Xy_g = t(sapply(1:M, function(i){t(X_g[i,])*y_g[i]}))
  
  return(
    list(XX=XX, yy=yy, Xy=Xy, XX_p=XX_p, yy_p=yy_p, Xy_p=Xy_p,
         XX_g=XX_g, yy_g=yy_g, Xy_g=Xy_g)
  )
}