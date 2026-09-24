h_simdata <- function(M, lmd) {
  n = rpois(M, lmd) + 1
  N = sum(n)
  
  Z_1 = bdiag(lapply(1:M, function(i){rep(1, 2*n[i])}))
  Z_2 = bdiag(lapply(1:N, function(i){rep(1, 2)}))
  Z = cbind(Z_1, Z_2)
  
  h_designmat <- function(i, n) {
    x = matrix(0, nrow=n*2, ncol=3)
    x[seq(from=2, to=2*n, by=2), i] = 1
    return(x)
  }
  
  X = cbind(rep(1, 2*N), do.call(rbind, lapply(1:M, function(i){h_designmat(i%%3 + 1, n[i])})))
  
  return(list(X=X, Z=Z, Z_1=Z_1, Z_2=Z_2, n=n, N=N, M=M))
}