Newey.West <- function( X, lags ) {
  S <- 0
  T <- dim(X)[[1]]
  for ( lag in lags:1 )
    S <- S + (lags + 1 - lag) / (lags+1) * t(X[(lag+1):T,]) %*% X[1:(T-lag),]
  1/T*(t(X)%*%X + S + t(S))
}