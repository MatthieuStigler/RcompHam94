\subsection{An Application of the Johansen Approach to the PPP data}
Section 20.3 reanalyzes the data used in Chapter 19 using the FIML approach. 
\begin{Scode}{}
data( ppp, package="RcompHam94" )
selection <- window( ppp, start=c(1973,1), end=c(1989,10) )
ppp.data <- cbind(
  pstar=100*log(selection[,"PC6IT"]/selection[[1,"PC6IT"]]),
  p=100*log(selection[,"PZUNEW"]/selection[[1,"PZUNEW"]]),
  ner=-100*log(selection[,"EXRITL"]/selection[[1,"EXRITL"]])
  )
ppp.data <- cbind( ppp.data, rer = ppp.data[,"p"] - ppp.data[,"ner"] - ppp.data[,"pstar"] )
y <- as.matrix(ppp.data)
\end{Scode}
First conduct the auxiliary regressions.  Given that the right hand sides consists of lagged values of the changes in y
for both [20.2.4] and [20.2.5], construct a regression with both lagged y and lagged changes of y as left hand side.
\begin{Scode}{}
delta.y <- diff(y)
lags <- 12
X <- embed(delta.y[-dim(delta.y)[[1]],],lags)
T <- dim(X)[[1]]
n <- dim(y)[[2]]
lhs <- cbind( delta.y[-1:(-lags),], y[c(-1:-lags,-(T+lags+1)),] )
aux.lm <- lm( lhs ~ 1 + X, list( lhs=lhs, X=X ) )
uv <- sapply(summary(aux.lm),FUN=function(x) { x$residuals })
u <- uv[,1:n]
v <- uv[,(n+1):(2*n)]
\end{Scode}
Now calculate the canonical correlations according to [20.2.6], [20.2.7], [20.2.8],
and calculate eigenvalues according to [20.2.9], and log likelihood as in [20.2.10].
Note that u is T rows by n columns
so that ut is the t-th row of matrix u, so only a single inner product, rather than sum of outer products, is needed.
\begin{Scode}{}
SigmaUU <- 1/T * t(u) %*% u
SigmaVV <- 1/T * t(v) %*% v
SigmaUV <- 1/T * t(u) %*% v
eigen.results <- eigen( solve(SigmaVV) %*% t(SigmaUV) %*% solve(SigmaUU) %*% SigmaUV)
lambda <- eigen.results$values
LRT <- -T*sum(log(1-lambda))
print(SigmaUU)
print(SigmaVV)
print(SigmaUV)
print(lambda)
print(T*log(1-lambda))
print(LRT)
\end{Scode}
Finally following page 648, calculate the first cointegrating vector normalized as in [20.3.9], and also
normalized to have unity for the first coefficient.
\begin{Scode}{}
ahat1 <- eigen.results$vectors[,1]
ahat1.tilde <- ahat1 / sqrt( t(ahat1) %*% SigmaVV %*% ahat1 )
ahat1.normal <- ahat1 / ahat1[[1]]
print(ahat1)
print(ahat1.tilde)
print(ahat1.normal)
\end{Scode}
\subsection{Likelihood Ratio Tests on the Cointegration Vector}
Page 649 shows how to conduct hypothesis tests on the cointegration vector.  The follow code implements
[20.3.10] - [20.3.14] and subsequent calculations.
\begin{Scode}{}
D = cbind( c(1, 0, 0), c(0, 0, 1) )
SigmaVV.tilde <- t(D) %*% SigmaVV %*% D
SigmaUV.tilde <- SigmaUV %*% D
eigen.results <- eigen( solve(SigmaVV.tilde) %*% t(SigmaUV.tilde) %*% solve(SigmaUU) %*% SigmaUV.tilde)
lambda.tilde <- eigen.results$values
h <- 1
LRT <- -T*sum(log(1-lambda[1:h])) + T*sum(log(1-lambda.tilde[1:h]))
ahat1.normal.tilde <- eigen.results$vectors[,1] / eigen.results$vectors[,1][[1]]
print(SigmaVV.tilde)
print(SigmaUV.tilde)
print(lambda.tilde)
print(T*log(1-lambda.tilde))
print(LRT)
print(ahat1.normal.tilde)
\end{Scode}
Page 650 shows a second example.
\begin{Scode}{}
h <- 1
D = c(1, -1,-1) %o% 1
SigmaVV.tilde <- t(D) %*% SigmaVV %*% D
SigmaUV.tilde <- SigmaUV %*% D
eigen.results <- eigen( solve(SigmaVV.tilde) %*% t(SigmaUV.tilde) %*% solve(SigmaUU) %*% SigmaUV.tilde)
lambda.tilde <- eigen.results$values
LRT <- -T*sum(log(1-lambda[1:h])) + T*sum(log(1-lambda.tilde[1:h]))
print(SigmaVV.tilde)
print(SigmaUV.tilde)
print(lambda.tilde)
print(T*log(1-lambda.tilde))
print(LRT)
\end{Scode}
