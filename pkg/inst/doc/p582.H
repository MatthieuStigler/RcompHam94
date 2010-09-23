\subsection{Testing Cointegration when the Cointegrating Vector is Known}
Section 19.2, beginning on page 582 describes cointegration testing of purchasing power parity between
Italian lire and US dollars.  The data used is 100 times log monthly price levels and spot nominal and real
exchange rates, normalized to a value of zero at the start of the series.  
\begin{Scode}{}
data( ppp, package="RcompHam94" )
selection <- window( ppp, start=c(1973,1), end=c(1989,10) )
ppp.data <- cbind(
  pstar=100*log(selection[,"PC6IT"]/selection[[1,"PC6IT"]]),
  p=100*log(selection[,"PZUNEW"]/selection[[1,"PZUNEW"]]),
  ner=-100*log(selection[,"EXRITL"]/selection[[1,"EXRITL"]])
  )
ppp.data <- cbind( ppp.data, rer = ppp.data[,"p"] - ppp.data[,"ner"] - ppp.data[,"pstar"] )
\end{Scode}
\begin{center}
\begin{Scode}{fig=TRUE, echo=FALSE}
plot( index(ppp.data), ppp.data[,"ner"], type="l",lty=3,xlab="Figure 19.2", ylab="", ylim=c(-150,250))
lines(index(ppp.data),ppp.data[,"p"],lty=2)
lines(index(ppp.data),ppp.data[,"pstar"],lty=1)
\end{Scode}
\end{center}
\begin{center}
\begin{Scode}{fig=TRUE, echo=FALSE}
plot(index(ppp.data), ppp.data[,"rer"], type="l",lty=1,xlab="Figure 19.3", ylab="")
\end{Scode}
\end{center}
To save time define a simple utility function to perform augmented Dickey Fuller analysis according to
the conventions in the text.
\begin{Scode}{}
do.DF <- function( series, lag )
{
  df.lms <- summary( dynlm(
    formula=as.formula(paste("y ~ L(y) + tt + L(d(y),1:",lag,") + 1",sep="")),
    data=cbind(y=series, tt=1:length(series))
      ))
  df.results <- Dickey.Fuller(
    T=length(df.lms$residuals),
    rho=df.lms$coefficients[["L(y)","Estimate"]],
    sigma.rho=df.lms$coefficients[["L(y)","Std. Error"]],
    zeta=df.lms$coefficients[paste("L(d(y), 1:",lag,")", 1:lag, sep = ""),"Estimate"] )
  F <- Wald.F.Test( R=cbind( rep(0,2), diag(2), rep(0,2) %o% rep(0,lag) ),
                      b=df.lms$coefficients[,"Estimate"],
                      r=c(1,0),
                      s2=df.lms$sigma^2,
                      XtX_1=df.lms$cov.unscaled )
  print( df.lms$coefficients )
  print( df.results )
  print(F)
}
\end{Scode}
Following the text, check each series with a Dickey Fuller test with a regression estimated with twelve lags.
\begin{Scode}{}
for ( series.name in c( "p", "pstar", "ner", "rer" ) )
  do.DF( series=ppp.data[,series.name], lag=12 )
\end{Scode}
Now check the real exchange rate with a Phillips Perron test

\begin{Scode}{}
pp.lms <- summary(dynlm( z ~ L(z) + 1, zooreg(cbind(z=as.vector(ppp.data[,"rer"]))) ))
PP.results <- Phillips.Perron(
  T=length(pp.lms$residuals),
  rho=pp.lms$coefficients[["L(z)","Estimate"]],
  sigma.rho=pp.lms$coefficients[["L(z)","Std. Error"]],
  s=pp.lms$sigma,
  lambda.hat.sq=as.numeric(Newey.West( pp.lms$residuals %o% 1, 12 )),
  gamma0=mean(pp.lms$residuals^2) )
print( pp.lms$coefficients )
print( PP.results)
\end{Scode}

Estimating the impulse response function gives a sense of the persistence of deviations from PPP.
\begin{Scode}{}
ar.results <- ar(ppp.data$rer, aic = FALSE, order.max = 13, method="ols", demean=TRUE)
tt <- seq(1,72)
start.innov <- rep(0,13)
et <- c(start.innov, 1, rep(0, length(tt) - 14))
arima.sim.output <- arima.sim( list(order=c(13,0,0), ar=ar.results$ar),
     n=length(tt), innov=et, n.start=length(start.innov), start.innov=start.innov )
irf <- as.vector(arima.sim.output)
\end{Scode}
\begin{Scode}{fig=TRUE, echo=FALSE}
plot( tt[-1:-length(start.innov)],irf[-1:-length(start.innov)], type = "l", xlab="Figure 19.4", ylab="")
lines( par("usr")[1:2], c(0,0) )
\end{Scode}
\subsection{Estimating the Cointegrating Vector}
Page 598 shows an example of the Phillips Ouliaris Hansen procedure for estimating a cointegrating vector.
\begin{Scode}{}
poh.cointegration.lm <- lm( p ~ 1 + ner + pstar, ppp.data )
poh.residual.lms <- summary( dynlm( u ~ 0 + L(u),
    zooreg(cbind(u=poh.cointegration.lm$residuals))
   ))
POH.results <- Phillips.Perron( T=length(poh.residual.lms$residuals),
  rho=poh.residual.lms$coefficients[["L(u)","Estimate"]], 
  sigma.rho=poh.residual.lms$coefficients[["L(u)","Std. Error"]], 
  s=poh.residual.lms$sigma, 
  lambda.hat.sq=as.numeric(Newey.West( poh.residual.lms$residuals %o% 1, 12 )), 
  gamma0=mean(poh.residual.lms$residuals^2) )
print( summary(poh.cointegration.lm)$coefficients )
print( poh.residual.lms$coefficients )
print( POH.results)
\end{Scode}
A second example performs a similar analysis on quarterly US consumption and income data from 1947Q1 to 1989Q3.
\begin{Scode}{}
data(coninc, package="RcompHam94")
coninc.data <- window(
	cbind( c = 100*log(coninc[,"GC82"]), y  = 100*log(coninc[,"GYD82"]) ),
	start=c(1947,1), end=c(1989,3) )
coninc.data <- cbind(coninc.data,tt=1:dim(coninc.data)[[1]])
\end{Scode}
\begin{Scode}{fig=TRUE, echo=FALSE}
plot( index(coninc.data), coninc.data[,"c"], type="l",lty=1,xlab="Figure 19.5", ylab="")
lines(index(coninc.data),coninc.data[,"y"],lty=2)
plot( index(coninc.data),coninc.data[,"c"] - coninc.data[,"y"], type="l",lty=1,xlab="Figure 19.6", ylab="")
\end{Scode}
Test individual series for unit root status using Dickey Fuller.
\begin{Scode}{}
for ( series.name in c(  "y", "c") )
  do.DF( series=coninc.data[,series.name], lag=6 )
\end{Scode}
Estimate cointegration vector, then check for unit root status of the residual using Phillips Perron.
\begin{Scode}{}  
poh.cointegration.lm <- lm( c ~ 1 + y, coninc.data )
poh.residual.lms <- summary( dynlm( u ~ 0 + L(u),
    zooreg(cbind(u=poh.cointegration.lm$residuals))
   ))
POH.results <- Phillips.Perron( T=length(poh.residual.lms$residuals),
  rho=poh.residual.lms$coefficients[["L(u)","Estimate"]],
  sigma.rho=poh.residual.lms$coefficients[["L(u)","Std. Error"]],
  s=poh.residual.lms$sigma,
  lambda.hat.sq=as.numeric(Newey.West( poh.residual.lms$residuals %o% 1, 6 )),
  gamma0=mean(poh.residual.lms$residuals^2) )
print( summary(poh.cointegration.lm)$coefficients )
print( poh.residual.lms$coefficients )
print( POH.results)
\end{Scode}
\subsection{Testing Hypotheses About the Cointegrating Vector}
Page 608-612 illustrate a technique that uses leads and lags to produce a stationary
vector for hypothesis testing.
The regression is estimated with both no trend and trend, and the corrected t-stat is calculated.
\begin{Scode}{}
no.trend.lm <- dynlm( c ~ 1 + y + L(d(y),-4:4), coninc.data )
trend.lm <- dynlm( c ~ 1 + y + tt + L(d(y),-4:4), coninc.data )
for ( model in list(no.trend.lm,trend.lm) )
{
  lags <- 2
  cms <- summary(model)
  T <- length(cms$residuals)
  cfs <- cms$coefficients
  t.rho <- (cfs[["y","Estimate"]]-1) / cfs[["y","Std. Error"]]
  rms <- summary(dynlm(
      as.formula(paste("u ~ 0 + L(u,1:",lags,")",sep="")),
      zooreg(cbind(u=as.vector(cms$residuals)))
    ))
  sigma1.hat.sq <- mean(rms$residuals^2)
  lambda.11 <- sigma1.hat.sq^.5 /  (1 - sum(rms$coefficients[ paste("L(u, 1:",lags,")",1:lags,sep=""),  "Estimate"]))
  t.a <- t.rho * cms$sigma / lambda.11
  print(cfs)
  print(rms$coefficients )
  print( T )
  print(cms$sigma)
  print(t.rho)
  print(sigma1.hat.sq)
  print(lambda.11)
  print(t.a)
}
\end{Scode}

