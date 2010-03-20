\subsection{Testing Cointegration when the Cointegrating Vector is Known}
Section 19.2, beginning on page 582 describes cointegration testing of purchasing power parity between
Italian lire and US dollars.  The data used is 100 times log monthly price levels and spot nominal and real
exchange rates, normalized to a value of zero at the start of the series.  
\begin{Scode}{}
data( ppp, package="Ham94" )
selection <- subset( ppp, Month >= "1973-01-01" & Month <= "1989-10-01" )
ppp.data <- data.frame(
  Month=selection$Month,
  pstar=100*log(selection$PC6IT/selection$PC6IT[[1]]),
  p=100*log(selection$PZUNEW/selection$PZUNEW[[1]]),
  ner=-100*log(selection$EXRITL/selection$EXRITL[[1]])
  )
ppp.data[["rer"]] <- ppp.data$p - ppp.data$ner - ppp.data$pstar
\end{Scode}
\begin{center}
\begin{Scode}{fig=TRUE, echo=FALSE}
plot( ppp.data$Month, ppp.data$ner, type="l",lty=3,xlab="Figure 19.2", ylab="", ylim=c(-150,250))
lines(ppp.data$Month,ppp.data$p,lty=2)
lines(ppp.data$Month,ppp.data$pstar,lty=1)
\end{Scode}
\end{center}
\begin{center}
\begin{Scode}{fig=TRUE, echo=FALSE}
plot(ppp.data$Month, ppp.data$rer, type="l",lty=1,xlab="Figure 19.3", ylab="")
\end{Scode}
\end{center}
To save time define a simple utility function to perform augmented Dickey Fuller analysis according to
the conventions in the text.
\begin{Scode}{}
do.DF <- function( series, lag )
{
  T <- length(series)
  df.lms <- summary( lm( yt ~ yt_1 + tt + delta.yt_ + 1,
    list(
      yt=series[-1:-(lag+1)],
      delta.yt_ = embed(diff(series[-T]),lag),
      yt_1=series[c(-1:-lag, -(T:T))],
      tt=(lag+2):T
      )))
  df.results <- Dickey.Fuller(
    T=length(series)- lag - 1,
    rho=df.lms$coefficients[["yt_1","Estimate"]],
    sigma.rho=df.lms$coefficients[["yt_1","Std. Error"]],
    zeta=df.lms$coefficients[paste("delta.yt_", 1:lag, sep = ""),"Estimate"] )
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
  do.DF( series=ppp.data[[series.name]], lag=12 )
\end{Scode}
Now check the real exchange rate with a Phillips Perron test

\begin{Scode}{}
pp.lms <- summary(lm( zt ~ zt_1 + 1,
            data.frame(zt=ppp.data$rer[-1], zt_1=ppp.data$rer[-length(ppp.data$rer)]) ))
PP.results <- Phillips.Perron(
  T=length(pp.lms$residuals),
  rho=pp.lms$coefficients[["zt_1","Estimate"]],
  sigma.rho=pp.lms$coefficients[["zt_1","Std. Error"]],
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
poh.residual.lms <- summary( lm( u ~ 0 + u_1,
  data.frame(u=poh.cointegration.lm$residuals[-1],
                u_1=poh.cointegration.lm$residuals[-length(poh.cointegration.lm$residuals)]) ))
POH.results <- Phillips.Perron( T=length(poh.residual.lms$residuals),
  rho=poh.residual.lms$coefficients[["u_1","Estimate"]], 
  sigma.rho=poh.residual.lms$coefficients[["u_1","Std. Error"]], 
  s=poh.residual.lms$sigma, 
  lambda.hat.sq=as.numeric(Newey.West( poh.residual.lms$residuals %o% 1, 12 )), 
  gamma0=mean(poh.residual.lms$residuals^2) )
print( summary(poh.cointegration.lm)$coefficients )
print( poh.residual.lms$coefficients )
print( POH.results)
\end{Scode}
A second example performs a similar analysis on quarterly US consumption and income data from 1947Q1 to 1989Q4.
\begin{Scode}{}
data(coninc, package="Ham94")
selection <- subset( coninc,  Quarter >= "1947-01-01" & Quarter <= "1989-07-01" )
coninc.data <- data.frame(Quarter=selection$Quarter,
  cons = 100*log(selection$GC82),
  inc  = 100*log(selection$GYD82))
\end{Scode}
\begin{Scode}{fig=TRUE, echo=FALSE}
plot( coninc.data$Quarter, coninc.data$cons, type="l",lty=1,xlab="Figure 19.5", ylab="")
lines(coninc.data$Quarter,coninc.data$inc,lty=2)
plot( coninc.data$Quarter,coninc.data$cons - coninc.data$inc, type="l",lty=1,xlab="Figure 19.6", ylab="")
\end{Scode}
Test individual series for unit root status using Dickey Fuller.
\begin{Scode}{}
for ( series.name in c(  "inc", "cons") )
  do.DF( series=coninc.data[[series.name]], lag=6 )
\end{Scode}
Estimate cointegration vector, then check for unit root status of the residual using Phillips Perron.
\begin{Scode}{}  
poh.cointegration.lm <- lm( cons ~ 1 + inc, coninc.data )
poh.residual.lms <- summary(lm( u ~ 0 + u_1,
  data.frame(u=poh.cointegration.lm$residuals[-1],
  u_1=poh.cointegration.lm$residuals[-length(poh.cointegration.lm$residuals)]) ))
POH.results <- Phillips.Perron( T=length(poh.residual.lms$residuals),
  rho=poh.residual.lms$coefficients[["u_1","Estimate"]],
  sigma.rho=poh.residual.lms$coefficients[["u_1","Std. Error"]],
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
\begin{Scode}{}
T <- length(coninc.data$Quarter)
lead.lag.data <- list(
  ct=coninc.data$cons[c(-1:-5,-((T-3):T))],
  yt=coninc.data$inc[c(-1:-5,-((T-3):T))],
  delta.yt=diff(coninc.data$inc[c(-1:-4,-((T-3):T))]),
  delta.yt_=embed(diff(coninc.data$inc[-((T-4):T)]),4),
  delta.yt.=embed(diff(coninc.data$inc[-1:-5])[(T-6):1],4)[(T-9):1,],
  tt=6:(T-4)
  )
\end{Scode}
The regression is estimated with both no trend and trend, and the corrected t-stat is calculated.
\begin{Scode}{}
no.trend.lm <- lm( ct ~ 1 + yt + delta.yt. + delta.yt + delta.yt_, lead.lag.data )
trend.lm <- lm( ct ~ 1 + yt + tt + delta.yt. + delta.yt + delta.yt_, lead.lag.data )
for ( model in list(no.trend.lm,trend.lm) )
{
  lags <- 2
  cms <- summary(model)
  T <- length(cms$residuals)
  cfs <- cms$coefficients
  t.rho <- (cfs[["yt","Estimate"]]-1) / cfs[["yt","Std. Error"]]
  rms <- summary(lm( u ~ 0 + u_,
      list(u=cms$residuals[-c(1:lags)],u_=embed(cms$residuals[-T],lags))
    ))
  sigma1.hat.sq <- mean(rms$residuals^2)
  lambda.11 <- sigma1.hat.sq^.5 /  (1 - sum(rms$coefficients[paste("u_",1:lags,sep=""),"Estimate"]))
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

