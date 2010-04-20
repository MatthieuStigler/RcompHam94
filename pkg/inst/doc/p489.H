\subsection{Preamble}
This section uses a few utility functions that follow procedures in the test for testing hypotheses about unit roots.
First is the Newey West estimator described by [10.5.10] and [10.5.15].
\begin{Scode}{}
print(Newey.West)
\end{Scode}
Next are the Dickey Fuller stats described in [17.4.7] and [17.4.9], with an optional correction for
serial correlation defined in [17.7.35] and [17.7.38].
\begin{Scode}{}
print(Dickey.Fuller)
\end{Scode}
The Phillips Perron stats are defined by [17.6.8] and [17.6.12]
\begin{Scode}{}
print(Phillips.Perron)
\end{Scode}
Finally the Wald form of an F test as defined by [8.1.32].
\begin{Scode}{}
print(Wald.F.Test)
\end{Scode}
\subsection{Dickey Fuller Tests for Unit Roots}
Page 489 describes the analysis of nominal
three month U.S. Treasury
yield data from
dataset gnptbill, shown below.
\begin{Scode}{}
data( gnptbill, package="Ham94" )
tbill.data <- data.frame(yt=gnptbill$TBILL[-1],yt_1=gnptbill$TBILL[-length(gnptbill$TBILL)])
\end{Scode}
\begin{center}
\begin{Scode}{fig=TRUE, echo=FALSE}
plot(gnptbill$Quarter,gnptbill$TBILL, type = "l", xlab="Figure 17.2 - Nominal Interest Rate", ylab="")
\end{Scode}
\end{center}
The regression model is shown in [17.4.13], and the results are shown below.
\begin{Scode}{}
case1.lms <- summary(lm( yt ~ 0 + yt_1 + 0, tbill.data))
case1.DF <- Dickey.Fuller( T=length(tbill.data$yt),
  rho=case1.lms$coefficients[["yt_1","Estimate"]],
  sigma.rho=case1.lms$coefficients[["yt_1","Std. Error"]] )
print( case1.lms$coefficients)
print( case1.DF )
\end{Scode}
A similar analysis is described on page 494 , but a constant is included in the regression model [17.4.37].
\begin{Scode}{}
case2.lms <- summary(lm( yt ~ 1 + yt_1, tbill.data))
case2.DF <- Dickey.Fuller( T=length(tbill.data$yt),
  rho=case2.lms$coefficients[["yt_1","Estimate"]],
  sigma.rho=case2.lms$coefficients[["yt_1","Std. Error"]] )
print( case2.lms$coefficients )
print( case2.DF )
\end{Scode}
Example 17.5 describes how to test the joint hypothesis that the trend coefficient is 0 and the autoregressive
coefficient is 1.
\begin{Scode}{}
F <- Wald.F.Test( R=diag(2),
                      b=case2.lms$coefficients[,"Estimate"],
                      r=c(0,1),
                      s2=case2.lms$sigma^2,
                      XtX_1=case2.lms$cov.unscaled )
print(F)
\end{Scode}
\subsection{Analyzing GNP data}
A similar analysis can be conducted on log real GNP data described beginning on page 501, shown below.
\begin{Scode}{}
logGNP <- 100*log(gnptbill$GNP)
gnp.data <- data.frame( tt=seq(1,length(gnptbill$GNP)-1),
    yt=logGNP[-1],
    yt_1=logGNP[-length(gnptbill$GNP)] )
\end{Scode}
\begin{Scode}{fig=TRUE, echo=FALSE}
plot(gnptbill$Quarter,gnptbill$GNP, type = "l", xlab="Figure 17.3 - Real GNP", ylab="")
\end{Scode}
The regression model here incorporates a time trend, based on the shape of the GDP graph
\begin{Scode}{}
case4.lms <- summary(lm( yt ~ 1 + yt_1 + tt, gnp.data ))
case4.DF <- Dickey.Fuller( T=length(gnp.data$yt),
  rho=case4.lms$coefficients[["yt_1","Estimate"]],
  sigma.rho=case4.lms$coefficients[["yt_1","Std. Error"]] )
print( case4.lms$coefficients )
print( case4.DF )
F <- Wald.F.Test( R=cbind( rep(0,2), diag(2) ),
                      b=case4.lms$coefficients[,"Estimate"],
                      r=c(1,0),
                      s2=case4.lms$sigma^2,
                      XtX_1=case4.lms$cov.unscaled )
print(F)
\end{Scode}

\subsection{Using Phillips Perron Tests}
Examples 17.6 and 17.7 reanalyze the case 2 and case 4 regressions above using the Phillips Perron tests
as shown on pages 511-513.
\begin{Scode}{}
case2.PP <- Phillips.Perron( T=length(case2.lms$residuals),
  rho=case2.lms$coefficients[["yt_1","Estimate"]],
  sigma.rho=case2.lms$coefficients[["yt_1","Std. Error"]],
  s=case2.lms$sigma,
  lambda.hat.sq=as.numeric(Newey.West( case2.lms$residuals %o% 1, 4 )),
  gamma0=mean(case2.lms$residuals^2) )
print( case2.lms$coefficients )
print( case2.PP)
case4.PP <- Phillips.Perron( T=length(case4.lms$residuals),
  rho=case4.lms$coefficients[["yt_1","Estimate"]],
  sigma.rho=case4.lms$coefficients[["yt_1","Std. Error"]],
  s=case4.lms$sigma,
  lambda.hat.sq=as.numeric(Newey.West( case4.lms$residuals %o% 1, 4 )),
  gamma0=mean(case4.lms$residuals^2) )
print( case4.lms$coefficients )
print( case4.PP)
\end{Scode}
\subsection{Augmented Dickey Fuller Tests}
Example 17.8 illustrates incorporates the use of lagged regressors to (putatively) eliminate serial
correlation in the residuals.  The function "embed" is useful for creating lagged regressors.
\begin{Scode}{}
tbill.data <- list(
  it=gnptbill$TBILL[-1:-5],
  delta.it_ = embed(diff(gnptbill$TBILL[-length(gnptbill$TBILL)]),4),
  it_1=gnptbill$TBILL[c(-1:-4, -(length(gnptbill$TBILL):length(gnptbill$TBILL)))]
  )
tbill.lms <- summary(lm( it ~ delta.it_ + 1 + it_1, tbill.data))
tbill.adf <- Dickey.Fuller(
  T=length(gnptbill$TBILL)-5,
  rho=tbill.lms$coefficients[["it_1","Estimate"]],
  sigma.rho=tbill.lms$coefficients[["it_1","Std. Error"]],
  zeta=tbill.lms$coefficients[paste("delta.it", 1:4, sep = "_"),"Estimate"] )
print( tbill.lms$coefficients)
print( tbill.adf )
\end{Scode}
The next test checks whether or not the farthest lag is different from zero, i.e. whether or not the right number
of lags are included in the equation.
\begin{Scode}{}
print( tbill.lms$coefficients[["delta.it_4","t value"]] )
\end{Scode}
Example 17.9 performs a similar analysis for the GNP data. 
\begin{Scode}{}
gnp.data <- list(
  yt=logGNP[-1:-5],
  delta.yt_ = embed(diff(logGNP[-length(logGNP)]),4),
  yt_1=logGNP[c(-1:-4, -(length(logGNP):length(logGNP)))],
  t=6:length(logGNP)
  )
gnp.lms <- summary(lm( yt ~ delta.yt_ + 1 + yt_1 + t, gnp.data))
gnp.adf <- Dickey.Fuller(
  T=length(logGNP)-5,
  rho=gnp.lms$coefficients[["yt_1","Estimate"]],
  sigma.rho=gnp.lms$coefficients[["yt_1","Std. Error"]],
  zeta=gnp.lms$coefficients[paste("delta.yt", 1:4, sep = "_"),"Estimate"] )
F <- Wald.F.Test( R=cbind( rep(0,2) %o% rep(0,5), diag(2) ),
                      b=gnp.lms$coefficients[,"Estimate"],
                      r=c(1,0),
                      s2=gnp.lms$sigma^2,
                      XtX_1=gnp.lms$cov.unscaled )
print( gnp.lms$coefficients )
print( gnp.adf )
print(F)
\end{Scode}
\subsection{Example 17.10 - Bayesian Test of Autoregressive Coefficient}
Page 532 describes a test on the autoregressive coefficient that weights prior probabilities.
\begin{Scode}{}
t.value <- (1 - gnp.lms$coefficients[["yt_1","Estimate"]]) / gnp.lms$coefficients[["yt_1","Std. Error"]]
print( t.value )
print( (1 - pt( t.value, 164 )) / 2 )
\end{Scode}
\subsection{Determining Lag Length}
Page 530 describes an iterative process to determine the correct lag length.  This is easily expressed
in terms of the structures used above.
\begin{Scode}{}
for ( lag in 10:1 )
{
  gnp.lm <- lm( yt ~ delta.yt_ + 1 + yt_1 + t,
    list(
      yt=logGNP[-1:-(lag+1)],
      delta.yt_ = embed(diff(logGNP[-length(logGNP)]),lag),
      yt_1=logGNP[c(-1:-lag, -(length(logGNP):length(logGNP)))],
      t=(lag+2):length(logGNP)
      ))
  if ( summary(gnp.lm)$coefficients[[paste("delta.yt",lag,sep="_"),"Pr(>|t|)"]] < .05 )
    break
}
print(lag)
\end{Scode}
\subsection{R Facilities for Testing Unit Roots}
TBD

