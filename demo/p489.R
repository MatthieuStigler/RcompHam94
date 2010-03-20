print(Newey.West)


print(Dickey.Fuller)


print(Phillips.Perron)


print(Wald.F.Test)


data( gnptbill, package="Ham94" )
tbill.data <- data.frame(yt=gnptbill$TBILL[-1],yt_1=gnptbill$TBILL[-length(gnptbill$TBILL)])


plot(gnptbill$Quarter,gnptbill$TBILL, type = "l", xlab="Figure 17.2 - Nominal Interest Rate", ylab="")


case1.lms <- summary(lm( yt ~ 0 + yt_1 + 0, tbill.data))
case1.DF <- Dickey.Fuller( T=length(tbill.data$yt),
  rho=case1.lms$coefficients[["yt_1","Estimate"]],
  sigma.rho=case1.lms$coefficients[["yt_1","Std. Error"]] )
print( case1.lms$coefficients)
print( case1.DF )


case2.lms <- summary(lm( yt ~ 1 + yt_1, tbill.data))
case2.DF <- Dickey.Fuller( T=length(tbill.data$yt),
  rho=case2.lms$coefficients[["yt_1","Estimate"]],
  sigma.rho=case2.lms$coefficients[["yt_1","Std. Error"]] )
print( case2.lms$coefficients )
print( case2.DF )


F <- Wald.F.Test( R=diag(2),
                      b=case2.lms$coefficients[,"Estimate"],
                      r=c(0,1),
                      s2=case2.lms$sigma^2,
                      XtX_1=case2.lms$cov.unscaled )
print(F)


logGNP <- 100*log(gnptbill$GNP)
gnp.data <- data.frame( tt=seq(1,length(gnptbill$GNP)-1),
    yt=logGNP[-1],
    yt_1=logGNP[-length(gnptbill$GNP)] )


plot(gnptbill$Quarter,gnptbill$GNP, type = "l", xlab="Figure 17.3 - Real GNP", ylab="")


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


print( tbill.lms$coefficients[["delta.it_4","t value"]] )


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


t.value <- (1 - gnp.lms$coefficients[["yt_1","Estimate"]]) / gnp.lms$coefficients[["yt_1","Std. Error"]]
print( t.value )
print( (1 - pt( t.value, 164 )) / 2 )


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


