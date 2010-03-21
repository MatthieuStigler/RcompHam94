###################################################
### chunk number 1: 
###################################################
arch.fitted.values <- function( THETA, YT )
{
  alpha <- THETA[  grep("alpha.*", names(THETA)) ]
  beta <- THETA[  grep("beta.*", names(THETA)) ]
  zeta <- THETA[  "zeta" ]
  m <- length(alpha)
  h <- rep( as.vector(zeta), m )  
  u <- YT$y - YT$x %*% beta  
  indices <- (m+1):length(YT$y)
  z <- array( 0, c(length(indices), 1 + m) )
  for ( tt in indices )
    h[tt] <- t( c(zeta,alpha) ) %*% c(1, u[ (tt-1):(tt-m) ]^2)    # [21.1.19]
  list( u=u, h=h )
}


###################################################
### chunk number 2: 
###################################################
arch.standard.errors <- function( THETA, YT )
{
  x <- YT$x
  y <- YT$y
  k <- dim(x)[[2]]
  alpha <- THETA[  grep("alpha.*", names(THETA)) ]
  zeta <- THETA[  "zeta" ]  
  m <- length(alpha)
  T <- length(y) - m
  a <- k + 1 + m
  fv <- arch.fitted.values( THETA, YT )
  h <- fv$h
  u2 <- fv$u^2

  S <- array( 0, c(a,a) )
  D <- array( 0, c(a,a) )
  for ( tt in (m+1):length(y) )
  {
    temp <- c(  t(alpha) %*% ( (u2[ (tt-1):(tt-m) ] %o% rep(1,k)) * x[ (tt-1):(tt-m), ]),
                    c(1, u[ (tt-1):(tt-m) ]^2) )      # [21.A.3]
    st <- (u2[tt] - h[tt]) /(2 * h[tt]^2) * temp + c( u2[tt]/h[tt] * x[tt,], rep(0, a-k) ) # [21.1.21]
    S <- S + 1/T * st %*% t(st)
    D <- D + 1/T * (
      1/(2*h[tt]^2) * temp %*% t(temp) +
      rbind(
        cbind( 1/h[tt] * x[tt,] %*% t(x[tt,]), array(0, c( k, a - k )) ),
        array(0, c( a - k, a) )
        )
      ) # [21.1.25]
  }
    
  diag(1/T*solve(D) %*% S %*% solve(D))^.5
}


###################################################
### chunk number 3: 
###################################################

arch.normal <- function(THETA, YT)
{
  fv <- arch.fitted.values( THETA, YT)
  m <- length(THETA[  grep("alpha.*", names(THETA)) ])
  h <- fv$h[-1:-m]
  u <- fv$u[-1:-m]
  -1/2*( length(h) *log(2*pi) - sum(log(h)) - sum(u^2/h) )
}
arch.scaled.t <- function(THETA,YT)
{
  fv <- arch.fitted.values( THETA, YT)
  m <- length(THETA[  grep("alpha.*", names(THETA)) ])
  h <- fv$h[-1:-m]
  u <- fv$u[-1:-m]
  nu <- THETA[  grep("nu", names(THETA)) ]
  result <- length(h) * log(gamma( (nu + 1)/2 ) / (sqrt(pi) * gamma(nu/2)) * (nu - 2)^-.5) -
    1/2 * sum(log(h)) -
    (nu + 1)/2 * sum( log(1 + u^2/(h * (nu - 2))) )
}


###################################################
### chunk number 4: 
###################################################
GMM.estimates <- function( YT, h, THETA, S)
{
  g <- function( YT, THETA )
  {
    apply( X=apply(X=YT,MARGIN=1, FUN=h,THETA=THETA), MARGIN=1, FUN=mean )
  }
  objective <- function( THETA, YT, W ) {
      g.value <- g( YT, THETA)
      as.numeric(t(g.value) %*% W %*% g.value)
  }
  r <- length(h(YT[1,], THETA))
  a <- length(THETA)
  stage.1.results <- optim( par=THETA, fn=objective, gr=NULL, YT=YT, W=diag(r))
  # S computed using [14.1.18]
  temp <- t(apply( X=YT, MARGIN=1, FUN=h, THETA=stage.1.results$par))
  ST <- S(temp)
  stage.2.results <- optim( par=stage.1.results$par, fn=objective, gr=NULL, YT=YT, W=solve(ST))
  list(stage.1.results=stage.1.results, stage.2.results=stage.2.results)
}


###################################################
### chunk number 5: 
###################################################
data(fedfunds, package="Ham94")
selection <- subset( fedfunds, Month >= "1955-01-01" & Month <= "2000-12-01" )
y <- selection$FFED


###################################################
### chunk number 6: 
###################################################
plot( selection$Month, y, type="l",lty=1,xlab="Figure 21.1 - US Fed Funds Rate", ylab="")


###################################################
### chunk number 7: 
###################################################
y.lm <- lm( y ~ 1 + y_1, list(y=y[-1],y_1=y[-length(y)]) )
u <- y.lm$residuals
u2.lm <- lm( u2 ~ 1 + u2_lag, list(u2=u[-1:-4]^2,u2_lag=embed(u[-length(u)]^2,4)) )
F34 <- Wald.F.Test( R=cbind( rep(0,2) %o% rep(0,3), diag(2) ),
                    b=u2.lm$coefficients,
                    r=c(0,0),
                    s2=summary(u2.lm)$sigma^2,
                    XtX_1=summary(u2.lm)$cov.unscaled )
F34.sig <- 1-pf(F34,2,length(u2.lm$residuals)-u2.lm$rank)
F234 <- Wald.F.Test( R=cbind( rep(0,3) %o% rep(0,2), diag(3) ),
                    b=u2.lm$coefficients,
                    r=c(0,0,0),
                    s2=summary(u2.lm)$sigma^2,
                    XtX_1=summary(u2.lm)$cov.unscaled )
F234.sig <- 1-pf(F234,3,length(u2.lm$residuals)-u2.lm$rank)
accept.arch <- pchisq(length(u2.lm$residuals)*summary(u2.lm)$r.squared,4)
print(F34)
print(F34.sig)
print(F234)
print(F234.sig)
print(accept.arch)


###################################################
### chunk number 8: 
###################################################
YT <- list(
  y=y[-1],
  x=cbind(rep(1,length(y)-1), y[-length(y)] ) )
THETA <- c( beta=y.lm$coefficients, zeta=var(y.lm$residuals), alpha=c(.1,.1) )
optimizer.results <- optim( par=THETA, fn=arch.normal, gr=NULL, YT=YT )
print(optimizer.results$par)
se <- arch.standard.errors( optimizer.results$par, YT )
print(se)


###################################################
### chunk number 9: 
###################################################
h <- function( wt, THETA )
{
  beta <- THETA[  grep("beta.*", names(THETA)) ]
  zeta <- THETA[  "zeta" ]
  alpha <- THETA[  grep("alpha.*", names(THETA)) ]
  m <- length(alpha)
  k <- length(beta)
  yt <- wt[ grep("yt.*", names(wt)) ]
  xt <- wt[ grep("xt.*", names(wt)) ]
  ylagt <- wt[ grep("ylagt.*", names(wt)) ]
  xlagt <- t( array( wt[ grep("xlagt.*", names(wt)) ], c(k, m) ) )
  ut <- yt - t(xt) %*% beta
  zt <- c(1, (ylagt - t(xlagt) %*% beta)^2 )
  c( ut * xt, (ut^2 - t(zt) %*% c(zeta, alpha)) * zt ) # page 664
}
S.estimator <- function( ht )
{
  1/dim(ht)[[1]] * t(ht) %*% ht
}

THETA <- c( beta=y.lm$coefficients, zeta=var(y.lm$residuals), alpha=c(.1,.1) )
m <- length(THETA[  grep("alpha.*", names(THETA)) ])
T <- length(YT$y) - m
w <- as.matrix(data.frame(
  yt=YT$y[-1:-m],
  xt=YT$x[-1:-m,],
  ylagt=embed(YT$y[-(T+m)],m),
  xlagt=embed(YT$x[-(T+m),],m)
  ))

estimates <- GMM.estimates( YT=w, h=h, THETA=THETA, S.estimator )
print(estimates$stage.1.results$par)
print(estimates$stage.2.results$par)


