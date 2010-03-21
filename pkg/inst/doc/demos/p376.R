###################################################
### chunk number 1: 
###################################################
kalman <- function( H, R, F, x, A, y, Q, xi.1.0, P.1.0 )
{
  T <- dim(x)[[2]]

  P.t.t_1 <- array( dim=c(dim(P.1.0),T+1) )
  P.t.t_1[,,1] <- P.1.0
  P.t.t <- array( dim=c(dim(P.1.0),T) )
  K.t <- array( dim=c(dim(H),T) )
  xi.t.t_1 <- array(dim=c(length(xi.1.0),T+1) )
  xi.t.t_1[,1] <- xi.1.0
  xi.t.t <- array(dim=c(length(xi.1.0),T) )
  L <- 0
  for ( tt in 1:T )
  {
    V <- solve( t(H) %*% P.t.t_1[,,tt] %*% H + R )
    K.t[,,tt] <- P.t.t_1[,,tt] %*% H %*% V            # [13.2.19] WITHOUT premultiplication by F
    P.t.t[,,tt] <- P.t.t_1[,,tt] - K.t[,,tt] %*%  t(H) %*% P.t.t_1[,,tt]       # [13.2.16]
    P.t.t_1[,,tt+1] <- F %*% P.t.t[,,tt] %*% t(F) + Q         # [13.2.21]
    w <- y[,tt] - t(A) %*% x[,tt] - t(H) %*% xi.t.t_1[,tt]
    xi.t.t[,tt] <- xi.t.t_1[,tt] + K.t[,,tt] %*% w        # [13.2.15]
    xi.t.t_1[,tt+1] <- F %*% xi.t.t[,tt]                # [13.2.17]
    L <- L - 1/2 * dim(y)[[1]] * log(2*pi) + 1/ 2 * log(det(V)) - 1/2 *t(w) %*% V %*% w        # [13.4.1]
  }

  xi.t.T <- array(dim=c(length(xi.1.0),T))
  xi.t.T[,T] <- xi.t.t[,T]
  P.t.T <- array(dim=c(dim(P.1.0),T))
  P.t.T[,,T] <- P.t.t[,,T]
  for ( tt in (T-1):1 )
  {
    Jt <- P.t.t[,,tt] %*% t(F) %*% solve(P.t.t_1[,,tt+1])   # [13.6.11]
    xi.t.T[,tt] <- xi.t.t[,tt] + Jt %*% (xi.t.T[,tt+1] - xi.t.t_1[,tt+1])    # [13.6.16]
    P.t.T[,,tt] <- P.t.t[,,tt] + Jt %*% (P.t.T[,,tt+1] - P.t.t_1[,,tt+1]) %*% t(Jt)   # [13.6.20]
  }
  list( xi.t.t=xi.t.t, xi.t.t_1=xi.t.t_1, P.t.t=P.t.t, P.t.t_1=P.t.t_1, K.t=K.t, log.likelihood=L,
        xi.t.T=xi.t.T, P.t.T=P.t.T)
}


###################################################
### chunk number 2: 
###################################################
data(coninc,package = "Ham94")
YGR <- diff( log(coninc$GYD82) )
CGR <- diff( log(coninc$GC82) )
y <- t(cbind(YGR-mean(YGR),CGR-mean(CGR)))


###################################################
### chunk number 3: 
###################################################
THETA <- c( "phic"=.9, "phi1"=.9, "phi2"=.9, "g1"=.5, "g2"=.5,
                  "sigc"=.05^.5, "sig11"=.05^.5, "sig22"=.05^.5,
                  "r11"=sd(YGR), "r22"=sd(CGR) )
theta.y.to.params <- function( THETA, y )
{
  params <- list(
    F=diag( THETA[c("phic", "phi1", "phi2")] ),
    Q=diag( THETA[c("sigc", "sig11", "sig22")]^2 ),
    H=rbind( THETA[c("g1","g2")], diag(2) ),
    R=diag(THETA[c("r11","r22")]^2),
    A=diag( c(0,0) ),
    x=c(1,1) %o% rep(1,dim(y)[[2]]),
    xi.1.0=c(0,0,0))
  c( params,
    list(P.1.0=array(
      solve( diag(length(params$xi.1.0)^2) - params$F %x% params$F, as.vector(params$Q) ),  # p 381 top
      c(length(params$xi.1.0),length(params$xi.1.0)))
    ))
}


###################################################
### chunk number 4: 
###################################################
objective <- function( THETA, y )
{
  params <- theta.y.to.params( THETA, y )
  kalman( params$H, params$R, params$F, params$x, params$A, y, params$Q, params$xi.1.0, params$P.1.0 )$log.likelihood
}
optimizer.results <- optim( par=THETA, fn=objective, gr=NULL, y=y, control=list(trace=0) )


###################################################
### chunk number 5: 
###################################################
params <- theta.y.to.params( optimizer.results$par, y)
smoothed.results <- kalman( params$H, params$R, params$F, params$x, params$A, y, params$Q, params$xi.1.0, params$P.1.0 )
smoothed.data <- smoothed.results$xi.t.T[1,]


###################################################
### chunk number 6: 
###################################################
plot(coninc$Quarter[-1],(y[1,])/sd(y[1,]),type="l",lty=1,ylab="CGR, YGR, and C")
lines(coninc$Quarter[-1],(y[2,])/sd(y[2,]),lty=2)
lines(coninc$Quarter[-1],(smoothed.data-mean(smoothed.data))/sd(smoothed.data),lty=3)


###################################################
### chunk number 7: 
###################################################
sigmasq <- 2
params <- list(
    F=array( c(0,1,0,0), c(2,2) ), # [13.3.4]
    Q=diag( c(sigmasq,0) ),         # [13.3.6]
    H= array(c(1, .8), c(2,1)),      # [13.3.10]
    R= array(0, c(1,1)),               # [13.3.12]
    A=array( .5, c(1,1) ),            # [13.3.8]
    x=1 %o% rep(1,5),               # [13.3.9]
    y=1 %o% c(1,seq(.5,4)),       # [13.3.7]
    xi.1.0=c(0,0))                      # [13.2.7]
params <-  c( params,
    list(P.1.0=array(
      solve( diag(length(params$xi.1.0)^2) - params$F %x% params$F, as.vector(params$Q) ),  # [10.2.18]
      c(length(params$xi.1.0),length(params$xi.1.0)))   #[13.2.8]
    ))
myResults <- kalman(params$H, params$R, params$F, params$x, params$A, params$y, params$Q, params$xi.1.0, params$P.1.0)


###################################################
### chunk number 8: 
###################################################
fkfResults <- FKF::fkf(
  a0=params$xi.1.0,
  P0=params$P.1.0,
  dt=rep(0,length(params$xi.1.0)) %o% 1,
  Tt=params$F %o% 1,
  HHt=params$Q %o% 1,
  ct=t(params$A) %*% params$x,
  Zt=t(params$H) %o% 1,
  GGt=params$R %o% 1,
  yt=params$y,
  check.input=TRUE)


###################################################
### chunk number 9: 
###################################################
print(myResults$xi.t.t)
print(fkfResults$att)


