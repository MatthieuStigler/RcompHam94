###################################################
### chunk number 1: 
###################################################
Y <- rt(500, 10)


###################################################
### chunk number 2: 
###################################################
objective <- function(nu, Y) { -sum( log(dt(Y,df=nu)) ) }   # top of p410
classical.results <- optimize( interval=c(1,30), f=objective, Y=Y)
mu2 <- mean(Y^2)        # [14.1.3]
nu <- 2*mu2/(mu2-1)   # [14.1.5]
print( classical.results )
print( nu )


###################################################
### chunk number 3: 
###################################################
compute.estimates <- function( Y, h, interval )
{
  g <- function( Y, THETA )
  {
    apply( X=apply(X=Y,MARGIN=1, FUN=h,THETA=THETA),
            MARGIN=1,
            FUN=mean )  # [14.1.10]
  }
  objective <- function( THETA, Y, W ) {
      g.value <- g( Y, THETA)
      t(g.value) %*% W %*% g.value #[14.1.7]
  }
  r <- length(h(Y[1,], interval[[1]]))  # dummy call to avoid finding the dimension of h
  a <- length(interval[[1]])
  T <- dim(Y)[[1]]
  stage.1.results <- optimize( interval = interval, f=objective, Y=Y, W=diag(r) )
  # S computed using [14.1.18]
  temp <- apply( X=Y, MARGIN=1, FUN=h, THETA=stage.1.results$objective)
  S <- 1/T * temp %*% t(temp)
  stage.2.results <- optimize( interval = interval, f=objective, Y=Y, W=solve(S))
  J.test <- 1 - pchisq(T*stage.2.results$objective,r-a)   #[14.1.27]
  list(stage.1.results=stage.1.results,
      stage.2.results=stage.2.results,
      overidentifying=J.test)
}


###################################################
### chunk number 4: 
###################################################
h <- function( Yt, THETA)
{
  nu <- THETA
  c(Yt^2 - nu/(nu-2), Yt^4 - 3*nu^2 / ((nu - 2) * (nu- 4)))
}
estimates <- compute.estimates( Y %o% 1, h, interval=c(5,30) )
print( estimates )


###################################################
### chunk number 5: 
###################################################
Yg <- rgamma(500,10) * sign( runif(500, -1, 1) )
hg <- function( Yt, THETA )
{
  k <- THETA
  nu <- k
  mu <- k
  sigma <- k
  skew <- 2 / sqrt(k)
  kurt <- 6 / k
  c(Yt^2 - sigma - mu^2, Yt^4 - (kurt*(sigma^2)+3) - 4*(skew*sigma^1.5)*mu - 6*sigma*mu^2 - mu^4)
}

gestimates <- compute.estimates( Yg %o% 1, hg, interval=c(5,30) )
print( gestimates )


