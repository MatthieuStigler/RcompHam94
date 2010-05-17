phi <- .8
T <- 20
w <- 1*cbind( 1:T == 6, 1:T >= 6 )
y <- array( dim=c(T,2) )
y[1:5,] <- 0
for ( j in 6:T )

  y[j,] <- phi * y[j-1,] + w[j,]


par( mfrow=c(2,2) )
for ( i in 1:2 )
{
  top <- max(w[,i], y[,])
  plot( 1:T,w[,i], type = "h", xlab="Value of w", ylab="", cex=.5, lwd=5,lend=1,ylim=c(0,top))
  plot( 1:T,y[,i], type = "h", xlab="Value of y", ylab="", cex=.5, lwd=5,lend=1,ylim=c(0,top))
}


