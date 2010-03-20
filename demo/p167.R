data(indprod,package = "Ham94")
selection <- subset( indprod, Month >= "1947-01-01" & Month <= "1989-11-01")
raw.data <- selection$IPMFG6
logdiff.data <- 100*diff(log(raw.data),lag=1)
yeardiff.data <- 100*diff(log(raw.data),lag=12)


s.Y.omega <- function( omega, gammas, params) # [6.2.5]
{
  1 / (2*pi) * (gammas[[1]] + 2*as.numeric( t( gammas[-1] ) %*% cos( 1:(length(gammas)-1)*omega ) ))
}
s.Y.omega.Bartlett <- function( omega, gammas, params) # [6.2.5]
{
  1 / (2*pi) * (gammas[[1]] + 2*as.numeric( t( (1 - 1:params/(params+1)) *gammas[2:(params+1)] )
                    %*% cos( 1:params*omega ) ))
}
generate.plot.data <- function( values, estimator, params )
{

  T <- length(values)
  acf.covariance <- acf(values, lag.max = T - 1, type = "covariance", plot = FALSE, demean = TRUE)
  sapply( 2 * pi / T * 1:((T - 1)/2), # [6.2.7]
          estimator, as.vector(acf.covariance$acf), params )
}
raw.s.Y.omega <- generate.plot.data( raw.data, s.Y.omega, NULL ) 
logdiff.s.Y.omega <- generate.plot.data( logdiff.data, s.Y.omega.Bartlett, 12 )
yeardiff.s.Y.omega <- generate.plot.data( yeardiff.data, s.Y.omega.Bartlett, 12 )


screens <- split.screen( figs=c(2,2) )
screen( n = screens[[1]], new = TRUE )
par( mar=c(4,2,1,2),cex=.75)
plot(selection$Month, raw.data, type = "l", xlab="Figure 6.3 - FRB IP Index, NSA", ylab="")
screen( n = screens[[2]], new = TRUE )
par( mar=c(4,2,1,2),cex=.75)
plot(1:length(raw.s.Y.omega),raw.s.Y.omega, type = "l", xlab="Figure 6.4 - Value of j", ylab="")
screen( n = screens[[3]], new = TRUE )
par( mar=c(4,2,1,2),cex=.75)
plot(1:length(logdiff.s.Y.omega),logdiff.s.Y.omega, type = "l", xlab="Figure 6.5 - Value of j", ylab="")
screen( n = screens[[4]], new = TRUE )
par( mar=c(4,2,1,2),cex=.75)
plot(1:length(yeardiff.s.Y.omega),yeardiff.s.Y.omega, type = "l", xlab="Figure 6.6 - Value of j", ylab="")
close.screen(all=TRUE)


