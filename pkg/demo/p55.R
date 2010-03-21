specifications <- list(
    list( label="f = 0", MA=vector(mode="numeric"), AR=vector(mode="numeric") ),
    list( label="f = .5", MA=vector(mode="numeric"), AR=c(.5) ),
    list( label="f = .9", MA=vector(mode="numeric"), AR=c(.9) )
  )
T <- 100
epsilon <- rnorm(T,0,1)


simulate.forward <- function( specification, epsilon )
{
  T <- length(epsilon)
  AR <- specification$AR
  MA <- specification$MA
  presample <- rep(0,max(length(AR),length(MA)))
  epsilon <- c(presample, epsilon)
  Y <- vector(mode="numeric", length=T+length(presample))
  Y[1:length(presample)] <- 0

  for ( i in (length(presample)+1):(T + length(presample)) )
  
    Y[i] <- epsilon[[i]] +
    
        ifelse(length(AR) > 0, t(AR) %*% Y[ (i-1):(i-length(AR))],0) +
        
        ifelse(length(MA) > 0, t(MA) %*% epsilon[(i-1):(i-length(MA))],0)

  Y[(length(presample)+1):(T+length(presample))]
}

for ( i in 1:length(specifications) )
  specifications[[i]]$Y <- simulate.forward( specifications[[i]], epsilon )



screens <- split.screen( figs=c(3,1) )
screen.index <- 1
for ( i in 1:length(specifications) )
{
  screen( n = screens[[i]], new = TRUE )
  specification <- specifications[[i]]
  par( mar=c(4,2,1,2), cex=.75)
  plot( 1:T, specification$Y, type = "l", xlab=specification$label, ylab="", font.lab=5)
  screen.index <- screen.index + 1
}
close.screen(all=TRUE)


for ( specification in specifications )
{
  AR <- specification$AR
  MA <- specification$MA
  shift <- max(length(AR),length(MA))
  Y <- arima.sim(
         model=list(order=c(length(AR),0,length(MA)), ar=AR,ma=MA),
         n=T, innov=epsilon[1:T],
         n.start=max(shift,1), start.innov=rep(0,max(shift,1)) )
  print(specification$Y[1:10])
  print(Y[1:10])
}


