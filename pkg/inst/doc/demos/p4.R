###################################################
### chunk number 1: 
###################################################
T <- 20
w <- 1*(1:T == 2)


###################################################
### chunk number 2: 
###################################################
phis <- c(.8,-.8,1.1,-1.1)
y <- array( dim=c(T,length(phis)) )
y[1,] <- rep(0, length(phis))
for ( j in 2:T )

  y[j,] <- phis * y[j-1,] + w[j]


###################################################
### chunk number 3: 
###################################################
print( y[2:T,1] )
print( phis[[1]]^seq(0,T-2) )


###################################################
### chunk number 4: 
###################################################
screens <- split.screen( figs=c(2,2) )
for ( i in 1:length(phis) )
{
  screen( n = screens[[i]], new = TRUE )
  plot( 1:T,y[,i], type = "h", xlab=paste("f = ",phis[[i]]), ylab="",
    font.lab = 5, cex=.5, lwd=5,lend=1,ylim=c(-max(y[,i]),max(y[,i])))
}
close.screen(all=TRUE)


