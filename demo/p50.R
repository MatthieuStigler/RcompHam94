T <- 20
specifications <- list(
    list( label="White Noise", MA=vector(mode="numeric"), AR=vector(mode="numeric") ),
    list( label="MA(1)", MA=c(.8), AR=vector(mode="numeric") ),
    list( label="MA(4)", MA=c(-.6,.5,-.5,.3), AR=vector(mode="numeric") ),
    list( label="AR(1) with 0.8", MA=vector(mode="numeric"), AR=c(.8) ),
    list( label="AR(1) with -0.8", MA=vector(mode="numeric"), AR=c(-.8) )
  )
sigmasq <- 1


specifications[[1]]$rho <- c( 1, rep(0, T-1))   # [3.2.2], [3.2.3]


for ( i in 2:3 )
{
  MA <- specifications[[i]]$MA
  q <- length(MA)
  gamma <- vector(mode="numeric", length=T)
  gamma[1] <- sigmasq * t(c(1,MA)) %*% c(1,MA) #[3.3.10]
  for ( j in 1:q )

    gamma[j+1] <- sigmasq * t(MA[j:q]) %*% c(1,MA)[1:(q-j+1)] #[3.3.12]
    
  gamma[(q+2):T] <- 0   #[3.3.12]
  specifications[[i]]$rho <- gamma/gamma[1]
}


for ( i in 4:5)
{
  AR <- specifications[[i]]$AR
  p <- length(AR)
  F <- rbind( AR, cbind(diag(p-1), rep(0, p-1)) )
  gamma <- vector(mode="numeric", length=T)
  gamma[1:p] <- sigmasq * solve(diag(p^2) - F %x% F)[1:p,1]
  for ( j in (p+1):T )
  
    gamma[[j]] <- t(gamma[(j-1):(j-p)]) %*% AR    #[3.4.36]
  specifications[[i]]$rho <- gamma/gamma[1]
}


screens <- split.screen( figs=c(3,2) )
screen.index <- 1
for ( i in 1:length(specifications) )
{
  screen( n = screens[[i]], new = TRUE )
  specification <- specifications[[i]]
  par( mar=c(4,2,1,2), cex=.75)
  plot( 1:T, specification$rho, type = "h", xlab=specification$label, ylab="", lwd=5,lend=1, ylim=c(-1,1))
  screen.index <- screen.index + 1
}
close.screen(all=TRUE)


g3 <- ARMAacf(ar = numeric(0), ma = specifications[[3]]$MA, lag.max = T, pacf = FALSE)
print(specifications[[3]]$rho)
print(g3)
g4<- ARMAacf(ar = specifications[[4]]$AR, ma = numeric(0), lag.max = T - 1, pacf = FALSE)
print(specifications[[4]]$rho)
print(g4)


