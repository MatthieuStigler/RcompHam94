\subsection{Dynamic Multipliers for Second Order Difference Equations}
Finally we use similar techniques to calculate the effects of an impulse on a second order system.
Here each column of phi represents the coefficients of a second order system.
\begin{Scode}{}
T <- 20
w <- 1*(1:20 == 3)
y <- array( dim=c(T,2) )
y[1:2,] <- 0
phi <- array( c(.6,.2, .5,-.8), c(2,2) )
for ( j in 3:T )

  y[j,] <- apply(X=phi * y[(j-1):(j-2),],MARGIN=2,FUN=sum) + w[j]
\end{Scode}
The results can be plotted reproducing figure 1.4.
\begin{center}
\begin{Scode}{fig=TRUE, echo=FALSE}
screens <- split.screen( figs=c(2,1) )
for ( i in 1:2 )
{
  screen( n = screens[i], new = TRUE )
  plot( 1:T,y[,i], type = "h",
        xlab=paste("f1 = ",phi[i,1], ", f2 = ", phi[i,2],sep=""), ylab="",
        font.lab = 5, cex=.5, lwd=5,lend=1,ylim=c(-1,1))
}
close.screen(all=TRUE)
\end{Scode}
\end{center}