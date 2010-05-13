\subsection{Dynamic Multipliers for First Order Difference Equations}
Page 3 describes calculations for dynamic multipliers for first order difference equations.  An example of these
calculations in action is given on page 4.  A simple method to calculate dynamic multipliers is to simulate
the difference equation calculating forward based on an initial shock at time t=1, assuming the value of y at time 0 is 0.
R indexes arrays starting at 1 instead of 0, so subscripts are one more than the convention used in the text, meaning that
the shock will be said to occur at time 2.
\begin{Scode}{}
T <- 20
w <- 1*(1:T == 2)
\end{Scode}
In the examples shown on page 4 there are actually four different equations being simulated,
so we will use a matrix, rather than a vector, to store the results.
\begin{Scode}{}
phis <- c(.8,-.8,1.1,-1.1)
y <- array( dim=c(T,length(phis)) )
y[1,] <- rep(0, length(phis))
for ( j in 2:T )

  y[j,] <- phis * y[j-1,] + w[j]
\end{Scode}
We can check this calculation against the closed form expression on page 3.
\begin{Scode}{}
print( y[2:T,1] )
print( phis[[1]]^seq(0,T-2) )
\end{Scode}
Finally we can plot the results using a histogram plot reproducing figure 1.1.
\begin{center}
\begin{Scode}{fig=TRUE,echo=FALSE}
par( mfrow=c(2,2) )
for ( i in 1:length(phis) )
  plot( 1:T,y[,i], type = "h", xlab=paste("f = ",phis[[i]]), ylab="",
    font.lab = 5, cex=.5, lwd=5,lend=1,ylim=c(-max(y[,i]),max(y[,i])))
\end{Scode}
\end{center}