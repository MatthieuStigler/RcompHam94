\subsection{Autocorrelations as a Function of the Moving Average Parameter}
Figure 3.2 is easily generated from the formula for autocorrelations of an MA(1) process.
\begin{Scode}{fig=TRUE}
theta <- (-300:300)*.01
corrs <- theta/(1+theta^2)
plot( theta, corrs, type="l")
grid(nx=2,ny=2)
\end{Scode}