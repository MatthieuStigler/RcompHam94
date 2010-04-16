\name{Dickey.Fuller}
\alias{Dickey.Fuller}

\title{Dickey Fuller Tests For Non-Stationary}
\description{
  Compute the (Augmented) Dickey Fuller rho and t statistics for non-stationarity.
}
\usage{
Dickey.Fuller(T, rho, sigma.rho, zeta=numeric(0) )
}

\arguments{
  \item{T}{Sample size}
  \item{rho}{Estimated autoregressive coefficient of a variable on its lagged value}
  \item{sigma.rho}{Standard error of the estimate of the autoregressive coefficient}
    \item{zeta}{Coefficients of the variable on its lagged changes}
}

\details{
\code{Dickey.Fuller} calculates the two statistics described in sections 17.5 and 17.8 of
\emph{Time Series Analysis}, specifically equation [17.7.35] and [17.7.36].
}

\value{
  A list containing the following items:
  \item{T, rho, sigma.rho}{Values provided to in the call to \code{Dickey.Fuller}}
  \item{rho.stat, t.stat}{Calculated Dickey Fuller Statistics}
}

\references{    
    Hamilton, J. D. (1994)
    \emph{Time Series Analysis}, Princeton University Press.
}

\examples{
data(gnptbill)
tbill.data <- data.frame(yt = gnptbill$TBILL[-1], yt_1 = gnptbill$TBILL[-length(gnptbill$TBILL)])
case1.lms <- summary(lm( yt ~ 0 + yt_1 + 0, tbill.data))
case1.DF <- Dickey.Fuller( T=length(tbill.data$yt),
  rho=case1.lms$coefficients[["yt_1","Estimate"]],
  sigma.rho=case1.lms$coefficients[["yt_1","Std. Error"]] )
}

\keyword{regression}
