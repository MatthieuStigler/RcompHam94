\name{Phillips.Perron}
\alias{Phillips.Perron}

\title{Phillips Perron Tests For Non-Stationary}
\description{
  Compute the Phillips Perron rho and t statistics for non-stationarity.
}
\usage{
Phillips.Perron( T, rho, sigma.rho, s, lambda.hat.sq, gamma0 )
}

\arguments{
  \item{T}{Sample size}
  \item{rho}{Estimated autoregressive coefficient of a variable on its lagged value}
  \item{sigma.rho}{Standard error of the estimate of the autoregressive coefficient}
  \item{s}{OLS tandard error of the residuals in the regression}
  \item{lambda.hat.sq}{Estimate of long run variance in residuals }
  \item{gamma0}{Estimate of contemporaneous variance of residuals}
}

\details{
\code{Phillips.Perron} calculates the two statistics described in sections 17.6of
\emph{Time Series Analysis}, specifically equation [17.6.8] and [17.6.12].
}

\value{
  A list containing the following items:
  \item{T, rho, sigma.rho, s, lambda.hat.sq, gamma0}{Values provided to in the call to \code{Phillips.Perron}}
  \item{rho.stat, t.stat}{Calculated Dickey Fuller Statistics}
}

\references{    
    Hamilton, J. D. (1994)
    \emph{Time Series Analysis}, Princeton University Press.
}


\keyword{regression}
