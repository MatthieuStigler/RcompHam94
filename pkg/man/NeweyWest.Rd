\name{Newey.West}
\alias{Newey.West}

\title{Newey West Estimator}
\description{
  Compute Newey West heteroskedasticity and autocorrelation consistent estimator of variance.
}
\usage{
Newey.West( X, lags )
}

\arguments{
  \item{X}{T rows by n variables matrix of observations}
  \item{lags}{Number of lags to include in calculation}
}

\details{
\code{Newey.West} calculates the Bartlett kernel based estimator described in section 10.5, specifically
equation [10.5.15].
}

\value{
  An \code{n x n} matrix of variance.
}

\references{    
    Hamilton, J. D. (1994)
    \emph{Time Series Analysis}, Princeton University Press.
}

\examples{
#Newey.West( annual.data.lm$residuals %o% 1, 12 )
}

\keyword{regression}
