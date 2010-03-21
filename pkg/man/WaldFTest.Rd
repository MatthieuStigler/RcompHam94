\name{Wald.F.Test}
\alias{Wald.F.Test}

\title{Wald F Test of a Linear Hypothesis}
\description{
  Compute F test of a linear hypothesis in an OLS regression.
}
\usage{
Wald.F.Test(R, b, r, s2, XtX_1)
}

\arguments{
  \item{R}{Linear combination of coefficients to be tested.}
  \item{b}{Estimated Coefficients.}
  \item{r}{Values that linear combination of coefficients are supposed to take on.}
  \item{s2}{The OLS estimate of variance }
  \item{XtX_1}{The matrix \eqn{(X'X)^{-1}} from the regression.  Available from the \code{summary.lm()} as
  attribute \code{cov.unscaled} }
}

\details{
\code{Wald.F.Test} implements [8.1.32] from page 206 in \emph{Time Series Analysis}.
}

\value{
  The value of the F statistic.
}

\references{    
    Hamilton, J. D. (1994)
    \emph{Time Series Analysis}, Princeton University Press.
}

\seealso{\code{\link[car:linear.hypothesis]{linear.hypothesis}}}

\examples{
## test in a regression on a constant and single regressor that the intercept is 0 and the slope is 1.
lm <- lm( y ~ 1 + x, list(y=rnorm(10) + 1:10,x=1:10))
F <- Wald.F.Test( R=diag(2),
                    b=lm$coefficients,
                    r=c(1,0),
                    s2=summary(lm)$sigma^2,
                    XtX_1=summary(lm)$cov.unscaled )
}

\keyword{regression}
