Phillips.Perron <- function( T, rho, sigma.rho, s, lambda.hat.sq, gamma0 )
{
  list(
    T=T,
    rho=rho,
    sigma.rho=sigma.rho,
    s.sq=s^2,
    lambda.hat.sq=lambda.hat.sq,
    gamma0=gamma0,
    rho.stat = T * ( rho - 1 )  - 1/2 * (T * sigma.rho / s ) ^ 2 * ( lambda.hat.sq - gamma0 ),
    t.stat = (gamma0/lambda.hat.sq)^.5 * (rho - 1)/sigma.rho -
      1/2 * (lambda.hat.sq - gamma0) * T * sigma.rho / s / (lambda.hat.sq ^.5)
  )
}