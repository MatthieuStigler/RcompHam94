Dickey.Fuller <- function( T, rho, sigma.rho, zeta=numeric(0) )
{
  list(
    T=T,
    rho=rho,
    sigma.rho=sigma.rho,
    zeta=zeta,
    rho.stat=T * (rho - 1) / (1 - sum(zeta)),
    t.stat=(rho - 1) / sigma.rho
  )
}