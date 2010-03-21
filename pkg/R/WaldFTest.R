Wald.F.Test <- function(R, b, r, s2, XtX_1)
{
  v <- R %*% b - r
  as.numeric(t(v) %*% solve( s2 * R %*% XtX_1 %*% t(R) ) %*% v / dim(R)[[1]])
}