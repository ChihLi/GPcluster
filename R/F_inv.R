F.mixed <- function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=s) )

# provide an initial bracket for the quantile. default is c(-3,3). 
F_inv <- function(p,w,u,s,br=c(-10,10))
{
  G <- function(x) F.mixed(x,w,u,s) - p
  return(uniroot(G,br)$root) 
}