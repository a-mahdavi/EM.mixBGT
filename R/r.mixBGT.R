r.mixBGT <- function(n , w, xi, s, b, nu, l, u){
	rgnorm <- function (n, mu = 0, alpha = 1, beta = 1) 
	{
    if (alpha <= 0 | beta <= 0) {
        cat("Not defined for negative values of alpha and/or beta.\n")
        return(rep(NaN, n))
    }
    lambda <- (1/alpha)^beta
    unifs <- runif(n)
    scales <- qgamma(unifs, shape = 1/beta, scale = 1/lambda)^(1/beta)
    return(scales * ((-1)^rbinom(n, 1, 0.5)) + mu)
	}
	rBGT=function(n,xi,s,b,nu,l,u){
	y=NULL;i=1
	while(i <=n){
	x=rgamma(1,shape=nu/b,rate=nu/b)^(-1/b)*rgnorm(1, mu = 0, alpha = 1, beta = b)*b^(1/b)
	x=xi+s*x
	if( l<x && x<u){
	y=c(y,x)
	i=i+1
	}
	}
	return(y)
	}
	g <- length(xi) ; y <-NULL ; Z <- rmultinom(n,size=1,prob=w); n <- rowSums(Z)
	for (j in 1:g)
	y <- c(y, rBGT(n[j], xi[j], s[j], b[j], nu[j], l, u ) )
		return(y)
		}







