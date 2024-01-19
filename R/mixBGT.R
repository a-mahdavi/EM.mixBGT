 mixBGT <- function(y, g=1, w=1, xi, s, b, nu, l, u,family="BGT", iter.max=50, tol=10^-6, get.init = TRUE, group=T){
 # g: number of components to fit
 # w: the vector of components probabilites.
 # xi: the vector of location parameters.
 # s: the variance parameters.
 # b: the vector of beta parameters.
 # nu: the vector of flatness parameters.
 #l: lower bound and u:upper bound.
 # family: distribution family to be used in fitting (default is the proposed BGT otherwise the BT model is fitted by fixing beta=2  ).
 # get.init: if TRUE, the initial values are generated.
 # iter.max: the maximum number of iterations of the EM algorithm. Default= 50.
 # tol: the covergence maximum error. Default= 10^-6.
 # group: if TRUE it returns the group id of each observation  
   begin <- proc.time()[3]
	dGT <- function(y, xi, s, b, nu){
	eta <- abs(y-xi)/s
	d <- log(b)-log(2)-log(nu)/b+lgamma(nu/b+1/b)-lgamma(1/b)-lgamma(nu/b)-log(s)-((nu+1)/b)*log(1+eta^b/nu)
	return(exp(d)) }
	dBGT <- function(y, xi, s, b, nu, l, u){
		d <- dGT(y, xi, s, b, nu)/integrate(dGT,lower=l,upper=u,xi=xi,s=s,b=b,nu=nu,stop.on.error=F)$value
		return(d) }
	dmixBGT <- function(y, w, xi, s, b, nu, l, u){
			d <- 0 ; g <- length(w)
			for ( j in 1:g)
			d <- d + w[j]*dBGT(y, xi[j], s[j], b[j], nu[j], l, u)
			return(d) }
	#Q <- function(y,s1,s2, xi, s, b, nu, l, u){
	#nu/b*log(nu/b)-log(gamma(nu/b))+(1-1/b)*log(b)-log(1/b)-log(s)-s1/b*(abs(y-xi)^b/s^b+nu)+
	#(nu+1)/b*s2-log(integrate(dGT,lower=l,upper=u,xi=xi,s=s,b=b,nu=nu)$value)}
  n <- length(y)   ;     dif <- 1   ;     count <- 0 
  if (get.init == TRUE) {
                 init <- kmeans(y, g,  algorithm="Hartigan-Wong")
                 w <- init$size/n 
                 xi <- as.vector(init$centers)
                 s <- sqrt(init$withinss/init$size) 
 			b=rep(2,g)+runif(g,-1,1) ; nu=runif(g,1,10)
			if(family!="BGT") 
			 b=rep(2,g) 
			  }
   LL <- 1 
  g0 <- function(x,xi,s,b,nu)
  (1+abs((x-xi)/s)^b/nu)^(-(nu+1)/b)
   g1 <- function(x,xi,s,b,nu)
   (nu+1)/(s^b*nu)*(1+abs((x-xi)/s)^b/nu)^(-(nu+1+b)/b)*sign(x-xi)*abs(x-xi)^(b-1)
  g2 <- function(x,xi,s,b,nu){
	eta <- abs(x-xi)/s
	(nu+1)/nu*eta^b*(1+eta^b/nu)^(-(nu+b+1)/b) }
  while ((dif > tol) && (count <= iter.max)) {
  z.hat  <- matrix(0,n,g)
  # E step
  for (j in 1:g){
  z.hat[,j] <- w[j]*dBGT(y,xi[j],s[j], b[j], nu[j], l, u)/dmixBGT(y, w, xi, s , b, nu, l, u)
  eta <- abs(y-xi[j])/s[j] 
  s1 <- (nu[j]+1)/(nu[j]+eta^b[j])
  #s2 <- digamma((nu[j]+1)/b[j])-log((nu[j]+eta^b[j])/b[j])
  # MCE steps
  w[j] <- sum(z.hat[,j])/n
  PB <- integrate(g0,lower=l,upper=u,xi=xi[j],s=s[j],b=b[j],nu=nu[j],stop.on.error=F)$value 
  K1 <- integrate(g1,lower=l,upper=u,xi=xi[j],s=s[j],b=b[j],nu=nu[j],stop.on.error=F)$value/PB 
  y.xi=y-xi[j] ; y.xi[which(y.xi==0)]=10^-6
  xi[j] <- sum(z.hat[,j]*(s1*abs(y.xi)^(b[j]-2)*y-s[j]^b[j]*K1))/sum(z.hat[,j]*(s1*abs(y.xi)^(b[j]-2)))
  y.xi <- y-xi[j]
  K2 <- integrate(g2,lower=l,upper=u,xi=xi[j],s=s[j],b=b[j],nu=nu[j],stop.on.error=F)$value/PB -1
  s[j] <- (sum(z.hat[,j]*s1*abs(y.xi)^b[j])/sum(z.hat[,j]*(1+K2)))^(1/b[j])
	if(family=="BGT"){
  cml <- optim(c(b[j],nu[j]),function(x){
		bj <- x[1] ; nuj=x[2]; b[j]=bj;nu[j]=nuj
		-sum(log(dmixBGT(y,w,xi,s,b,nu, l, u)))
		},method="L-BFGS-B",lower=c(.01,.01),upper=c(50,100))$par
	b[j] <- cml[1] ;   nu[j] <- cml[2] 
	} else 
  nu[j] <- optim(nu[j],function(x){
		nuj <- x[1] ; nu[j]=nuj
		-sum(log(dmixBGT(y,w,xi,s,b,nu, l, u)))
		},method="L-BFGS-B",lower=.01,upper=100)$par
	}
  LL.new <- sum(log(dmixBGT(y,w,xi,s,b,nu,l,u))) # log-likelihood function
  count <- count +1 
  dif <- abs(LL.new/LL-1)
	LL <- LL.new
  }
  if(family=="BGT"){
  aic <- -2 * LL.new + 2 * (4*g+g-1)
  bic <- -2 * LL.new + log(n) * (4*g+g-1)
  edc <- -2 * LL.new + 0.2*sqrt(n) *(4*g+g-1)
	 }
 else {
  aic <- -2 * LL.new + 2 * (3*g+g-1)
  bic <- -2 * LL.new + log(n) * (3*g+g-1) 
  edc <- -2 * LL.new + 0.2*sqrt(n) *(3*g+g-1)
	}
  end <- proc.time()[3]
  time <- end-begin
  obj.out <- list(w=w, xi=xi, s=s, b=b, nu=nu , loglik=LL.new, aic=aic, bic=bic, edc=edc, iter=count,elapsed=as.numeric(time),group = apply(z.hat, 
                  1, which.max))
  if (group==FALSE)
  obj.out <- obj.out[names(obj.out)!="group"]
  obj.out
  }

