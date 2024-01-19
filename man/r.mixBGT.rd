\name{r.mixBGT}
\alias{r.mixBGT}
\title{r.mixBGT function}
\usage{
r.mixBGT(n , w, xi, s, b, nu, l, u)
}
\description{
 Generating random samples from mixtures of SMSMN distributions 
 n: number of samples.
 w: the vector of probability of each component.
 xi: the of vector of location parameters.
 s: the vector of variances.
 b: the of vector of beta parameters.
 nu: the of vector of flatness parameters.
 l: lower bound and u:upper bound.
	}
\examples{
 # Example:
 # Simulating 1000 samples from two component BGTM distribution:
   y<-r.mixBGT(n=1000,w=c(.5,.5), xi=c(-.9,.9), s=c(.5,1), b=c(5,.5) , nu=c(5,5),l=-1,u=1)




