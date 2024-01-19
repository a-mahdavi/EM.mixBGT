\name{mixBGT}
\alias{mixBGT}
\title{mixBGT}
\usage{
mixBGT(y, g=1, w=1, xi, s, b, nu, l, u, family="BGT", iter.max=50, tol=10^-6, get.init = TRUE, group=T)
}
\description{
 Fit the mixures of BGT distributions using EM-algorithm
  g: number of components to fit
  w: the vector of probability of each component.
  xi: the of vector of location parameters.
  s: the vector of variances.
  b: the of vector of beta parameters.
  nu: the of vector of flatness parameters.
  l: lower bound and u:upper bound.
  family: distribution family to be used in fitting (default is the proposed BGT otherwise the BT model is fitted by fixing beta=2  ).
  get.init: if TRUE, the initial values are generated.
  iter.max: the maximum number of iterations of the EM algorithm. Default= 50.
  tol: the covergence maximum error. Default= 10^-6.
  group: if TRUE it returns the group id of each observation
 }
\examples{

 #  Example:
 # Simulating 1000 samples from two component BGTM distribution:
   y<-r.mixBGT(n=1000, w=c(.5,.5), xi=c(-.9,.9), s=c(.5,1), b=c(5,.5) , nu=c(5,5),l=-1,u=1)

 # EM output by given initial values: 
   mixBGT(y,w=c(.3,.7),xi=c(-.7,.5),s=c(.5,1),b=c(1,5),nu=c(3,7),get.init=F,g=2,l=-1,u=1)

  # EM output by getting initial values using K-mean algorithm: 
    mixBGT(y,get.init=T,g=2,l=-1,u=1)

