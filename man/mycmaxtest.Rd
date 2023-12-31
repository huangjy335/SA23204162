% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SA23204162R.R
\name{mycmaxtest}
\alias{mycmaxtest}
\title{A permutation test}
\usage{
mycmaxtest(x, y, R)
}
\arguments{
\item{x}{a vector of type 1 samples}

\item{y}{a vector of type 2 samples}

\item{R}{permutation times}
}
\value{
a p-value
}
\description{
A permutation test for equal variance based on the maximum number of extreme points that applies when sample sizes are not necessarily equal
}
\examples{
\dontrun{
n1<-20
n2<-30
mu1<-mu2<-0
s1<-s2<-1
set.seed(1234)
x<-rnorm(n1,mu1,s1)
y<-rnorm(n2,mu2,s2)
R<-1e3
mycmaxtest(x,y,R)
}
}
