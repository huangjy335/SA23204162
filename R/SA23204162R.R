intsample<-function(x,size,pro){
  p <-pro
  cp <-cumsum(p)
  m <-size
  U = runif(m)
  r <-findInterval(U,cp)+1
  return(r)
}
#' @title Reproduce the function Sample
#' @description Use the inverse transformation method to reproduce the function Sample when replace=TRUE
#' @param x a vector of one or more elements from which to choose, or a positive integer
#' @param size a non-negative integer giving the number of items to choose
#' @param pro description a vector of probability weights for obtaining the elements of the vector being sampled
#' @return a random sample of size \code{size}
#' @examples
#' \dontrun{
#' set.seed(1234)
#' x1<-mysample(1:3,1000,c(0.2,0.3,0.5))
#' x1<-as.vector(table(x1))
#' x1/sum(x1)/c(0.2,0.3,0.5)
#' }
#' @export
mysample<-function(x,size,pro){
  if(length(x)==1&&is.numeric(x)&&x>=1){
    if(missing(size))
      size<-x
    intsample(x,size,pro)
  }
  else{
    if(missing(size))
      size<-length(x)
    if(missing(pro))
      pro<-rep(1,length(x))/length(x)
    x[intsample(length(x),size,pro)]
  }
}

#' @title Generate random sample from the standard Laplace distribution
#' @description Use the inverse transform method to generate a random sample of size \code{n} from the standard Laplace distribution
#' @param n the sample size
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' set.seed(1234)
#' y<-invlap(1000)
#' print(y)
#' }
#' @export
invlap<-function(n){
  x<-numeric(n)
  u<-runif(n)
  for (i in 1:n) {
    if(u[i]>=0.5)
      x[i]=-log(2-2*u[i])
    else
      x[i]=log(2*u[i]) 
  }
  return(x)
}

#' @title Generate random sample from the beta(a,b) distribution 
#' @description A function to generate a random sample of size \code{n} from the Beta(\code{a},\code{b}) distribution by the acceptance-rejection method
#' @param n the sample size
#' @param a the Beta param1
#' @param b the Beta param2
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' set.seed(1234)
#' y<-mybeta(1000,3,2)
#' print(y)
#' }
#' @export
mybeta<-function(n,a,b){
  n<-1e4
  j<-0
  y<-numeric(n)
  while(j<n) {
    u<-runif(1)
    x<-runif(1) 
    if(x^(a-1)*(1-x)^(b-1)>u){
      j<-j+1
      y[j] <- x
    }
  }
  return(y)
}

#' @title Generate random sample from the the distribution which pdf is rescaled Epanechnikov kernel
#' @description A function to generate a random sample of size \code{n} from the the pdf is rescaled Epanechnikov kernel
#' @param n the sample size
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' set.seed(1234)
#' y<-Epan(1e4)
#' print(y)
#' }
#' @export
Epan<-function(n){
  x<-numeric(n)
  u1<-runif(n,-1,1)
  u2<-runif(n,-1,1)
  u3<-runif(n,-1,1)
  for (i in 1:n) {
    if(abs(u3[i])>=abs(u2[i])&&abs(u3[i])>=abs(u1[i]))
      x[i]=u2[i]
    else
      x[i]=u3[i] 
  }
  return(x)
}

#' @import boot
#' @import bootstrap
#' @import DAAG
#' @importFrom stats runif rnorm
NULL

cv.test<-function(x,y){
  n<-length(x)
  m<-length(y)
  z<-c(x,y)
  N<-n+m
  Fn<-numeric(N)
  Gm<-numeric(N)
  for (i in 1:N) {
    Fn[i]<-mean(as.integer(z[i]<=x))
    Gm[i]<-mean(as.integer(z[i]<=y))
  }
  return(((n*m)/(N^2))*sum((Fn-Gm)^2))
}
#' @title Two-sample Cramer-von Mises test
#' @description The two-sample Cramer-von Mises test for equal distributions as a permutation test
#' @param x a vector of type 1 samples
#' @param y a vector of type 2 samples
#' @param R permutation times
#' @return a p-value
#' @examples
#' \dontrun{
#' attach(chickwts)
#' x<-sort(as.vector(weight[feed=="soybean"]))
#' y<-sort(as.vector(weight[feed=="linseed"]))
#' detach(chickwts)
#' R<-1e4
#' mycvtest(x,y,R)
#' }
#' @export
mycvtest<-function(x,y,R){
  z<-c(x,y)
  K<-1:length(z)
  n<-length(x)
  set.seed(1234)
  cv<-numeric(R)
  cv0<-cv.test(x,y)
  for (i in 1:R) {
    k<-sample(K,size=n,replace=FALSE)
    x1<-z[k]
    y1<-z[-k] 
    cv[i]<-cv.test(x1,y1)
  }
  mean(c(cv0,cv)>=cv0)
}

cmax.test<-function(x, y) {
  X<-x-mean(x)
  Y<-y-mean(y)
  outx<-sum(X>max(Y))+sum(X<min(Y))
  outy<-sum(Y>max(X))+sum(Y<min(X))
  return(max(c(outx, outy)))
}
#' @title A permutation test
#' @description  A permutation test for equal variance based on the maximum number of extreme points that applies when sample sizes are not necessarily equal
#' @param x a vector of type 1 samples
#' @param y a vector of type 2 samples
#' @param R permutation times
#' @return a p-value
#' @examples
#' \dontrun{
#' n1<-20
#' n2<-30
#' mu1<-mu2<-0
#' s1<-s2<-1
#' set.seed(1234)
#' x<-rnorm(n1,mu1,s1)
#' y<-rnorm(n2,mu2,s2)
#' R<-1e3
#' mycmaxtest(x,y,R)
#' }
#' @export
mycmaxtest<-function(x,y,R){
  z<-c(x,y)
  K<-1:length(z)
  n<-length(x)
  co0<-cmax.test(x,y)
  co<-numeric(R)
  for (i in 1:R) {
    k<-sample(K,n,replace=FALSE)
    x1<-z[k]
    y1<-z[-k] 
    co[i]<-cmax.test(x1,y1)
  }
  mean(c(co0,co)>=co0)
}

#' @title A random walk Metropolis sampler
#' @description A random walk Metropolis sampler for generating the standard Laplace distribution
#' @param sigma the sd of the proposal distribution which is normal distribution
#' @return a MC chain
#' @examples
#' \dontrun{
#' set.seed(1234)
#' sigma<-c(0.05,0.5,5,50)
#' result<-lapply(sigma,Rwalk)
#' }
#' @export
Rwalk<-function(sigma){
  m<-1e4
  k=0
  u<-runif(m)
  x<-numeric(m)
  x[1]<-0
  for (i in 2:m) {
    y<-rnorm(1,x[i-1],sigma) 
    if(u[i]<=(0.5*exp(-abs(y)))/(0.5*exp(-abs(x[i-1])))) x[i]<-y 
    else{
      x[i]<-x[i-1]
      k=k+1
    }
  }
  return(list(x,(m-k)/m))
}

#' @title A illustration of Gibbs sampler
#' @description A illustration of Gibbs sampler
#' @param N the length of whole chain
#' @param burn the length of pre-burn 
#' @return a MC chain
#' @examples
#' \dontrun{
#' set.seed(1234)
#' result<-gbR(1e4,1e3)
#' }
#' @export
gbR<-function(N,burn){
  X<-matrix(0,N,2)
  rho<-0.9
  mu1<-mu2<-0
  sigma1<-sigma2<-1
  va1<-va2<-sqrt(1-rho^2)
  X[1,]<-c(mu1,mu2)
  for(i in 2:N){
    x2<-X[i-1,2]
    m1<-mu1+rho*(x2-mu2)*sigma1/sigma2
    X[i,1]<-rnorm(1,m1,va1)
    x1<-X[i,1]
    m2<-mu2+rho*(x1-mu1)*sigma2/sigma1
    X[i,2]<-rnorm(1,m2,va2)
  }
  x<-X[burn:N,1]
  y<-X[burn:N,2]
  return(list(x,y))
}

#' @title MLE and EM compare
#' @description A illustration of the MLE and EM method
#' @param lamda the true value
#' @return the compare of  EM and MLE
#' @examples
#' \dontrun{
#' set.seed(12345)
#' U<-c(11,8,27,13,16,0,23,10,24,2)
#' V<-U+1
#' sol1<-uniroot(f_MLE,c(0,10))
#' lamda1<-round(sol1$root,6)
#' lamdan<-0.08
#' for (i in 1:100) {
#'   lamdan<-lamdan
#'   sum2<-numeric(10)
#'   for (i in 1:10) {
#'     fenzi1=(U[i]-(V[i]+1/lamdan)*exp(lamdan*(U[i]-V[i]))+1/lamdan)
#'     fenmu1=1-exp(-lamdan*(V[i]-U[i]))
#'     sum2[i]=fenzi1/fenmu1
#'   }
#'   lamdan<-1/mean(sum2)
#' }
#' lamda2<-round(lamdan,6)
#' c(lamda1,lamda2)
#' }
#' @export 
f_MLE<-function(lamda){
  U<-c(11,8,27,13,16,0,23,10,24,2)
  V<-U+1
  sum<-0
  for (i in 1:10) {
    fenzi<-V[i]*exp(-lamda*V[i])-U[i]*exp(-lamda*U[i])
    fenmu<-exp(-lamda*U[i])-exp(-lamda*V[i])
    sum=sum+fenzi/fenmu
  }
  return(sum)
}

#' @title Use microbenchmark compare the performance of R and Rcpp
#' @name microben
#' @description Use R package microbenchmark to compare the performance of R function \code{GibbsR} and C function \code{GibbsC}.
#' @examples
#' \dontrun{
#' tm1 <- microbenchmark(
#'   rnR = GibbsR(1e3,10,2,3,10),
#'   rnC = GibbsC(1e3,10,2,3,10)
#' )
#' summary(tm1)[,c(1,3,5,6)]
#' }
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats rbinom rbeta
#' @useDynLib SA23204162
NULL

#' @title A illustration of Gibbs sampler using R
#' @description A Gibbs sampler using R
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @param a the first parameter of Gibbs sampler proposition distribution of rbeta of y
#' @param b the second parameter of Gibbs sampler proposition distribution of rbeta of y
#' @param n the parameter of Gibbs sampler proposition distribution of rbinom of x
#' @return a random sample of size \code{N}
#' @examples
#' \dontrun{
#' rnR <- GibbsR(1e3,10,2,3,10)
#' par(mfrow=c(2,1));
#' plot(rnR[,1],type='l')
#' plot(rnR[,2],type='l')
#' }
#' @export
GibbsR <- function(N,thin,a,b,n){
  mat <- matrix(nrow = N, ncol = 2)
  x <-rbinom(1,prob=0.5,size=n)
  y <-rbeta(1,x+a,n-x+b)
  for (i in 1:N) {
    for (j in 1:thin) {
      x <- rbinom(1,prob = y, size = n)
      y <- rbeta(1,x+ a, n - x+ b)
    }
    mat[i, ] <- c(x, y)
  }
  mat
}