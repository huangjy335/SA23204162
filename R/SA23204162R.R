
intsample<-function(x,size,pro){
  p <-pro
  cp <-cumsum(p)
  m <-size
  U = runif(m)
  r <-findInterval(U,cp)+1
  return(r)
}

#' @title The hw1 Question1
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

#' @title The hw1 Question2 
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

#' @title The hw1 Question3 
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

#' @title The hw1 Question4 
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

#' @title The hw2 Question1 
#' @name hw2Q1
#' @description What value should take to minimize the asymptotic variance of puffen?
#' @examples
#' \dontrun{
#' rho<-c(0.5,0.8,1)
#' set.seed(1234)
#' va_pi<-rep(0,length(rho))
#' K<-100
#' pi_hat<-rep(0,K)
#' for (i in 1:3) {
#'  for (k in 1:100) {
#'    d<-1
#'    l<-rho[i]*d
#'  n_times<-1e6
#'    x_i<-runif(n_times,0,d/2)
#'    y_i<-runif(n_times,0,pi/2)
#'    pi_hat[k]<-2*l/d/mean(l/2*sin(y_i)>x_i)
#'  }
#'  va_pi[i]<-var(pi_hat)
#' }
#' paste("the variance in different rho(0.5,0.8,1.0) is:",va_pi[1],";",va_pi[2],";",va_pi[3])
#' }
#' @importFrom stats runif rnorm
NULL

#' @title The hw2 Question2 
#' @name hw2Q2
#' @description A illustration of antithetic variate approach of Monte Carlo integration
#' @examples
#' \dontrun{
#' result1<-1/2*(1/2*(exp(2)-1)-(exp(1)-1)^2)
#' result2<- -3/4*exp(2)+5/2*exp(1)-5/4
#' print((result1-result2)/result1)
#' }
NULL

#' @title The hw2 Question3
#' @name hw2Q3
#' @description A illustration of compare of the antithetic variate approach and the simple Monte Carlo method
#' @examples
#' \dontrun{
#' n_times<-1e5
#' set.seed(1234)
#' mcsimple<-rep(0,1000)
#' mcanti<-rep(0,1000)
#' for (k in 1:1000) {
#'  mcsimple[k]<-mean(exp(runif(n_times)))
#'  U_k<-runif(n_times/2)
#'  one_minus_U_k<-1-U_k
#'  mcanti[k]<-mean((exp(U_k)+exp(one_minus_U_k))/2)
#'  }
#' paste(round(mean(mcsimple),5),round(mean(mcanti),5),round(exp(1)-1,5))
#' paste(round((var(mcsimple)-var(mcanti))/var(mcsimple),6))
#' }
NULL

#' @title The hw3
#' @name hw3
#' @description A illustration of importance sampling
#' @examples
#' \dontrun{
#' m<-1e6
#' set.seed(123)
#' est<-sd<-numeric(2)
#' g<-function(x)x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)
#' true_value<-integrate(g,1,Inf)$value
#' x1<-sqrt(rchisq(m,1))+1
#' f_1<-2*dnorm(x1,1,1)
#' f_g_1<-g(x1)/f_1
#' est[1]<-mean(f_g_1)
#' sd[1]<-sd(f_g_1)
#' x2<-rgamma(m,3/2,2)+1
#' f_2<-dgamma(x2-1,3/2,2)
#' f_g_2<-g(x2)/f_2
#' est[2]<-mean(f_g_2)
#' sd[2]<-sd(f_g_2)
#' }
NULL

#' @title The hw4 Question1
#' @name hw4Q1
#' @description A illustration compare the p-value adjust methods
#' @examples
#' \dontrun{
#' M<-m<-1000
#' l1<-0.95*m
#' alpha<-0.1
#' FWER_bon<-FWER_bh<-FDR_bon<-FDR_bh<-TPR_bon<-TPR_bh<-numeric(M)
#' set.seed(123)
#' for (i in 1:M) {
#'   p1<-runif(l1)
#'   p2<-rbeta(0.05*m,0.1,1)
#'   p<-c(p1,p2)
#'   p_bon<-p.adjust(p,method='bonferroni')
#'   p_bh<-p.adjust(p,method='fdr')
#'   FWER_bon[i]<-1-all(p_bon[1:l1]>=alpha)
#'   FWER_bh[i]<-1-all(p_bh[1:l1]>=alpha)
#'   FDR_bon[i]<-sum(p_bon[1:l1]<=alpha)/sum(p_bon<=alpha)
#'   FDR_bh[i]<-sum(p_bh[1:l1]<=alpha)/sum(p_bh<=alpha)
#'   TPR_bon[i]<-mean(p_bon[l1:m]<=alpha)
#'   TPR_bh[i]<-mean(p_bh[l1:m]<=alpha)
#' }
#' FW_bon<-mean(FWER_bon)
#' FW_bh<-mean(FWER_bh)
#' FD_bon<-mean(FDR_bon)
#' FD_bh<-mean(FDR_bh)
#' TPRbon<-mean(TPR_bon)
#' TPRbh<-mean(TPR_bh)
#' re<-data.frame(FWER=round(c(FW_bon,FW_bh),4),
#'                FDR=round(c(FD_bon,FD_bh),4),
#'                TPR=round(c(TPRbon,TPRbh),4),row.names=c("bonf","b-h"))
#' re 
#' }
NULL

#' @title The hw4 Question2
#' @name hw4Q2
#' @description A illustration of bootstrap estimate of bias and sd
#' @examples
#' \dontrun{
#' lamda<-2
#' n<-c(5,10,20)
#' B<-1e3
#' m<-1e3
#' set.seed(123)
#' lamdaf<-function(x,i) 1/mean(x[i])
#' bias_the<-sd_the<-bias_boot<-sd_boot<-numeric(length(n))
#' bias_bootm<-sd_bootm<-numeric(m)
#' for (i in 1:3) {
#'   bias_the[i]<-lamda/(n[i]-1)
#'   sd_the[i]<-lamda*n[i]/((n[i]-1)*sqrt(n[i]-2))
#'   for (j in 1:m) {
#'     data<-rexp(n[i],lamda)
#'     om<-boot(data,lamdaf,B)
#'     bias_bootm[j]<-mean(om$t)-om$t0
#'     sd_bootm[j]<-sd(om$t)
#'   }
#'   bias_boot[i]<-mean(bias_bootm)
#'   sd_boot[i]<-mean(sd_bootm)
#' }
#' re<-data.frame(bias_boot=round(c(bias_boot[1],bias_boot[2],bias_boot[3]),4),
#'                bias_theory=round(c(bias_the[1],bias_the[2],bias_the[3]),4),
#'                sd_boot=round(c(sd_boot[1],sd_boot[2],sd_boot[3]),4),
#'                sd_theory=round(c(sd_the[1],sd_the[2],sd_the[3]),4),
#'                row.names=c("n=5","n=10","n=20"))
#' re 
#' }
#' @import boot
NULL

#' @title The hw4 Question3
#' @name hw4Q3
#' @description A illustration of the bootstrap t confidence interval estimate for correlation statistic
#' @examples
#' \dontrun{
#' set.seed(123)
#' m<-10
#' B<-1e2
#' n<-nrow(law)
#' LCL<-UCL<-numeric(m)
#' boot.cor<-function(x,i) cor(x[i,1],x[i,2])
#' corr<-boot.cor(law,1:n)
#' corr_b<-t<-numeric(B)
#' index<-numeric(n)
#' for (i in 1:m) {
#'   for (b in 1:B) {
#'     index<-sample(1:n,replace=TRUE)
#'     corr_b[b]<-boot.cor(law,index)
#'     ob<-boot(law[index,],boot.cor,B)
#'     t[b]<-(corr_b[b]-corr)/sd(ob$t)
#'   }
#'   LCL[i]<-corr-quantile(t,0.975)*sd(corr_b)
#'   UCL[i]<-corr-quantile(t,0.025)*sd(corr_b)
#' }
#' cat("the 95% CI is:","[",round(mean(LCL),4),",",round(mean(UCL),4),"]")
#' }
#' @import boot
#' @import bootstrap
NULL

#' @title The hw5 Question2
#' @name hw5Q2
#' @description A illustration of the jackknife estimates of bias and standard error 
#' @examples
#' \dontrun{
#' attach(scor)
#' data2<-as.matrix(scor)
#' n<-nrow(data2)
#' theta_jacm<-numeric(n)
#' lambda<-eigen(cov(data2))$values
#' theta_hat<-max(lambda)/sum(lambda)
#' for (i in 1:n) {
#'   lambda_n<-eigen(cov(data2[-i,]))$values
#'   theta_jacm[i]<-max(lambda_n)/sum(lambda_n)
#' }
#' bias_jacm<-(n-1)*(mean(theta_jacm)-theta_hat)
#' se_jacm<-sqrt((n-1)/n*sum((theta_jacm-mean(theta_jacm))^2))
#' cat('bias','=',round(bias_jacm,4),'sd','=',round(se_jacm,4))
#' }
#' @import bootstrap
NULL

#' @title The hw5 Question3
#' @name hw5Q3
#' @description A illustration of the leave-two-out cross validation to compare the models
#' @examples
#' \dontrun{
#' attach(ironslag)
#' n<-length(magnetic)
#' N<-n*(n-1)/2
#' m<-1
#' err1<-err2<-err3<-err4<-numeric(N)
#' for (i in 1:(n-1)) {
#'   for (j in (i+1):n){
#'     k<-c(i,j)
#'     y<-magnetic[-k]
#'     x<-chemical[-k]
#'     Mod1<-lm(y~x)
#'     yhat1<-Mod1$coef[1]+Mod1$coef[2]*chemical[k]
#'     err1[m]<-sum((magnetic[k]-yhat1)^2)
#'     Mod2<-lm(y~x+I(x^2))
#'     yhat2<-Mod2$coef[1]+Mod2$coef[2]*chemical[k]+Mod2$coef[3]*chemical[k]^2
#'     err2[m]<-sum((magnetic[k]-yhat2)^2)
#'     Mod3<-lm(log(y)~x)
#'     logyhat3<-Mod3$coef[1]+Mod3$coef[2]*chemical[k]
#'     yhat3<-exp(logyhat3)
#'     err3[m]<-sum((magnetic[k]-yhat3)^2)
#'     Mod4<-lm(log(y)~log(x))
#'     logyhat4<-Mod4$coef[1]+Mod4$coef[2]*log(chemical[k])
#'     yhat4<-exp(logyhat4)
#'     err4[m]<-sum((magnetic[k]-yhat4)^2)
#'     m<-m+1
#'   }
#' }
#' cat('the errors are:',round(c(mean(err1),mean(err2),mean(err3),mean(err4)),4))
#' }
#' @import DAAG
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
#' @title The hw6 Question2
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
#' @title The hw6 Question3
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

#' @title The hw7 Question2
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

#' @title The hw7 Question3
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

#' @title The hw8 Question1
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

#' @title The hw8 Question2
#' @name hw8Q2
#' @description A illustration of the simplex method
#' @examples
#' \dontrun{
#' solve.game <- function(A) {
#'   min.A <- min(A)
#'   A <- A - min.A 
#'   max.A <- max(A)
#'   A <- A / max(A)
#'   m <- nrow(A)
#'   n <- ncol(A)
#'   it <- n^3
#'   a <- c(rep(0, m), 1) 
#'   A1 <- -cbind(t(A), rep(-1, n)) 
#'   b1 <- rep(0, n)
#'   A3 <- t(as.matrix(c(rep(1, m), 0))) 
#'   b3 <- 1
#'   sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
#'                 maxi=TRUE, n.iter=it)
#'   a <- c(rep(0, n), 1) 
#'   A1 <- cbind(A, rep(-1, m))
#'   b1 <- rep(0, m)
#'   A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
#'   b3 <- 1
#'   sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
#'                 maxi=FALSE, n.iter=it)
#'   soln <- list("A" = A * max.A + min.A,
#'                "x" = sx$soln[1:m],
#'                "y" = sy$soln[1:n],
#'                "v" = sx$soln[m+1] * max.A + min.A)
#'   soln
#' }
#' A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
#'                2,0,0,0,-3,-3,4,0,0,
#'                2,0,0,3,0,0,0,-4,-4,
#'                -3,0,-3,0,4,0,0,5,0,
#'                0,3,0,-4,0,-4,0,5,0,
#'                0,3,0,0,4,0,-5,0,-5,
#'                -4,-4,0,0,0,5,0,0,6,
#'                0,0,4,-5,-5,0,0,0,6,
#'                0,0,4,0,0,5,-6,-6,0), 9, 9)
#' B<-A+2
#' sol <- solve.game(B)
#' sol$v
#' round(cbind(sol$x, sol$y), 7)
#' fractions(round(sol$x , 7))
#' }
NULL

#' @title The hw9 Question6
#' @name hw9Q6
#' @description Use R package \code{microbenchmark} to compare the performance of R function \code{GibbsR} and C function \code{GibbsC}.
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