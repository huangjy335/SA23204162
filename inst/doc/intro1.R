## ----eval=FALSE---------------------------------------------------------------
#  intsample<-function(x,size,pro){
#    p <-pro
#    cp <-cumsum(p)
#    m <-size
#    U = runif(m)
#    r <-findInterval(U,cp)+1
#    return(r)
#  }
#  mysample<-function(x,size,pro){
#    if(length(x)==1&&is.numeric(x)&&x>=1){
#      if(missing(size))
#        size<-x
#      intsample(x,size,pro)
#    }
#    else{
#      if(missing(size))
#        size<-length(x)
#      if(missing(pro))
#        pro<-rep(1,length(x))/length(x)
#      x[intsample(length(x),size,pro)]
#    }
#  }
#  invlap<-function(n){
#    x<-numeric(n)
#    u<-runif(n)
#    for (i in 1:n) {
#      if(u[i]>=0.5)
#        x[i]=-log(2-2*u[i])
#      else
#        x[i]=log(2*u[i])
#    }
#    return(x)
#  }
#  mybeta<-function(n,a,b){
#    n<-1e4
#    j<-0
#    y<-numeric(n)
#    while(j<n) {
#      u<-runif(1)
#      x<-runif(1)
#    if(x^(a-1)*(1-x)^(b-1)>u){
#        j<-j+1
#        y[j] <- x
#    }
#      }
#    return(y)
#  }
#  Epan<-function(n){
#    x<-numeric(n)
#    u1<-runif(n,-1,1)
#    u2<-runif(n,-1,1)
#    u3<-runif(n,-1,1)
#    for (i in 1:n) {
#      if(abs(u3[i])>=abs(u2[i])&&abs(u3[i])>=abs(u1[i]))
#        x[i]=u2[i]
#      else
#        x[i]=u3[i]
#    }
#    return(x)
#  }

## ----eval=TRUE----------------------------------------------------------------
set.seed(1234)
library(SA23204162)
x1<-mysample(1:3,1000,c(0.2,0.3,0.5))
x2<-invlap(1000)
x3<-mybeta(1000,3,2)
x4<-Epan(1e4)

## ----eval=TRUE----------------------------------------------------------------
n_times<-1e5
set.seed(1234)
mcsimple<-numeric(1e3)
mcanti<-numeric(1e3)
for (k in 1:1000) {
  mcsimple[k]<-mean(exp(runif(n_times)))
  U_k<-runif(n_times/2)
  one_minus_U_k<-1-U_k
  mcanti[k]<-mean((exp(U_k)+exp(one_minus_U_k))/2)
}
paste("the results are:",round(mean(mcsimple),5),round(mean(mcanti),5),round(exp(1)-1,5))
paste("the reduction in variance is:",round((var(mcsimple)-var(mcanti))/var(mcsimple),6))

## ----eval=TRUE----------------------------------------------------------------
m<-1e6
set.seed(123)
est<-sd<-numeric(2)
g<-function(x)x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)
true_value<-integrate(g,1,Inf)$value
x1<-sqrt(rchisq(m,1))+1
f_1<-2*dnorm(x1,1,1)
f_g_1<-g(x1)/f_1
est[1]<-mean(f_g_1)
sd[1]<-sd(f_g_1)
x2<-rgamma(m,3/2,2)+1
f_2<-dgamma(x2-1,3/2,2)
f_g_2<-g(x2)/f_2
est[2]<-mean(f_g_2)
sd[2]<-sd(f_g_2)
round(c(true_value,est[1],est[2]),5)
round(c(sd[1],sd[2]),5)

## ----eval=TRUE----------------------------------------------------------------
n<-20
N<-50000
set.seed(123)
sqrtn<-sqrt(n)
t_nminus1l<-qt(0.025,n-1)
t_nminus1u<-qt(0.975,n-1)
LCL<-UCL<-numeric(N)
for (i in 1:N) {
  x<-rchisq(n,2)
  LCL[i]<-mean(x)+t_nminus1l*sd(x)/sqrtn
  UCL[i]<-mean(x)+t_nminus1u*sd(x)/sqrtn
}
mean(LCL<2 & UCL>2)

## ----eval=TRUE----------------------------------------------------------------
n<-40
N<-1e4
p1<-p2<-p3<-numeric(N)
miu1_true<-1
miu2_true<-1
miu3_true<-1
set.seed(123)
for (i in 1:N) {
  x1<-rchisq(n,1)
  x2<-runif(n,0,2)
  x3<-rexp(n,1)
  p1[i]<-t.test(x1,alternative="two.sided",mu=miu1_true)$p.value
  p2[i]<-t.test(x2,alternative="two.sided",mu=miu2_true)$p.value
  p3[i]<-t.test(x3,alternative="two.sided",mu=miu3_true)$p.value
}
mean(p1<=0.05)
mean(p2<=0.05)
mean(p3<=0.05)

## ----eval=TRUE----------------------------------------------------------------
lamda<-2
n<-c(5,10,20)
B<-1e3
m<-1e3
library(boot)
set.seed(123)
lamdaf<-function(x,i) 1/mean(x[i])
bias_the<-sd_the<-bias_boot<-sd_boot<-numeric(length(n))
bias_bootm<-sd_bootm<-numeric(m)
for (i in 1:3) {
  bias_the[i]<-lamda/(n[i]-1)
  sd_the[i]<-lamda*n[i]/((n[i]-1)*sqrt(n[i]-2))
  for (j in 1:m) {
   data<-rexp(n[i],lamda)
   om<-boot(data,lamdaf,B)
   bias_bootm[j]<-mean(om$t)-om$t0
   sd_bootm[j]<-sd(om$t)
  }
  bias_boot[i]<-mean(bias_bootm)
  sd_boot[i]<-mean(sd_bootm)
}
re<-data.frame(bias_boot=round(c(bias_boot[1],bias_boot[2],bias_boot[3]),4),
               bias_theory=round(c(bias_the[1],bias_the[2],bias_the[3]),4),
               sd_boot=round(c(sd_boot[1],sd_boot[2],sd_boot[3]),4),
               sd_theory=round(c(sd_the[1],sd_the[2],sd_the[3]),4),
               row.names=c("n=5","n=10","n=20"))
re

## ----eval=TRUE----------------------------------------------------------------
library(bootstrap)
library(boot)
set.seed(123)
m<-10
B<-1e2
n<-nrow(law)
LCL<-UCL<-numeric(m)
boot.cor<-function(x,i) cor(x[i,1],x[i,2])
corr<-boot.cor(law,1:n)
corr_b<-t<-numeric(B)
index<-numeric(n)
for (i in 1:m) {
  for (b in 1:B) {
    index<-sample(1:n,replace=TRUE)
    corr_b[b]<-boot.cor(law,index)
    ob<-boot(law[index,],boot.cor,B)
    t[b]<-(corr_b[b]-corr)/sd(ob$t)
  }
  LCL[i]<-corr-quantile(t,0.975)*sd(corr_b)
  UCL[i]<-corr-quantile(t,0.025)*sd(corr_b)
}
cat("the 95% CI is:","[",round(mean(LCL),4),",",round(mean(UCL),4),"]")

## ----eval=TRUE----------------------------------------------------------------
library(boot)
data1<-aircondit[1]
B<-1e3
m<-1e3
set.seed(123)
lamdaf<-function(x,i) mean(as.matrix(x[i,]))
ci_norm<-ci_basic<-ci_perc<-ci_bca<-matrix(NA,m,2)
for(i in 1:m){
 re<-boot(data1,statistic=lamdaf,B)
 ci<-boot.ci(re,type=c("norm","basic","perc","bca"))
 ci_norm[i,]<-ci$norm[2:3];ci_basic[i,]<-ci$basic[4:5]
 ci_perc[i,]<-ci$percent[4:5];ci_bca[i,]<-ci$bca[4:5]
}
cat('norm =',"[",round(mean(ci_norm[,1]),4) ,",",round(mean(ci_norm[,2]),4),"]",
'basic =','[',round(mean(ci_basic[,1]),4),",",round(mean(ci_basic[,2]),4),']',
'perc =','[',round(mean(ci_perc[,1]),4),",",round(mean(ci_perc[,2]),4),']',
'BCa =','[',round(mean(ci_bca[,1]),4),",",round(mean(ci_bca[,2]),4),']'
)

## ----eval=TRUE----------------------------------------------------------------
library(bootstrap)
data2<-as.matrix(scor)
n<-nrow(data2)
theta_jacm<-numeric(n)
lambda<-eigen(cov(data2))$values
theta_hat<-max(lambda)/sum(lambda)
for (i in 1:n) {
  lambda_n<-eigen(cov(data2[-i,]))$values
  theta_jacm[i]<-max(lambda_n)/sum(lambda_n)
}
bias_jacm<-(n-1)*(mean(theta_jacm)-theta_hat)
se_jacm<-sqrt((n-1)/n*sum((theta_jacm-mean(theta_jacm))^2))
cat('bias','=',round(bias_jacm,4),'sd','=',round(se_jacm,4))

## ----eval=TRUE----------------------------------------------------------------
library(DAAG)
attach(ironslag)
n<-length(magnetic)
N<-n*(n-1)/2
m<-1
err1<-err2<-err3<-err4<-numeric(N)
for (i in 1:(n-1)) {
  for (j in (i+1):n){
   k<-c(i,j)
   y<-magnetic[-k]
   x<-chemical[-k]
   Mod1<-lm(y~x)
   yhat1<-Mod1$coef[1]+Mod1$coef[2]*chemical[k]
   err1[m]<-sum((magnetic[k]-yhat1)^2)
   Mod2<-lm(y~x+I(x^2))
   yhat2<-Mod2$coef[1]+Mod2$coef[2]*chemical[k]+Mod2$coef[3]*chemical[k]^2
   err2[m]<-sum((magnetic[k]-yhat2)^2)
   Mod3<-lm(log(y)~x)
   logyhat3<-Mod3$coef[1]+Mod3$coef[2]*chemical[k]
   yhat3<-exp(logyhat3)
   err3[m]<-sum((magnetic[k]-yhat3)^2)
   Mod4<-lm(log(y)~log(x))
   logyhat4<-Mod4$coef[1]+Mod4$coef[2]*log(chemical[k])
   yhat4<-exp(logyhat4)
   err4[m]<-sum((magnetic[k]-yhat4)^2)
   m<-m+1
  }
}
cat('the errors are:',round(c(mean(err1),mean(err2),mean(err3),mean(err4)),4))

## ----eval=FALSE---------------------------------------------------------------
#  cv.test<-function(x,y){
#    n<-length(x)
#    m<-length(y)
#    z<-c(x,y)
#    N<-n+m
#    Fn<-numeric(N)
#    Gm<-numeric(N)
#    for (i in 1:N) {
#      Fn[i]<-mean(as.integer(z[i]<=x))
#      Gm[i]<-mean(as.integer(z[i]<=y))
#    }
#    return(((n*m)/(N^2))*sum((Fn-Gm)^2))
#  }
#  mycvtest<-function(x,y,R){
#    z<-c(x,y)
#    K<-1:length(z)
#    n<-length(x)
#    set.seed(1234)
#    cv<-numeric(R)
#    cv0<-cv.test(x,y)
#    for (i in 1:R) {
#      k<-sample(K,size=n,replace=FALSE)
#      x1<-z[k]
#      y1<-z[-k]
#      cv[i]<-cv.test(x1,y1)
#    }
#    mean(c(cv0,cv)>=cv0)
#  }
#  cmax.test<-function(x, y) {
#    X<-x-mean(x)
#    Y<-y-mean(y)
#    outx<-sum(X>max(Y))+sum(X<min(Y))
#    outy<-sum(Y>max(X))+sum(Y<min(X))
#    return(max(c(outx, outy)))
#  }
#  mycmaxtest<-function(x,y,R){
#    z<-c(x,y)
#    K<-1:length(z)
#    n<-length(x)
#    co0<-cmax.test(x,y)
#    co<-numeric(R)
#    for (i in 1:R) {
#      k<-sample(K,n,replace=FALSE)
#      x1<-z[k]
#      y1<-z[-k]
#      co[i]<-cmax.test(x1,y1)
#    }
#    mean(c(co0,co)>=co0)
#  }

## ----eval=TRUE----------------------------------------------------------------
set.seed(1234)
library(SA23204162)
attach(chickwts)
x1<-sort(as.vector(weight[feed=="soybean"]))
y1<-sort(as.vector(weight[feed=="linseed"]))
detach(chickwts)
R1<-1e4
mycvtest(x1,y1,R1)

## ----eval=TRUE----------------------------------------------------------------
n1<-20
n2<-30
mu1<-mu2<-0
s1<-s2<-1
x2<-rnorm(n1,mu1,s1)
y2<-rnorm(n2,mu2,s2)
R2<-1e3
mycmaxtest(x2,y2,R2)

## ----eval=FALSE---------------------------------------------------------------
#  Rwalk<-function(sigma){
#    m<-1e4
#    k=0
#    u<-runif(m)
#    x<-numeric(m)
#    x[1]<-0
#    for (i in 2:m) {
#      y<-rnorm(1,x[i-1],sigma)
#      if(u[i]<=(0.5*exp(-abs(y)))/(0.5*exp(-abs(x[i-1])))) x[i]<-y
#      else{
#        x[i]<-x[i-1]
#        k=k+1
#      }
#    }
#    return(list(x,(m-k)/m))
#  }
#  gbR<-function(N,burn){
#    X<-matrix(0,N,2)
#    rho<-0.9
#    mu1<-mu2<-0
#    sigma1<-sigma2<-1
#    va1<-va2<-sqrt(1-rho^2)
#    X[1,]<-c(mu1,mu2)
#   for(i in 2:N){
#   x2<-X[i-1,2]
#   m1<-mu1+rho*(x2-mu2)*sigma1/sigma2
#   X[i,1]<-rnorm(1,m1,va1)
#   x1<-X[i,1]
#   m2<-mu2+rho*(x1-mu1)*sigma2/sigma1
#   X[i,2]<-rnorm(1,m2,va2)
#  }
#  x<-X[burn:N,1]
#  y<-X[burn:N,2]
#  return(list(x,y))
#  }

## ----eval=TRUE----------------------------------------------------------------
set.seed(1234)
library(SA23204162)
sigma<-c(0.05,0.5,5,50)
result<-lapply(sigma,Rwalk)
index<-5000:7000
plot(index,result[[1]][[1]][index],type="l",ylab="x")
plot(index,result[[2]][[1]][index],type="l",col=2,ylab="x")
plot(index,result[[3]][[1]][index],type="l",col=3,ylab="x")
plot(index,result[[4]][[1]][index],type="l",col=4,ylab="x")
result2<-gbR(1e4,1e3)

## ----eval=FALSE---------------------------------------------------------------
#  f_MLE<-function(lamda){
#    U<-c(11,8,27,13,16,0,23,10,24,2)
#    V<-U+1
#    sum<-0
#    for (i in 1:10) {
#      fenzi<-V[i]*exp(-lamda*V[i])-U[i]*exp(-lamda*U[i])
#      fenmu<-exp(-lamda*U[i])-exp(-lamda*V[i])
#      sum=sum+fenzi/fenmu
#    }
#    return(sum)
#  }

## ----eval=TRUE----------------------------------------------------------------
set.seed(12345)
library(SA23204162)
U<-c(11,8,27,13,16,0,23,10,24,2)
V<-U+1
sol1<-uniroot(f_MLE,c(0,10))
lamda1<-round(sol1$root,6)
lamdan<-0.08
for (i in 1:100) {
  lamdan<-lamdan
  sum2<-numeric(10)
  for (i in 1:10) {
    fenzi1=(U[i]-(V[i]+1/lamdan)*exp(lamdan*(U[i]-V[i]))+1/lamdan)
    fenmu1=1-exp(-lamdan*(V[i]-U[i]))
    sum2[i]=fenzi1/fenmu1
  }
  lamdan<-1/mean(sum2)
}
lamda2<-round(lamdan,6)
c(lamda1,lamda2)

## ----eval=FALSE---------------------------------------------------------------
#  GibbsR <- function(N,thin,a,b,n){
#    mat <- matrix(nrow = N, ncol = 2)
#    x <-rbinom(1,prob=0.5,size=n)
#    y <-rbeta(1,x+a,n-x+b)
#    for (i in 1:N) {
#      for (j in 1:thin) {
#        x <- rbinom(1,prob = y, size = n)
#        y <- rbeta(1,x+ a, n - x+ b)
#      }
#      mat[i, ] <- c(x, y)
#    }
#    mat
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix GibbsC(int N, int thin, int a, int b,int n) {
#     NumericMatrix mat(N, 2);
#     double x = rbinom(1,0.5,n)[0];
#     double y = rbeta(1,x+a,n-x+b)[0];
#     for(int i = 0; i < N; i++) {
#       for(int j = 0; j < thin; j++) {
#         x = rbinom(1,y,n)[0];
#         y = rbeta(1,x+a,n-x+b)[0];
#       }
#       mat(i, 0) = x;
#       mat(i, 1) = y;
#     }
#     return(mat);
#  }

## ----eval=TRUE----------------------------------------------------------------
library(SA23204162)
library(microbenchmark)
tm1 <- microbenchmark::microbenchmark(
rnR = GibbsR(1e3,10,2,3,10),
rnC = GibbsC(1e3,10,2,3,10)
 )
 print(summary(tm1)[,c(1,3,5,6)])

