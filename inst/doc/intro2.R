## ----eval=TRUE----------------------------------------------------------------
library(MASS)
set.seed(1234)
#variable dimension is 40
pd=36
beta0=c(1,-1,1.5,-2,rep(0,pd))
C=diag(length(beta0))
for (i in 1:length(beta0)) {
  for (j in 1:length(beta0)) {
    if(i==j)
        C[i,j]=0.5
    else
        C[i,j]=0.2
    }
}
#repeat observe 5 times
Xi=mvrnorm(5,rep(0,length(beta0)),C)
e=mvrnorm(5,c(0),diag(1))
#sample size is 3000
n=3000
UX=list(Xi)
yi=Xi%*%beta0+e
Uy=list(yi)
for(i in 1:n){
    Xi=mvrnorm(n=5,rep(0,length(beta0)), C)
    UX[[i]]=Xi
    e=mvrnorm(n=5,c(0),diag(1))
    Uy[[i]]=Xi%*%beta0+e
}

## ----eval=FALSE---------------------------------------------------------------
#  gethn<-function(DXi,betan)
#  {
#    return(DXi%*%betan)
#  }
#  sq<-function(epsiloni,Ai)
#  {return(solve(sqrt(Ai))%*%epsiloni%*%
#            t(epsiloni)%*%solve(sqrt(Ai)))}
#  getGn<-function(DXi,R_bar,Ai)
#  {return(t(DXi)%*%sqrt(Ai)%*%solve(R_bar)
#          %*%sqrt(Ai)%*%DXi)}
#  getgn<-function(Xi,R_bar,Ai)
#  {return(t(Xi)%*%sqrt(Ai)%*%
#            solve(R_bar)%*%sqrt(Ai)%*%Xi)
#  }
#  ASER<-function(UX,Uy,noo){
#    #doo:the width of confidence interval, which influnce the estimate accuracy
#    doo=0.2
#    #the following four parameters are suggested by Chen et al.
#    gama=1
#    dlta=0.55
#    theta=0.75
#    a=0.05
#    DX<-Dy<-list()
#    for (i in 1:noo) {
#      DX[[i]]=UX[[i]]
#      UX[[i]]<-NULL
#      Dy[[i]]=Uy[[i]]
#      Uy[[i]]<-NULL
#    }
#    Ri=diag(dim(DX[[1]])[1])
#    #the following two parameters can control the shape of the final confidence inteval
#      v_n=6
#      an2=1
#      d=doo
#      betan=rep(0,dim(DX[[1]])[2])
#      h_n=lapply(DX,gethn,betan)
#      epsilon_n=h_n
#  f<-function(betan,method.args=list(DX,Dy))
#  {
#    Sn=0
#    h<-lapply(DX, function(x) x%*%betan)
#    for (i in 1:length(DX)) {
#      Sn=Sn+t(t(DX[[i]])%*%(Dy[[i]]-h[[i]]))
#    }
#    c(Sn)
#  }
#      pk=-solve(jacobian(func=f,x=betan,method.args=list(DX,Dy)))%*%f(betan,method.args=list(DX,Dy))
#      betan=betan+pk
#      p0_hat=length(betan)
#      beta_hat=betan
#      tt=1
#      while(v_n>d^2*tt/an2)
#      {
#        On=diag(length(beta_hat))
#        for(i in 1:p0_hat)
#        {
#          for (j in i:length(beta_hat))
#        {
#            On_hat=diag(length(beta_hat))
#        if(abs(beta_hat[j])>0.001&&i!=j&&abs(beta_hat[i])<=0.001)
#          {
#          On_hat[i,i]=0
#          On_hat[i,j]=1
#          On_hat[j,j]=0
#          On_hat[j,i]=1
#          On=On%*%On_hat
#        }
#        }
#        }
#        beta_hat1=On%*%beta_hat
#        #variable select
#        detl=-Inf
#        q=1
#        R_bar=0
#        Gn=0
#        gn=0
#        Ai<-diag(dim(DX[[1]])[1])
#        h_n=lapply(DX,gethn,betan)
#        for (i in 1:length(DX)) {
#          epsilon_n[[i]]<-Dy[[i]]-h_n[[i]]
#        }
#        for (i in 1:length(epsilon_n)) {
#          R_bar=R_bar+lapply(epsilon_n,sq,Ai)[[i]]
#        }
#        for (i in 1:length(DX)) {
#          Gn=Gn+lapply(DX,getGn,R_bar,Ai)[[i]][1:p0_hat,1:p0_hat]
#        }
#        gn=lapply(UX,getgn,R_bar,Ai)
#        for (j in 1:length(UX)) {
#          detj=det(Gn+(gn[[j]][1:p0_hat,1:p0_hat]))
#          if(detj>detl){
#            detl=detj
#            q=j
#          }
#        }
#        DX[length(DX)+1]=UX[q]
#        Dy[length(Dy)+1]=Uy[q]
#        UX[q]=NULL
#        Uy[q]=NULL
#        #compute betan
#        pk=-solve(jacobian(func=f,x=betan,method.args=list(DX,Dy)))%*%f(betan,method.args=list(DX,Dy))
#        betan=betan+pk
#        #compute beta_hat
#        S=0
#        for(i in 1:length(DX)){
#          S=S+DX[[i]]%*%t(DX[[i]])
#        }
#        ev<-eigen(S)
#        lamedabarn<-max(ev$val)
#        lameda_n<-min(ev$val)
#        Lr=sqrt(lamedabarn*log(lamedabarn))*log(log(lamedabarn))^(0.5+a)/lameda_n
#        lb=Lr^(-1*theta)
#        In=matrix(rep(0,length(betan)*length(betan)),length(betan))
#        #QIC
#        mse<-function(Dy,DX,betan){
#          h<-lapply(DX, function(x) x%*%betan)
#          b<-as.list(1:length(DX))
#          return(sum(unlist
#                     (lapply(b,function(x) t(Uy[[x]]-h[[x]])%*%(Uy[[x]]-h[[x]])))))
#        }
#        epsilon=log(mse(Dy,DX,betan)/length(Dy))+2*p0_hat/length(Dy)
#        for (i in 1:length(betan)) {
#          for (j in 1:length(betan)) {
#            if(sqrt(Lr)*lb*(abs(betan[j])^(-1*gama))<epsilon&&i==j)
#              In[i,j]=1
#          }
#        }
#        p0_hat=sum(In)
#        beta_hat=In%*%betan
#        #compute v_n
#        Hn=0
#        Mn=0
#        for (i in 1:length(DX)) {
#          hi=DX[[i]]%*%beta_hat
#          epsiloni<-Dy[[i]]-hi
#          Hn=Hn+t(DX[[i]])%*%DX[[i]]
#          Mn=Mn+t(DX[[i]])%*%epsiloni%*%t(epsiloni)%*%DX[[i]]
#        }
#        ei<-eigen(tt*In%*%solve(Hn%*%solve(Mn)%*%Hn)%*%In)
#        v_n=max(ei$val)
#        #compute Chi square quantile:an
#        #freedom:p0_hat，confidence level:1-a
#        an2<-qchisq(1-a,p0_hat)
#        tt=tt+1
#      }
#  result<-list(beta_hat,tt+noo)
#  return(result)
#  }

## ----eval=TRUE----------------------------------------------------------------
library(numDeriv)
library(SA23204162)
result1<-ASER(UX,Uy,50)
knitr::kable(data.frame("beta_estimate"=result1[[1]],"beta_true"=beta0))

## ----eval=FALSE---------------------------------------------------------------
#  MQLE<-function(UX,Uy,noo){
#    #doo:the width of confidence interval, which influnce the estimate accuracy
#    doo=0.2
#    tt=1
#    #the following four parameters are suggested by Chen et al.
#    gama=1
#    dlta=0.55
#    theta=0.75
#    a=0.05
#    DX<-Dy<-list()
#    for (i in 1:noo) {
#      DX[[i]]=UX[[i]]
#      UX[[i]]<-NULL
#      Dy[[i]]=Uy[[i]]
#      Uy[[i]]<-NULL
#    }
#    Ri=diag(dim(DX[[1]])[1])
#    #the following two parameters can control the shape of the final confidence inteval
#      q=1
#      v_n=6
#      an2=1
#      d=doo
#      betan=rep(0,dim(DX[[1]])[2])
#      h_n=lapply(DX,gethn,betan)
#      epsilon_n=h_n
#  f<-function(betan,method.args=list(DX,Dy))
#  {
#    Sn=0
#    h<-lapply(DX, function(x) x%*%betan)
#    for (i in 1:length(DX)) {
#      Sn=Sn+t(t(DX[[i]])%*%(Dy[[i]]-h[[i]]))
#    }
#    c(Sn)
#  }
#      pk=-solve(jacobian(func=f,x=betan,method.args=list(DX,Dy)))%*%f(betan,method.args=list(DX,Dy))
#      betan=betan+pk
#      p0_hat=length(betan)
#      beta_hat=betan
#      while(v_n>d^2*tt/an2)
#      {
#        DX[length(DX)+1]=UX[q]
#        Dy[length(Dy)+1]=Uy[q]
#        UX[q]=NULL
#        Uy[q]=NULL
#        #compute betan
#        pk=-solve(jacobian(func=f,x=betan,method.args=list(DX,Dy)))%*%f(betan,method.args=list(DX,Dy))
#        betan=betan+pk
#        Hn=0
#        Mn=0
#        for (i in 1:length(DX)) {
#          hi=DX[[i]]%*%beta_hat
#          epsiloni<-Dy[[i]]-hi
#          Hn=Hn+t(DX[[i]])%*%DX[[i]]
#          Mn=Mn+t(DX[[i]])%*%epsiloni%*%t(epsiloni)%*%DX[[i]]
#        }
#        p0_hat=0
#        In=matrix(rep(0,length(betan)*length(betan)),length(betan))
#        for (i in 1:length(betan)) {
#          for (j in 1:length(betan)) {
#            if(abs(betan[j])>0.1&&i==j){
#              In[i,j]=1
#              p0_hat=p0_hat+1
#            }
#          }
#        }
#        ei<-eigen(tt*In%*%solve(Hn%*%solve(Mn)%*%Hn)%*%In)
#        v_n=max(ei$val)
#        #compute Chi square quantile:an
#        #freedom:p0_hat，confidence level:1-a
#        an2<-qchisq(1-a,p0_hat)
#        tt=tt+1
#      }
#  result<-list(beta_hat,tt+noo)
#  return(result)
#  }

## ----eval=TRUE----------------------------------------------------------------
result1<-ASER(UX,Uy,noo=50)
result2<-MQLE(UX,Uy,noo=50)
knitr::kable(data.frame("MQLE_beta"=round(result2[[1]],3),
           "ASE_beta"=round(result1[[1]],3),
           "true_beta"=beta0))

## ----eval=TRUE----------------------------------------------------------------
knitr::kable(data.frame("MQLE_samplesize"=round(result2[[2]],0),
           "ASE_samplesize"=round(result1[[2]],0)))

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix inverseMatrix(NumericMatrix mat) {
#    int n = mat.nrow();
#    NumericMatrix identity(n, n);
#    for (int i = 0; i < n; i++) {
#      identity(i, i) = 1;
#    }
#    for (int i = 0; i < n; i++) {
#      double pivot = mat(i, i);
#      for (int j = 0; j < n; j++) {
#        mat(i, j) /= pivot;
#        identity(i, j) /= pivot;
#      }
#      for (int k = 0; k < n; k++) {
#        if (k != i) {
#          double factor = mat(k, i);
#  
#          for (int j = 0; j < n; j++) {
#            mat(k, j) -= factor * mat(i, j);
#            identity(k, j) -= factor * identity(i, j);
#          }
#        }
#      }
#    }
#    return identity;
#  }

## ----eval=FALSE---------------------------------------------------------------
#  solv<-function(mat){
#      n <- nrow(mat)
#      identity <- diag(n)
#      for (i in 1:n) {
#        pivot <- mat[i, i]
#        mat[i, ] <- mat[i, ] / pivot
#        identity[i, ] <- identity[i, ] / pivot
#        for (k in 1:n) {
#          if (k != i) {
#            factor <- mat[k, i]
#            mat[k, ] <- mat[k, ] - factor * mat[i, ]
#            identity[k, ] <- identity[k, ] - factor * identity[i, ]
#          }
#        }
#      }
#      return(identity)
#    }

## ----eval=TRUE----------------------------------------------------------------
library(SA23204162)
library(microbenchmark)
set.seed(1234)
l<-10
A=diag(l)
for (i in 1:l) {
  for (j in 1:l) {
    if(i==j)
        A[i,j]=0.5
    else
        A[i,j]=0.2
    }
}
tm2 <- microbenchmark::microbenchmark(
rnR = solv(A),
rnC = inverseMatrix(A)
 )
print(summary(tm2)[,c(1,3,5,6)])

