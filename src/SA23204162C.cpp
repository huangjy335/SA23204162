#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param thin the number of between-sample random numbers
//' @param a the first parameter of Gibbs sampler proposition distribution of rbeta of y
//' @param b the second parameter of Gibbs sampler proposition distribution of rbeta of y
//' @param n the parameter of Gibbs sampler proposition distribution of rbinom of x
//' @return a random sample of size \code{N}
//' @examples
//' \dontrun{
//' rnC <- GibbsC(1e3,10,2,3,10)
//' par(mfrow=c(2,1));
//' plot(rnC[,1],type='l')
//' plot(rnC[,2],type='l')
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix GibbsC(int N, int thin, int a, int b,int n) {
   NumericMatrix mat(N, 2);
   double x = rbinom(1,0.5,n)[0]; 
   double y = rbeta(1,x+a,n-x+b)[0];
   for(int i = 0; i < N; i++) {
     for(int j = 0; j < thin; j++) {
       x = rbinom(1,y,n)[0];
       y = rbeta(1,x+a,n-x+b)[0];
     }
     mat(i, 0) = x;
     mat(i, 1) = y;
   }
   return(mat);
}
