#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat update_U_cpp(arma::mat Z, arma::mat Omega, arma::mat Sigma,
                       arma::mat eta, arma::mat Lambda, int n, int Q) {
  arma::mat umat(n,Q);
  arma::mat Sigma_inv = inv(Sigma);
  for(int i=0;i<n;++i){
    arma::vec Omegai=Omega.row(i).t();
    arma::mat psi_i=diagmat(Omegai)+Sigma_inv;
    arma::mat V=inv(psi_i);
    arma::vec M=V*(Omegai % (Z.row(i).t()-Lambda.t()*eta.row(i).t()));
    umat.row(i)=arma::mvnrnd(M,V).t();
  }
  return umat;
}

// [[Rcpp::export]]
arma::vec update_R_cpp(arma::vec response_types, arma::mat Y, arma::mat eta,
                       arma::mat Lambda, arma::mat U, int Q, int n, arma::vec r){
 for(int q=0;q<Q;++q){
  if(response_types(q)==3){
   arma::vec thetaq=eta*Lambda.col(q)+U.col(q);
   arma::vec expthetaq=exp(thetaq);
   arma::vec pp=expthetaq/(1+expthetaq);
   pp.transform([](double v){return std::min(v,0.9999);});
   arma::vec l(n,fill::zeros),Yq=Y.col(q);
   for(int i=0;i<n;++i){
    int Yiq=Yq(i);
    arma::vec pq=r(q)/(r(q)+(arma::linspace(1,Yiq,Yiq))-1);
    arma::vec l1(Yiq,fill::zeros);
    for(int j=0;j<Yiq;++j) l1(j)=R::rbinom(1,Rf_fround(pq(j),6));
    l(i)=sum(l1);
   }
   r(q)=R::rgamma(10+sum(l),1/(1-sum(log(1-pp))));
  }
 }
 return r;
}
