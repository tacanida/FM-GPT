#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Eigen::MatrixXd update_eta_cpp2(Eigen::MatrixXd Omega,Eigen::MatrixXd Lambda,Eigen::MatrixXd Z,int p,Eigen::MatrixXd U,Eigen::MatrixXd Sigma,Eigen::MatrixXd B,int n,Eigen::MatrixXd X,Eigen::MatrixXd rand_mat){
Eigen::MatrixXd emat(n,p),SSigma=Sigma.inverse();
for(int i=0;i<n;++i){
Eigen::MatrixXd O=Omega.row(i).asDiagonal();
Eigen::MatrixXd V=(Lambda*O*Lambda.transpose()+SSigma).inverse();
Eigen::MatrixXd M=V*(Lambda*O*(Z.row(i)-U.row(i)).transpose()+SSigma*B.transpose()*X.row(i).transpose());
Eigen::VectorXd mean=M.col(0);
emat.row(i)=(mean+V.llt().matrixL()*rand_mat.row(i).transpose()).transpose();}
return emat;
}

// [[Rcpp::export]]
Eigen::VectorXd gamma_updateC2(Eigen::MatrixXd X,Eigen::MatrixXd Y,Eigen::MatrixXd b,
Rcpp::IntegerVector group,Rcpp::NumericVector pi1,Eigen::MatrixXd amat,Eigen::MatrixXd gmat,Eigen::MatrixXd omat,const Eigen::Map<Eigen::MatrixXd> Sigma,int ncolY,int ncolX){
Eigen::VectorXd gamma(ncolX);
for(int x=0;x<ncolX;++x){
Eigen::MatrixXd gmats=Eigen::MatrixXd::Zero(ncolX,ncolY);gmats.row(x).setConstant(1);
Eigen::MatrixXd gmatns=gmat;gmatns.row(x).setConstant(0);
Eigen::MatrixXd Zs=amat.array()*gmats.array()*omat.array();
Eigen::MatrixXd Zns=amat.array()*gmatns.array()*omat.array();
Eigen::MatrixXd Yns=Y-X*(Zns.array()*b.array()).matrix();
Eigen::MatrixXd XB=X*(Zs.array()*b.array()).matrix();
double c1=-0.5*((((Yns-XB)*Sigma).array()*(Yns-XB).array()).sum());
double c2=-0.5*((Yns*Sigma).array()*Yns.array()).sum();
double prob=pi1[x]/(pi1[x]+(1-pi1[x])*exp(c2-c1));
gamma[x]=R::rbinom(1,prob);gmat.row(x).setConstant(gamma[x]);}
return gamma;}

// [[Rcpp::export]]
Eigen::VectorXd omega_updateC2(Eigen::MatrixXd X,Eigen::MatrixXd Y,Eigen::MatrixXd b,Rcpp::IntegerVector group,double pi2,Eigen::MatrixXd amat,Eigen::MatrixXd gmat,Eigen::MatrixXd omat,const Eigen::Map<Eigen::MatrixXd> Sigma,int ncolY,int ncolX){
Eigen::VectorXd omega(ncolY*ncolX);int t=0;
for(int y=0;y<ncolY;++y)for(int x=0;x<ncolX;++x){
Eigen::MatrixXd om1=Eigen::MatrixXd::Zero(ncolX,ncolY);om1(x,y)=1;
Eigen::MatrixXd om0=omat;om0(x,y)=0;
Eigen::MatrixXd Ynt=Y-X*((amat.array()*gmat.array()*om0.array()*b.array()).matrix());
Eigen::MatrixXd XB=X*((amat.array()*gmat.array()*om1.array()*b.array()).matrix());
double c1=-0.5*((((Ynt-XB)*Sigma).array()*(Ynt-XB).array()).sum());
double c2=-0.5*((Ynt*Sigma).array()*Ynt.array()).sum();
double prob=pi2/(pi2+(1-pi2)*exp(c2-c1));
omega[t]=R::rbinom(1,prob);omat(x,y)=omega[t];t++;}
return omega;}
