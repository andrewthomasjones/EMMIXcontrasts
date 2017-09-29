//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//'@importFrom Rcpp sourceCpp
//'@useDynLib EMMIXwire


#include <vector>
#include <string>
#include <algorithm>
#include <RcppArmadillo.h>


//using namespace arma;
using namespace Rcpp;
//using namespace boost::math;


typedef std::vector<double> stdvec;

Rcpp::NumericVector export_vec(arma::vec y)
{
  Rcpp::NumericVector tmp = Rcpp::wrap(y);
  tmp.attr("dim") = R_NilValue;
  return tmp;
}

Rcpp::NumericVector export_mat(arma::mat y)
{
  Rcpp::NumericMatrix tmp = Rcpp::wrap(y) ;
  tmp.attr("dim") = R_NilValue;
  return(tmp) ;
}


Rcpp::NumericVector export_uvec(arma::uvec y)
{
  Rcpp::NumericVector tmp = Rcpp::wrap(y);
  tmp.attr("dim") = R_NilValue;
  return tmp;
}


arma::vec mvn_norm(arma::mat dat, arma::vec mu_i, arma::mat sigma_i){
  int n = dat.n_rows;
  int m = mu_i.n_elem;
  arma::vec dens = arma::zeros<arma::vec>(n);
  arma::vec mh = arma::zeros<arma::vec>(n);
  
  double front = std::pow(2*arma::datum::pi,(m/-2.0))*std::pow(arma::det(sigma_i),-0.5);
  arma::mat sigma_inv = inv(sigma_i);
  arma::mat diff = dat.each_row() - mu_i.t();
  arma::mat diff_t = diff.t();
  
  
  for(int i=0; i<n; i++){
    arma::mat test = diff.row(i) *sigma_inv* diff_t.col(i);
  }

  dens = arma::exp(-0.5*mh);
  dens *= front;
  
  return(dens);
}

// [[Rcpp::export]]
Rcpp::List estep(const arma::mat& dat, int n, int m, double g, arma::vec pro, arma::mat mu, arma::cube sigma){
  
  arma::vec s_tau = arma::zeros<arma::vec>(g);
  arma::vec pi = arma::ones<arma::vec>(g);
  arma::mat tau(g,n, arma::fill::zeros);
  arma::mat tau2(g,n, arma::fill::zeros);
  
  //arma::mat tau = mtauk;
  double LL = 0;
  
  for(int i=0;i<g;i++){
    arma::vec mu_i = mu.col(i);
    arma::mat sigma_i = sigma.slice(i);
    arma::vec dens = mvn_norm(dat, mu_i, sigma_i);
    tau.row(i)  = pro.at(i)*dens.t();
  }
  
  tau2 = tau.each_row() / sum(tau,0);
  
  
  s_tau = sum(tau,1);
  pi = s_tau/n;
  
  
  
  Rcpp::List ret = List::create(
    Named("tau")= export_mat(tau2),
    Named("loglik")= LL,
    Named("pro")= export_vec(pi)
  );
  
  
  return(ret);
}


