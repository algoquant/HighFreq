////////////////////////////
// Test functions to be included in HighFreq.cpp
////////////////////////////

// Compile this file in R by running this command:
// Rcpp::sourceCpp(file="/Users/jerzy/Develop/HighFreq/sandbox/test_fun.cpp")

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
using namespace arma;
// Use STL
using namespace std;


// Calculate the alpha and beta regression coefficients
// without using the Moore-Penrose inverse arma::pinv().
// [[Rcpp::export]]
arma::mat calc_betas(const arma::mat& response, 
                     const arma::mat& predictor) {
  
  // arma::uword nrows = predictor.n_rows;
  // arma::uword ncols = predictor.n_cols;
  // Response and Predictor means
  arma::mat respmeans = arma::mean(response, 0);
  arma::mat predmeans = arma::mean(predictor, 0);
  arma::mat predz = predictor.each_row() - predmeans;
  
  // Initialize the variables
  arma::mat covarespred = arma::trans(response.each_row()-respmeans)*predz;
  arma::mat covarpred = arma::trans(predz)*predz;
  arma::mat betas = covarespred*arma::inv(covarpred);
  arma::mat alphas = respmeans - arma::dot(betas, predmeans);
  
  return arma::join_rows(alphas, betas);
  
}  // end calc_betas


// Calculate a limited number of eigen values.
// arma::eigs_sym() works only for sparse matrices which are 
// not standard R matrices.
// [[Rcpp::export]]
Rcpp::List calc_eigensp(const arma::sp_mat& matrixv, const arma::uword& neigen) {
  
  arma::mat eigenvec;
  arma::vec eigenval;
  arma::eigs_sym(eigenval, eigenvec, matrixv, neigen);

  // Reverse the order of elements from largest eigenvalue to smallest, similar to R
  return Rcpp::List::create(Rcpp::Named("values") = arma::flipud(eigenval),
                            Rcpp::Named("vectors") = arma::fliplr(eigenvec));
  
}  // end calc_eigensp


// [[Rcpp::export]]
void push_covar(const arma::rowvec& newdata, // Row of new asset returns
                arma::mat& covmat,  // Covariance matrix
                arma::rowvec& meanv, // Trailing means of the returns
                const double& lambda) { // Decay factor to multiply the past values
  
  double lambda1 = 1-lambda;
  
  // Update the means of the returns
  meanv = lambda*meanv + lambda1*newdata;
  // Calculate the de-meaned returns
  arma::rowvec datav = (newdata - meanv);
  
  // Update the covariance of the returns
  covmat = lambda*covmat + lambda1*arma::trans(datav)*datav;
  
}  // end push_covar


// Slower version below - because push_covar() calls arma::trans()
// [[Rcpp::export]]
arma::mat run_covar(const arma::mat& tseries, double lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::uword ncols = tseries.n_cols;
  arma::mat vars(nrows, ncols);
  arma::mat covar(nrows, 1);
  arma::rowvec meanv = tseries.row(0);
  arma::mat covmat = arma::trans(meanv)*meanv;
  
  // Copy the covariance data
  vars.row(0) = arma::trans(covmat.diag());
  covar.row(0) = covmat(0, 1);
  
  // Perform loop over the rows
  for (arma::uword it = 1; it < nrows; it++) {
    // Update the covariance matrix
    push_covar(tseries.row(it), covmat, meanv, lambda);
    // Copy the covariance data
    vars.row(it) = arma::trans(covmat.diag());
    covar.row(it) = covmat(0, 1);
  }  // end for
  
  return arma::join_rows(covar, vars);
  
  // Faster code below - because it doesn't call arma::trans()
  // It doesn't call push_covar()
  // arma::uword nrows = tseries.n_rows;
  // arma::uword ncols = tseries.n_cols;
  // arma::mat meanm = tseries.row(0);
  // arma::mat meand(1, ncols);
  // arma::mat vars(nrows, ncols);
  // arma::mat covar(nrows, 1);
  // double lambda1 = 1-lambda;
  // 
  // // Perform loop over the rows
  // vars.row(0) = arma::square(tseries.row(0));
  // covar.row(0) = tseries(0, 0)*tseries(0, 1);
  // for (arma::uword it = 1; it < nrows; it++) {
  //   // Calculate the mean as the weighted sum
  //   meanm = lambda*meanm + lambda1*tseries.row(it);
  //   meand = tseries.row(it) - meanm;
  //   // Calculate the covariance as the weighted sum of products of returns
  //   vars.row(it) = lambda*vars.row(it-1) + lambda1*arma::square(meand);
  //   covar.row(it) = lambda*covar.row(it-1) + lambda1*(meand(0)*meand(1));
  // }  // end for
  // 
  // return arma::join_rows(covar, vars);
  // 
  
}  // end run_covar


// [[Rcpp::export]]
void push_eigen(const arma::rowvec& newdata, // Row of new asset returns
                arma::mat& covmat,  // Covariance matrix
                arma::vec& eigenval, // Eigen values
                arma::mat& eigenvec, // Eigen vectors
                arma::rowvec& reteigen, // Row of eigen portfolio returns
                arma::rowvec& meanv, // Trailing means of the returns
                const double& lambda) { // Decay factor to multiply the past values
  
  // Calculate the eigen portfolio returns - the products of the previous eigen vectors times the scaled returns
  reteigen = (newdata/arma::trans(arma::sqrt(covmat.diag())))*eigenvec;
  // Update the covariance matrix
  push_covar(newdata, covmat, meanv, lambda);
  // Calculate the eigen decomposition
  arma::eig_sym(eigenval, eigenvec, covmat);
  
}  // end push_eigen


// [[Rcpp::export]]
arma::mat run_eigen(const arma::mat& tseries, 
                    const double& lambda) {
  
  arma::uword nrows = tseries.n_rows;
  arma::uword ncols = tseries.n_cols;
  arma::mat eigenvec(nrows, ncols);
  arma::vec eigenval(nrows);
  arma::rowvec meanv = tseries.row(0);
  arma::mat retmat(nrows, ncols);
  arma::rowvec reteigen(nrows);
  arma::mat covmat = arma::trans(meanv)*meanv;
  
  // Perform loop over the rows
  for (arma::uword it = 0; it < nrows; it++) {
    // Update the covariance matrix
    push_eigen(tseries.row(it), covmat, eigenval, eigenvec, reteigen, meanv, lambda);
    retmat.row(it) = reteigen;
  }  // end for

  return retmat;
  
}  // end run_eigen


// [[Rcpp::export]]
void run_sgao(arma::colvec& eigenval, // eigen values
             arma::mat& eigenvec, // eigen vectors
             const arma::colvec& newdata, // row of new data
             arma::colvec& meanv, // mean of the data
             arma::colvec& volv, // volatility of the data
             const double& lambda, // decay factor to multiply the past mean and volatility
             const double& gamma) { // gain factor to multiply the past eigenelements
  
  // Update the mean and volatility of the data
  meanv = lambda*meanv + (1-lambda)*newdata;
  // volv = lambda*volv + (1-lambda)*arma::square(newdata-meanv);
  // arma::colvec datav = (newdata-meanv)/volv;
  arma::colvec datav = (newdata-meanv);
  // std::cout << "datav: " << std::endl << arma::trans(datav) << std::endl;
  
  // Calculate the phis - the products of the eigen vectors times the new data
  arma::mat vecd = arma::trans(eigenvec)*datav;
  // std::cout << "vecd: " << std::endl << arma::trans(vecd) << std::endl;
  // Calculate the updated eigen values
  arma::colvec neigenval = (1-gamma)*eigenval + gamma*arma::square(vecd);
  // std::cout << "eigenval: " << std::endl << eigenval << std::endl;
  // std::cout << "neigenval: " << std::endl << arma::trans(neigenval) << std::endl;
  
  // Calculate diagonal matrix of gain factors
  arma::uword ncols = eigenvec.n_cols;
  arma::mat dgamma = arma::mat(ncols, ncols, arma::fill::zeros);
  dgamma.diag().fill(gamma);
  
  // Perform QR decomposition of Q to get the eigen vectors W
  arma::mat W, R, Q = eigenvec;
  Q += datav*arma::trans(vecd)*dgamma;
  // std::cout << "Q: " << std::endl << Q << std::endl;
  arma::qr_econ(W, R, Q);
  
  arma::uvec sorti = arma::sort_index(neigenval, "descend");
  // std::cout << "sorti: " << std::endl << sorti << std::endl;
  neigenval = neigenval(sorti);
  eigenval = neigenval;
  W = W.cols(sorti);
  eigenvec = W;
  
}  // end run_sgao





// [[Rcpp::export]]
void sgapca_exCn(arma::colvec& eigenval, 
                 arma::mat& Qs, 
                 const arma::colvec& x,
                 arma::colvec& meanv,
                 const double& n,
                 const double& gamma) {
  
  arma::mat W, R, Q = Qs;
  double lambda = 1/(n+1);
  meanv = (1-lambda)*meanv + lambda*x;
  arma::colvec x1 = x - meanv;
  
  arma::mat ux = arma::trans(Qs)*x1;
  // std::cout << "ux: " << std::endl << ux << std::endl;
  // arma::colvec gamy = gamma*ux;
  arma::uword ncols = Qs.n_cols;
  arma::mat dgamma = arma::mat(ncols, ncols, arma::fill::zeros);
  dgamma.diag().fill(gamma);
  // arma::mat dgamma = arma::diagmat(gamma);
  Q += x1*arma::trans(ux)*dgamma;
  // std::cout << "Q: " << std::endl << Q << std::endl;
  arma::qr_econ(W, R, Q);
  arma::colvec neigenval = (1 - gamma)*eigenval + gamma*arma::square(ux);
  // std::cout << "eigenval: " << std::endl << eigenval << std::endl;
  // std::cout << "neigenval: " << std::endl << neigenval << std::endl;
  // return Rcpp::wrap(W);
  arma::uvec sorti = arma::sort_index(neigenval, "descend");
  
  neigenval = neigenval(sorti);
  eigenval = neigenval;
  W = W.cols(sorti);
  Qs = W;
  
  // return Rcpp::List::create(Rcpp::Named("values") = neigenval,
  //                           Rcpp::Named("vectors") = W);
  // Rcpp::Named("ux") = arma::sort(ux, "descend"));
  
}  // end sgapca_exCn




// [[Rcpp::export]]
Rcpp::NumericMatrix sgapca_nnC (const arma::mat& Q, const arma::colvec& x, 
                                const arma::colvec& y, const arma::colvec& gamma) {
  
  int m = Q.n_rows, n = Q.n_cols;
  arma::colvec gamy = gamma % y;
  arma::colvec b = y(0)*Q.col(0);
  arma::mat A(m,n);
  A.col(0) = Q.col(0) - gamy(0)*b;	
  for (int i=1; i<n; i++) {
    b += y(i-1)*Q.col(i-1) + y(i)*Q.col(i);
    A.col(i) = Q.col(i) - gamy(i)*b;
  }		
  A += x*arma::trans(gamy);
  return Rcpp::wrap(A);
  
}  // end calc_eigensp



// [[Rcpp::export]]
Rcpp::NumericMatrix ghapca_C (const arma::mat& Q, const arma::colvec& x, 
                              const arma::colvec& y, const arma::colvec& gamma) {
  
  int m = Q.n_rows, n = Q.n_cols;
  arma::colvec gamy = gamma % y;
  arma::colvec b = y(0)*Q.col(0);
  arma::mat A(m,n);
  A.col(0) = Q.col(0) - gamy(0)*b;	
  for (int i=1; i<n; i++) {
    b += y(i)*Q.col(i);
    A.col(i) = Q.col(i) - gamy(i)*b;
  }		
  A += x*arma::trans(gamy);
  return Rcpp::wrap(A);
  
}  // end calc_eigensp


// [[Rcpp::export]]
Rcpp::List ccipca_C (arma::colvec lambda, arma::mat U, 
                     arma::colvec x, int n, int q, double l, double tol) {
  
  int i, d = x.n_elem, k = lambda.n_elem;
  if (q != k) {
    U.resize(d,q);
    lambda.resize(q);
  }
  arma::colvec v(d);  
  double f = (1.0+l)/(1.0+n);
  double nrm;
  
  for (i=0; i<q; i++) {  
    
    nrm = arma::norm(x);
    if (nrm < tol) {
      lambda.tail(q-i) = (1-f)*lambda.tail(q-i);
      break;
    }			
    
    if (i == n) {
      lambda[i] = nrm;			
      U.col(i) = x/nrm;
      break;
    }
    
    v = ((1-f)*arma::as_scalar(lambda[i]))*U.col(i) + 
      (f*arma::dot(U.col(i),x))*x;
    nrm = arma::norm(v);
    if (nrm < tol) {
      lambda[i] = 0.0;
      break; 
    } 	
    lambda[i] = nrm; 
    U.col(i) = v/nrm;
    x -= arma::dot(U.col(i), x)*U.col(i);		
    
  }
  
  // 	typedef	std::vector<double> stdvec;
  // 	stdvec lambdav = arma::conv_to< stdvec >::from(lambda); 
  // 	return Rcpp::List::create(Rcpp::Named("values") = lambdav, 
  // 			Rcpp::Named("vectors") = U);
  return Rcpp::List::create(Rcpp::Named("values") = lambda, 
                            Rcpp::Named("vectors") = U);
  
}  // end calc_eigensp


