#include <pedxing.h>

//' @export
// [[Rcpp::export]]
Rcpp::S4 rlca(unsigned int n, arma::vec kappa, arma::vec beta, arma::vec Z,
              arma::vec t0, arma::vec I, arma::vec x0, double dtovertau,
              double tau, unsigned int maxiter, bool nonLinear, bool random) 
{
  lca * obj = new lca(kappa, beta, Z, t0, I, x0, dtovertau, tau, maxiter, 
                      nonLinear, random);
  unsigned nacc = I.n_elem;
  arma::uvec counter(n);
  arma::mat output(n, 2); 

  arma::cube act_cube(maxiter, nacc, n);
  arma::cube lea_cube(maxiter, nacc, n);
  arma::cube inh_cube(maxiter, nacc, n);
  arma::cube dx_cube(maxiter, nacc, n);
  obj->r(n, output, act_cube, lea_cube, inh_cube, dx_cube, counter);
  delete obj;
  
  Rcpp::S4 out ("LCA");
  out.slot("choice_RT")  = output;
  out.slot("activation") = act_cube;
  out.slot("inhibition") = inh_cube;
  out.slot("leakage")    = lea_cube;
  out.slot("dx")         = dx_cube;
  out.slot("counter")    = counter;
  return out;
}

