#ifndef PDA_HPP
#define PDA_HPP

#include <RcppArmadillo.h>

arma::vec spdf(arma::vec x, arma::vec RT, int n, double h_in, bool debug);

#endif