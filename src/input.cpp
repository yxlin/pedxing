#include <pedxing.h>

arma::mat rInput_internal(unsigned n, arma::vec par, arma::vec I, arma::vec x0,
                        unsigned maxiter, bool nonlinear, bool random, bool debug);
  
//' @export
// [[Rcpp::export]]
Rcpp::S4 rInput(unsigned n, arma::vec par, arma::vec I, arma::vec x0,
              unsigned maxiter, bool nonlinear, bool random, bool debug) {
  // Prepare the process ---------
  if (par.n_elem != 4) Rcpp::stop("par must has 4 elements.");
  unsigned nacc = I.n_elem;
  double Z = par[0];
  double t0 = par[1];
  double dtovertau = par[2];
  double dt = par[3];
  double sdv = std::sqrt(dtovertau);
  
  arma::uvec counter(n);
  arma::mat output(2, n); 
  arma::cube act_cube(nacc, maxiter+1, n);
  arma::cube dx_cube(nacc, maxiter+1, n);

  arma::vec dx(nacc);
  arma::mat act_mat(nacc, maxiter+1);
  arma::mat dx_mat(nacc, maxiter+1);
  arma::vec activation(nacc);
  arma::vec choice_RT(2);
  bool undone; 
  unsigned counterj;
  
  for (size_t j=0; j<n; j++) {
    choice_RT.zeros(); 
    undone = true;
    counterj = 0;
    dx.zeros();
    activation.zeros();
    activation = x0;    // make activation same size as x0
    act_mat.zeros();
    dx_mat.zeros();
    act_mat.col(0) = activation; 

    do {
      counterj++;
      for (size_t i=0; i<nacc; i++) {
         dx[i] = dtovertau * I[i];
         activation[i] += dx[i];
         if (random) { activation[i] += Rf_rnorm(0, sdv); }
        
         if (activation[i] >= Z) {
            choice_RT[0] = ( (double)counterj * dt ) - 0.5*dt + t0;
            choice_RT[1] = i+1;
            undone = false;
            if (debug) Rcpp::Rcout << "Activation state: " << activation[i] << std::endl;
         }
         if (nonlinear && activation[i] < 0) activation[i] = 0;
      }
      if(debug) Rcpp::Rcout << "counterj " << counterj << std::endl;
      
      act_mat.col(counterj) = activation;
      dx_mat.col(counterj) = dx;
   } while (counterj < maxiter && undone);
  

   act_cube.slice(j) = act_mat;
   dx_cube.slice(j)  = dx_mat;
   output.col(j)     = choice_RT;
   counter(j)        = counterj;
  }
  
  Rcpp::S4 out ("Input");
  out.slot("choice_RT")  = output;
  out.slot("activation") = act_cube;
  out.slot("dx")         = dx_cube;
  out.slot("counter")    = counter;
  return out;
}


//' @export
// [[Rcpp::export]]
arma::mat rInput_internal(unsigned n, arma::vec par, arma::vec I, arma::vec x0,
              unsigned maxiter, bool nonlinear, bool random) {
  // Prepare the process ---------
  if (par.n_elem != 4) Rcpp::stop("par must has 4 elements.");
  unsigned nacc = I.n_elem;
  double Z = par[0];
  double t0 = par[1];
  double dtovertau = par[2];
  double dt = par[3];
  double sdv = std::sqrt(dtovertau);
  
  arma::mat output(2, n); 
  arma::vec dx(nacc);
  arma::vec activation(nacc);
  arma::vec choice_RT(2);
  bool undone; 
  unsigned counterj;
  
  for (size_t j=0; j<n; j++) {
    choice_RT.zeros(); 
    undone = true;
    counterj = 0;
    dx.zeros();
    activation.zeros();
    activation = x0;    // make activation same size as x0

    do {
      counterj++;
      for (size_t i=0; i<nacc; i++) {
        dx[i] = dtovertau * I[i];
        activation[i] += dx[i];
        if (random) { activation[i] += (sdv * arma::randn()); }
        
        if (activation[i] >= Z) {
          choice_RT[0] = ( (double)counterj * dt ) - 0.5*dt + t0;
          choice_RT[1] = i;
          undone = false;
        }
        if (nonlinear && activation[i] < 0) activation[i] = 0;
      }

    } while (counterj < maxiter && undone);
    
    output.col(j) = choice_RT;
  }
  
  return output;
}

//' @export
// [[Rcpp::export]]
arma::field<arma::vec> dInput(arma::vec RT, arma::vec R,
    unsigned nsim, arma::vec par, arma::vec I, arma::vec x0,
              unsigned maxiter, bool nonlinear, bool random, bool debug) {
 unsigned nacc = I.n_elem;
 arma::mat samples = rInput_internal(nsim, par, I, x0, maxiter, nonlinear, 
                                     random);
 arma::vec sRT = samples.row(0).t();
 arma::vec sR  = samples.row(1).t();
 
 if (debug) Rcpp::Rcout<< "n sim " << sRT.n_elem << std::endl;
 
 arma::field<arma::uvec> index_data(nacc), index_sam(nacc);
 arma::field<arma::vec> RTi(nacc), sRTi(nacc), out(nacc);
 
 for(size_t i=0; i<nacc; i++)
 {
   index_data[i] = arma::find(R == i);
   index_sam[i]  = arma::find(sR == i);
   
   RTi[i]   = RT.elem(index_data[i]);
   // RTi[i].print("RT from data");
   sRTi[i]  = sRT.elem(index_sam[i]);
   
   if (debug) {
     Rcpp::Rcout << "Simulated RT for accumulator " << i << " has " << 
       sRTi[i].n_elem << " elements" << std::endl;
   }
   
   if (sRTi[i].n_elem == 0) {
     out[i] = RTi[i].n_elem * 1e-10;
   } else {
     out[i] = spdf(RTi[i], sRTi[i], nsim, 0, debug);    
   }
   
   
 }
  return out;
}



