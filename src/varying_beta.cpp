#include <pedxing.h>
arma::mat rvB_internal(unsigned n, arma::vec par, arma::vec I, arma::vec x0,
                        unsigned maxiter, bool nonlinear, bool random, bool debug);
  
//' @export
// [[Rcpp::export]]
Rcpp::S4 rvB(unsigned n, arma::vec par, arma::vec I, arma::vec x0,
               unsigned maxiter, bool nonlinear, bool random, bool debug) {
  
  // Prepare the process ---------
  if (par.n_elem != 5) Rcpp::stop("par must has 5 elements.");
  unsigned nacc = I.n_elem;
  double Z = par[0];
  double t0 = par[1];
  double dtovertau = par[2];
  double dt = par[3];
  
  double K = par[4];
  double B = par[5];
  
  double sdv = std::sqrt(dtovertau);
  
  arma::uvec counter(n);
  arma::mat output(2, n); 
  arma::cube act_cube(nacc, maxiter+1, n);
  arma::cube dx_cube(nacc, maxiter+1, n);
  arma::cube leak_cube(nacc, maxiter+1, n);
  arma::cube inhibition_cube(nacc, maxiter+1, n);
  arma::cube inhibition_rec_cube(nacc, maxiter+1, n);
  
  arma::mat act_mat(nacc, maxiter+1);
  arma::mat dx_mat(nacc, maxiter+1);
  arma::mat leak_mat(nacc, maxiter+1);
  arma::mat inhibition_mat(nacc, maxiter+1);
  arma::mat inhibition_rec_mat(nacc, maxiter+1);
  
  arma::vec dx(nacc);
  arma::vec activation(nacc);
  arma::vec leak(nacc);
  arma::vec inhibition(nacc);
  arma::vec inhibition_rec(nacc);
  
  arma::vec choice_RT(2);
  bool undone; 
  unsigned counterj;
  double total_inhibition;
  for (size_t j = 0; j < n; j++) {
    choice_RT.zeros(); 
    undone = true;
    counterj = 0;
    dx.zeros();
    activation.zeros();
    inhibition.zeros();
    inhibition_rec.zeros();
    
    activation = x0;    // make activation same size as x0
    act_mat.zeros();
    dx_mat.zeros();
    leak_mat.zeros();
    inhibition_mat.zeros();
    inhibition_rec_mat.zeros();
    
    act_mat.col(0) = activation; 

    if(debug) Rcpp::Rcout << "counterj ";
    do {
      counterj++;
      leak = dtovertau*K*activation;
      inhibition = dtovertau*B*activation;  
      total_inhibition = arma::accu(inhibition);
        
      for (size_t i=0; i<nacc; i++) {
         inhibition_rec[i] = (total_inhibition - inhibition[i]);
         dx[i] = dtovertau * I[i] - leak[i] - inhibition_rec[i];
         activation[i] += dx[i];
         
         if (random) { activation[i] += (sdv * arma::randn()); }
         
         
         if (activation[i] >= Z) {
            choice_RT[0] = ( (double)counterj * dt ) + t0;
            choice_RT[1] = i+1;
            undone = false;
            if (debug) Rcpp::Rcout << "Activation state: " << activation[i] << std::endl;
         }
         if (nonlinear && activation[i] < 0) activation[i] = 0;
      }
      if(debug) Rcpp::Rcout << " " << counterj << " ";
      
      act_mat.col(counterj) = activation;
      dx_mat.col(counterj) = dx;
      leak_mat.col(counterj) = leak;
      inhibition_mat.col(counterj) = inhibition;
      inhibition_rec_mat.col(counterj) = inhibition_rec;
      
   } while (counterj < maxiter && undone);
  
   act_cube.slice(j) = act_mat;
   dx_cube.slice(j)  = dx_mat;
   leak_cube.slice(j)= leak_mat;
   inhibition_cube.slice(j)= inhibition_mat;
   inhibition_rec_cube.slice(j)= inhibition_rec_mat;
   output.col(j)     = choice_RT;
   counter(j)        = counterj;
  }
  
  Rcpp::S4 out ("Inhibition");
  out.slot("choice_RT")  = output;
  out.slot("activation") = act_cube;
  out.slot("dx")         = dx_cube;
  out.slot("leak")       = leak_cube;
  out.slot("inhibition")     = inhibition_cube;
  out.slot("inhibition_rec") = inhibition_rec_cube;
  out.slot("counter")    = counter;
  return out;
}

//' @export
// [[Rcpp::export]]
arma::mat rvB_internal(unsigned n, arma::vec par, arma::vec I, arma::vec x0,
              unsigned maxiter, bool nonlinear, bool random) {
  if (par.n_elem != 6) Rcpp::stop("par must has 6 elements.");
  unsigned nacc = I.n_elem;
  double Z = par[0];
  double t0 = par[1];
  double dtovertau = par[2];
  double dt = par[3];
  
  double K = par[4];
  double B = par[5];
  
  double sdv = std::sqrt(dtovertau);
  
  arma::mat output(2, n); 

  arma::vec dx(nacc);
  arma::vec activation(nacc);
  arma::vec leak(nacc);
  arma::vec inhibition(nacc);
  arma::vec inhibition_rec(nacc);
  
  arma::vec choice_RT(2);
  bool undone; 
  unsigned counterj;
  double total_inhibition;
  for (size_t j = 0; j < n; j++) {
    choice_RT.zeros(); 
    undone = true;
    counterj = 0;
    dx.zeros();
    activation.zeros();
    inhibition.zeros();
    inhibition_rec.zeros();
    
    activation = x0;    // make activation same size as x0

    do {
      counterj++;
      leak = dtovertau*K*activation;
      inhibition = dtovertau*B*activation;  
      total_inhibition = arma::accu(inhibition);
      
      for (size_t i=0; i<nacc; i++) {
        inhibition_rec[i] = (total_inhibition - inhibition[i]);
        dx[i] = dtovertau * I[i] - leak[i] - inhibition_rec[i] ;
        activation[i] += dx[i];
        
        if (random) { activation[i] += (sdv * arma::randn()); }
        
        if (activation[i] >= Z) {
          choice_RT[0] = ( (double)counterj * dt )  + t0;
          choice_RT[1] = i;
          undone = false;
        }
        if (nonlinear && activation[i] < 0) activation[i] = 0;
      }

    } while ((counterj < maxiter) && undone);
    
    output.col(j) = choice_RT;
  }

  return output;
}

//' @export
// [[Rcpp::export]]
arma::field<arma::vec> dvB(arma::vec RT, arma::vec R,
    unsigned nsim, arma::vec par, arma::vec I, arma::vec x0,
    unsigned maxiter, bool nonlinear, bool random, bool debug) {
  
 unsigned nacc = I.n_elem;
 arma::mat samples = rvB_internal(nsim, par, I, x0, maxiter, nonlinear, random);

 arma::vec sRT = samples.row(0).t();
 arma::vec sR  = samples.row(1).t();
 
 if (debug) Rcpp::Rcout<< "nsim = " << sRT.n_elem << std::endl;
 
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
     // if (debug) {
     //   Rcpp::Rcout << "0 element at accumulator " << i << std::endl;
     // }
     
   } else {
     out[i] = spdf(RTi[i], sRTi[i], nsim, 0, debug);    
   }
   
   
 }
  return out;
}



