#ifndef LCA_HPP
#define LCA_HPP

#include <RcppArmadillo.h>

class lca {
public:
  arma::vec m_kappa, m_beta, m_Z, m_t0, m_I, m_x0;
  
  double m_dtovertau;
  double m_tau;
  unsigned int m_maxiter;
  bool m_nonLinear;
  
  unsigned int m_nAcc0, m_nAcc1;
  double m_dt;    // LCA pseudo distribution.
  double m_sdv;
  
  bool m_random; 

  lca (arma::vec kappa, arma::vec beta, arma::vec Z, arma::vec t0,
       arma::vec I, arma::vec x0, double dtovertau, double tau, 
       unsigned int maxiter, bool nonLinear, bool is_random) :
    m_kappa(kappa), m_beta(beta), m_Z(Z), m_t0(t0), m_I(I), m_x0(x0), 
    m_dtovertau(dtovertau), m_tau(tau), m_maxiter(maxiter), 
    m_nonLinear(nonLinear), m_random(is_random) 
  {
      // Rcout << "LCA constructor\n"; 
      m_nAcc0 = m_I.n_elem;
      // m_nAcc1 = m_x0.n_elem;
      // if (m_nAcc0 != m_nAcc1) Rcpp::stop("The number of accumulators unequal");
      unsigned int nkappa = m_kappa.n_elem;
      unsigned int nbeta  = m_beta.n_elem;
      unsigned int nZ     = m_Z.n_elem;
      unsigned int nt0    = m_t0.n_elem;
      unsigned int nI     = m_I.n_elem;
      unsigned int nx0    = m_x0.n_elem;
      if (nkappa == 1) m_kappa = arma::repmat(m_kappa, m_nAcc0, 1);
      if (nbeta == 1)  m_beta  = arma::repmat(m_beta,  m_nAcc0, 1);
      if (nZ == 1)     m_Z     = arma::repmat(m_Z,     m_nAcc0, 1);
      if (nt0 == 1)    m_t0    = arma::repmat(m_t0,    m_nAcc0, 1);
      if (nI == 1)     m_I     = arma::repmat(m_I,     m_nAcc0, 1);
      if (nx0 == 1)    m_x0    = arma::repmat(m_x0,    m_nAcc0, 1);    

      /*
      for (size_t i=0; i<m_nAcc0; i++)
      {
        if (m_t0[i] > 1.0) 
        {
          Rcpp::Rcout << "In accumulator " << i << 
            ", you provide a non-decision time, which is larger than 1." << 
            " I assume you meant " << m_t0[i] << " miliseconds\n";
          m_t0[i] = m_t0[i]/1e3;
        }
        if (m_t0[i] < 0 || m_dt < 0) Rcpp::stop("time step and nondecision time must be greater than 0.");
      }
      */
      
      m_sdv = std::sqrt(dtovertau);
      m_dt  = dtovertau*tau;

  }
  // rlca. 
  
  ~lca() 
  { 
    // Rcout << "LCA destructor\n"; 
  };

  void r_ (unsigned int & n, arma::mat & output)
  {
    // output n x 2
    arma::vec tmp(2); 
    tmp.fill(NA_REAL); // x0 values are linked with tmp?
    // tmp.print("tmp in r");

    for (size_t i=0; i<n; i++)
    {
       rt_(tmp);
       output(i, 0) = tmp[0];
       output(i, 1) = tmp[1];
    }
  }
  
  void r (unsigned int & n, arma::mat& output, arma::cube& act, 
          arma::cube& lea, arma::cube& inh, arma::cube& dx, arma::uvec& counter)
  {
    // output n x 2
    arma::vec tmp(2); 
    tmp.fill(NA_REAL); // x0 values are linked with tmp?
    
    for (size_t i=0; i < n; i++)
    {
      counter(i) = 0;
      rt(tmp, act.slice(i), lea.slice(i), inh.slice(i), dx.slice(i), counter(i));
      output(i, 0) = tmp[0];
      output(i, 1) = tmp[1];
    }
  }
  
  void print(const std::string & x) const
  {
    Rcpp::Rcout << x << "dt " << m_dt << "; ";
    m_kappa.print("kappa");
    m_beta.print("beta");
    m_Z.print("Z");
    m_t0.print("t0");
    m_I.print("I");
    m_x0.print("x0");
  }
  
  bool ValidateParams (bool print)
  {
    using namespace Rcpp;
    bool valid = true;
    for (size_t i=0; i<m_nAcc0; i++)
    {
      if (m_Z[i] < 0) {
        valid = false;
        if (print) Rcout << "invalid parameter Z = " << m_Z[i] << 
          " in accumulator " << i << std::endl; 
      } else if (m_Z[i] <= m_x0[i]) {
        valid = false;
        if (print) Rcout << "invalid parameters Z, " << m_Z[i] << " <= x0, " << 
          m_x0[i] << " in accumulator " << i << std::endl; 
      } else if (m_Z[i] <= (m_I[i]*m_dt)) {
        valid = false;
        if (print) Rcout << "invalid parameters Z, " << m_Z[i] << " <= I*dt, " << 
          m_I[i]*m_dt << " in accumulator " << i << std::endl; 
      }
      
      if (m_t0[i] < 0)
      { 
        valid = false; 
        if (print) Rcout << "invalid parameter t0 = " << m_t0[i] << 
          " in accumulator " << i << std::endl; 
      }

    }
    // 
    return valid;
  }

private:
  bool winner;
  unsigned int counter;
  double total_inhibition;
  void rt_(arma::vec& output)
  {
    winner  = false;
    counter = 0;
    
    arma::vec activation(m_nAcc0);   
    arma::vec inhibition(m_nAcc0);
    arma::vec inhibition_received(m_nAcc0);
    arma::vec leak(m_nAcc0);
    arma::vec dx(m_nAcc0);
    inhibition.zeros(); 
    inhibition_received.zeros(); 
    leak.zeros();
    dx.zeros();
    
    for (size_t i=0; i<m_nAcc0; i++) {
      activation[i] = m_x0[i];
      inhibition[i] = 0;
      inhibition_received[i] = 0;
      leak[i] = 0;
      dx[i] = 0;
    }
    

    do {
      total_inhibition = 0;
      for(size_t i = 0; i < m_nAcc0; i++) {
        inhibition_received[i] = total_inhibition - inhibition[i];
        dx[i] = m_dtovertau * (m_I[i] - leak[i] - inhibition_received[i]);
        activation[i] += dx[i]; 
        if (m_random) { activation[i] += Rf_rnorm(0, m_sdv); }
        
        if(activation[i] >= m_Z[i]) { 
          output[0] = ((double)(counter+1) * m_dt) - 0.5*m_dt + m_t0[i];
          output[1] = i+1; 
          winner = true; 
        }
        
        if (m_nonLinear && activation[i] < 0) activation[i] = 0;
      }
      counter += 1;
      
    } while (counter < m_maxiter && winner==false);
  }
  
  void rt(arma::vec& output, arma::mat& act_mat, arma::mat& lea_mat, 
          arma::mat& inh_mat, arma::mat& dx_mat, unsigned int& counter)
  {
    winner  = false;
    // counter = 0;
    
    // initialize accumulator array to record activation value
    arma::vec activation(m_nAcc0);   
    arma::vec inhibition(m_nAcc0);
    arma::vec inhibition_received(m_nAcc0);
    arma::vec leak(m_nAcc0);
    arma::vec dx(m_nAcc0);
    inhibition.zeros(); 
    inhibition_received.zeros(); 
    leak.zeros();
    dx.zeros();
    
    for (size_t i=0; i<m_nAcc0; i++) { activation[i] = m_x0[i]; }

    act_mat.row(0) = activation.t();
    lea_mat.row(0).zeros(); 
    inh_mat.row(0).zeros(); 
    dx_mat.row(0).zeros(); 
    
    do {
      total_inhibition = 0;
      for (size_t i=0; i<m_nAcc0; i++) {
        inhibition[i] = m_dtovertau * m_beta[i]*activation[i];  
        leak[i] = m_dtovertau * m_kappa[i]*activation[i];
        total_inhibition += inhibition[i];
      }
      for(size_t i = 0; i < m_nAcc0; i++) {
        inhibition_received[i] = total_inhibition - inhibition[i];
        dx[i] = m_dtovertau * (m_I[i] - leak[i] - inhibition_received[i]);
        activation[i] += dx[i]; 
        if (m_random) { activation[i] += Rf_rnorm(0, m_sdv); }

        if(activation[i] >= m_Z[i]) { 
          output[0] = ((double)(counter+1) * m_dt) - 0.5*m_dt + m_t0[i];
          output[1] = i+1; 
          winner = true; 
        }
        
        if (m_nonLinear && activation[i] < 0) activation[i] = 0;
      }
      act_mat.row(counter) = activation.t();
      dx_mat.row(counter) = dx.t();
      inh_mat.row(counter) = inhibition_received.t();
      lea_mat.row(counter) = leak.t();

      counter += 1;
      
    } while (counter < m_maxiter && winner==false);
    
    
  }
};

// arma::vec dlca_(arma::vec RT, 
//            arma::vec R, 
//            arma::vec kappa, 
//            arma::vec beta, 
//            arma::vec Z,
//            arma::vec s, 
//            arma::vec t0,
//            arma::vec I, 
//            arma::vec x0, 
//            double dt, 
//            unsigned int maxiter, 
//            bool nonLinear, 
//            unsigned int nsim);
#endif
