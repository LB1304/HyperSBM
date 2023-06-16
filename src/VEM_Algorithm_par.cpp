#include <RcppArmadillo.h>
#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace RcppParallel;



// [[Rcpp::export]]
arma::rowvec next_combination (int n, arma::rowvec first_comb) {
  int m = first_comb.n_cols;
  int iter = m-1;
  while (iter >= 0 && first_comb[iter] == n - m + (iter+1)) {
    iter--;
  }
  if (iter < 0) {
    arma::rowvec out(1);
    out.fill(arma::datum::nan);
    return out;
  }
  first_comb[iter]++;
  for (int i = iter+1; i<m; i++) {
    first_comb[i] = first_comb[iter]+(i-iter);
  }
  return first_comb;
}

// [[Rcpp::export]]
arma::rowvec next_combination_rep (int n, arma::rowvec first_comb) {
  arma::uvec ind_n_reached = arma::find(first_comb == n);
  
  if (any(ind_n_reached == 0)) {
    arma::rowvec out(1);
    out.fill(arma::datum::nan);
    return out;
  } else if (ind_n_reached.is_empty()) {
    first_comb.tail(1) += 1;
    return first_comb;
  } else {
    unsigned int first_ind = ind_n_reached(0);
    first_comb(first_ind-1) += 1;
    double new_val = first_comb(first_ind-1);
    
    first_comb(ind_n_reached) *= 0;
    first_comb(ind_n_reached) += new_val;
    return first_comb;
  }
}

// [[Rcpp::export]]
arma::rowvec next_permutation (arma::rowvec first_perm) {
  int m = first_perm.n_cols;
  if (m == 1) {
    arma::rowvec out(1);
    out.fill(arma::datum::nan);
    return out;
  }
  
  int i = m-1;
  while (first_perm[i-1] >= first_perm[i]) {
    if (--i == 0) {
      arma::rowvec out(1);
      out.fill(arma::datum::nan);
      return out;
    }
  }
  
  int j = m-1;
  while (j > i && first_perm[j] <= first_perm[i-1]) {
    j--;
  }
  
  first_perm.swap_cols(i-1, j);
  arma::rowvec first_perm_head = first_perm.head(i);
  arma::rowvec first_perm_tail = arma::reverse(first_perm.tail(m-i));
  
  return arma::join_horiz(first_perm_head, first_perm_tail);
}

// [[Rcpp::export]]
int bin_coeff(const int n, const int m) {
  if (m == 0) {
    return 1;
  }
  
  int step1 = n - m + 1;
  int step0;
  for (int i = 1; i < m; ++i) {
    step1 = (step0 = step1) * (n - m + 1 + i) / (i + 1);
  }
  
  return step1;
}

// [[Rcpp::export]]
int find_comb (int n, arma::uvec comb) {
  int m = comb.n_rows;
  arma::vec bc_m(m);
  for (int i = 0; i < m; i++) {
    bc_m[i] = bin_coeff(n-comb[m-i-1], i+1);
  }
  int pos = bin_coeff(n, m) - arma::accu(bc_m);
  return pos;
}

// [[Rcpp::export]]
int find_comb_rep (int n, arma::uvec comb) {
  int m = comb.n_rows;
  arma::vec bc_m(m);
  for (int i = 0; i < m; i++) {
    bc_m[i] = bin_coeff(n+i-comb[m-i-1], i+1);
  }
  int pos = bin_coeff(n+m-1, m) - arma::accu(bc_m);
  return pos;
}





// [[Rcpp::export]]
Rcpp::List compute_Y (int n, int M, Rcpp::List list_edges) {
  Rcpp::List Y(M-1);
  
  for (int m = 0; m < M-1; m++) {
    if (list_edges[m] != R_NilValue) {
      arma::imat edges_m = list_edges[m];
      int num_edges_m = edges_m.n_rows;
      arma::uvec Y_m = arma::zeros<arma::uvec>(num_edges_m);
      for (int e = 0; e < num_edges_m; e++) {
        arma::uvec edge_me = arma::conv_to<arma::uvec>::from(edges_m.row(e));
        int ind = find_comb(n, edge_me);
        Y_m(e) = ind-1;
      }
      Y[m] = Y_m;
    }
  }
  
  return Y;
}


struct B_full : public Worker {
  
  // Input matrix to read from
  const RMatrix<double> all_latents_m_par;
  
  // Output matrix to write to
  RVector<double> B_m_par;
  
  // Additional arguments
  int num_all_latents_m;
  int n;
  int m;
  arma::mat tau;
  arma::uvec Y_m;
  
  
  // Initialize arguments
  B_full(const NumericMatrix all_latents_m_par, NumericVector B_m_par, int num_all_latents_m, int n, int m, arma::mat tau, arma::uvec Y_m)
    : all_latents_m_par(all_latents_m_par), B_m_par(B_m_par), num_all_latents_m(num_all_latents_m), n(n), m(m), tau(tau), Y_m(Y_m) {}
  
  // convert RVector/RMatrix into arma type for Rcpp function
  // and the follwing arma data will be shared in parallel computing
  arma::mat convert() {
    RMatrix<double> all_latents_m_aux = all_latents_m_par;
    arma::mat all_latents_m_conv(all_latents_m_aux.begin(), num_all_latents_m, m+2, false);
    return all_latents_m_conv;
  }
  
  // Function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++) {
      arma::mat all_latents_m_conv = convert();
      double num = 0, den = 0;
      arma::rowvec perm = arma::conv_to<arma::rowvec>::from(all_latents_m_conv.row(j));
      do {
        arma::uvec latent_m = arma::conv_to<arma::uvec>::from(perm);
        int count = 0;
        arma::rowvec comb = arma::linspace<arma::rowvec>(1, m+2, m+2);
        do {
          arma::uvec edge_m = arma::conv_to<arma::uvec>::from(comb);
          arma::vec tau_elem = tau.elem((latent_m - 1) * tau.n_rows + (edge_m - 1));
          double new_Tau = arma::prod(tau_elem);
          den +=new_Tau;
          if (any((Y_m - count) == 0)) {
            num +=new_Tau;
          }
          count +=1;
          
          arma::rowvec new_comb = next_combination(n, comb);
          comb = new_comb;
        } while (!comb.has_nan());
        arma::rowvec new_perm = next_permutation(perm);
        perm = new_perm;
      } while (!perm.has_nan());
      // double B_m_j = num/den;
      // B_m_par[j] = B_m_j;
      double B_m_j = 0;
      if (den != 0) {
        B_m_j = num/den;
      }
      B_m_par[j] = B_m_j;
      
    }
  }
};

struct B_aff : public Worker {
  
  // Input matrix to read from
  const RVector<double> B_input;
  // Output matrix to write to
  RMatrix<double> B_output;
  // Additional arguments
  int num_all_latents_m;
  int n;
  int m;
  arma::mat tau;
  arma::uvec Y_m;
  arma::imat all_latents_m;
  arma::uvec intra_or_inter;
  double num_intra, den_intra, num_inter, den_inter;
  
  // Initialize arguments
  B_aff(const NumericVector B_input, NumericMatrix B_output, int num_all_latents_m, int n, int m, arma::mat tau, arma::uvec Y_m, arma::imat all_latents_m, arma::uvec intra_or_inter, double num_intra, double den_intra, double num_inter, double den_inter)
    : B_input(B_input), B_output(B_output), num_all_latents_m(num_all_latents_m), n(n), m(m), tau(tau), Y_m(Y_m), all_latents_m(all_latents_m), intra_or_inter(intra_or_inter), num_intra(num_intra), den_intra(den_intra), num_inter(num_inter), den_inter(den_inter) {}
  
  // Function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++) {
      arma::rowvec perm = arma::conv_to<arma::rowvec>::from(all_latents_m.row(j));
      arma::rowvec unique_perm = arma::unique(perm);

      if (unique_perm.n_elem == 1) {         // q1 = ... = qm (intra)

        do {
          arma::uvec latent_m = arma::conv_to<arma::uvec>::from(perm);
          int count = 0;
          arma::rowvec comb = arma::linspace<arma::rowvec>(1, m+2, m+2);
          do {
            arma::uvec edge_m = arma::conv_to<arma::uvec>::from(comb);
            arma::vec tau_elem = tau.elem((latent_m - 1) * tau.n_rows + (edge_m - 1));
            double new_Tau = arma::prod(tau_elem);
            den_intra +=new_Tau;
            if (any((Y_m - count) == 0)) {
              num_intra +=new_Tau;
            }
            count +=1;
            arma::rowvec new_comb = next_combination(n, comb);
            comb = new_comb;
          } while (!comb.has_nan());
          arma::rowvec new_perm = next_permutation(perm);
          perm = new_perm;
        } while (!perm.has_nan());

      } else {            // there exist at least q_j different
        
        intra_or_inter[j] = 1;
        do {
          arma::uvec latent_m = arma::conv_to<arma::uvec>::from(perm);
          int count = 0;
          arma::rowvec comb = arma::linspace<arma::rowvec>(1, m+2, m+2);
          do {
            arma::uvec edge_m = arma::conv_to<arma::uvec>::from(comb);
            arma::vec tau_elem = tau.elem((latent_m - 1) * tau.n_rows + (edge_m - 1));
            double new_Tau = arma::prod(tau_elem);
            den_inter +=new_Tau;
            if (any((Y_m - count) == 0)) {
              num_inter +=new_Tau;
            }
            count +=1;
            arma::rowvec new_comb = next_combination(n, comb);
            comb = new_comb;
          } while (!comb.has_nan());
          arma::rowvec new_perm = next_permutation(perm);
          perm = new_perm;
        } while (!perm.has_nan());

      }
    }
    
    for (int xx = 0; xx < num_all_latents_m; xx++) {
      if (intra_or_inter[xx] == 0) {
        B_output(0, xx) = num_intra;
        B_output(1, xx) = den_intra;
      } else {
        B_output(0, xx) = num_inter;
        B_output(1, xx) = den_inter;
      }
    }
    
  }
};

struct B_aff_m : public Worker {
  
  // Input matrix to read from
  const RMatrix<double> all_latents_m_par;
  // Output matrix to write to
  RVector<double> B_m_par;
  // Additional arguments
  int num_all_latents_m;
  int n;
  int m;
  arma::mat tau;
  arma::uvec Y_m;
  arma::uvec intra_or_inter;
  double num_intra, den_intra, num_inter, den_inter;
  
  // Initialize arguments
  B_aff_m(const NumericMatrix all_latents_m_par, NumericVector B_m_par, int num_all_latents_m, int n, int m, arma::mat tau, arma::uvec Y_m, arma::uvec intra_or_inter, double num_intra, double den_intra, double num_inter, double den_inter)
    : all_latents_m_par(all_latents_m_par), B_m_par(B_m_par), num_all_latents_m(num_all_latents_m), n(n), m(m), tau(tau), Y_m(Y_m), intra_or_inter(intra_or_inter), num_intra(num_intra), den_intra(den_intra), num_inter(num_inter), den_inter(den_inter) {}
  
  // Convert RVector/RMatrix into arma type
  arma::mat convert() {
    RMatrix<double> all_latents_m_aux = all_latents_m_par;
    arma::mat all_latents_m_conv(all_latents_m_aux.begin(), num_all_latents_m, m+2, false);
    return all_latents_m_conv;
  }
  
  // Function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++) {
      arma::mat all_latents_m_conv = convert();
      arma::rowvec perm = arma::conv_to<arma::rowvec>::from(all_latents_m_conv.row(j));
      arma::rowvec unique_perm = arma::unique(perm);
      
      if (unique_perm.n_elem == 1) {         // q1 = ... = qm (intra)
        
        do {
          arma::uvec latent_m = arma::conv_to<arma::uvec>::from(perm);
          int count = 0;
          arma::rowvec comb = arma::linspace<arma::rowvec>(1, m+2, m+2);
          do {
            arma::uvec edge_m = arma::conv_to<arma::uvec>::from(comb);
            arma::vec tau_elem = tau.elem((latent_m - 1) * tau.n_rows + (edge_m - 1));
            double new_Tau = arma::prod(tau_elem);
            den_intra +=new_Tau;
            if (any((Y_m - count) == 0)) {
              num_intra +=new_Tau;
            }
            count +=1;
            arma::rowvec new_comb = next_combination(n, comb);
            comb = new_comb;
          } while (!comb.has_nan());
          arma::rowvec new_perm = next_permutation(perm);
          perm = new_perm;
        } while (!perm.has_nan());
        
      } else {            // there exist at least q_j different
        
        intra_or_inter[j] = 1;
        do {
          arma::uvec latent_m = arma::conv_to<arma::uvec>::from(perm);
          int count = 0;
          arma::rowvec comb = arma::linspace<arma::rowvec>(1, m+2, m+2);
          do {
            arma::uvec edge_m = arma::conv_to<arma::uvec>::from(comb);
            arma::vec tau_elem = tau.elem((latent_m - 1) * tau.n_rows + (edge_m - 1));
            double new_Tau = arma::prod(tau_elem);
            den_inter +=new_Tau;
            if (any((Y_m - count) == 0)) {
              num_inter +=new_Tau;
            }
            count +=1;
            arma::rowvec new_comb = next_combination(n, comb);
            comb = new_comb;
          } while (!comb.has_nan());
          arma::rowvec new_perm = next_permutation(perm);
          perm = new_perm;
        } while (!perm.has_nan());
        
      }
    }
    
    double B_m_intra = num_intra/den_intra;
    double B_m_inter = num_inter/den_inter;
    for (int xx = 0; xx < num_all_latents_m; xx++) {
      if (intra_or_inter[xx] == 0) {
        B_m_par[xx] = B_m_intra;
      } else {
        B_m_par[xx] = B_m_inter;
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::List compute_B (int n, int M, arma::mat tau, Rcpp::List Y, Rcpp::List all_latents, int model, int n_threads) {
  Rcpp::List B(M-1);
  
  // FULL MODEL
  if (model == 0) {
    for (int m = 0; m < M-1; m++) {
      if (Y[m] != R_NilValue) {
        arma::uvec Y_m = Y[m];
        arma::imat all_latents_m = all_latents[m];
        int num_all_latents_m = all_latents_m.n_rows;
        
        // Converto da arma::mat a Rcpp::NumericMatrix
        NumericMatrix transformed_input = as<NumericMatrix>(wrap(all_latents_m));
        // Allocate the matrix we will return
        NumericVector output(num_all_latents_m);
        // Create the worker
        B_full B_full(transformed_input, output, num_all_latents_m, n, m, tau, Y_m);
        // Call the worker with parallelFor
        parallelFor(0, num_all_latents_m, B_full, 1, n_threads);
        // Converto da Rcpp::NumericMatrix a arma::mat
        arma::rowvec B_m = arma::rowvec(output.begin(), output.length(), false);
        
        B[m] = B_m;
      }
    }
  } 
  
  // AFFILIATION MODEL
  else if (model == 1) {
    arma::mat alpha_mat(2, M-1);
    arma::mat beta_mat(2, M-1);
    for (int m = 0; m < M-1; m++) {
      if (Y[m] != R_NilValue) {
        arma::uvec Y_m = Y[m];
        arma::imat all_latents_m = all_latents[m];
        int num_all_latents_m = all_latents_m.n_rows;
        arma::uvec intra_or_inter = arma::zeros<arma::uvec>(num_all_latents_m);
        double num_intra = 0, den_intra = 0;
        double num_inter = 0, den_inter = 0;
        
        // Converto da arma::mat a Rcpp::NumericMatrix
        NumericVector input = {0};
        // Allocate the matrix we will return
        NumericMatrix output(2, num_all_latents_m);
        // Create the worker
        B_aff B_aff(input, output, num_all_latents_m, n, m, tau, Y_m, all_latents_m, intra_or_inter, num_intra, den_intra, num_inter, den_inter);
        // Call the worker with parallelFor
        parallelFor(0, num_all_latents_m, B_aff, 1, n_threads);
        // Converto da Rcpp::NumericMatrix a arma::mat
        arma::mat alpha_beta = arma::mat(output.begin(), output.nrow(), output.ncol(), false);
        
        B[m] = alpha_beta;
        alpha_mat.col(m) = alpha_beta.col(0);
        beta_mat.col(m) = alpha_beta.col(1);
      }
    }
    
    arma::colvec alpha_vec = sum(alpha_mat, 1);
    arma::colvec beta_vec = sum(beta_mat, 1);
    double alpha = alpha_vec(0)/alpha_vec(1);
    double beta = beta_vec(0)/beta_vec(1);
    
    for (int m = 0; m < M-1; m++) {
      if (Y[m] != R_NilValue) {
        arma::mat alpha_beta = B[m];
        arma::rowvec alpha_beta_0 = alpha_beta.row(0);
        arma::uvec ind = find(alpha_beta_0 == alpha_beta_0(0));
        arma::rowvec B_m = arma::ones<arma::rowvec>(alpha_beta.n_cols);
        B_m *= beta;
        B_m(ind) *= alpha;
        B[m] = B_m;
      }
    }
  }
  
  // M-AFFILIATION MODEL
  else if (model == 2) {
    for (int m = 0; m < M-1; m++) {
      if (Y[m] != R_NilValue) {
        arma::uvec Y_m = Y[m];
        arma::imat all_latents_m = all_latents[m];
        int num_all_latents_m = all_latents_m.n_rows;
        arma::uvec intra_or_inter = arma::zeros<arma::uvec>(num_all_latents_m);
        double num_intra = 0, den_intra = 0;
        double num_inter = 0, den_inter = 0;
        
        // Converto da arma::mat a Rcpp::NumericMatrix
        NumericMatrix transformed_input = as<NumericMatrix>(wrap(all_latents_m));
        // Allocate the matrix we will return
        NumericVector output(num_all_latents_m);
        // Create the worker
        B_aff_m B_aff_m(transformed_input, output, num_all_latents_m, n, m, tau, Y_m, intra_or_inter, num_intra, den_intra, num_inter, den_inter);
        // Call the worker with parallelFor
        parallelFor(0, num_all_latents_m, B_aff_m, 1, n_threads);
        // Converto da Rcpp::NumericMatrix a arma::mat
        arma::rowvec B_m = arma::rowvec(output.begin(), output.length(), false);
        
        B[m] = B_m;
      }
    }
  }
  
  return B;
}


// [[Rcpp::export]]
arma::rowvec compute_pi (arma::mat tau) {
  int n = tau.n_rows;
  arma::rowvec pi = sum(tau, 0)/n;
  
  return pi;
}

// [[Rcpp::export]]
double compute_ELBO (int M, arma::mat tau, arma::rowvec pi, Rcpp::List Y, Rcpp::List B) {
  int n = tau.n_rows;
  int Q = tau.n_cols;
  
  arma::mat log_tau = log(tau);
  arma::uvec ind_log_tau = arma::find(tau == 0);
  log_tau.elem(ind_log_tau) = arma::zeros<arma::vec>(ind_log_tau.n_elem);
  arma::rowvec log_pi = log(pi);
  arma::uvec ind_log_pi = arma::find(pi == 0);
  log_pi.elem(ind_log_pi) = arma::zeros<arma::vec>(ind_log_pi.n_elem);
  
  int nrow = tau.n_rows;
  int ncol = tau.n_cols;
  arma::mat log_Pi(nrow, ncol);
  for (int i = 0; i < nrow; i++) {
    log_Pi.row(i) = log_pi;
  }
  double ELBO = accu(tau % log_Pi - tau % log_tau);
  
  for (int m = 0; m < M-1; m++) {
    if (Y[m] != R_NilValue) {
      arma::uvec Y_m = Y[m];
      arma::rowvec B_m = B[m];
      
      arma::rowvec edge = arma::linspace<arma::rowvec>(1, m+2, m+2);
      do {
        arma::uvec edge_conv = arma::conv_to<arma::uvec>::from(edge);
        int pos_edge = find_comb(n, edge_conv)-1;
        int Y_edge = 0;
        if (any(Y_m - pos_edge == 0)) {
          Y_edge = 1;
        }
        
        arma::rowvec latent = arma::ones<arma::rowvec>(m+2);
        do {
          arma::uvec latent_conv = arma::conv_to<arma::uvec>::from(latent);
          int pos_latent = find_comb_rep(Q, latent_conv)-1;
          double B_latent = B_m(pos_latent);
          double log_B_latent1 = 0.0, log_B_latent2 = 0.0;
          if (B_latent != 0 && B_latent != 1) {
            log_B_latent1 = log(B_latent);
            log_B_latent2 = log(1-B_latent);
          }
          
          double YlogB = Y_edge*log_B_latent1 + (1-Y_edge)*log_B_latent2;
          
          double tau_prod = 0;
          arma::rowvec latent_perm = latent;
          do {
            arma::uvec latent_perm_as_ind = arma::conv_to<arma::uvec>::from(latent_perm);
            arma::uvec edge_as_ind = arma::conv_to<arma::uvec>::from(edge);
            arma::vec tau_elems = tau.elem((latent_perm_as_ind - 1) * tau.n_rows + (edge_as_ind - 1));
            double new_Tau = arma::prod(tau_elems);
            tau_prod += new_Tau;
            
            arma::rowvec new_latent_perm = next_permutation(latent_perm);
            latent_perm = new_latent_perm;
          } while (!latent_perm.has_nan());
          
          double aux = tau_prod * YlogB;
          ELBO += aux;
          
          arma::rowvec new_latent = next_combination_rep(Q, latent);
          latent = new_latent;
        } while (!latent.has_nan());
        
        arma::rowvec new_edge = next_combination(n, edge);
        edge = new_edge;
      } while (!edge.has_nan());
      
    }
  }
  
  return ELBO;
}



struct compute_tau_parallel : public Worker {
  
  // Input matrix to read from
  const RMatrix<double> tau_old;
  
  // Output matrix to write to
  RMatrix<double> tau_new;
  
  // Additional arguments
  int n;
  int m;
  int Q;
  arma::rowvec pi;
  arma::rowvec B_m;
  arma::uvec Y_m;
  
  
  // Initialize arguments
  compute_tau_parallel(const NumericMatrix tau_old, NumericMatrix tau_new, int n, int m, int Q, arma::rowvec pi, arma::rowvec B_m, arma::uvec Y_m)
    : tau_old(tau_old), tau_new(tau_new), n(n), m(m), Q(Q), pi(pi), B_m(B_m), Y_m(Y_m) {}
  
  // convert RVector/RMatrix into arma type for Rcpp function
  // and the follwing arma data will be shared in parallel computing
  arma::mat convert() {
    RMatrix<double> tau_old_aux = tau_old;
    arma::mat tau_old_conv(tau_old_aux.begin(), n, Q, false);
    return tau_old_conv;
  }
  
  // Function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (int q = 0; q < Q; q++) {
        arma::mat tau_old_conv = convert();
        arma::rowvec edge = arma::linspace<arma::rowvec>(1, m+1, m+1);
        do {
          arma::rowvec new_edge = next_combination(n-1, edge);
          edge.elem(find(edge > i)) += 1;
          arma::rowvec augm_e(1);
          augm_e.fill(i+1);
          arma::uvec edge_augm = arma::conv_to<arma::uvec>::from(arma::sort(arma::join_rows(edge, augm_e)));
          int pos_edge_augm = find_comb(n, edge_augm)-1;
          int Y_edge_augm = 0;
          if (any(Y_m - pos_edge_augm == 0)) {
            Y_edge_augm = 1;
          }
          
          arma::rowvec latent = arma::ones<arma::rowvec>(m+1);
          do {
            arma::rowvec augm_l(1);
            augm_l.fill(q+1);
            arma::uvec latent_augm = arma::conv_to<arma::uvec>::from(arma::sort(arma::join_rows(latent, augm_l)));
            int pos_latent_augm = find_comb_rep(Q, latent_augm)-1;
            double B_latent_augm = B_m(pos_latent_augm);
            
            double log_BY_1 = 0.0, log_BY_2 = 0.0;   //MODIFICA
            if (B_latent_augm != 0 && B_latent_augm != 1) {   //MODIFICA
              log_BY_1 = log(B_latent_augm);   //MODIFICA
              log_BY_2 = log(1-B_latent_augm);   //MODIFICA
            }   //MODIFICA
            double BY = Y_edge_augm * log_BY_1 + (1 - Y_edge_augm) * log_BY_2;   //MODIFICA
            
            double tau_prod = 0;
            arma::rowvec latent_perm = latent;
            do {
              arma::uvec latent_perm_as_ind = arma::conv_to<arma::uvec>::from(latent_perm);
              arma::uvec edge_as_ind = arma::conv_to<arma::uvec>::from(edge);
              arma::vec tau_elems = tau_old_conv.elem((latent_perm_as_ind - 1) * tau_old_conv.n_rows + (edge_as_ind - 1));
              double new_Tau = arma::prod(tau_elems);
              tau_prod += new_Tau;
              
              arma::rowvec new_latent_perm = next_permutation(latent_perm);
              latent_perm = new_latent_perm;
            } while (!latent_perm.has_nan());
            
            tau_new(i, q) += (BY * tau_prod);   //MODIFICA
            
            arma::rowvec new_latent = next_combination_rep(Q, latent);
            latent = new_latent;
          } while (!latent.has_nan());
          
          edge = new_edge;
        } while (!edge.has_nan());
      }
    }
  }
};

// [[Rcpp::export]]
arma::mat compute_tau_par (int M, arma::mat tau_prec, arma::rowvec pi, Rcpp::List Y, Rcpp::List B, int n_threads) {
  int n = tau_prec.n_rows;
  int Q = tau_prec.n_cols;
  
  arma::mat tau(n, Q);   //MODIFICA
  for (int i = 0; i < n; i++) {   //MODIFICA
    for (int q = 0; q < Q; q++) {   //MODIFICA
      tau(i, q) = log(pi(q));  //MODIFICA
    }   //MODIFICA
  }   //MODIFICA
  
  for (int m = 0; m < M-1; m++) {
    if (Y[m] != R_NilValue) {
      arma::uvec Y_m = Y[m];
      arma::rowvec B_m = B[m];
      
      // Converto da arma::mat a Rcpp::NumericMatrix
      NumericMatrix transformed_input = as<NumericMatrix>(wrap(tau_prec));
      
      // Allocate the matrix we will return
      NumericMatrix output(n, Q);
      
      // Create the worker
      compute_tau_parallel compute_tau_parallel(transformed_input, output, n, m, Q, pi, B_m, Y_m);
      
      // Call the worker with parallelFor
      parallelFor(0, n, compute_tau_parallel, 1, n_threads);
      
      // Converto da Rcpp::NumericMatrix a arma::mat
      arma::mat tau_aux = arma::mat(output.begin(), output.nrow(), output.ncol(), false);
      tau += tau_aux;
    }
  }
  
  for (int i = 0; i < n; i++) {   //MODIFICA
    double tau_i_max = arma::max(tau.row(i));
    for (int q = 0; q < Q; q++) {   //MODIFICA
      tau(i, q) = exp(tau(i, q) - tau_i_max);    //MODIFICA
    }   //MODIFICA
    double rsum = arma::sum(tau.row(i));   //MODIFICA
    arma::mat tau_sum = arma::ones<arma::mat>(n, Q);   //MODIFICA
    tau_sum.row(i).fill(1/rsum);   //MODIFICA
    tau %= tau_sum;   //MODIFICA
  }   //MODIFICA
  
  return tau;
}


// [[Rcpp::export]]
Rcpp::List VE_step_par (int M, arma::mat tau_old, arma::rowvec pi, Rcpp::List Y, Rcpp::List B, double tol, int maxit, int n_threads) {
  int n = tau_old.n_rows;
  int Q = tau_old.n_cols;
  arma::mat tau(n, Q);
  
  bool stop = false;
  int it = 0;
  while (!stop) {
    it += 1;
    tau = compute_tau_par(M, tau_old, pi, Y, B, n_threads);
    arma::mat diff_taus = abs(tau - tau_old);
    bool converged = (diff_taus.max() < tol);
    bool maxit_reached = (it > maxit-1);
    stop = converged + maxit_reached;
    tau_old = tau;
  }
  
  return Rcpp::List::create(Rcpp::Named("tau") = tau,
                            Rcpp::Named("it") = it);
}

// [[Rcpp::export]]
Rcpp::List VEM_par (int M, int Q, int n, arma::mat tau_old, Rcpp::List all_latents, Rcpp::List Y, double tol, int maxit_VEM, int maxit_FP, int model, int start, bool print, int n_threads) {
  Rcpp::List B(M-1);
  arma::rowvec pi(Q);
  arma::mat tau(n, Q);
  double J;
  arma::vec J_vec(maxit_VEM);
  
  if (print) {
    Rcout << "|------------|-------------|-------------|-------------|-------------|-------------|-------------|\n";
    Rcout << "|    model   |      Q      |    start    |     step    |      J      | discrepancy |   FP_steps  |\n";
    Rcout << "|------------|-------------|-------------|-------------|-------------|-------------|-------------|\n";
  }
  
  bool stop = false;
  double J_old = R_NegInf;
  int it = 0;
  while (!stop) {
    it +=1;
    
    // 1. M-Step
    B = compute_B(n, M, tau_old, Y, all_latents, model, n_threads);
    pi = compute_pi(tau_old);
    
    // 2. VE-Step
    Rcpp::List VE = VE_step_par(M, tau_old, pi, Y, B, tol, maxit_FP, n_threads);
    arma::mat tau_aux = VE["tau"];
    tau = tau_aux;
    int it_FP = VE["it"];
    
    // 3. Compute ELBO
    J = compute_ELBO(M, tau, pi, Y, B);
    J_vec(it-1) = J;
    
    // 4. Print
    if (print && it == 1) {
      Rprintf("|%12d|%13d|%13d|%13d|%13f|             |%13d|\n", model, Q, start, it, J, it_FP);
    }
    if (print && it > 1) {
      Rprintf("|%12d|%13d|%13d|%13d|%13f|%13f|%13d|\n", model, Q, start, it, J, abs(J - J_old)/abs(J_old), it_FP);
    }
    
    // 5. Update convergence conditions
    bool converged = (abs(J - J_old)/abs(J_old) < tol);
    bool maxit_reached = (it > maxit_VEM-1);
    bool fp_converged = (it_FP == 1);
    stop = maxit_reached + (converged * fp_converged);
    
    // 6. Update parameters
    tau_old = tau;
    J_old = J;
  }
  
  if (print) {
    Rcout << "|------------|-------------|-------------|-------------|-------------|-------------|-------------|\n" ;
  }
  
  return Rcpp::List::create(Rcpp::Named("Q") = Q,
                            Rcpp::Named("B") = B,
                            Rcpp::Named("pi") = pi,
                            Rcpp::Named("tau") = tau,
                            Rcpp::Named("J") = J,
                            Rcpp::Named("J_vec") = J_vec.head(it-1),
                            Rcpp::Named("it") = it);
}

// [[Rcpp::export]]
double compute_LogLik (int n, int M, int Q, arma::mat tau, arma::rowvec pi, Rcpp::List Y, Rcpp::List B) {
  // Creo z come modifica "hard" di tau (0,1)
  arma::mat z(n, Q);
  for (int i = 0; i < n; i++) {
    arma::rowvec tau_i = tau.row(i);
    arma::uword z_i = tau_i.index_max();
    z(i, z_i) = 1;
  }
  
  arma::rowvec z_row = arma::sum(z, 0);
  arma::rowvec log_pi = log(pi);
  arma::uvec ind_log_pi = arma::find(pi == 0);
  log_pi.elem(ind_log_pi) = arma::zeros<arma::vec>(ind_log_pi.n_elem);
  double LogLik = accu(z_row % log_pi);
  
  for (int m = 0; m < M-1; m++) {
    if (Y[m] != R_NilValue) {
      arma::uvec Y_m = Y[m];
      arma::rowvec B_m = B[m];
      
      arma::rowvec edge = arma::linspace<arma::rowvec>(1, m+2, m+2);
      do {
        arma::uvec edge_conv = arma::conv_to<arma::uvec>::from(edge);
        int pos_edge = find_comb(n, edge_conv)-1;
        int Y_edge = 0;
        if (any(Y_m - pos_edge == 0)) {
          Y_edge = 1;
        }
        
        arma::rowvec latent = arma::ones<arma::rowvec>(m+2);
        do {
          arma::uvec latent_conv = arma::conv_to<arma::uvec>::from(latent);
          int pos_latent = find_comb_rep(Q, latent_conv)-1;
          double B_latent = B_m(pos_latent);
          double log_B_latent1 = 0.0, log_B_latent2 = 0.0;
          if (B_latent != 0 && B_latent != 1) {
            log_B_latent1 = log(B_latent);
            log_B_latent2 = log(1-B_latent);
          }
          
          double YlogB = Y_edge*log_B_latent1 + (1-Y_edge)*log_B_latent2;
          
          double z_prod = 0;
          arma::rowvec latent_perm = latent;
          do {
            arma::uvec latent_perm_as_ind = arma::conv_to<arma::uvec>::from(latent_perm);
            arma::uvec edge_as_ind = arma::conv_to<arma::uvec>::from(edge);
            arma::vec z_elems = z.elem((latent_perm_as_ind - 1) * z.n_rows + (edge_as_ind - 1));
            double new_Z = arma::prod(z_elems);
            z_prod += new_Z;
            
            arma::rowvec new_latent_perm = next_permutation(latent_perm);
            latent_perm = new_latent_perm;
          } while (!latent_perm.has_nan());
          
          double aux = z_prod * YlogB;
          LogLik += aux;
          
          arma::rowvec new_latent = next_combination_rep(Q, latent);
          latent = new_latent;
        } while (!latent.has_nan());
        
        arma::rowvec new_edge = next_combination(n, edge);
        edge = new_edge;
      } while (!edge.has_nan());
      
    }
  }
  
  return LogLik;
}









