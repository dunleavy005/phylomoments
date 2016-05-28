
#include <RcppArmadilloExtensions/sample.h>


using namespace Rcpp;
using namespace RcppArmadillo;
using namespace arma;




// [[Rcpp::export]]
List moments_ctmcjumps(const arma::vec& edge_lengths, const arma::mat& rate_mat,
const arma::imat& label_mat, const int& num_edges, const int& num_states) {
  
  // diagonalizing rate matrix
  cx_vec cx_rate_eigval;
  cx_mat cx_rate_eigvec;
  eig_gen(cx_rate_eigval, cx_rate_eigvec, rate_mat);
  vec rate_eigval = conv_to<vec>::from(cx_rate_eigval);
  mat rate_eigvec = conv_to<mat>::from(cx_rate_eigvec);
  
  // calculating probability matrices
  cube moment_zeroth(num_states, num_states, num_edges, fill::zeros);
  for (int i = 0; i < num_edges; i++) moment_zeroth.slice(i) = expmat(rate_mat * edge_lengths(i));
  
  // calculating 1st restricted factorial moment matrices - labeled markov jumps
  cube spectr_mat(num_states, num_states, num_states, fill::zeros);
  mat E_mat = zeros<mat>(num_states, num_states);
  mat rate_mat_L = rate_mat % label_mat;
  for (int i = 0; i < num_states; i++) {
    E_mat(i,i) = 1;
    spectr_mat.slice(i) = rate_eigvec * E_mat * inv(rate_eigvec);
    E_mat.zeros();
  }
  
  cube moment_first(num_states, num_states, num_edges, fill::zeros);
  double aux_int;
  for (int t = 0; t < num_edges; t++) {
    for (int i = 0; i < num_states; i++) {
      for (int j = 0; j < num_states; j++) {
        
        if(fabs(rate_eigval(i) - rate_eigval(j)) < pow(10,-5)) {
          
          aux_int = edge_lengths(t) * exp(rate_eigval(i) * edge_lengths(t));
          
        } else {
          
          aux_int = (exp(rate_eigval(i) * edge_lengths(t)) - exp(rate_eigval(j) * edge_lengths(t))) / (rate_eigval(i) - rate_eigval(j));
          
        }
        
        moment_first.slice(t) += spectr_mat.slice(i) * rate_mat_L * spectr_mat.slice(j) * aux_int;
      }
    }
  }
  
  // calculating 2nd restricted factorial moment matrices - labeled markov jumps
  cube moment_second(num_states, num_states, num_edges, fill::zeros);
  for (int t = 0; t < num_edges; t++) {
    for (int i = 0; i < num_states; i++) {
      for (int j = 0; j < num_states; j++) {
        for (int k = 0; k < num_states; k++) {
          
          if(fabs(rate_eigval(i) - rate_eigval(j)) < pow(10,-5) && fabs(rate_eigval(j) - rate_eigval(k)) < pow(10,-5)) {
            
            aux_int = (pow(edge_lengths(t),2)/2) * exp(rate_eigval(i) * edge_lengths(t));
            
          } else if(fabs(rate_eigval(i) - rate_eigval(j)) < pow(10,-5) && fabs(rate_eigval(j) - rate_eigval(k)) >= pow(10,-5)) {
            
            aux_int = ((edge_lengths(t) * exp(rate_eigval(j) * edge_lengths(t))) / (rate_eigval(j) - rate_eigval(k))) -
            ((exp(rate_eigval(j) * edge_lengths(t)) - exp(rate_eigval(k) * edge_lengths(t))) / (pow(rate_eigval(j) - rate_eigval(k),2)));
            
          } else if(fabs(rate_eigval(i) - rate_eigval(j)) >= pow(10,-5) && fabs(rate_eigval(i) - rate_eigval(k)) < pow(10,-5)) {
            
            aux_int = (1 / (rate_eigval(i) - rate_eigval(j))) * ((edge_lengths(t) * exp(rate_eigval(i) * edge_lengths(t))) -
            ((exp(rate_eigval(j) * edge_lengths(t)) - exp(rate_eigval(k) * edge_lengths(t))) / (rate_eigval(j) - rate_eigval(k))));
            
          } else if(fabs(rate_eigval(i) - rate_eigval(j)) >= pow(10,-5) && fabs(rate_eigval(j) - rate_eigval(k)) < pow(10,-5)) {
            
            aux_int = (1 / (rate_eigval(i) - rate_eigval(j))) * (((exp(rate_eigval(i) * edge_lengths(t)) - exp(rate_eigval(k) * edge_lengths(t))) /
            (rate_eigval(i) - rate_eigval(k))) - (edge_lengths(t) * exp(rate_eigval(j) * edge_lengths(t))));
            
          } else {
            
            aux_int = (1 / (rate_eigval(i) - rate_eigval(j))) * (((exp(rate_eigval(i) * edge_lengths(t)) - exp(rate_eigval(k) * edge_lengths(t))) /
            (rate_eigval(i) - rate_eigval(k))) - ((exp(rate_eigval(j) * edge_lengths(t)) - exp(rate_eigval(k) * edge_lengths(t))) / (rate_eigval(j) - rate_eigval(k))));
            
          }
          
          moment_second.slice(t) += 2 * spectr_mat.slice(i) * rate_mat_L * spectr_mat.slice(j) * rate_mat_L * spectr_mat.slice(k) * aux_int;
        }
      }
    }
  }
  
  return List::create(
    _["zeroth"] = moment_zeroth,
    _["first"] = moment_first,
    _["second"] = moment_second
    );
}





// [[Rcpp::export]]
NumericMatrix ctmc_sim(const double& t, const arma::mat& rate_mat, const int& init_state) {
  
  // initializing variables
  IntegerVector states = seq_len(rate_mat.n_rows) - 1;    
  int new_state = init_state;
  double time = 0;
  IntegerVector state_vec = IntegerVector::create(init_state);
  NumericVector time_vec = NumericVector::create(time);
  
  while (time <= t) {
    
    // simulating time until moving to a new state
    double new_time = rexp(1, -rate_mat(new_state, new_state))[0];
    
    // sampling a new state
    vec state_prob = rate_mat.row(new_state).t() / -rate_mat(new_state, new_state);
    state_prob(new_state) = 0;
    new_state = sample(states, 1, true, wrap(state_prob))[0];
    
    // caching results
    time += new_time;
    time_vec.push_back(time);
    state_vec.push_back(new_state);
  }
  
  // last cached entry must be at time "t"
  time_vec(time_vec.size() - 1) = t;
  state_vec(state_vec.size() - 1) = state_vec(state_vec.size() - 2);
  
  // storing output
  NumericMatrix outp(2, time_vec.size());
  outp(0, _) = state_vec;
  outp(1, _) = time_vec;
  
  return outp;
}





List stoch_mapping_lists(const arma::imat& edge_mat, const arma::ivec& edge_set, const int& num_edges,
const int& num_states, const int& num_term_nodes, const List& edge_moments, const arma::ivec& tip_data) {
  
  // storing probability matrices & restricted factorial moment matrices - labeled markov jumps
  NumericVector tmp_vec;
  
  tmp_vec = edge_moments["zeroth"];
  cube prob_mat(tmp_vec.begin(), num_states, num_states, num_edges, false);
  tmp_vec = edge_moments["first"];
  cube moments_first(tmp_vec.begin(), num_states, num_states, num_edges, false);
  tmp_vec = edge_moments["second"];
  cube moments_second(tmp_vec.begin(), num_states, num_states, num_edges, false);
  
  // preparing for bottom-up tree traversal
  mat downward_lik(2*num_term_nodes - 1, num_states, fill::zeros);
  mat direction_lik(num_edges, num_states, fill::zeros);
  mat V_list_first(num_edges, num_states, fill::zeros);
  mat V_list_second(num_edges, num_states, fill::zeros);
  mat W_list(num_edges, num_states, fill::zeros);
  
  // initializing downward partial likelihoods - tips
  if(tip_data.n_elem == 0) {
    // are we looking for prior moments?
    downward_lik.rows(0, num_term_nodes - 1).ones();
  } else {
    // or posterior moments?
    for (int i = 0; i < num_term_nodes; i++) downward_lik(i, tip_data(i)) = 1;
  }
  
  for (int i = 0; i < num_edges; i += 2) {
    
    // calculating downward partial likelihoods - internal nodes & directional partial likelihoods - branches
    direction_lik.row(i) = trans(prob_mat.slice(i) * downward_lik.row(edge_mat(i,1) - 1).t());
    direction_lik.row(i+1) = trans(prob_mat.slice(i+1) * downward_lik.row(edge_mat(i+1,1) - 1).t());
    downward_lik.row(edge_mat(i,0) - 1) = direction_lik.row(i) % direction_lik.row(i+1);
    
    // calculating the V/W lists - branches
    
    for (int j = i; j < i+2; j++) {
      
      // is the branch internal?
      if(edge_mat(j,1) > num_term_nodes) {
        
        uvec child_ind = find(edge_mat.col(0) == edge_mat(j,1));
        
        // calculating (internal) non-branch parts of V/W lists
        V_list_first.row(j) = trans(prob_mat.slice(j) * trans((V_list_first.row(child_ind(0)) %
        direction_lik.row(child_ind(1))) + (V_list_first.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
        
        V_list_second.row(j) = trans(prob_mat.slice(j) * trans((V_list_second.row(child_ind(0)) %
        direction_lik.row(child_ind(1))) + (V_list_second.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
        
        W_list.row(j) = trans(prob_mat.slice(j) * trans(2 * (V_list_first.row(child_ind(0)) %
        V_list_first.row(child_ind(1))) + (W_list.row(child_ind(0)) % direction_lik.row(child_ind(1))) +
        (W_list.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
        
        // is the (internal) branch in the edge set?
        if(any(edge_set == j+1)) {
          
          W_list.row(j) += trans(2 * moments_first.slice(j) * trans((V_list_first.row(child_ind(0)) % direction_lik.row(child_ind(1))) +
          (V_list_first.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
          
        }
        
      }
      
      // is the branch in the edge set?
      if(any(edge_set == j+1)) {
        
        V_list_first.row(j) += trans(moments_first.slice(j) * downward_lik.row(edge_mat(j,1) - 1).t());
        
        V_list_second.row(j) += trans(moments_second.slice(j) * downward_lik.row(edge_mat(j,1) - 1).t());
        
      }
    }
  }
  
  return List::create(
    _["downward_lik"] = downward_lik,
    _["direction_lik"] = direction_lik,
    _["V_list_first"] = V_list_first,
    _["V_list_second"] = V_list_second,
    _["W_list"] = W_list
    );
}





List joint_stoch_mapping_lists(const arma::imat& edge_mat, const arma::ivec& edge_set1, const arma::ivec& edge_set2,
const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments, const arma::ivec& tip_data) {
  
  // storing probability matrices & restricted factorial moment matrices - labeled markov jumps
  NumericVector tmp_vec;
  
  tmp_vec = edge_moments["zeroth"];
  cube prob_mat(tmp_vec.begin(), num_states, num_states, num_edges, false);
  tmp_vec = edge_moments["first"];
  cube moments_first(tmp_vec.begin(), num_states, num_states, num_edges, false);
  tmp_vec = edge_moments["second"];
  cube moments_second(tmp_vec.begin(), num_states, num_states, num_edges, false);
  
  // preparing for bottom-up tree traversal
  mat downward_lik(2*num_term_nodes - 1, num_states, fill::zeros);
  mat direction_lik(num_edges, num_states, fill::zeros);
  mat V_list1_first(num_edges, num_states, fill::zeros);
  mat V_list2_first(num_edges, num_states, fill::zeros);
  mat V_list12_first(num_edges, num_states, fill::zeros);
  mat V_list1_second(num_edges, num_states, fill::zeros);
  mat V_list2_second(num_edges, num_states, fill::zeros);
  mat V_list12_second(num_edges, num_states, fill::zeros);
  mat W_list1(num_edges, num_states, fill::zeros);
  mat W_list2(num_edges, num_states, fill::zeros);
  mat W_list12(num_edges, num_states, fill::zeros);
  
  // initializing downward partial likelihoods - tips
  if(tip_data.n_elem == 0) {
    // are we looking for prior moments?
    downward_lik.rows(0, num_term_nodes - 1).ones();
  } else {
    // or posterior moments?
    for (int i = 0; i < num_term_nodes; i++) downward_lik(i, tip_data(i)) = 1;
  }
  
  for (int i = 0; i < num_edges; i += 2) {
    
    // calculating downward partial likelihoods - internal nodes & directional partial likelihoods - branches
    direction_lik.row(i) = trans(prob_mat.slice(i) * downward_lik.row(edge_mat(i,1) - 1).t());
    direction_lik.row(i+1) = trans(prob_mat.slice(i+1) * downward_lik.row(edge_mat(i+1,1) - 1).t());
    downward_lik.row(edge_mat(i,0) - 1) = direction_lik.row(i) % direction_lik.row(i+1);
    
    // calculating the V/W lists - branches
    
    for (int j = i; j < i+2; j++) {
      
      // is the branch internal?
      if(edge_mat(j,1) > num_term_nodes) {
        
        uvec child_ind = find(edge_mat.col(0) == edge_mat(j,1));
        
        // calculating (internal) non-branch parts of V/W lists
        V_list1_first.row(j) = trans(prob_mat.slice(j) * trans((V_list1_first.row(child_ind(0)) %
        direction_lik.row(child_ind(1))) + (V_list1_first.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
        
        V_list2_first.row(j) = trans(prob_mat.slice(j) * trans((V_list2_first.row(child_ind(0)) %
        direction_lik.row(child_ind(1))) + (V_list2_first.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
        
        V_list12_first.row(j) = trans(prob_mat.slice(j) * trans((V_list12_first.row(child_ind(0)) %
        direction_lik.row(child_ind(1))) + (V_list12_first.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
        
        V_list1_second.row(j) = trans(prob_mat.slice(j) * trans((V_list1_second.row(child_ind(0)) %
        direction_lik.row(child_ind(1))) + (V_list1_second.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
        
        V_list2_second.row(j) = trans(prob_mat.slice(j) * trans((V_list2_second.row(child_ind(0)) %
        direction_lik.row(child_ind(1))) + (V_list2_second.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
        
        V_list12_second.row(j) = trans(prob_mat.slice(j) * trans((V_list12_second.row(child_ind(0)) %
        direction_lik.row(child_ind(1))) + (V_list12_second.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
        
        W_list1.row(j) = trans(prob_mat.slice(j) * trans(2 * (V_list1_first.row(child_ind(0)) %
        V_list1_first.row(child_ind(1))) + (W_list1.row(child_ind(0)) % direction_lik.row(child_ind(1))) +
        (W_list1.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
        
        W_list2.row(j) = trans(prob_mat.slice(j) * trans(2 * (V_list2_first.row(child_ind(0)) %
        V_list2_first.row(child_ind(1))) + (W_list2.row(child_ind(0)) % direction_lik.row(child_ind(1))) +
        (W_list2.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
        
        W_list12.row(j) = trans(prob_mat.slice(j) * trans((V_list1_first.row(child_ind(0)) % V_list2_first.row(child_ind(1))) +
        (V_list2_first.row(child_ind(0)) % V_list1_first.row(child_ind(1))) + (W_list12.row(child_ind(0)) %
        direction_lik.row(child_ind(1))) + (W_list12.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
        
        // is the (internal) branch in the 1st edge set?
        if(any(edge_set1 == j+1)) {
          
          W_list1.row(j) += trans(2 * moments_first.slice(j) * trans((V_list1_first.row(child_ind(0)) % direction_lik.row(child_ind(1))) +
          (V_list1_first.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
          
          W_list12.row(j) += trans(moments_first.slice(j) * trans((V_list2_first.row(child_ind(0)) % direction_lik.row(child_ind(1))) +
          (V_list2_first.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
          
        }
        
        // is the (internal) branch in the 2nd edge set?
        if(any(edge_set2 == j+1)) {
          
          W_list2.row(j) += trans(2 * moments_first.slice(j) * trans((V_list2_first.row(child_ind(0)) % direction_lik.row(child_ind(1))) +
          (V_list2_first.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
          
          W_list12.row(j) += trans(moments_first.slice(j) * trans((V_list1_first.row(child_ind(0)) % direction_lik.row(child_ind(1))) +
          (V_list1_first.row(child_ind(1)) % direction_lik.row(child_ind(0)))));
          
        }
        
      }
      
      // is the branch in the 1st edge set?
      if(any(edge_set1 == j+1)) {
        
        V_list1_first.row(j) += trans(moments_first.slice(j) * downward_lik.row(edge_mat(j,1) - 1).t());
        
        V_list1_second.row(j) += trans(moments_second.slice(j) * downward_lik.row(edge_mat(j,1) - 1).t());
        
      }
      
      // is the branch in the 2nd edge set?
      if(any(edge_set2 == j+1)) {
        
        V_list2_first.row(j) += trans(moments_first.slice(j) * downward_lik.row(edge_mat(j,1) - 1).t());
        
        V_list2_second.row(j) += trans(moments_second.slice(j) * downward_lik.row(edge_mat(j,1) - 1).t());
        
      }
      
      // is the branch in both edge sets?
      if(any(edge_set1 == j+1) && any(edge_set2 == j+1)) {
        
        V_list12_first.row(j) += trans(moments_first.slice(j) * downward_lik.row(edge_mat(j,1) - 1).t());
        
        V_list12_second.row(j) += trans(moments_second.slice(j) * downward_lik.row(edge_mat(j,1) - 1).t());
        
      }
    }
  }
  
  return List::create(
    _["downward_lik"] = downward_lik,
    _["direction_lik"] = direction_lik,
    _["V_list1_first"] = V_list1_first,
    _["V_list2_first"] = V_list2_first,
    _["V_list12_first"] = V_list12_first,
    _["V_list1_second"] = V_list1_second,
    _["V_list2_second"] = V_list2_second,
    _["V_list12_second"] = V_list12_second,
    _["W_list1"] = W_list1,
    _["W_list2"] = W_list2,
    _["W_list12"] = W_list12
    );
}






mat phylo_likelihood_list(const imat& edge_mat, const int& num_edges, const int& num_states,
const int& num_term_nodes, const cube& prob_mat, const ivec& tip_data) {
  
  // preparing for bottom-up tree traversal
  mat downward_lik(2*num_term_nodes - 1, num_states, fill::zeros);
  
  // initializing downward partial likelihoods - tips
  for (int i = 0; i < num_term_nodes; i++) downward_lik(i, tip_data(i)) = 1;
  
  for (int i = 0; i < num_edges; i += 2) {
    
    // calculating downward partial likelihoods - internal nodes
    rowvec left_lik = trans(prob_mat.slice(i) * downward_lik.row(edge_mat(i,1) - 1).t());
    rowvec right_lik = trans(prob_mat.slice(i+1) * downward_lik.row(edge_mat(i+1,1) - 1).t());
    downward_lik.row(edge_mat(i,0) - 1) = left_lik % right_lik;
    
  }
  
  return downward_lik;
}






// [[Rcpp::export]]
double phylo_log_likelihood(const arma::imat& edge_mat, const arma::vec& edge_lengths, const arma::mat& rate_mat,
const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const arma::imat& seq_data) {
  
  // initializing variables
  int seq_ncols = seq_data.n_cols;
  double log_lik = 0;
  
  // calculating probability matrices
  cube prob_mat(num_states, num_states, num_edges, fill::zeros);
  for (int i = 0; i < num_edges; i++) prob_mat.slice(i) = expmat(rate_mat * edge_lengths(i));
  
  for (int i = 0; i < seq_ncols; i++) {
    
    // calculating downward partial likelihoods
    mat downward_lik_i = phylo_likelihood_list(edge_mat, num_edges, num_states, num_term_nodes, prob_mat, seq_data.col(i));
    
    log_lik += log(as_scalar(downward_lik_i.row(num_term_nodes) * root_dist));
    
  }
  
  return log_lik;
}






// [[Rcpp::export]]
arma::imat int_states_sim(const arma::imat& edge_mat, const arma::vec& edge_lengths, const arma::mat& rate_mat, const arma::vec& root_dist,
const int& num_edges, const int& num_states, const int& num_term_nodes, const arma::ivec& tip_data, const int& N) {
  
  // initializing variables
  int num_to_sample = num_term_nodes - 1;
  IntegerVector states = seq_len(num_states) - 1;
  imat node_states(2*num_term_nodes - 1, N, fill::zeros);
  node_states.rows(0, num_term_nodes - 1).each_col() = tip_data;
  
  // calculating probability matrices
  cube prob_mat(num_states, num_states, num_edges, fill::zeros);
  for (int i = 0; i < num_edges; i++) prob_mat.slice(i) = expmat(rate_mat * edge_lengths(i));
  
  // storing downward partial likelihoods
  mat downward_lik = phylo_likelihood_list(edge_mat, num_edges, num_states, num_term_nodes, prob_mat, tip_data);
  
  // sampling root node states
  vec state_prob = (downward_lik.row(num_term_nodes).t() % root_dist) / sum(downward_lik.row(num_term_nodes).t() % root_dist);
  
  IntegerVector sample_vec = sample(states, N, true, wrap(state_prob));
  node_states.row(num_term_nodes) = irowvec(sample_vec.begin(), sample_vec.size(), false);
  
  num_to_sample -= 1;
  
  // finding the next internal nodes to sample at
  uvec j;
  j << 1;
  ivec next_int_nodes = edge_mat(find(edge_mat.col(0) == (num_term_nodes + 1)), j);
  
  // making sure the vector of internal nodes actually contains internal nodes
  ivec int_nodes = next_int_nodes(find(next_int_nodes >= (num_term_nodes + 1)));
  
  // within this loop, we'll sample states for all remaining internal nodes
  while (num_to_sample > 0) {
    
    int int_nodes_len = int_nodes.n_elem;
    ivec next_int_nodes = zeros<ivec>(2*int_nodes_len);
    
    for (int i = 0; i < int_nodes_len; i++) {
      
      // finding the parent node and branch associated with the current internal node
      uvec current_ind;
      current_ind << int_nodes(i) - 1;
      uvec edge_ind = find(edge_mat.col(1) == (current_ind(0) + 1));
      int parent_ind = edge_mat(edge_ind(0), 0) - 1;
      
      for (int m = 0; m < num_states; m++) {
        
        uvec col_ind = find(node_states.row(parent_ind) == m);
        
        vec state_prob = trans(downward_lik.row(current_ind(0)) % prob_mat.slice(edge_ind(0)).row(m)) /
        sum(downward_lik.row(current_ind(0)) % prob_mat.slice(edge_ind(0)).row(m));
        
        sample_vec = sample(states, col_ind.n_elem, true, wrap(state_prob));
        node_states(current_ind, col_ind) = irowvec(sample_vec.begin(), sample_vec.size(), false);
      }
      
      num_to_sample -= 1;
      
      // finding the next internal nodes to sample at
      next_int_nodes.subvec(2*i, 2*i + 1) = edge_mat(find(edge_mat.col(0) == (current_ind(0) + 1)), j);
    }
    
    // making sure the vector of internal nodes actually contains internal nodes
    int_nodes = next_int_nodes(find(next_int_nodes >= (num_term_nodes + 1)));
  }
  
  return node_states;
}






// [[Rcpp::export]]
NumericVector prior_moments_phylojumps(const arma::imat& edge_mat, const arma::ivec& edge_set, const arma::vec& root_dist,
const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments) {
  
  // calculating/storing stochastic mapping lists for prior mean/variance calculations
  ivec fake_tips = zeros<ivec>(0);
  List phyl_list = stoch_mapping_lists(edge_mat, edge_set, num_edges, num_states, num_term_nodes, edge_moments, fake_tips);
  
  mat direction_lik = phyl_list["direction_lik"];
  mat V_list_first = phyl_list["V_list_first"];
  mat V_list_second = phyl_list["V_list_second"];
  mat W_list = phyl_list["W_list"];
  
  // calculating prior moments
  uvec root_edges_ind, root_edges_ind_rev;
  root_edges_ind << num_edges - 1 << num_edges - 2;
  root_edges_ind_rev << num_edges - 2 << num_edges - 1;
  
  double prior_mean = sum( (V_list_first.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) * root_dist );
  
  double prior_var = sum( ( ((V_list_second.rows(root_edges_ind) + V_list_first.rows(root_edges_ind)) %
  direction_lik.rows(root_edges_ind_rev)) + (W_list.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) +
  (V_list_first.rows(root_edges_ind) % V_list_first.rows(root_edges_ind_rev)) ) * root_dist ) - pow(prior_mean, 2);
  
  return NumericVector::create(
    _["mean"] = prior_mean,
    _["var"] = prior_var
    );
}






// [[Rcpp::export]]
NumericVector post_moments_phylojumps(const arma::imat& edge_mat, const arma::ivec& edge_set, const arma::vec& root_dist,
const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments, const arma::imat& seq_data) {
  
  int seq_ncols = seq_data.n_cols;
  uvec root_edges_ind, root_edges_ind_rev;
  root_edges_ind << num_edges - 1 << num_edges - 2;
  root_edges_ind_rev << num_edges - 2 << num_edges - 1;
  double post_mean = 0, post_var = 0;
  
  for (int i = 0; i < seq_ncols; i++) {
    
    // calculating/storing stochastic mapping lists for posterior mean/variance calculations  
    List phyl_list = stoch_mapping_lists(edge_mat, edge_set, num_edges, num_states, num_term_nodes, edge_moments, seq_data.col(i));
    
    mat downward_lik = phyl_list["downward_lik"];
    mat direction_lik = phyl_list["direction_lik"];
    mat V_list_first = phyl_list["V_list_first"];
    mat V_list_second = phyl_list["V_list_second"];
    mat W_list = phyl_list["W_list"];
    
    double post_mean_i = sum( (V_list_first.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) * root_dist ) /
    sum(downward_lik.row(num_term_nodes) * root_dist);
    
    post_mean += post_mean_i;
    
    double post_var_i = sum( ( ((V_list_second.rows(root_edges_ind) + V_list_first.rows(root_edges_ind)) %
    direction_lik.rows(root_edges_ind_rev)) + (W_list.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) +
    (V_list_first.rows(root_edges_ind) % V_list_first.rows(root_edges_ind_rev)) ) * root_dist ) /
    sum(downward_lik.row(num_term_nodes) * root_dist) - pow(post_mean_i, 2);
    
    post_var += post_var_i;
  }
  
  return NumericVector::create(
    _["mean"] = post_mean,
    _["var"] = post_var
    );
}






// [[Rcpp::export]]
arma::ivec find_subtree_edges(const int& node, const arma::imat& edge_mat) {
  
  // finding branch above "node"
  int branch_above = as_scalar(find(edge_mat.col(1) == node)) + 1;
  ivec outp;
  outp << branch_above;
  
  // finding "node" children (if there are any)
  uvec node_children_ind = find(edge_mat.col(0) == node);
  
  // recursion base case
  if(node_children_ind.is_empty()) return outp;
  
  // else, append "branch_above" to results obtained from recursion on the "node" children
  uvec j;
  j << 1;
  ivec node_children = edge_mat(node_children_ind, j);
  
  return join_cols(
    join_cols(find_subtree_edges(node_children(0), edge_mat),
    find_subtree_edges(node_children(1), edge_mat)), outp
    );
}






// [[Rcpp::export]]
arma::imat tips_sim(const arma::imat& edge_mat, const arma::vec& edge_lengths, const arma::mat& rate_mat,
const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const int& N) {
  
  // initializing variables
  IntegerVector states = seq_len(num_states) - 1;
  imat node_states(2*num_term_nodes - 1, N, fill::zeros);
  
  // calculating probability matrices
  cube prob_mat(num_states, num_states, num_edges, fill::zeros);
  for (int i = 0; i < num_edges; i++) prob_mat.slice(i) = expmat(rate_mat * edge_lengths(i));
  
  // sampling root node states
  IntegerVector sample_vec = sample(states, N, true, wrap(root_dist));
  node_states.row(num_term_nodes) = irowvec(sample_vec.begin(), sample_vec.size(), false);
  
  // sampling states for remaining nodes
  for (int i = (num_edges - 1); i > 0; i -= 2) {
    
    int parent_ind = edge_mat(i,0) - 1;
    uvec lchild_ind, rchild_ind;
    lchild_ind << edge_mat(i-1,1) - 1;
    rchild_ind << edge_mat(i,1) - 1;
    
    for (int m = 0; m < num_states; m++) {
      
      uvec col_ind = find(node_states.row(parent_ind) == m);
      
      sample_vec = sample(states, col_ind.n_elem, true, wrap(prob_mat.slice(i-1).row(m).t()));
      node_states(lchild_ind, col_ind) = irowvec(sample_vec.begin(), sample_vec.size(), false);
      
      sample_vec = sample(states, col_ind.n_elem, true, wrap(prob_mat.slice(i).row(m).t()));
      node_states(rchild_ind, col_ind) = irowvec(sample_vec.begin(), sample_vec.size(), false);
    }
  }
  
  return node_states.rows(0, num_term_nodes - 1);
}






// [[Rcpp::export]]
arma::imat phylojumps_sim(const arma::imat& edge_mat, const arma::vec& edge_lengths, const arma::mat& rate_mat, const arma::vec& root_dist,
const int& num_edges, const int& num_states, const int& num_term_nodes, const arma::imat& seq_data, const int& N) {
  
  // initializing variables
  IntegerVector states = seq_len(num_states) - 1;
  int seq_ncols = seq_data.n_cols;
  imat simjumps(N, seq_ncols, fill::zeros);
  
  for (int i = 0; i < seq_ncols; i++) {
    
    ivec tip_data = seq_data.col(i);
    imat node_states = int_states_sim(edge_mat, edge_lengths, rate_mat, root_dist, num_edges, num_states, num_term_nodes, tip_data, N);
    
    for (int j = 0; j < num_edges; j++) {
      
      int parent_ind = edge_mat(j,0) - 1;
      int child_ind = edge_mat(j,1) - 1;
      int k = 0;
      double edgelen = edge_lengths(j);
      ivec edge_jumps = zeros<ivec>(N);
      
      while (k < N) {
        
        NumericVector ctmc_states;
        
        if(node_states(parent_ind,k) != node_states(child_ind,k) && edge_jumps(k) == 0) {
          
          double parent_rate = -rate_mat(node_states(parent_ind,k), node_states(parent_ind,k));
          
          // sampling first waiting time conditional on at least 1 jump
          double wait_time = -log(1 - runif(1)[0] * (1 - exp(-parent_rate * edgelen))) / parent_rate;
          
          // sampling first state change
          vec state_prob = rate_mat.row(node_states(parent_ind,k)).t() / parent_rate;
          state_prob(node_states(parent_ind,k)) = 0;
          int new_state = sample(states, 1, true, wrap(state_prob))[0];
          
          // sampling rest of trajectory
          NumericMatrix ctmc_samp = ctmc_sim(edgelen - wait_time, rate_mat, new_state);
          ctmc_states = ctmc_samp(0, _);
          ctmc_states.push_front(node_states(parent_ind,k));
          
        } else {
          
          NumericMatrix ctmc_samp = ctmc_sim(edgelen, rate_mat, node_states(parent_ind,k));
          ctmc_states = ctmc_samp(0, _);
          
        }
        
        if(ctmc_states(ctmc_states.size() - 1) == node_states(child_ind,k)) {
          
          edge_jumps(k) = ctmc_states.size() - 2;
          k++;
          
        }
        
      }
      
      simjumps.col(i) += edge_jumps;
    }
  }
  
  return simjumps;
}






// [[Rcpp::export]]
NumericVector postmean_moments_phylojumps(const arma::imat& edge_mat, const arma::vec& edge_lengths, const arma::mat& rate_mat, const arma::ivec& edge_set,
const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments, const int& N) {
  
  // calculating prior moments
  NumericVector prior_moments = prior_moments_phylojumps(edge_mat, edge_set, root_dist, num_edges, num_states, num_term_nodes, edge_moments);
  
  double prior_mean = prior_moments["mean"];
  double prior_var = prior_moments["var"];
  
  // getting a monte carlo estimate of the mean posterior variance
  imat seq_sim = tips_sim(edge_mat, edge_lengths, rate_mat, root_dist, num_edges, num_states, num_term_nodes, N);
  
  double mean_post_var = post_moments_phylojumps(edge_mat, edge_set, root_dist, num_edges,
  num_states, num_term_nodes, edge_moments, seq_sim)["var"] / N;
  
  return NumericVector::create(
    _["mean"] = prior_mean,
    _["var"] = prior_var - mean_post_var
    );
}






// [[Rcpp::export]]
NumericVector joint_prior_moments_phylojumps(const arma::imat& edge_mat, const arma::ivec& edge_set1, const arma::ivec& edge_set2,
const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments) {
  
  // calculating/storing stochastic mapping lists for prior mean/variance/covariance calculations
  ivec fake_tips = zeros<ivec>(0);
  List phyl_list = joint_stoch_mapping_lists(edge_mat, edge_set1, edge_set2, num_edges, num_states, num_term_nodes, edge_moments, fake_tips);
  
  mat direction_lik = phyl_list["direction_lik"];
  mat V_list1_first = phyl_list["V_list1_first"];
  mat V_list2_first = phyl_list["V_list2_first"];
  mat V_list12_first = phyl_list["V_list12_first"];
  mat V_list1_second = phyl_list["V_list1_second"];
  mat V_list2_second = phyl_list["V_list2_second"];
  mat V_list12_second = phyl_list["V_list12_second"];
  mat W_list1 = phyl_list["W_list1"];
  mat W_list2 = phyl_list["W_list2"];
  mat W_list12 = phyl_list["W_list12"];
  
  // calculating prior moments
  uvec root_edges_ind, root_edges_ind_rev;
  root_edges_ind << num_edges - 1 << num_edges - 2;
  root_edges_ind_rev << num_edges - 2 << num_edges - 1;
  
  double prior_mean1 = sum( (V_list1_first.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) * root_dist );
  
  double prior_mean2 = sum( (V_list2_first.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) * root_dist );
  
  double prior_var1 = sum( ( ((V_list1_second.rows(root_edges_ind) + V_list1_first.rows(root_edges_ind)) %
  direction_lik.rows(root_edges_ind_rev)) + (W_list1.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) +
  (V_list1_first.rows(root_edges_ind) % V_list1_first.rows(root_edges_ind_rev)) ) * root_dist ) - pow(prior_mean1, 2);
  
  double prior_var2 = sum( ( ((V_list2_second.rows(root_edges_ind) + V_list2_first.rows(root_edges_ind)) %
  direction_lik.rows(root_edges_ind_rev)) + (W_list2.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) +
  (V_list2_first.rows(root_edges_ind) % V_list2_first.rows(root_edges_ind_rev)) ) * root_dist ) - pow(prior_mean2, 2);
  
  double prior_cov = sum( ( ((V_list12_second.rows(root_edges_ind) + V_list12_first.rows(root_edges_ind)) %
  direction_lik.rows(root_edges_ind_rev)) + (W_list12.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) +
  (V_list1_first.rows(root_edges_ind) % V_list2_first.rows(root_edges_ind_rev)) ) * root_dist ) - prior_mean1 * prior_mean2;
  
  return NumericVector::create(
    _["mean1"] = prior_mean1,
    _["mean2"] = prior_mean2,
    _["var1"] = prior_var1,
    _["var2"] = prior_var2,
    _["cov"] = prior_cov
    );
}






// [[Rcpp::export]]
NumericVector joint_post_moments_phylojumps(const arma::imat& edge_mat, const arma::ivec& edge_set1, const arma::ivec& edge_set2, const arma::vec& root_dist,
const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments, const arma::imat& seq_data) {
  
  int seq_ncols = seq_data.n_cols;
  uvec root_edges_ind, root_edges_ind_rev;
  root_edges_ind << num_edges - 1 << num_edges - 2;
  root_edges_ind_rev << num_edges - 2 << num_edges - 1;
  double post_mean1 = 0, post_mean2 = 0, post_var1 = 0, post_var2 = 0, post_cov = 0;
  
  for (int i = 0; i < seq_ncols; i++) {
    
    // calculating/storing stochastic mapping lists for posterior mean/variance/covariance calculations  
    List phyl_list = joint_stoch_mapping_lists(edge_mat, edge_set1, edge_set2, num_edges, num_states, num_term_nodes, edge_moments, seq_data.col(i));
    
    mat downward_lik = phyl_list["downward_lik"];
    mat direction_lik = phyl_list["direction_lik"];
    mat V_list1_first = phyl_list["V_list1_first"];
    mat V_list2_first = phyl_list["V_list2_first"];
    mat V_list12_first = phyl_list["V_list12_first"];
    mat V_list1_second = phyl_list["V_list1_second"];
    mat V_list2_second = phyl_list["V_list2_second"];
    mat V_list12_second = phyl_list["V_list12_second"];
    mat W_list1 = phyl_list["W_list1"];
    mat W_list2 = phyl_list["W_list2"];
    mat W_list12 = phyl_list["W_list12"];
    
    double post_mean1_i = sum( (V_list1_first.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) * root_dist ) /
    sum(downward_lik.row(num_term_nodes) * root_dist);
    
    post_mean1 += post_mean1_i;
    
    double post_mean2_i = sum( (V_list2_first.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) * root_dist ) /
    sum(downward_lik.row(num_term_nodes) * root_dist);
    
    post_mean2 += post_mean2_i;
    
    double post_var1_i = sum( ( ((V_list1_second.rows(root_edges_ind) + V_list1_first.rows(root_edges_ind)) %
    direction_lik.rows(root_edges_ind_rev)) + (W_list1.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) +
    (V_list1_first.rows(root_edges_ind) % V_list1_first.rows(root_edges_ind_rev)) ) * root_dist ) /
    sum(downward_lik.row(num_term_nodes) * root_dist) - pow(post_mean1_i, 2);
    
    post_var1 += post_var1_i;
    
    double post_var2_i = sum( ( ((V_list2_second.rows(root_edges_ind) + V_list2_first.rows(root_edges_ind)) %
    direction_lik.rows(root_edges_ind_rev)) + (W_list2.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) +
    (V_list2_first.rows(root_edges_ind) % V_list2_first.rows(root_edges_ind_rev)) ) * root_dist ) /
    sum(downward_lik.row(num_term_nodes) * root_dist) - pow(post_mean2_i, 2);
    
    post_var2 += post_var2_i;
    
    double post_cov_i = sum( ( ((V_list12_second.rows(root_edges_ind) + V_list12_first.rows(root_edges_ind)) %
    direction_lik.rows(root_edges_ind_rev)) + (W_list12.rows(root_edges_ind) % direction_lik.rows(root_edges_ind_rev)) +
    (V_list1_first.rows(root_edges_ind) % V_list2_first.rows(root_edges_ind_rev)) ) * root_dist ) /
    sum(downward_lik.row(num_term_nodes) * root_dist) - post_mean1_i * post_mean2_i;
    
    post_cov += post_cov_i;
  }
  
  return NumericVector::create(
    _["mean1"] = post_mean1,
    _["mean2"] = post_mean2,
    _["var1"] = post_var1,
    _["var2"] = post_var2,
    _["cov"] = post_cov
    );
}






// [[Rcpp::export]]
NumericVector joint_postmean_moments_phylojumps(const arma::imat& edge_mat, const arma::vec& edge_lengths, const arma::mat& rate_mat, const arma::ivec& edge_set1,
const arma::ivec& edge_set2, const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments, const int& N) {
  
  // calculating prior moments
  NumericVector prior_moments = joint_prior_moments_phylojumps(edge_mat, edge_set1, edge_set2,
  root_dist, num_edges, num_states, num_term_nodes, edge_moments);
  
  double prior_mean1 = prior_moments["mean1"];
  double prior_mean2 = prior_moments["mean2"];
  double prior_var1 = prior_moments["var1"];
  double prior_var2 = prior_moments["var2"];
  double prior_cov = prior_moments["cov"];
  
  // getting monte carlo estimates of the mean posterior variances/covariance
  imat seq_sim = tips_sim(edge_mat, edge_lengths, rate_mat, root_dist, num_edges, num_states, num_term_nodes, N);
  
  NumericVector post_moments = joint_post_moments_phylojumps(edge_mat, edge_set1, edge_set2,
  root_dist, num_edges, num_states, num_term_nodes, edge_moments, seq_sim);
  
  double mean_post_var1 = post_moments["var1"] / N;
  double mean_post_var2 = post_moments["var2"] / N;
  double mean_post_cov = post_moments["cov"] / N;
  
  return NumericVector::create(
    _["mean1"] = prior_mean1,
    _["mean2"] = prior_mean2,
    _["var1"] = prior_var1 - mean_post_var1,
    _["var2"] = prior_var2 - mean_post_var2,
    _["cov"] = prior_cov - mean_post_cov
    );
}


