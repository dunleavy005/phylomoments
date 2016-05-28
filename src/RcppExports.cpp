// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// moments_ctmcjumps
List moments_ctmcjumps(const arma::vec& edge_lengths, const arma::mat& rate_mat, const arma::imat& label_mat, const int& num_edges, const int& num_states);
RcppExport SEXP phylomoments_moments_ctmcjumps(SEXP edge_lengthsSEXP, SEXP rate_matSEXP, SEXP label_matSEXP, SEXP num_edgesSEXP, SEXP num_statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec& >::type edge_lengths(edge_lengthsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rate_mat(rate_matSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type label_mat(label_matSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_edges(num_edgesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_states(num_statesSEXP);
    __result = Rcpp::wrap(moments_ctmcjumps(edge_lengths, rate_mat, label_mat, num_edges, num_states));
    return __result;
END_RCPP
}
// ctmc_sim
NumericMatrix ctmc_sim(const double& t, const arma::mat& rate_mat, const int& init_state);
RcppExport SEXP phylomoments_ctmc_sim(SEXP tSEXP, SEXP rate_matSEXP, SEXP init_stateSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const double& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rate_mat(rate_matSEXP);
    Rcpp::traits::input_parameter< const int& >::type init_state(init_stateSEXP);
    __result = Rcpp::wrap(ctmc_sim(t, rate_mat, init_state));
    return __result;
END_RCPP
}
// phylo_log_likelihood
double phylo_log_likelihood(const arma::imat& edge_mat, const arma::vec& edge_lengths, const arma::mat& rate_mat, const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const arma::imat& seq_data);
RcppExport SEXP phylomoments_phylo_log_likelihood(SEXP edge_matSEXP, SEXP edge_lengthsSEXP, SEXP rate_matSEXP, SEXP root_distSEXP, SEXP num_edgesSEXP, SEXP num_statesSEXP, SEXP num_term_nodesSEXP, SEXP seq_dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::imat& >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type edge_lengths(edge_lengthsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rate_mat(rate_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type root_dist(root_distSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_edges(num_edgesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_term_nodes(num_term_nodesSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type seq_data(seq_dataSEXP);
    __result = Rcpp::wrap(phylo_log_likelihood(edge_mat, edge_lengths, rate_mat, root_dist, num_edges, num_states, num_term_nodes, seq_data));
    return __result;
END_RCPP
}
// int_states_sim
arma::imat int_states_sim(const arma::imat& edge_mat, const arma::vec& edge_lengths, const arma::mat& rate_mat, const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const arma::ivec& tip_data, const int& N);
RcppExport SEXP phylomoments_int_states_sim(SEXP edge_matSEXP, SEXP edge_lengthsSEXP, SEXP rate_matSEXP, SEXP root_distSEXP, SEXP num_edgesSEXP, SEXP num_statesSEXP, SEXP num_term_nodesSEXP, SEXP tip_dataSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::imat& >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type edge_lengths(edge_lengthsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rate_mat(rate_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type root_dist(root_distSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_edges(num_edgesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_term_nodes(num_term_nodesSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type tip_data(tip_dataSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    __result = Rcpp::wrap(int_states_sim(edge_mat, edge_lengths, rate_mat, root_dist, num_edges, num_states, num_term_nodes, tip_data, N));
    return __result;
END_RCPP
}
// prior_moments_phylojumps
NumericVector prior_moments_phylojumps(const arma::imat& edge_mat, const arma::ivec& edge_set, const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments);
RcppExport SEXP phylomoments_prior_moments_phylojumps(SEXP edge_matSEXP, SEXP edge_setSEXP, SEXP root_distSEXP, SEXP num_edgesSEXP, SEXP num_statesSEXP, SEXP num_term_nodesSEXP, SEXP edge_momentsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::imat& >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type edge_set(edge_setSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type root_dist(root_distSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_edges(num_edgesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_term_nodes(num_term_nodesSEXP);
    Rcpp::traits::input_parameter< const List& >::type edge_moments(edge_momentsSEXP);
    __result = Rcpp::wrap(prior_moments_phylojumps(edge_mat, edge_set, root_dist, num_edges, num_states, num_term_nodes, edge_moments));
    return __result;
END_RCPP
}
// post_moments_phylojumps
NumericVector post_moments_phylojumps(const arma::imat& edge_mat, const arma::ivec& edge_set, const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments, const arma::imat& seq_data);
RcppExport SEXP phylomoments_post_moments_phylojumps(SEXP edge_matSEXP, SEXP edge_setSEXP, SEXP root_distSEXP, SEXP num_edgesSEXP, SEXP num_statesSEXP, SEXP num_term_nodesSEXP, SEXP edge_momentsSEXP, SEXP seq_dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::imat& >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type edge_set(edge_setSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type root_dist(root_distSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_edges(num_edgesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_term_nodes(num_term_nodesSEXP);
    Rcpp::traits::input_parameter< const List& >::type edge_moments(edge_momentsSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type seq_data(seq_dataSEXP);
    __result = Rcpp::wrap(post_moments_phylojumps(edge_mat, edge_set, root_dist, num_edges, num_states, num_term_nodes, edge_moments, seq_data));
    return __result;
END_RCPP
}
// find_subtree_edges
arma::ivec find_subtree_edges(const int& node, const arma::imat& edge_mat);
RcppExport SEXP phylomoments_find_subtree_edges(SEXP nodeSEXP, SEXP edge_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const int& >::type node(nodeSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type edge_mat(edge_matSEXP);
    __result = Rcpp::wrap(find_subtree_edges(node, edge_mat));
    return __result;
END_RCPP
}
// tips_sim
arma::imat tips_sim(const arma::imat& edge_mat, const arma::vec& edge_lengths, const arma::mat& rate_mat, const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const int& N);
RcppExport SEXP phylomoments_tips_sim(SEXP edge_matSEXP, SEXP edge_lengthsSEXP, SEXP rate_matSEXP, SEXP root_distSEXP, SEXP num_edgesSEXP, SEXP num_statesSEXP, SEXP num_term_nodesSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::imat& >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type edge_lengths(edge_lengthsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rate_mat(rate_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type root_dist(root_distSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_edges(num_edgesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_term_nodes(num_term_nodesSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    __result = Rcpp::wrap(tips_sim(edge_mat, edge_lengths, rate_mat, root_dist, num_edges, num_states, num_term_nodes, N));
    return __result;
END_RCPP
}
// phylojumps_sim
arma::imat phylojumps_sim(const arma::imat& edge_mat, const arma::vec& edge_lengths, const arma::mat& rate_mat, const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const arma::imat& seq_data, const int& N);
RcppExport SEXP phylomoments_phylojumps_sim(SEXP edge_matSEXP, SEXP edge_lengthsSEXP, SEXP rate_matSEXP, SEXP root_distSEXP, SEXP num_edgesSEXP, SEXP num_statesSEXP, SEXP num_term_nodesSEXP, SEXP seq_dataSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::imat& >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type edge_lengths(edge_lengthsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rate_mat(rate_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type root_dist(root_distSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_edges(num_edgesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_term_nodes(num_term_nodesSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type seq_data(seq_dataSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    __result = Rcpp::wrap(phylojumps_sim(edge_mat, edge_lengths, rate_mat, root_dist, num_edges, num_states, num_term_nodes, seq_data, N));
    return __result;
END_RCPP
}
// postmean_moments_phylojumps
NumericVector postmean_moments_phylojumps(const arma::imat& edge_mat, const arma::vec& edge_lengths, const arma::mat& rate_mat, const arma::ivec& edge_set, const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments, const int& N);
RcppExport SEXP phylomoments_postmean_moments_phylojumps(SEXP edge_matSEXP, SEXP edge_lengthsSEXP, SEXP rate_matSEXP, SEXP edge_setSEXP, SEXP root_distSEXP, SEXP num_edgesSEXP, SEXP num_statesSEXP, SEXP num_term_nodesSEXP, SEXP edge_momentsSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::imat& >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type edge_lengths(edge_lengthsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rate_mat(rate_matSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type edge_set(edge_setSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type root_dist(root_distSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_edges(num_edgesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_term_nodes(num_term_nodesSEXP);
    Rcpp::traits::input_parameter< const List& >::type edge_moments(edge_momentsSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    __result = Rcpp::wrap(postmean_moments_phylojumps(edge_mat, edge_lengths, rate_mat, edge_set, root_dist, num_edges, num_states, num_term_nodes, edge_moments, N));
    return __result;
END_RCPP
}
// joint_prior_moments_phylojumps
NumericVector joint_prior_moments_phylojumps(const arma::imat& edge_mat, const arma::ivec& edge_set1, const arma::ivec& edge_set2, const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments);
RcppExport SEXP phylomoments_joint_prior_moments_phylojumps(SEXP edge_matSEXP, SEXP edge_set1SEXP, SEXP edge_set2SEXP, SEXP root_distSEXP, SEXP num_edgesSEXP, SEXP num_statesSEXP, SEXP num_term_nodesSEXP, SEXP edge_momentsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::imat& >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type edge_set1(edge_set1SEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type edge_set2(edge_set2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type root_dist(root_distSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_edges(num_edgesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_term_nodes(num_term_nodesSEXP);
    Rcpp::traits::input_parameter< const List& >::type edge_moments(edge_momentsSEXP);
    __result = Rcpp::wrap(joint_prior_moments_phylojumps(edge_mat, edge_set1, edge_set2, root_dist, num_edges, num_states, num_term_nodes, edge_moments));
    return __result;
END_RCPP
}
// joint_post_moments_phylojumps
NumericVector joint_post_moments_phylojumps(const arma::imat& edge_mat, const arma::ivec& edge_set1, const arma::ivec& edge_set2, const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments, const arma::imat& seq_data);
RcppExport SEXP phylomoments_joint_post_moments_phylojumps(SEXP edge_matSEXP, SEXP edge_set1SEXP, SEXP edge_set2SEXP, SEXP root_distSEXP, SEXP num_edgesSEXP, SEXP num_statesSEXP, SEXP num_term_nodesSEXP, SEXP edge_momentsSEXP, SEXP seq_dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::imat& >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type edge_set1(edge_set1SEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type edge_set2(edge_set2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type root_dist(root_distSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_edges(num_edgesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_term_nodes(num_term_nodesSEXP);
    Rcpp::traits::input_parameter< const List& >::type edge_moments(edge_momentsSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type seq_data(seq_dataSEXP);
    __result = Rcpp::wrap(joint_post_moments_phylojumps(edge_mat, edge_set1, edge_set2, root_dist, num_edges, num_states, num_term_nodes, edge_moments, seq_data));
    return __result;
END_RCPP
}
// joint_postmean_moments_phylojumps
NumericVector joint_postmean_moments_phylojumps(const arma::imat& edge_mat, const arma::vec& edge_lengths, const arma::mat& rate_mat, const arma::ivec& edge_set1, const arma::ivec& edge_set2, const arma::vec& root_dist, const int& num_edges, const int& num_states, const int& num_term_nodes, const List& edge_moments, const int& N);
RcppExport SEXP phylomoments_joint_postmean_moments_phylojumps(SEXP edge_matSEXP, SEXP edge_lengthsSEXP, SEXP rate_matSEXP, SEXP edge_set1SEXP, SEXP edge_set2SEXP, SEXP root_distSEXP, SEXP num_edgesSEXP, SEXP num_statesSEXP, SEXP num_term_nodesSEXP, SEXP edge_momentsSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::imat& >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type edge_lengths(edge_lengthsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rate_mat(rate_matSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type edge_set1(edge_set1SEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type edge_set2(edge_set2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type root_dist(root_distSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_edges(num_edgesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_term_nodes(num_term_nodesSEXP);
    Rcpp::traits::input_parameter< const List& >::type edge_moments(edge_momentsSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    __result = Rcpp::wrap(joint_postmean_moments_phylojumps(edge_mat, edge_lengths, rate_mat, edge_set1, edge_set2, root_dist, num_edges, num_states, num_term_nodes, edge_moments, N));
    return __result;
END_RCPP
}
