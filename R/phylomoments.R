



#' phylomoments: Efficient computations of phylogenetic stochastic mapping 
#' moments.
#' 
#' This package implements a simulation-free dynamic programming algorithm for 
#' calculating higher-order moments of stochastic mapping summaries on a 
#' phylogeny. Utility functions related to simulation-based stochastic mapping 
#' and continuous-time Markov chains (CTMCs) are also provided.
#' 
#' @section : Please see the vignettes for some examples of using the functions 
#'   provided in this package.  They can be accessed using 
#'   \code{browseVignettes("phylomoments")}.
#'   
#' @author Amrit Dhar \email{adhar@@uw.edu}
#'   
#' @references OUR PAPER!!!
#'   
#'   R. Nielsen (2002) \dQuote{Mapping mutations on phylogenies,} 
#'   \emph{Systematic Biology,} 51(5):729-739.
#'   
#'   V.N. Minin and M.A. Suchard (2008) \dQuote{Counting labeled transitions in 
#'   continuous-time Markov models of evolution,} \emph{Journal of Mathematical 
#'   Biology,} 56(3):391-412.
#'   
#' @docType package
#' @name phylomoments
#'   
#' @useDynLib phylomoments
#' @importFrom Rcpp evalCpp
NULL




scale.rate.mat = function(rate.mat, root.dist) rate.mat / sum(-diag(rate.mat) * root.dist)






#' Restricted factorial moments of labeled substitution counts for a reversible 
#' CTMC.
#' 
#' This function calculates the 0th, 1st, and 2nd restricted factorial moment 
#' matrices of labeled substitution counts for a reversible continuous-time 
#' Markov chain (CTMC) model.
#' 
#' The 0th restricted factorial moment matrix of labeled CTMC substitution 
#' counts is defined as the CTMC transition probability matrix.
#' 
#' @usage moments.ctmcjumps(t, rate.mat, label.mat,
#'                   root.dist = NULL, scale = FALSE)
#' 
#' @param t A nonnegative numeric vector of CTMC time interval lengths.
#' @param rate.mat A reversible CTMC rate matrix.
#' @param label.mat A 0-1 matrix of the same size as \code{rate.mat}.  This 
#'   matrix labels the substitutions of interest.
#' @param root.dist A numeric vector defining the CTMC stationary distribution. 
#'   This only needs to be specified when \code{scale = TRUE}.
#' @param scale A logical indicating whether to scale the time dimension.  If 
#'   \code{TRUE}, then time is specified in terms of the expected number of CTMC
#'   substitutions.
#'   
#' @return A list with named elements \code{"zeroth"}, \code{"first"}, and 
#'   \code{"second"}.  Each list element is a three-dimensional array with 
#'   dimensions \code{c(nrow(rate.mat), ncol(rate.mat), length(t))}.  The 
#'   \code{"zeroth"}, \code{"first"}, and \code{"second"} elements of this list 
#'   contain the 0th, 1st, and 2nd restricted factorial moment matrices of 
#'   labeled CTMC substitution counts, respectively, for the time interval 
#'   lengths defined by \code{t}.
#'   
#' @references V.N. Minin and M.A. Suchard (2008) \dQuote{Counting labeled 
#'   transitions in continuous-time Markov models of evolution,} \emph{Journal 
#'   of Mathematical Biology,} 56(3):391-412.
#'   
#' @seealso \code{\link{ctmc.sim}}
#'   
#' @export moments.ctmcjumps
moments.ctmcjumps = function(t, rate.mat, label.mat, root.dist = NULL, scale = FALSE) {
  
  if(any(dim(rate.mat) != dim(label.mat)))
    stop("Dimensions of 'rate.mat' and 'label.mat' don't match", call. = FALSE)
  
  if(scale == TRUE & is.null(root.dist))
    stop("'root.dist' must be specified to scale 'rate.mat'", call. = FALSE)
  
  if(scale == TRUE & any(dim(rate.mat) != length(root.dist)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'root.dist'", call. = FALSE)
  
  if(scale == TRUE) rate.mat = scale.rate.mat(rate.mat, root.dist)
  num.edges = length(t)
  num.states = nrow(rate.mat)
  
  moments_ctmcjumps(t, rate.mat, label.mat, num.edges, num.states)
}







#' Simulation of a CTMC path.
#' 
#' Generates a realization of a continuous-time Markov chain (CTMC).
#' 
#' @usage ctmc.sim(t, rate.mat, states = c("a","c","g","t"), init.state)
#' 
#' @param t A nonnegative numeric scalar representing the CTMC time interval 
#'   length.
#' @param rate.mat A CTMC rate matrix.
#' @param states A character or integer vector denoting the CTMC state space. 
#'   Defaults to the set of DNA nucleotides \code{c("a","c","g","t")}.
#' @param init.state The initial state of the CTMC.  This state must be an 
#'   element of \code{states}.
#'   
#' @return A data frame with two named columns \code{"state"} and \code{"time"}.
#'   The \code{"state"} column stores the observed states of the CTMC path and 
#'   the \code{"time"} column specifies the times at which each of the observed 
#'   states is entered.
#'   
#' @seealso \code{\link{moments.ctmcjumps}}
#'   
#' @export ctmc.sim
ctmc.sim = function(t, rate.mat, states = c("a","c","g","t"), init.state) {
  
  if(!(init.state %in% states))
    stop("'init.state' is not an element of 'states'", call. = FALSE)
  
  if(any(dim(rate.mat) != length(states)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'states'", call. = FALSE)
  
  init.state = match(init.state, states) - 1
  
  outp = data.frame(t(ctmc_sim(t, rate.mat, init.state)))
  outp[,1] = states[outp[,1] + 1]
  colnames(outp) = c("state","time")
  
  outp
}






#' Simulation of internal node states on a phylogenetic tree.
#' 
#' Samples internal node states on a phylogeny conditional on the observed tip 
#' states.
#' 
#' This simulation function follows the procedure outlined in Nielsen (2002).
#' 
#' More information on phylogenetic tree objects of class \code{"phylo"} can be 
#' found at \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
#' 
#' @usage int.states.sim(tree, rate.mat, root.dist, scale = FALSE,
#'                states = c("a","c","g","t"), tip.data, N)
#' 
#' @param tree A tree object of class \code{"phylo"} (\strong{ape} format).
#' @param rate.mat A continuous-time Markov chain (CTMC) rate matrix.
#' @param root.dist A numeric vector defining the root distribution of the 
#'   evolutionary process on the tree. If this process starts at stationarity 
#'   (as is commonly assumed), then the root distribution is equal to the CTMC 
#'   stationary distribution.
#' @param scale A logical indicating whether to scale the time dimension.  If 
#'   \code{TRUE}, then time is specified in terms of the expected number of CTMC
#'   substitutions per site.
#' @param states A character or integer vector denoting the CTMC state space. 
#'   Defaults to the set of DNA nucleotides \code{c("a","c","g","t")}.
#' @param tip.data A character or integer vector of observed tip states.  These 
#'   tip states must be contained in the vector \code{states}.
#' @param N An integer specifying the number of samples to draw.
#'   
#' @return A character or integer matrix with dimensions 
#'   \code{c(2*length(tree$tip.label)-1, N)}, where each column contains an 
#'   independent sample of internal node states along with the observed tip 
#'   states.
#'   
#' @references R. Nielsen (2002) \dQuote{Mapping mutations on phylogenies,} 
#'   \emph{Systematic Biology,} 51(5):729-739.
#'   
#' @seealso \code{\link{tips.sim}}, \code{\link{phylojumps.sim}}
#'   
#' @export int.states.sim
int.states.sim = function(tree, rate.mat, root.dist, scale = FALSE, states = c("a","c","g","t"), tip.data, N) {
  
  if(!inherits(tree, "phylo"))
    stop("Tree object must be of class 'phylo'", call. = FALSE)
  
  if(!all(tip.data %in% states))
    stop("Elements of 'tip.data' are not consistent with 'states'", call. = FALSE)
  
  if(length(root.dist) != length(states))
    stop("Dimensions of 'root.dist' and 'states' don't match", call. = FALSE)
  
  if(any(dim(rate.mat) != length(states)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'states'", call. = FALSE)
  
  edge.mat = tree$edge
  edge.lengths = tree$edge.length
  if(scale == TRUE) rate.mat = scale.rate.mat(rate.mat, root.dist)
  num.edges = nrow(edge.mat)
  num.states = nrow(rate.mat)
  num.term.nodes = length(tree$tip.label)
  tip.data = match(tip.data, states) - 1
  
  if(length(tip.data) != num.term.nodes)
    stop("'tip.data' doesn't have the correct number of elements", call. = FALSE)
  
  outp = int_states_sim(edge.mat, edge.lengths, rate.mat, root.dist,
                        num.edges, num.states, num.term.nodes, tip.data, N)
  
  matrix(states[outp + 1], ncol = N, byrow = FALSE)
}






#' Simulation of terminal node states on a phylogenetic tree.
#' 
#' Simulates tip states on a phylogeny according to a continuous-time Markov 
#' model of evolution.
#' 
#' More information on phylogenetic tree objects of class \code{"phylo"} can be 
#' found at \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
#' 
#' @usage tips.sim(tree, rate.mat, root.dist, scale = FALSE,
#'          states = c("a","c","g","t"), N)
#' 
#' @param tree A tree object of class \code{"phylo"} (\strong{ape} format).
#' @param rate.mat A continuous-time Markov chain (CTMC) rate matrix.
#' @param root.dist A numeric vector defining the root distribution of the 
#'   evolutionary process on the tree. If this process starts at stationarity 
#'   (as is commonly assumed), then the root distribution is equal to the CTMC 
#'   stationary distribution.
#' @param scale A logical indicating whether to scale the time dimension.  If 
#'   \code{TRUE}, then time is specified in terms of the expected number of CTMC
#'   substitutions per site.
#' @param states A character or integer vector denoting the CTMC state space. 
#'   Defaults to the set of DNA nucleotides \code{c("a","c","g","t")}.
#' @param N An integer specifying the number of samples to draw.
#'   
#' @return A character or integer matrix with dimensions 
#'   \code{c(length(tree$tip.label), N)}, where each column contains an 
#'   independent sample of terminal node states.
#'   
#' @seealso \code{\link{int.states.sim}}, \code{\link{phylojumps.sim}}
#'   
#' @export tips.sim
tips.sim = function(tree, rate.mat, root.dist, scale = FALSE, states = c("a","c","g","t"), N) {
  
  if(!inherits(tree, "phylo"))
    stop("Tree object must be of class 'phylo'", call. = FALSE)
  
  if(attr(tree, "order") != "pruningwise")
    stop("Edge matrix must be in 'pruningwise' order", call. = FALSE)
  
  if(length(root.dist) != length(states))
    stop("Dimensions of 'root.dist' and 'states' don't match", call. = FALSE)
  
  if(any(dim(rate.mat) != length(states)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'states'", call. = FALSE)
  
  edge.mat = tree$edge
  edge.lengths = tree$edge.length
  if(scale == TRUE) rate.mat = scale.rate.mat(rate.mat, root.dist)
  num.edges = nrow(edge.mat)
  num.states = nrow(rate.mat)
  num.term.nodes = length(tree$tip.label)
  
  outp = tips_sim(edge.mat, edge.lengths, rate.mat, root.dist,
                  num.edges, num.states, num.term.nodes, N)
  
  matrix(states[outp + 1], ncol = N, byrow = FALSE)
}






#' Simulation of substitution counts on a phylogenetic tree.
#' 
#' Simulates the number of substitutions on a phylogeny conditional on the 
#' observed tip states.
#' 
#' This simulation function follows the procedure outlined in Nielsen (2002).
#' 
#' More information on phylogenetic tree objects of class \code{"phylo"} can be 
#' found at \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
#' 
#' @usage phylojumps.sim(tree, rate.mat, root.dist, scale = FALSE,
#'                states = c("a","c","g","t"), seq.data, N)
#' 
#' @param tree A tree object of class \code{"phylo"} (\strong{ape} format).
#' @param rate.mat A continuous-time Markov chain (CTMC) rate matrix.
#' @param root.dist A numeric vector defining the root distribution of the 
#'   evolutionary process on the tree. If this process starts at stationarity 
#'   (as is commonly assumed), then the root distribution is equal to the CTMC 
#'   stationary distribution.
#' @param scale A logical indicating whether to scale the time dimension.  If 
#'   \code{TRUE}, then time is specified in terms of the expected number of CTMC
#'   substitutions per site.
#' @param states A character or integer vector denoting the CTMC state space. 
#'   Defaults to the set of DNA nucleotides \code{c("a","c","g","t")}.
#' @param seq.data A character or integer matrix representing the observed 
#'   sequence alignment with rows (columns) that correspond to the sequences 
#'   (sites) under study.  All elements of \code{seq.data} must be contained in 
#'   the vector \code{states}.
#' @param N An integer specifying the number of samples to draw.
#'   
#' @return An integer matrix with dimensions \code{c(N, ncol(seq.data))}, where 
#'   the \emph{i}'th column contains \code{N} independent samples of 
#'   substitution counts for the \emph{i}'th site in \code{seq.data}.
#'   
#' @references R. Nielsen (2002) \dQuote{Mapping mutations on phylogenies,} 
#'   \emph{Systematic Biology,} 51(5):729-739.
#'   
#' @seealso \code{\link{int.states.sim}}, \code{\link{tips.sim}}
#'   
#' @export phylojumps.sim
phylojumps.sim = function(tree, rate.mat, root.dist, scale = FALSE, states = c("a","c","g","t"), seq.data, N) {
  
  if(!inherits(tree, "phylo"))
    stop("Tree object must be of class 'phylo'", call. = FALSE)
  
  if(attr(tree, "order") != "pruningwise")
    stop("Edge matrix must be in 'pruningwise' order", call. = FALSE)
  
  if(!all(seq.data %in% states))
    stop("Elements of 'seq.data' are not consistent with 'states'", call. = FALSE)
  
  if(length(root.dist) != length(states))
    stop("Dimensions of 'root.dist' and 'states' don't match", call. = FALSE)
  
  if(any(dim(rate.mat) != length(states)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'states'", call. = FALSE)
  
  edge.mat = tree$edge
  edge.lengths = tree$edge.length
  if(scale == TRUE) rate.mat = scale.rate.mat(rate.mat, root.dist)
  num.edges = nrow(edge.mat)
  num.states = nrow(rate.mat)
  num.term.nodes = length(tree$tip.label)
  seq.data = matrix(match(seq.data, states) - 1, ncol = ncol(seq.data), byrow = FALSE)
  
  if(nrow(seq.data) != num.term.nodes)
    stop("'seq.data' doesn't have the correct number of rows", call. = FALSE)
  
  phylojumps_sim(edge.mat, edge.lengths, rate.mat, root.dist,
                 num.edges, num.states, num.term.nodes, seq.data, N)
}






#' Prior moments of labeled substitution counts on a phylogenetic tree.
#' 
#' This function computes the prior mean and variance of labeled substitution 
#' counts over a prespecified set of branches on a phylogeny using a 
#' simulation-free dynamic programming algorithm.
#' 
#' More information on phylogenetic tree objects of class \code{"phylo"} can be 
#' found at \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
#' 
#' @usage prior.moments.phylojumps(tree, rate.mat, label.mat,
#'                          edge.set, root.dist, scale = FALSE)
#' 
#' @param tree A tree object of class \code{"phylo"} (\strong{ape} format).
#' @param rate.mat A reversible continuous-time Markov chain (CTMC) rate matrix.
#' @param label.mat A 0-1 matrix of the same size as \code{rate.mat}.  This 
#'   matrix labels the substitutions of interest.
#' @param edge.set An integer vector of branch indices specifying the branch set
#'   over which the prior moments will be calculated.
#' @param root.dist A numeric vector defining the root distribution of the 
#'   evolutionary process on the tree. If this process starts at stationarity 
#'   (as is commonly assumed), then the root distribution is equal to the CTMC 
#'   stationary distribution.
#' @param scale A logical indicating whether to scale the time dimension.  If 
#'   \code{TRUE}, then time is specified in terms of the expected number of CTMC
#'   substitutions per site.
#'   
#' @return A numeric vector holding the prior mean and variance of labeled 
#'   substitution counts.
#'   
#' @references OUR PAPER!!!
#' 
#' @seealso \code{\link{post.moments.phylojumps}}, 
#'   \code{\link{postmean.moments.phylojumps}}, 
#'   \code{\link{joint.prior.moments.phylojumps}}
#'   
#' @export prior.moments.phylojumps
prior.moments.phylojumps = function(tree, rate.mat, label.mat, edge.set, root.dist, scale = FALSE) {
  
  if(!inherits(tree, "phylo"))
    stop("Tree object must be of class 'phylo'", call. = FALSE)
  
  if(attr(tree, "order") != "pruningwise")
    stop("Edge matrix must be in 'pruningwise' order", call. = FALSE)
  
  if(any(dim(rate.mat) != dim(label.mat)))
    stop("Dimensions of 'rate.mat' and 'label.mat' don't match", call. = FALSE)
  
  if(any(dim(rate.mat) != length(root.dist)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'root.dist'", call. = FALSE)
  
  edge.mat = tree$edge
  edge.lengths = tree$edge.length
  if(scale == TRUE) rate.mat = scale.rate.mat(rate.mat, root.dist)
  num.edges = nrow(edge.mat)
  num.states = nrow(rate.mat)
  num.term.nodes = length(tree$tip.label)
  edge.moments = moments_ctmcjumps(edge.lengths, rate.mat, label.mat, num.edges, num.states)
  
  if(any(edge.set > num.edges | edge.set < 1))
    stop("Nonexistent edge(s) passed to 'edge.set'", call. = FALSE)
  
  prior_moments_phylojumps(edge.mat, edge.set, root.dist, num.edges,
                           num.states, num.term.nodes, edge.moments)
}






#' Posterior moments of labeled substitution counts on a phylogenetic tree.
#' 
#' This function computes the posterior mean and variance of labeled 
#' substitution counts over a prespecified set of branches on a phylogeny using 
#' a simulation-free dynamic programming algorithm.
#' 
#' More information on phylogenetic tree objects of class \code{"phylo"} can be 
#' found at \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
#' 
#' @usage post.moments.phylojumps(tree, rate.mat, label.mat,
#'                         edge.set, root.dist, scale = FALSE,
#'                         states = c("a","c","g","t"), seq.data)
#' 
#' @param tree A tree object of class \code{"phylo"} (\strong{ape} format).
#' @param rate.mat A reversible continuous-time Markov chain (CTMC) rate matrix.
#' @param label.mat A 0-1 matrix of the same size as \code{rate.mat}.  This 
#'   matrix labels the substitutions of interest.
#' @param edge.set An integer vector of branch indices specifying the branch set
#'   over which the posterior moments will be calculated.
#' @param root.dist A numeric vector defining the root distribution of the 
#'   evolutionary process on the tree. If this process starts at stationarity 
#'   (as is commonly assumed), then the root distribution is equal to the CTMC 
#'   stationary distribution.
#' @param scale A logical indicating whether to scale the time dimension.  If 
#'   \code{TRUE}, then time is specified in terms of the expected number of CTMC
#'   substitutions per site.
#' @param states A character or integer vector denoting the CTMC state space. 
#'   Defaults to the set of DNA nucleotides \code{c("a","c","g","t")}.
#' @param seq.data A character or integer matrix representing the observed 
#'   sequence alignment with rows (columns) that correspond to the sequences 
#'   (sites) under study.  All elements of \code{seq.data} must be contained in 
#'   the vector \code{states}.
#'   
#' @return A numeric vector holding the posterior mean and variance of labeled 
#'   substitution counts summed over all sites in \code{seq.data}.
#'   
#' @references OUR PAPER!!!
#' 
#' @seealso \code{\link{prior.moments.phylojumps}}, 
#'   \code{\link{postmean.moments.phylojumps}}, 
#'   \code{\link{joint.post.moments.phylojumps}}
#'   
#' @export post.moments.phylojumps
post.moments.phylojumps = function(tree, rate.mat, label.mat, edge.set, root.dist, scale = FALSE, states = c("a","c","g","t"), seq.data) {
  
  if(!inherits(tree, "phylo"))
    stop("Tree object must be of class 'phylo'", call. = FALSE)
  
  if(attr(tree, "order") != "pruningwise")
    stop("Edge matrix must be in 'pruningwise' order", call. = FALSE)
  
  if(!all(seq.data %in% states))
    stop("Elements of 'seq.data' are not consistent with 'states'", call. = FALSE)
  
  if(any(dim(rate.mat) != dim(label.mat)))
    stop("Dimensions of 'rate.mat' and 'label.mat' don't match", call. = FALSE)
  
  if(any(dim(rate.mat) != length(root.dist)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'root.dist'", call. = FALSE)
  
  if(any(dim(rate.mat) != length(states)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'states'", call. = FALSE)
  
  edge.mat = tree$edge
  edge.lengths = tree$edge.length
  if(scale == TRUE) rate.mat = scale.rate.mat(rate.mat, root.dist)
  num.edges = nrow(edge.mat)
  num.states = nrow(rate.mat)
  num.term.nodes = length(tree$tip.label)
  edge.moments = moments_ctmcjumps(edge.lengths, rate.mat, label.mat, num.edges, num.states)
  seq.data = matrix(match(seq.data, states) - 1, ncol = ncol(seq.data), byrow = FALSE)
  
  if(any(edge.set > num.edges | edge.set < 1))
    stop("Nonexistent edge(s) passed to 'edge.set'", call. = FALSE)
  
  if(nrow(seq.data) != num.term.nodes)
    stop("'seq.data' doesn't have the correct number of rows", call. = FALSE)
  
  post_moments_phylojumps(edge.mat, edge.set, root.dist, num.edges,
                          num.states, num.term.nodes, edge.moments, seq.data)
}






#' Moments of the posterior mean of labeled substitution counts on a 
#' phylogenetic tree.
#' 
#' This function calculates the mean and variance of the posterior mean of 
#' labeled substitution counts over a prespecified set of branches on a 
#' phylogeny using a simulation-free dynamic programming algorithm and Monte 
#' Carlo sampling.
#' 
#' The mean obtained in this function is computed exactly, while the variance is
#' partially approximated using Monte Carlo simulations.
#' 
#' More information on phylogenetic tree objects of class \code{"phylo"} can be 
#' found at \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
#' 
#' @usage postmean.moments.phylojumps(tree, rate.mat, label.mat,
#'                             edge.set, root.dist,
#'                             scale = FALSE, N = 100000)
#' 
#' @param tree A tree object of class \code{"phylo"} (\strong{ape} format).
#' @param rate.mat A reversible continuous-time Markov chain (CTMC) rate matrix.
#' @param label.mat A 0-1 matrix of the same size as \code{rate.mat}.  This 
#'   matrix labels the substitutions of interest.
#' @param edge.set An integer vector of branch indices specifying the branch set
#'   over which the moments will be calculated.
#' @param root.dist A numeric vector defining the root distribution of the 
#'   evolutionary process on the tree. If this process starts at stationarity 
#'   (as is commonly assumed), then the root distribution is equal to the CTMC 
#'   stationary distribution.
#' @param scale A logical indicating whether to scale the time dimension.  If 
#'   \code{TRUE}, then time is specified in terms of the expected number of CTMC
#'   substitutions per site.
#' @param N An integer specifying the number of Monte Carlo samples to draw.
#'   
#' @return A numeric vector holding the mean and variance of the posterior mean 
#'   of labeled substitution counts.
#'   
#' @references OUR PAPER!!!
#' 
#' @seealso \code{\link{prior.moments.phylojumps}}, 
#'   \code{\link{post.moments.phylojumps}}, 
#'   \code{\link{joint.postmean.moments.phylojumps}}
#'   
#' @export postmean.moments.phylojumps
postmean.moments.phylojumps = function(tree, rate.mat, label.mat, edge.set, root.dist, scale = FALSE, N = 100000) {
  
  if(!inherits(tree, "phylo"))
    stop("Tree object must be of class 'phylo'", call. = FALSE)
  
  if(attr(tree, "order") != "pruningwise")
    stop("Edge matrix must be in 'pruningwise' order", call. = FALSE)
  
  if(any(dim(rate.mat) != dim(label.mat)))
    stop("Dimensions of 'rate.mat' and 'label.mat' don't match", call. = FALSE)
  
  if(any(dim(rate.mat) != length(root.dist)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'root.dist'", call. = FALSE)
  
  edge.mat = tree$edge
  edge.lengths = tree$edge.length
  if(scale == TRUE) rate.mat = scale.rate.mat(rate.mat, root.dist)
  num.edges = nrow(edge.mat)
  num.states = nrow(rate.mat)
  num.term.nodes = length(tree$tip.label)
  edge.moments = moments_ctmcjumps(edge.lengths, rate.mat, label.mat, num.edges, num.states)
  
  if(any(edge.set > num.edges | edge.set < 1))
    stop("Nonexistent edge(s) passed to 'edge.set'", call. = FALSE)
  
  postmean_moments_phylojumps(edge.mat, edge.lengths, rate.mat, edge.set, root.dist,
                              num.edges, num.states, num.term.nodes, edge.moments, N)
}






#' Joint prior moments of labeled substitution counts on a phylogenetic tree.
#' 
#' This function computes the prior means, variances, and covariance of labeled 
#' substitution counts over two prespecified sets of branches on a phylogeny 
#' using a simulation-free dynamic programming algorithm.
#' 
#' More information on phylogenetic tree objects of class \code{"phylo"} can be 
#' found at \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
#' 
#' @usage joint.prior.moments.phylojumps(tree, rate.mat, label.mat,
#'                                edge.set1 = NULL, edge.set2 = NULL,
#'                                subtree = NULL, root.dist, scale = FALSE)
#' 
#' @param tree A tree object of class \code{"phylo"} (\strong{ape} format).
#' @param rate.mat A reversible continuous-time Markov chain (CTMC) rate matrix.
#' @param label.mat A 0-1 matrix of the same size as \code{rate.mat}.  This 
#'   matrix labels the substitutions of interest.
#' @param edge.set1,edge.set2 Two integer vectors of branch indices specifying 
#'   the branch sets over which the joint prior moments will be calculated.
#' @param subtree An integer denoting a node index in the tree.  If provided, 
#'   then \code{edge.set1} will consist of the branches in the subtree beneath 
#'   the specified node (including the branch above \code{subtree}), while 
#'   \code{edge.set2} will contain the branches in the complementary supertree.
#' @param root.dist A numeric vector defining the root distribution of the 
#'   evolutionary process on the tree. If this process starts at stationarity 
#'   (as is commonly assumed), then the root distribution is equal to the CTMC 
#'   stationary distribution.
#' @param scale A logical indicating whether to scale the time dimension.  If 
#'   \code{TRUE}, then time is specified in terms of the expected number of CTMC
#'   substitutions per site.
#'   
#' @return A numeric vector holding the prior means, variances, and covariance 
#'   of labeled substitution counts.
#'   
#' @references OUR PAPER!!!
#' 
#' @seealso \code{\link{joint.post.moments.phylojumps}}, 
#'   \code{\link{joint.postmean.moments.phylojumps}}, 
#'   \code{\link{prior.moments.phylojumps}}
#'   
#' @export joint.prior.moments.phylojumps
joint.prior.moments.phylojumps = function(tree, rate.mat, label.mat, edge.set1 = NULL, edge.set2 = NULL,
                                          subtree = NULL, root.dist, scale = FALSE) {
  
  if(!inherits(tree, "phylo"))
    stop("Tree object must be of class 'phylo'", call. = FALSE)
  
  if(attr(tree, "order") != "pruningwise")
    stop("Edge matrix must be in 'pruningwise' order", call. = FALSE)
  
  if(any(dim(rate.mat) != dim(label.mat)))
    stop("Dimensions of 'rate.mat' and 'label.mat' don't match", call. = FALSE)
  
  if(any(dim(rate.mat) != length(root.dist)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'root.dist'", call. = FALSE)
  
  if(is.null(subtree) & (is.null(edge.set1) | is.null(edge.set2)))
    stop("Either 'subtree' or ('edge.set1','edge.set2') must be specified", call. = FALSE)
  
  edge.mat = tree$edge
  edge.lengths = tree$edge.length
  if(scale == TRUE) rate.mat = scale.rate.mat(rate.mat, root.dist)
  num.edges = nrow(edge.mat)
  num.states = nrow(rate.mat)
  num.term.nodes = length(tree$tip.label)
  edge.moments = moments_ctmcjumps(edge.lengths, rate.mat, label.mat, num.edges, num.states)
  
  ## constructing the two edge sets if "subtree" is specified
  if(!is.null(subtree)) {
    
    if(subtree > (2*num.term.nodes - 1) | subtree < 1)
      stop("Node passed to 'subtree' doesn't exist", call. = FALSE)
    
    if(subtree == edge.mat[num.edges,1])
      stop("Root node can't be an argument to 'subtree'", call. = FALSE)
    
    edge.set1 = find_subtree_edges(subtree, edge.mat)
    edge.set2 = setdiff(1:num.edges, edge.set1)
  }
  
  if(any(edge.set1 > num.edges | edge.set1 < 1))
    stop("Nonexistent edge(s) passed to 'edge.set1'", call. = FALSE)
  
  if(any(edge.set2 > num.edges | edge.set2 < 1))
    stop("Nonexistent edge(s) passed to 'edge.set2'", call. = FALSE)
  
  joint_prior_moments_phylojumps(edge.mat, edge.set1, edge.set2, root.dist,
                                 num.edges, num.states, num.term.nodes, edge.moments)
}






#' Joint posterior moments of labeled substitution counts on a phylogenetic 
#' tree.
#' 
#' This function computes the posterior means, variances, and covariance of 
#' labeled substitution counts over two prespecified sets of branches on a 
#' phylogeny using a simulation-free dynamic programming algorithm.
#' 
#' More information on phylogenetic tree objects of class \code{"phylo"} can be 
#' found at \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
#' 
#' @usage joint.post.moments.phylojumps(tree, rate.mat, label.mat,
#'                               edge.set1 = NULL, edge.set2 = NULL,
#'                               subtree = NULL, root.dist, scale = FALSE,
#'                               states = c("a","c","g","t"), seq.data)
#' 
#' @param tree A tree object of class \code{"phylo"} (\strong{ape} format).
#' @param rate.mat A reversible continuous-time Markov chain (CTMC) rate matrix.
#' @param label.mat A 0-1 matrix of the same size as \code{rate.mat}.  This 
#'   matrix labels the substitutions of interest.
#' @param edge.set1,edge.set2 Two integer vectors of branch indices specifying 
#'   the branch sets over which the joint posterior moments will be calculated.
#' @param subtree An integer denoting a node index in the tree.  If provided, 
#'   then \code{edge.set1} will consist of the branches in the subtree beneath 
#'   the specified node (including the branch above \code{subtree}), while 
#'   \code{edge.set2} will contain the branches in the complementary supertree.
#' @param root.dist A numeric vector defining the root distribution of the 
#'   evolutionary process on the tree. If this process starts at stationarity 
#'   (as is commonly assumed), then the root distribution is equal to the CTMC 
#'   stationary distribution.
#' @param scale A logical indicating whether to scale the time dimension.  If 
#'   \code{TRUE}, then time is specified in terms of the expected number of CTMC
#'   substitutions per site.
#' @param states A character or integer vector denoting the CTMC state space. 
#'   Defaults to the set of DNA nucleotides \code{c("a","c","g","t")}.
#' @param seq.data A character or integer matrix representing the observed 
#'   sequence alignment with rows (columns) that correspond to the sequences 
#'   (sites) under study.  All elements of \code{seq.data} must be contained in 
#'   the vector \code{states}.
#'   
#' @return A numeric vector holding the posterior means, variances, and 
#'   covariance of labeled substitution counts summed over all sites in 
#'   \code{seq.data}.
#'   
#' @references OUR PAPER!!!
#' 
#' @seealso \code{\link{joint.prior.moments.phylojumps}}, 
#'   \code{\link{joint.postmean.moments.phylojumps}}, 
#'   \code{\link{post.moments.phylojumps}}
#'   
#' @export joint.post.moments.phylojumps
joint.post.moments.phylojumps = function(tree, rate.mat, label.mat, edge.set1 = NULL, edge.set2 = NULL,
                                         subtree = NULL, root.dist, scale = FALSE, states = c("a","c","g","t"), seq.data) {
  
  if(!inherits(tree, "phylo"))
    stop("Tree object must be of class 'phylo'", call. = FALSE)
  
  if(attr(tree, "order") != "pruningwise")
    stop("Edge matrix must be in 'pruningwise' order", call. = FALSE)
  
  if(!all(seq.data %in% states))
    stop("Elements of 'seq.data' are not consistent with 'states'", call. = FALSE)
  
  if(any(dim(rate.mat) != dim(label.mat)))
    stop("Dimensions of 'rate.mat' and 'label.mat' don't match", call. = FALSE)
  
  if(any(dim(rate.mat) != length(root.dist)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'root.dist'", call. = FALSE)
  
  if(any(dim(rate.mat) != length(states)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'states'", call. = FALSE)
  
  if(is.null(subtree) & (is.null(edge.set1) | is.null(edge.set2)))
    stop("Either 'subtree' or ('edge.set1','edge.set2') must be specified", call. = FALSE)
  
  edge.mat = tree$edge
  edge.lengths = tree$edge.length
  if(scale == TRUE) rate.mat = scale.rate.mat(rate.mat, root.dist)
  num.edges = nrow(edge.mat)
  num.states = nrow(rate.mat)
  num.term.nodes = length(tree$tip.label)
  edge.moments = moments_ctmcjumps(edge.lengths, rate.mat, label.mat, num.edges, num.states)
  seq.data = matrix(match(seq.data, states) - 1, ncol = ncol(seq.data), byrow = FALSE)
  
  ## constructing the two edge sets if "subtree" is specified
  if(!is.null(subtree)) {
    
    if(subtree > (2*num.term.nodes - 1) | subtree < 1)
      stop("Node passed to 'subtree' doesn't exist", call. = FALSE)
    
    if(subtree == edge.mat[num.edges,1])
      stop("Root node can't be an argument to 'subtree'", call. = FALSE)
    
    edge.set1 = find_subtree_edges(subtree, edge.mat)
    edge.set2 = setdiff(1:num.edges, edge.set1)
  }
  
  if(any(edge.set1 > num.edges | edge.set1 < 1))
    stop("Nonexistent edge(s) passed to 'edge.set1'", call. = FALSE)
  
  if(any(edge.set2 > num.edges | edge.set2 < 1))
    stop("Nonexistent edge(s) passed to 'edge.set2'", call. = FALSE)
  
  if(nrow(seq.data) != num.term.nodes)
    stop("'seq.data' doesn't have the correct number of rows", call. = FALSE)
  
  joint_post_moments_phylojumps(edge.mat, edge.set1, edge.set2, root.dist, num.edges,
                                num.states, num.term.nodes, edge.moments, seq.data)
}







#' Joint moments of the posterior mean of labeled substitution counts on a 
#' phylogenetic tree.
#' 
#' This function calculates the means, variances, and covariance of the 
#' posterior mean of labeled substitution counts over two prespecified sets of 
#' branches on a phylogeny using a simulation-free dynamic programming algorithm
#' and Monte Carlo sampling.
#' 
#' The means obtained in this function are computed exactly, while the variances
#' and covariance are partially approximated using Monte Carlo simulations.
#' 
#' More information on phylogenetic tree objects of class \code{"phylo"} can be 
#' found at \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
#' 
#' @usage joint.postmean.moments.phylojumps(tree, rate.mat, label.mat,
#'                                   edge.set1 = NULL, edge.set2 = NULL, 
#'                                   subtree = NULL, root.dist,
#'                                   scale = FALSE, N = 100000)
#'                                   
#' @param tree A tree object of class \code{"phylo"} (\strong{ape} format).
#' @param rate.mat A reversible continuous-time Markov chain (CTMC) rate matrix.
#' @param label.mat A 0-1 matrix of the same size as \code{rate.mat}.  This 
#'   matrix labels the substitutions of interest.
#' @param edge.set1,edge.set2 Two integer vectors of branch indices specifying 
#'   the branch sets over which the joint posterior moments will be calculated.
#' @param subtree An integer denoting a node index in the tree.  If provided, 
#'   then \code{edge.set1} will consist of the branches in the subtree beneath 
#'   the specified node (including the branch above \code{subtree}), while 
#'   \code{edge.set2} will contain the branches in the complementary supertree.
#' @param root.dist A numeric vector defining the root distribution of the 
#'   evolutionary process on the tree. If this process starts at stationarity 
#'   (as is commonly assumed), then the root distribution is equal to the CTMC 
#'   stationary distribution.
#' @param scale A logical indicating whether to scale the time dimension.  If 
#'   \code{TRUE}, then time is specified in terms of the expected number of CTMC
#'   substitutions per site.
#' @param N An integer specifying the number of Monte Carlo samples to draw.
#' 
#' @return A numeric vector holding the means, variances, and covariance of the 
#' posterior mean of labeled substitution counts.
#' 
#' @references OUR PAPER!!!
#' 
#' @seealso \code{\link{joint.prior.moments.phylojumps}}, 
#'   \code{\link{joint.post.moments.phylojumps}}, 
#'   \code{\link{postmean.moments.phylojumps}}
#'   
#' @export joint.postmean.moments.phylojumps
joint.postmean.moments.phylojumps = function(tree, rate.mat, label.mat, edge.set1 = NULL, edge.set2 = NULL,
                                             subtree = NULL, root.dist, scale = FALSE, N = 100000) {
  
  if(!inherits(tree, "phylo"))
    stop("Tree object must be of class 'phylo'", call. = FALSE)
  
  if(attr(tree, "order") != "pruningwise")
    stop("Edge matrix must be in 'pruningwise' order", call. = FALSE)
  
  if(any(dim(rate.mat) != dim(label.mat)))
    stop("Dimensions of 'rate.mat' and 'label.mat' don't match", call. = FALSE)
  
  if(any(dim(rate.mat) != length(root.dist)))
    stop("Dimensions of 'rate.mat' aren't compatible with 'root.dist'", call. = FALSE)
  
  if(is.null(subtree) & (is.null(edge.set1) | is.null(edge.set2)))
    stop("Either 'subtree' or ('edge.set1','edge.set2') must be specified", call. = FALSE)
  
  edge.mat = tree$edge
  edge.lengths = tree$edge.length
  if(scale == TRUE) rate.mat = scale.rate.mat(rate.mat, root.dist)
  num.edges = nrow(edge.mat)
  num.states = nrow(rate.mat)
  num.term.nodes = length(tree$tip.label)
  edge.moments = moments_ctmcjumps(edge.lengths, rate.mat, label.mat, num.edges, num.states)
  
  ## constructing the two edge sets if "subtree" is specified
  if(!is.null(subtree)) {
    
    if(subtree > (2*num.term.nodes - 1) | subtree < 1)
      stop("Node passed to 'subtree' doesn't exist", call. = FALSE)
    
    if(subtree == edge.mat[num.edges,1])
      stop("Root node can't be an argument to 'subtree'", call. = FALSE)
    
    edge.set1 = find_subtree_edges(subtree, edge.mat)
    edge.set2 = setdiff(1:num.edges, edge.set1)
  }
  
  if(any(edge.set1 > num.edges | edge.set1 < 1))
    stop("Nonexistent edge(s) passed to 'edge.set1'", call. = FALSE)
  
  if(any(edge.set2 > num.edges | edge.set2 < 1))
    stop("Nonexistent edge(s) passed to 'edge.set2'", call. = FALSE)
  
  joint_postmean_moments_phylojumps(edge.mat, edge.lengths, rate.mat, edge.set1, edge.set2,
                                    root.dist, num.edges, num.states, num.term.nodes, edge.moments, N)
}


