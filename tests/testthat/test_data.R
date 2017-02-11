
set.seed(0)

## example phylogenetic tree
my.tree = ape::rtree(10, rooted = TRUE)
my.tree = reorder(my.tree, order = "pruningwise")
my.tree$node.label = as.character(1:(length(my.tree$tip.label) - 1))

## CTMC states
states = c("a", "b", "c")

## F81 rate matrix & root distribution
root.dist = c(0.2, 0.3, 0.5)
rate.mat = matrix(rep(root.dist, times = 3), nrow = 3, ncol = 3, byrow = TRUE)
diag(rate.mat) = 0
diag(rate.mat) = -apply(rate.mat, 1, sum)

## PHAST phylogenetic tree
phast.tree = rphast::tm(ape::write.tree(my.tree), subst.mod = "F81",
                        rate.matrix = rate.mat, backgd = root.dist, alphabet = "abc")

## label matrix (all states) & edge set (all edges)
label.mat = matrix(1, nrow = 3, ncol = 3) - diag(1, 3)
edge.set = 1:nrow(my.tree$edge)

## tip data
tip.data = matrix(sample(states, size = 10, replace = TRUE))
phast.data = rphast::msa(apply(tip.data, 1, paste, collapse = ""),
                         alphabet = "abc", names = my.tree$tip.label)

## initializing useful tree variables
edge.mat = my.tree$edge
edge.lengths = my.tree$edge.length

## (relative) sum of squares function
SS = function(x, y) sum(((x - y) / x)^2)
