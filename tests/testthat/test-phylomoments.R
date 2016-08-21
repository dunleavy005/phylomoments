context("All tests")
source("test_data.R")


test_that("moments.ctmcjumps computes the correct moments", {
  true.moments = moments.ctmcjumps(1, rate.mat, label.mat)
  approx.moments = unlist(sapply(states, function(init.state) {
    raw.samp = replicate(10000, ctmc.sim(1, rate.mat, states, init.state))
    agg.samp = apply(raw.samp, 2, function(col) {
      c("njumps" = length(col$state) - 2,
        "end.state" = col$state[length(col$state)])
    })
    moments = tapply(as.numeric(agg.samp["njumps",]), agg.samp["end.state",], function(subset) {
      c("zeroth" = length(subset) / 10000,
        "first" = sum(subset) / 10000,
        "second" = sum(subset^2 - subset) / 10000)
    })
  }))
  
  true.approx.dist = function(string) {
    SS(true.moments[[string]][,,1], matrix(approx.moments[names(approx.moments) == string],
                                           nrow = 3, ncol = 3, byrow = TRUE))
  }
  
  expect_true(true.approx.dist("zeroth") <= 1e-2)
  expect_true(true.approx.dist("first") <= 1e-2)
  expect_true(true.approx.dist("second") <= 1e-1)
})

test_that("(prior/post).moments.phylojumps compute the correct moments", {
  exp.moments = rphast::phyloP.sph(phast.tree, phast.data, fit.model = FALSE)
  my.prior.moments = prior.moments.phylojumps(my.tree, rate.mat, label.mat, edge.set, root.dist)
  my.post.moments = post.moments.phylojumps(my.tree, rate.mat, label.mat, edge.set,
                                            root.dist, states = states, seq.data = tip.data)
  
  expect_equal(my.prior.moments[["mean"]], exp.moments[["prior.mean"]], tol = 1e-5)
  ## PHAST does not output the correct prior variance
  expect_equal(my.post.moments[["mean"]], exp.moments[["posterior.mean"]], tol = 1e-5)
  expect_equal(my.post.moments[["var"]], exp.moments[["posterior.var"]], tol = 1e-5)
  expect_equal_to_reference(my.prior.moments[["mean"]], "prior_mean.rds")
  expect_equal_to_reference(my.prior.moments[["var"]], "prior_var.rds")
  expect_equal_to_reference(my.post.moments[["mean"]], "post_mean.rds")
  expect_equal_to_reference(my.post.moments[["var"]], "post_var.rds")
  
  edge.moments = moments.ctmcjumps(edge.lengths, rate.mat, label.mat)
  int.states = int.states.sim(my.tree, rate.mat, root.dist,
                              states = states, tip.data = tip.data, N = 20000)
  
  total.moments = sapply(1:ncol(int.states), function(i) {
    node.states = match(int.states[,i], states)
    total.var.mat = ((edge.moments[["second"]] + edge.moments[["first"]]) / edge.moments[["zeroth"]]) -
      (edge.moments[["first"]] / edge.moments[["zeroth"]])^2
    total.exp.mat = edge.moments[["first"]] / edge.moments[["zeroth"]]
    total.var = 0
    total.exp = 0
    
    for (j in edge.set) {
      parent = edge.mat[j,1]
      child = edge.mat[j,2]
      total.var = total.var + total.var.mat[node.states[parent], node.states[child], j]
      total.exp = total.exp + total.exp.mat[node.states[parent], node.states[child], j]
    }
    
    c("total.var" = total.var, "total.exp" = total.exp)
  })
  
  expect_true(SS(my.post.moments[["mean"]], mean(total.moments["total.exp",])) <= 1e-5)
  expect_true(SS(my.post.moments[["var"]], mean(total.moments["total.var",]) + var(total.moments["total.exp",])) <= 1e-5)
  
  njumps = phylojumps.sim(my.tree, rate.mat, root.dist, states = states, seq.data = tip.data, N = 50000)
  expect_true(SS(my.post.moments[["mean"]], mean(njumps)) <= 1e-5)
  expect_true(SS(my.post.moments[["var"]], var(njumps)) <= 1e-4)
})

test_that("postmean.moments.phylojumps computes the correct moments", {
  my.postmean.moments = postmean.moments.phylojumps(my.tree, rate.mat, label.mat, edge.set, root.dist, N = 50000)
  seq.sim = tips.sim(my.tree, rate.mat, root.dist, states = states, N = 20000)
  post.means = sapply(1:ncol(seq.sim), function(i) {
    post.moments.phylojumps(my.tree, rate.mat, label.mat, edge.set, root.dist,
                            states = states, seq.data = matrix(seq.sim[,i]))[["mean"]]
  })
  
  expect_true(SS(my.postmean.moments[["mean"]], mean(post.means)) <= 1e-5)
  expect_true(SS(my.postmean.moments[["var"]], var(post.means)) <= 1e-2)
  expect_equal_to_reference(my.postmean.moments[["mean"]], "prior_mean.rds")
  expect_equal_to_reference(my.postmean.moments[["var"]], "postmean_var.rds")
})

test_that("joint.(prior/post/postmean).moments.phylojumps compute the correct moments", {
  sub.node = 14
  sub.edges = edge.set[apply(edge.mat, 1, function(row) any(row == sub.node))]
  sup.edges = setdiff(edge.set, sub.edges)
  exp.moments = rphast::phyloP.sph(phast.tree, phast.data, subtree = "t10-t6", fit.model = FALSE)
  my.joint.prior.moments = joint.prior.moments.phylojumps(my.tree, rate.mat, label.mat, sub.edges,
                                                          sup.edges, root.dist = root.dist)
  my.joint.post.moments = joint.post.moments.phylojumps(my.tree, rate.mat, label.mat, sub.edges, sup.edges,
                                                        root.dist = root.dist, states = states, seq.data = tip.data)
  my.joint.postmean.moments = joint.postmean.moments.phylojumps(my.tree, rate.mat, label.mat, sub.edges,
                                                                sup.edges, root.dist = root.dist, N = 50000)
  
  expect_identical(my.joint.prior.moments, joint.prior.moments.phylojumps(my.tree, rate.mat, label.mat,
                                                                          subtree = sub.node, root.dist = root.dist))
  expect_identical(my.joint.post.moments, joint.post.moments.phylojumps(my.tree, rate.mat, label.mat, subtree = sub.node,
                                                                        root.dist = root.dist, states = states, seq.data = tip.data))
  expect_true(SS(my.joint.postmean.moments, joint.postmean.moments.phylojumps(my.tree, rate.mat, label.mat, subtree = sub.node,
                                                                              root.dist = root.dist, N = 50000)) <= 1e-2)
  expect_equal(my.joint.prior.moments[["mean1"]], exp.moments[["prior.subtree.mean"]], tol = 1e-5)
  expect_equal(my.joint.prior.moments[["mean2"]], exp.moments[["prior.supertree.mean"]], tol = 1e-5)
  expect_equal(my.joint.prior.moments[["var1"]], exp.moments[["prior.subtree.var"]], tol = 1e-5)
  expect_equal(my.joint.prior.moments[["var2"]], exp.moments[["prior.supertree.var"]], tol = 1e-5)
  expect_equal(my.joint.post.moments[["mean1"]], exp.moments[["post.subtree.mean"]], tol = 1e-5)
  expect_equal(my.joint.post.moments[["mean2"]], exp.moments[["post.supertree.mean"]], tol = 1e-5)
  expect_equal(my.joint.post.moments[["var1"]], exp.moments[["post.subtree.var"]], tol = 1e-5)
  expect_equal(my.joint.post.moments[["var2"]], exp.moments[["post.supertree.var"]], tol = 1e-5)
  expect_equal(my.joint.postmean.moments[["mean1"]], exp.moments[["prior.subtree.mean"]], tol = 1e-5)
  expect_equal(my.joint.postmean.moments[["mean2"]], exp.moments[["prior.supertree.mean"]], tol = 1e-5)
  
  seq.sim = tips.sim(my.tree, rate.mat, root.dist, states = states, N = 20000)
  post.means = sapply(1:ncol(seq.sim), function(i) {
    joint.post.moments.phylojumps(my.tree, rate.mat, label.mat, sub.edges, sup.edges, root.dist = root.dist,
                                  states = states, seq.data = matrix(seq.sim[,i]))[c("mean1", "mean2")]
  })
  
  expect_true(SS(my.joint.postmean.moments[["mean1"]], mean(post.means["mean1",])) <= 1e-5)
  expect_true(SS(my.joint.postmean.moments[["mean2"]], mean(post.means["mean2",])) <= 1e-5)
  expect_true(SS(my.joint.postmean.moments[["var1"]], var(post.means["mean1",])) <= 1e-2)
  expect_true(SS(my.joint.postmean.moments[["var2"]], var(post.means["mean2",])) <= 1e-2)
  expect_true(SS(my.joint.postmean.moments[["cov"]], cov(post.means["mean1",], post.means["mean2",])) <= 1e-2)
  
  exp.prior.var = sum(my.joint.prior.moments[c("var1", "var2")]) + 2 * my.joint.prior.moments[["cov"]]
  exp.post.var = sum(my.joint.post.moments[c("var1", "var2")]) + 2 * my.joint.post.moments[["cov"]]
  exp.postmean.var = sum(my.joint.postmean.moments[c("var1", "var2")]) + 2 * my.joint.postmean.moments[["cov"]]
  
  expect_equal_to_reference(exp.prior.var, "prior_var.rds")
  expect_equal_to_reference(exp.post.var, "post_var.rds")
  expect_equal_to_reference(exp.postmean.var, "postmean_var.rds", tol = 1e-3)
})

