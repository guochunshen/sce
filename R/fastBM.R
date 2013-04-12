#'
#'Fast Brownian simulation
#'
#'
#'@param tree  is a phylogenetic tree in "phylo" format.
#'@param a a value for ancestral state at the root node.
#'@param mu	an optional value for the mean of random normal changes along branches of the tree - can be used to simulate a trend if mu!=0.
#'@param sig2	instantaneous variance of the BM process.
#'@param bounds	a vector with the lower and upper bounds (respectively) for bounded Brownian simulation - by default simulation is unbounded.
#'@param internal	logical value indicating whether or not to return states for internal nodes.
#'@param nsim	number of simulations.
#'
#'
#'@details
#'it is an only copy of fastBM from phytools package. please see details in that package.
#'
#'




fastBM<-function (tree, a = 0, mu = 0, sig2 = 1, bounds = c(-Inf, Inf), 
          internal = FALSE, nsim = 1) 
{
  if (class(tree) != "phylo") 
    stop("tree object must be of class 'phylo.'")
  if (bounds[2] < bounds[1]) {
    warning("bounds[2] must be > bounds[1]. Simulating without bounds.")
    bounds <- c(-Inf, Inf)
  }
  if (bounds[1] == -Inf && bounds[2] == Inf) 
    no.bounds = TRUE
  else no.bounds = FALSE
  if (a < bounds[1] || a > bounds[2]) {
    warning("a must be bounds[1]<a<bounds[2]. Setting a to midpoint of bounds.")
    a <- bounds[1] + (bounds[2] - bounds[1])/2
  }
  if (sig2 < 0) {
    warning("sig2 must be > 0.  Setting sig2 to 1.0.")
    sig2 = 1
  }
  reflect <- function(yy, bounds) {
    while (yy < bounds[1] || yy > bounds[2]) {
      if (yy < bounds[1]) 
        yy <- 2 * bounds[1] - yy
      if (yy > bounds[2]) 
        yy <- 2 * bounds[2] - yy
    }
    return(yy)
  }
  n <- length(tree$tip)
  x <- matrix(data = rnorm(n = length(tree$edge.length) * nsim, 
                           mean = rep(mu * tree$edge.length, nsim), sd = rep(sqrt(sig2 * 
                                                                                    tree$edge.length), nsim)), length(tree$edge.length), 
              nsim)
  y <- array(0, dim = c(nrow(tree$edge), ncol(tree$edge), nsim))
  for (i in 1:nrow(x)) {
    if (tree$edge[i, 1] == (n + 1)) 
      y[i, 1, ] <- a
    else y[i, 1, ] <- y[match(tree$edge[i, 1], tree$edge[, 
                                                         2]), 2, ]
    y[i, 2, ] <- y[i, 1, ] + x[i, ]
    if (!no.bounds) 
      y[i, 2, ] <- apply(as.matrix(y[i, 2, ]), 1, function(yy) reflect(yy, 
                                                                       bounds))
  }
  rm(x)
  x <- matrix(data = rbind(y[1, 1, ], as.matrix(y[, 2, ])), 
              length(tree$edge.length) + 1, nsim)
  rownames(x) <- c(n + 1, tree$edge[, 2])
  x <- as.matrix(x[as.character(1:(n + tree$Nnode)), ])
  rownames(x)[1:n] <- tree$tip.label
  if (internal == TRUE) 
    return(x[1:nrow(x), ])
  else return(x[1:length(tree$tip.label), ])
}