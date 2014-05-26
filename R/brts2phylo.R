brts2phylo <- function(times,root=FALSE,tip.label=NULL)
{
# This code is taken from the package TESS by Sebastian Hoehna, where the function is called tess.create.phylo
# It takes a set of branching times and adds a random topology.
  times = sort(times)
  n <- as.integer(length(times))+1
  if ( root ) {
    n <- n-1
  }
  nbr <- 2*n - 2 

  # create the data types for edges and edge-lengths
  edge <- matrix(NA, nbr, 2)
  edge.length <- numeric(nbr)
  
  h <- numeric(2*n - 1) # initialized with 0's
  pool <- 1:n
  # VERY VERY IMPORTANT: the root MUST have index n+1 !!!
  nextnode <- 2L*n - 1L
  if ( n > 1) {
    for (i in 1:(n - 1)) {
      # sample two nodes that have no parent yet
      y <- sample(pool, size = 2)
      # compute the edge indices (we just order the edges from 1 to 2n-2)
      ind <- (i - 1)*2 + 1:2
      # set the source node of the new edges (i.e. the new internal node)
      edge[ind, 1] <- nextnode
      # set the destination of the new edges (i.e. the two sampled nodes)
      edge[ind, 2] <- y
      # compute the edge length from the difference between the node heights (child <-> parent)
      edge.length[ind] <- times[i] - h[y]
      # store the node height of this new internal node
      # we cannot use times because then we would get into trouble with the indices and we would need to check for tip nodes ...
      h[nextnode] <- times[i]
      # reset the pool of available nodes to merge
      pool <- c(pool[! pool %in% y], nextnode)
      # increase the node index counter
      nextnode <- nextnode - 1L
    }
  }

  phy <- list(edge = edge, edge.length = edge.length)
  if (is.null(tip.label))
    tip.label <- paste("t", 1:n, sep = "")
  phy$tip.label <- sample(tip.label)
  phy$Nnode <- n - 1L

  if ( root ) {
    phy$root.edge <- times[n] - times[n-1]
    phy$root <- times[n] - times[n-1]
  }

  class(phy) <- "phylo"
  
  phy <- reorder(phy)
  ## to avoid crossings when converting with as.hclust:
  phy$edge[phy$edge[, 2] <= n, 2] <- 1:n

  return(phy)
}

