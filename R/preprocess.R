# preprocess
#' Run pass-to-rank on a weighted graph.
#'
#' It extracts (non-zero) edge weight vector \eqn{W} from a graph and replaces it with \eqn{2*R / (|E|+1)} where \eqn{R} is the rank of \eqn{W} and \eqn{|E|} is the number of edges. This does 'no-op' for an unweighted graph.
#'
#' @param g a graph in \code{igraph} format or an n x 2 edge list or an n x n adjacency matrix
#'
#' @author Youngser Park <youngser@jhu.edu>
#' @export
#' @import igraph
ptr <- function(g)
{
  if (class(g) != "igraph") {
    if (!is.matrix(g)) stop("the input has to be either an igraph object or a matrix!")
    else {
      if (ncol(g)==2) g <- igraph::graph_from_edgelist(g)
      else if (nrow(g)==ncol(g)) g <- igraph::graph_from_adjacency_matrix(g, weighted = TRUE)
      else stop("the input matrix is not a graph format!")
    }
  }

  if (igraph::is.weighted(g)) {
    W <- E(g)$weight
  } else { # no-op!
    W <- rep(1,igraph::ecount(g))
  }

  E(g)$weight <- rank(W)*2 / (igraph::ecount(g)+1)
  return(g)
}
#' find largest connected component in a graph
#' It extracts (non-zero) largest connected subgraph .
#'
#' @param graph a graph in \code{igraph} format
#'
#' @author Guodong Chen <gchen35@jhu.edu>
#' @export
#' @import igraph
giant.component <- function(graph, ...) {
  #require(igraph)
  cl <- igraph::clusters(graph, ...)
  igraph::induced.subgraph(graph, which(cl$membership == which.max(cl$csize)))
}
#' remove edges which has zero weights for all graphs (if any) and find jointly largest connected component in graphs. Finally it removes all self-loops,
#' It extracts (non-zero) igraph list \eqn{gip} and removes all edges with zero edge weights and return a list of jointly largest connected component in graphs without self-loops .
#'
#' @param gip a list of graphs in \code{igraph} format
#'
#' @author Guodong Chen <gchen35@jhu.edu>
#' @export
#' @import igraph
#' @return A list of graphs in igraph format
jlcc <- function(gip){
  l.length <- length(gip)
  for (i in 1:l.length) {
    gip[[i]] <- igraph::delete.edges(gip[[i]], which(E(gip[[i]])$weight==0))
  }
  #find joint largest connected component
  df1 <- igraph::as_data_frame(gip[[1]])[,1:2]
  df2 <- igraph::as_data_frame(gip[[2]])[,1:2]
  g <- igraph::graph.intersection(gip[[1]],gip[[2]],keep.all.vertices = TRUE)
  if(l.length>2){
    for (i in 3:l.length) {
      g <- igraph::graph.intersection(g,gip[[i]],keep.all.vertices = TRUE)
    }
  }
  lcc <- giant.component(g)
  glist <- gip
  gip <- list()
  for (i in 1:l.length) {
    gip[[i]] <- igraph::induced_subgraph(glist[[i]], V(lcc)$name);
    gip[[i]] <- igraph::permute.vertices(gip[[i]], match(V(gip[[i]])$name, V(gip[[1]])$name));
  }
  for (i in 1:l.length) {
    gip[[i]] <- igraph::simplify(gip[[i]])
  }
  return(gip)
}



#' create a planted clique
#' It creates planted clique for a specific graph for a list of graphs.
#'
#' @param gip a list of graphs in \code{igraph} format.
#' @param p is the index of graph to be inserted a planted clique
#' @param art.anomaly.v is the vertex index in igraph.vs format to be planted clique.
#' @author Guodong Chen <gchen35@jhu.edu>
#' @export
#' @import igraph
#' @return A list containing a planted clique size as size of art.anomaly.v at p-th graph

#
pltclique <- function(gip, p, art.anomaly.v){
  middle.max.inx <- which(V(gip[[1]])%in% art.anomaly.v )
  art.gip <- gip
  adj <- as.matrix(igraph::get.adjacency(art.gip[[p]], attr = attrweigth))
  adj[middle.max.inx,middle.max.inx] <- (adj[middle.max.inx,middle.max.inx] + 1)
  art.gip[[p]] <- igraph::graph.adjacency(adj,mode = "undirected", weighted = TRUE,diag = FALSE)
  art.anomaly.v <- V(art.gip[[p]])[middle.max.inx]
  glist <- list()
  l.length <- length(gip)
  for (i in 1:l.length) {
    glist[[i]] <- ptr(art.gip[[i]])
  }
  return(glist)
}


#get degree distribution for graphs
#'
#' It extracts (non-zero) degree change \eqn{deg.change} matrix n by m-1 from a list of graphs.
#'
#' @param gip a list of graphs in \code{igraph} format.
#' @author Guodong Chen <gchen35@jhu.edu>
#' @export
#' @import igraph
#' @return A matrix of size n x t-1, with each element to be degree changes
getdegchange <- function(gip){
  l.length <- length(gip)
  num.edge <- matrix(0,length(V(gip[[1]])),l.length )
  for (i in 1:l.length) {
    num.edge[,i] <- igraph::degree(gip[[i]])
  }
  deg.change <- matrix(0,length(V(gip[[1]])),l.length-1)
  for (i in 1:(l.length-1)) {
    deg.change[,i] <- num.edge[,i+1]- num.edge[,i]
  }
  return(deg.change)
  # plot(density(deg.change))
  # plot(density(abs(deg.change)))
}


#' get weighted degree change for a list of weighted graphs
#' It extracts (non-zero) weighted degree change matrix n by m-1 \eqn{deg.change} from a list of graphs.
#'
#' @param gip a list of graphs in \code{igraph} format.
#' @author Guodong Chen <gchen35@jhu.edu>
#' @export
#' @import igraph
#' @return A matrix of size n x t-1, with each element to be weight degree changes
getweightchange <- function(gip){
  l.length <- length(gip)
  num.edge.weight <- matrix(0,length(V(gip[[1]])), l.length )
  for (i in 1:l.length) {
    num.edge.weight[,i] <- igraph::graph.strength(gip[[i]],loops = FALSE)
  }
  deg.change.weight <- matrix(0,length(V(gip[[1]])),l.length-1)
  for (i in 1:(l.length-1)) {
    deg.change.weight[,i] <- num.edge.weight[,i+1]- num.edge.weight[,i]
  }
  return(deg.change.weight)
}
