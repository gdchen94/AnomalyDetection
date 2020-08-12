#' Sample RDPG graph with latent position
#' @param X latent position matrix of \eqn{d} columns and \eqn{n} rows.
#' @return A un-directed hollow symmeric unweighted graph generated from bernoulli \eqn{EA=P=XX^{T}}
rdpg.sample <- function(X) {
  P <- X %*% t(X)
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  U[col(U) > row(U)] <- runif(n*(n-1)/2)
  U <- (U + t(U))
  diag(U) <- runif(n)
  A <- (U < P) + 0 ;
  diag(A) <- 0
  return(graph.adjacency(A,"undirected"))
}
#' Sample a time series of RDPG graph (length tmax > 17) with same 1-1 matched vertices unweighted
#' hollow symmetric undirected graphs, the latent positions i.i.d uniform.
#' Some vertices in 16-th and 17-th graphs are given perturbations so there exists anomalies at 16:17.
#' @param n number of vertices
#' @param nperturb number of perturbed vertices
#' @param cperturb number of perturbation. Larger cperturb means more obvious anomalies.
#' @param rmin,rmax parameter for uniform[rmin, rmax].
#' @param tmax number of graphs must be greater than 17.
#' @return A list containing a vector tnorm of length tmax-1, with the
#' latent position difference for graphs, and a matrix pdist
#'  with latent position estmate difference for vertices (size n x tmax-1 ),
#'  and a list of length tmax of undirectected hollow symmeric unweighted graphs
#' @examples
#' #' @examples
#' # Sample a time series of RDPG graph (length tmax > 17) with same 1-1 matched vertices unweighted
#' # hollow symmetric undirected graphs, the latent positions i.i.d uniform.
#' # Some vertices in 16-th and 17-th graphs are given perturbations so there exists anomalies at 16:17.
#' n <- 100 #number of vertices
#' nperturb <- 20 #number of perturbed vertices
#' cperturb <- .12 #number of perturbation, larger cperturb means more obvious anomalies.
#' rmin <- .2 # parameter for uniform[rmin, rmax].
#' rmax <- .8 # parameter for uniform[rmin, rmax].
#' tmax <- 22 # number of graphs must be greater than 17.
#' #Generate data or load the data you want
#' glist <- generate.tsg(n, nperturb, cperturb=NULL, rmin, rmax, tmax)$glist

generate.tsg <- function(n, nperturb, cperturb=NULL, rmin, rmax, tmax){
  X1 <- runif(n, rmin, rmax)
  pert1 <- ifelse(is.null(cperturb), min(X1), cperturb)
  pert2 <- ifelse(is.null(cperturb), 1-max(X1), cperturb)
  Xlist <- rep(list(X1), tmax)
  Xlist[[10+6]] <- Xlist[[1]] + c(rep(c(pert2, -pert1), each=nperturb/2), rep(0, n-nperturb))
  Xlist[[10+7]] <- Xlist[[1]] + c(rep(c(-pert1, pert2), each=nperturb/2), rep(0, n-nperturb))

  norm <- sapply(1:(tmax-1), function(x) norm(Xlist[[x]]-Xlist[[x+1]], "2"))
  pdist <- sapply(1:(tmax-1), function(x) pdistXY(Xlist[[x]], Xlist[[x+1]]))
  glist <- lapply(Xlist, function(x) rdpg.sample(matrix(x,ncol=1)))
  return(list(tnorm=norm, pdist=pdist, glist=glist))
}
