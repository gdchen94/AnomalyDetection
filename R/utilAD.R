# library(igraph)
# library(dplyr)
# library(ggplot2)
# library(latex2exp)
# library(qcc)
# library(gtools)
#'
#' Given a decreasingly sorted vector, return the given number of elbows
#'
#' @param dat a input vector (e.g. a vector of standard deviations), or a input feature matrix
#' @param n the number of returned elbows
#' @param threshold either FALSE or a number. If threshold is a number, then all the elements in d that are not larger than the threshold will be ignored.
#' @param plot logical. When T, it depicts a scree plot with highlighted elbows
#' @param main a string of the plot title
#'
#' @return a vector of length \eqn{n}
#'
#' @references Zhu, Mu and Ghodsi, Ali (2006), Automatic dimensionality selection from
#' the scree plot via the use of profile likelihood, Computational
#' Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006.
#'
#' @export
#' @importFrom graphics plot points
#' @importFrom stats sd dnorm
#'
getElbows <- function(dat, n = 3, threshold = FALSE, plot = TRUE, main="", ...) {
  ## Given a decreasingly sorted vector, return the given number of elbows
  ##
  ## Args:
  ##   dat: a input vector (e.g. a vector of standard deviations), or a input feature matrix.
  ##   n: the number of returned elbows.
  ##   threshold: either FALSE or a number. If threshold is a number, then all
  ##   the elements in d that are not larger than the threshold will be ignored.
  ##   plot: logical. When T, it depicts a scree plot with highlighted elbows.
  ##
  ## Return:
  ##   q: a vector of length n.
  ##
  ## Reference:
  ##   Zhu, Mu and Ghodsi, Ali (2006), "Automatic dimensionality selection from
  ##   the scree plot via the use of profile likelihood", Computational
  ##   Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006.

  #  if (is.unsorted(-d))

  if (is.matrix(dat)) {
    d <- sort(apply(dat,2,sd), decreasing=TRUE)
  } else {
    d <- sort(dat,decreasing=TRUE)
  }

  if (!is.logical(threshold))
    d <- d[d > threshold]

  p <- length(d)
  if (p == 0)
    stop(paste("d must have elements that are larger than the threshold ",
               threshold), "!", sep="")

  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])              # = NaN when q = p
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
      (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
      sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }

  q <- which.max(lq)
  # print(n)
  # print(q)
  # print(p)
  if (n > 1 && q < (p-1)) {
    q <- c(q, q + getElbows(d[(q+1):p], n-1, plot=FALSE))
  }

  if (plot==TRUE) {
    if (is.matrix(dat)) {
      sdv <- d # apply(dat,2,sd)
      plot(sdv,type="b",xlab="dim",ylab="stdev",main=main,...)
      points(q,sdv[q],col=2,pch=19)
    } else {
      plot(dat, type="b",main=main,...)
      points(q,dat[q],col=2,pch=19)
    }
  }

  return(q)
}

#' Function to perform diagonal augmentation for graph adjacency matrix.
#' @param A a hollow adjacency matrix
#' @return a non-hollow adjacency matrix
#' @export
diagAug <- function(A) {
  diag(A) <- Matrix::rowSums(A) / (nrow(A)-1)
  return(A)
}
#' Calculate \eqn{l_2} distance between latent positions for vertices.
#' @param X,Y are matrices with n rows and d columns
#' @return A vector of size \eqn{n}, with element being \eqn{l_2} distance between rows of X and Y.

#' @export
pdistXY <- function(X,Y) {
  D <- Matrix::rowSums((as.matrix(X) - as.matrix(Y))^2)
  return(sqrt(D))
}

#'
#' Function to perform graph adjacency spectral embedding (ASE)
#' @param A adjacency matrix
#' @param d number of embedding dimension. If NULL, dimension is chosen automatically.
#' @param approx whether to find a few approximate singular values and corresponding singular vectors of a matrix using irlba package (TRUE/FALSE).
#' @param diagaug whether to do diagonal augmentation (TRUE/FALSE)
#' @param d.max maximum number of embedding dimensions to try when d is not provided. Default is round(log(nrow(A))).
#' @param elbow number of elbow selected in Zhu & Ghodsi method for the scree plot of each individual graph singular values. Default is \eqn{2}.
#' @references Zhu, Mu and Ghodsi, Ali (2006), Automatic dimensionality selection from
#' the scree plot via the use of profile likelihood, Computational
#' Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006.
#' @return A matrix with n rows and d columns containing the estimated latent positions
#' @export
#' @references Sussman, D.L., Tang, M., Fishkind, D.E., Priebe, C.E.  A
#' Consistent Adjacency Spectral Embedding for Stochastic Blockmodel Graphs,
#' \emph{Journal of the American Statistical Association}, Vol. 107(499), 2012
ase <- function(A, d=NULL, d.max=round(log(nrow(A))),diagaug=TRUE, approx=TRUE, elbow=2) {
  n <- nrow(A)
  #require(Matrix)
  # diagaug
  if (diagaug) {
    diag(A) <- Matrix::rowSums(A) / (n-1)
  }

  if (approx) {
    #require(irlba)
    if (is.null(d)){
      A.svd <- irlba::irlba(A, d.max, maxit=10000, tol=1e-10)
      result = getElbows(A.svd$d, plot = TRUE)
      if(length(result)==1){
        d = result
        print(paste0("do not have ",elbow," elbows, output first elbow instead"))
      }else{
        d = result[elbow]
      }
    }
    A.svd <- irlba::irlba(A, d, maxit=10000, tol=1e-10)
  } else {
    #require(rARPACK)
    if (is.null(d)){
      A.svd <- svd(A,d.max)
      result = getElbows(A.svd$d, plot = TRUE)
      if(length(result)==1){
        d = result
        print(paste0("do not have ",elbow," elbows, output first elbow instead"))
      }else{
        d = result[elbow]
      }
    }
    A.svd <- svd(A,d)
  }

  Xhat <- as.matrix(A.svd$u[,1:d]) %*% diag(sqrt(A.svd$d[1:d]), nrow=d, ncol=d)
  return(list(eval=A.svd$d, Xhat=Matrix(Xhat)))
}

#' Function to build OMNI matrix
#' @param Alist a list (length M) of n x n adjacency matrices or igraph objects
#' @param diagaug whether to do diagonal augmentation (TRUE/FALSE)
#' @return Omnibus matrix of Mn x Mn

#' @import igraph
buildOmni <- function(Alist, diagaug=FALSE) {
  # require(igraph)
  # require(Matrix)

  if (class(Alist[[1]]) == "igraph") {
    Alist <- lapply(Alist, igraph::get.adjacency)
  }


  if (diagaug) {
    Alist <- lapply(Alist, diagAug)
  }

  m <- length(Alist)
  nvec <- sapply(Alist, nrow)
  nsum <- c(0, cumsum(nvec))

  omni <- (bdiag(Alist))
  for (i in 1:(m-1)) {
    irng <- (nsum[i]+1):nsum[i+1]
    for (j in (i+1):m) {
      jrng <- (nsum[j]+1):nsum[j+1]
      omni[irng,jrng] <- (Alist[[i]] + Alist[[j]]) / 2
      omni[jrng,irng] <- omni[irng,jrng]
    }
  }
  return(omni)
}
#' Function to build OMNI matrix with two matrix
#' @param Alist a list (length \eqn{2}) of n x n adjacency matrices or igraph objects
#' @param diagaug whether to do diagonal augmentation (TRUE/FALSE)
#' @return Omnibus matrix of 2n x 2n

fast2buildOmni <- function(Alist, diagaug=FALSE, attrweight=NULL ) {
  # require(igraph)
  # require(Matrix)

  if (class(Alist[[1]]) == "igraph") {
    Alist <- lapply(Alist, function(x) igraph::get.adjacency(x, attr = attrweight ) )
  }


  if (diagaug) {
    Alist <- lapply(Alist, diagAug)
  }

  omni <- cbind(Alist[[1]], (Alist[[1]]+Alist[[2]])/2 )
  omni <- rbind(omni, cbind((Alist[[1]]+Alist[[2]])/2, Alist[[2]]))
  return(omni)
}
#' Function to calculate \eqn{l_2} distance between adjacent latent position estmate for graphs and vertices with OMNI
#' @param n number of vertices in each graph
#' @param Z a matrix of size 2n x d as the latent estimate for the omnibus matrix contructed by adjacent graphs
#' @return A list containing a numeric value tnorm, with the \eqn{l_2} norm of
#' latent position estmate difference for adjacent graphs, and a vector pdist
#'  with \eqn{l_2} distance between latent position estmates for each vertex (size n).

fast2doOmni <- function(n, Z) {
  norm <- sapply(Z, function(x) norm(x$Xhat[1:n,] - x$Xhat[-(1:n),], "2"))
  pdist <- sapply(Z, function(x) pdistXY(x$Xhat[1:n,],x$Xhat[-(1:n),]))
  return(list(tnorm=norm, pdist=pdist))
}
#' Function to obtain latent position estmates for MASE (Multiple Adjancency Spectral Embedding) with individual ASE estimates for corresponding graphs
#' @param A a list (length M-1) of n x n adjacency matrices
#' @param latpos.list a list (length M-1) of individual ASE estimates for the corresponding adjacency matrices
#' @param dsvd dimension for joint embedding. If NULL then dimension is chosen automatically as the second elbow selected in Zhu & Ghodsi method for the scree plot of the singular values of the concatenated spectral embeddings of invidual ASE estimates.
#' @return A list with length (M-1) consists of matrices with size n x dsvd each.

mase.latent <- function(A, latpos.list, dsvd=NULL){
  m <- length(A)
  A <- lapply(A, function(x) (x) )
  jrdpg <- mase(A, latpos.list,dsvd, show.scree.results = TRUE)
  d2 <- dim(jrdpg$R[[1]])[1]
  B.svd <- lapply(jrdpg$R, function (x) svd(x) )
  Xhat <- lapply(B.svd, function(x) (jrdpg$V) %*% as.matrix(x$u) %*% diag(sqrt(x$d), nrow=d2, ncol=d2) %*% diag(sign(x$d)) %*% t(as.matrix(x$u) )  )

  return(list(Xhat=Xhat))
}

#' Function to calculate \eqn{l_2} distance between adjacent latent position estmate for graphs and vertices with MASE
#' @param glist a list (length M) of n x n adjacency matrices or igraph objects
#' @param latpos.list a list (length M) of ASE estimates for each graphs
#' @param nmase number of graphs to do joint embedding of MASE. It can only be 2 or M.
#' @param dsvd dsvd is number of dimension to do joint svd for MASE. If NULL then dimension is chosen automatically as the second elbow selected in Zhu & Ghodsi method for the scree plot of the singular values of the concatenated spectral embeddings of invidual ASE estimates.
#' @return A list containing a vector tnorm of length M-1, with the
#' latent position estmate difference for graphs, and a matrix pdist
#'  with latent position estmate difference for vertices (size n x M-1 ).
doMase <- function(glist, latpos.list, nmase=2, dsvd=NULL,attrweight=NULL)
{
  n <- igraph::vcount(glist[[1]])
  tmax <- length(glist)

  if (nmase == 2) {
    Xhat <- lapply(1:(tmax-1), function(x) mase.latent(lapply(glist[x:(x+1)], function(y) igraph::get.adjacency(y,attr=attrweight) ), latpos.list[x:(x+1)] , dsvd))
    norm <- sapply(1:(tmax-1), function(x) norm((Xhat[[x]]$Xhat)[[1]] - (Xhat[[x]]$Xhat)[[2]], "2"))
    pdist <- sapply(1:(tmax-1), function(x) pdistXY((Xhat[[x]]$Xhat)[[1]] , (Xhat[[x]]$Xhat)[[2]]))
    dsvd <- sapply(Xhat, function(x) x$d)
  } else {
    adj <- lapply(glist, function(y) igraph::get.adjacency(y, attr = attrweight) )
    Xhat <- mase.latent(adj, latpos.list,dsvd)
    norm <- sapply(1:(tmax-1), function(x) norm(Xhat$Xhat[[x]]-Xhat$Xhat[[x+1]], "2"))
    pdist <- sapply(1:(tmax-1), function(x) pdistXY(Xhat$Xhat[[x]], Xhat$Xhat[[x+1]]))
    dsvd <- Xhat$d
  }
  return(list(tnorm=norm, pdist=pdist,d=dsvd))
}
#' Function to perform joint embedding part of MASE
#'
#' @param Adj_list a list of adjacency matrices with the same size n x n
#' @param d number of joint embedding dimensions. If NA, dimension is chosen automatically
#' @param latpos.list Individual ASE estimate for the adjacent graphs
#' @param elbow_mase number of elbow selected in Zhu & Ghodsi method for the scree plot of the singular values of the concatenated spectral embeddings of MASE.
#' @param show.scree.results when TRUE, the histogram of the estimated d for each graph, and the scree plot of the singular values of  the graph is shown if d is not specified.
#'
#' @return A list containing a matrix V of size n x d, with the
#' estimated invariant subspace, and a list R with the individual score parameters for each graph (size d x d each).
#'
#' @references
#' @export

mase <- function(Adj_list, latpos.list, dsvd = NULL,
                 elbow_mase = 2,
                 show.scree.results = FALSE) {
  V_all  <- Reduce(cbind, latpos.list)
  #require(rARPACK)
  jointsvd <- svd(V_all)
  if(is.null(dsvd)) {
    # if(show.scree.results) {
    #    hist(sapply(latpos.list, ncol), main = "Estimated d for each graph")
    # }

    result = getElbows(jointsvd$d, plot = show.scree.results)#[elbow_mase]
    if(length(result)==1){
      dsvd = result
    }else{
      dsvd = result[elbow_mase]
    }
  }
  V = jointsvd$u[, 1:dsvd, drop = FALSE]
  R <- project_networks(Adj_list, V)
  return(list(V = V, sigma = jointsvd$d, R = R))
}

#'
#' Function to estimated the score matrices of a list of graphs given the common invariant subspace V
#'
#' @param Adj_list list of adjacency matrices, of size n x n
#' @param V common invariant subspace. A matrix of size n x d.
#' @return A list containing the score matrices
#'
project_networks <- function(Adj_list, V) {
  # require(Matrix)
  lapply(Adj_list, function(A) Matrix::crossprod(Matrix::crossprod(A, V), V))
}
#' @param plot.LCL a Boolean variable to decide whether to show the anomalies lower than lower limits (LCL \eqn{\mu^{t}-3\sigma^t}) .
#' Return a list of cases beyond limits and violating runs
#' @references see R package qcc
shewhart.rules <- function(object, limits = object$limits, run.length = qcc.options("run.length"),plot.LCL=FALSE)
{
  # Return a list of cases beyond limits and violating runs
  bl <- beyond.limits(object, limits = limits, plot.LCL=plot.LCL )
  vr <- qcc::violating.runs(object, run.length = run.length)
  list(beyond.limits = bl, violating.runs = vr)
}
#' @param plot.LCL a Boolean variable to decide whether to show the anomalies lower than lower limits (LCL \eqn{\mu^{t}-3\sigma^t}) .
#' Return cases beyond limits
#' @references see R package qcc
beyond.limits <- function(object, limits = object$limits, plot.LCL=FALSE)
{
  statistics <- c(object$statistics, object$newstats)
  lcl <- limits[,1]
  ucl <- limits[,2]
  index.above.ucl <- seq(along = statistics)[statistics > ucl]
  if(plot.LCL){
    index.below.lcl <- seq(along = statistics)[statistics < lcl]
    return(c(index.below.lcl,index.above.ucl))
  }else{
    return(c(index.above.ucl))
  }

}
#' Function to plot Shewhart chart for GraphAD
#' @param x a qcc object
#' @param l length of previous graphs to estimate moving averages and moving standard deviation.
#' @param title a string of the plot title
#' @param plot.LCL a Boolean variable to decide whether to show the anomalies lower than lower limits (LCL \eqn{\mu^{t}-3\sigma^t}) .
#' @return A control chart ggplot for GraphAD
#'
plot.qcc <- function(x, l=l, title, plot.LCL=FALSE){
  object <- x  # Argh.  Really want to use 'object' anyway
  if ((missing(object)) | (!inherits(object, "qcc")))
    stop("an object of class `qcc' is required")
  # collect info from object
  type <- object$type
  std.dev <- object$std.dev
  data.name <- object$data.name
  center <- object$center
  stats <- object$statistics
  limits <- object$limits
  lcl <- limits[,1]
  ucl <- limits[,2]
  newstats <- object$newstats
  newdata.name <- object$newdata.name
  violations <- object$violations
  statistics <- c(stats, newstats)
  indices <- 1:length(statistics)

  #library(ggplot2)
  l.length <- length(indices)
  tmax <- (l.length+l-1)
  tvec <- 1:(l.length+l)
  m2 <- names(gtools::running(tvec, width=2))
  minx <- m2[l:(l.length+l-1)]
  vioinx <- rep(0,length(indices))
  runinx <- rep(0,length(indices))
  vioinx[violations$beyond.limits] <- 1
  runinx[violations$violating.runs] <- 1
  df <- data.frame(time=indices, y= statistics,lcl=lcl,ucl=ucl,center=center,lim=vioinx,run=runinx)
  if(plot.LCL){
    p <- ggplot2::ggplot(df,aes(x=time, y=y))+
      geom_step(aes(x=time, y=ucl), linetype="dashed")+
      geom_step(aes(x=time, y=lcl), linetype="dashed")+
      geom_step(aes(x=time, y=center), linetype="solid") +
      geom_point(data = df, alpha=1, color="black")+ylab(TeX("$y^{(t)}$"))+
      geom_line(aes(x=time, y=y))+
      geom_point(data = df %>% filter(run==1), alpha=1, color="yellow")+
      geom_point(data = df %>% filter(lim==1), alpha=1, color="red") +
      scale_x_discrete(name ="time points",
                       limits=minx)+theme_bw()+
      annotate("text", label = "UCL",
               x =  tmax-l+1.1, y = rev(ucl)[1]+.3)+
      annotate("text", label = "CL",
               x =  tmax-l+1.1, y = rev(center)[1]+.1 )+
      annotate("text", label = "LCL",
               x =  tmax-l+1.1, y = rev(lcl)[1]-.3 )+theme_classic(base_size = 18)+
      theme(axis.text.x = element_text(angle=45),legend.position = "none",plot.title = element_text(hjust = 0.5,size=20, face='bold'),plot.subtitle = element_text(hjust = 0.5)) +ggtitle(title)
  }else{
    p <- ggplot2::ggplot(df,aes(x=time, y=y))+
      geom_step(aes(x=time, y=ucl), linetype="dashed")+
      geom_step(aes(x=time, y=center), linetype="solid") +
      geom_point(data = df, alpha=1, color="black")+ylab(TeX("$y^{(t)}$"))+
      geom_line(aes(x=time, y=y))+
      geom_point(data = df %>% filter(run==1), alpha=1, color="yellow")+
      geom_point(data = df %>% filter(lim==1), alpha=1, color="red") +
      scale_x_discrete(name ="time points",
                       limits=minx)+theme_bw()+
      annotate("text", label = "UCL",
               x =  tmax-l+1.1, y = rev(ucl)[1]+.3)+
      annotate("text", label = "CL",
               x =  tmax-l+1.1, y = rev(center)[1]+.1 )+theme_classic(base_size = 18)+
      theme(axis.text.x = element_text(angle=45),legend.position = "none",plot.title = element_text(hjust = 0.5,size=20, face='bold'),plot.subtitle = element_text(hjust = 0.5)) +ggtitle(title)
  }
  return(p)
  invisible()
}

#' Function to plot Shewhart chart for VertexAD
#' @param x a list of qcc object
#' @param plot.LCL a Boolean variable to decide whether to show the anomalies lower than lower limits (LCL \eqn{\mu_i^{t}-3\sigma_i^t}) .
#' @param l length of previous graphs to estimate moving averages and moving standard deviation.
#' @param title a string of the plot title
#' @return A control chart ggplot for VertexAD
plot.qcc.vertex <- function(x, l=l, title, plot.LCL=FALSE){
  l.length <- length(x)
  df <- c() #c(2,6)
  for (i in seq(l.length)) {
    object <- x[[i]]  # Argh.  Really want to use 'object' anyway
    if ((missing(object)) | (!inherits(object, "qcc")))
      stop("an object of class `qcc' is required")

    # collect info from object
    type <- object$type
    std.dev <- object$std.dev
    data.name <- object$data.name
    center <- object$center
    stats <- object$statistics
    limits <- object$limits
    lcl <- limits[,1]
    ucl <- limits[,2]
    newstats <- object$newstats
    newdata.name <- object$newdata.name
    violations <- object$violations
    statistics <- c(stats, newstats)
    indices <- 1:length(statistics)

    #library(ggplot2)
    tmax <- l.length+1
    tvec <- 1:(l.length+l)
    m2 <- names(gtools::running(tvec, width=2))
    minx <- m2[l:(l.length+l-1)]
    vioinx <- rep(0,length(indices))
    runinx <- rep(0,length(indices))
    vioinx[violations$beyond.limits] <- 1
    runinx[violations$violating.runs] <- 1
    idf <- data.frame(vertex=indices, y= statistics, timepoints=factor(minx[i], levels = minx),lcl=lcl,ucl=ucl,center=center,lim=vioinx,run=runinx)
    df <- rbind(df,idf)

  }
  if(plot.LCL){

    p <- ggplot2::ggplot(df,aes(x=vertex, y=y))+
      geom_point(data = df, alpha=.1, color="black",shape=20,size=.5)+
      geom_point(data = df %>% filter(lim==1), alpha=.3, color="red",shape=20,size=1)+
      geom_hline(aes(yintercept=center), linetype="solid")+facet_wrap(~timepoints,nrow = 1,switch = "x", scales = "fixed") +
      geom_hline(aes(yintercept=ucl), linetype="dashed")+
      geom_hline(aes(yintercept=lcl), linetype="dashed")+
      labs(x="vertex")+theme_bw()+ylab(TeX("$y_{i}^{(t)}$"))+
      theme(axis.text.x = element_text(angle=45),legend.position = "none",plot.title = element_text(hjust = 0.5,size=10, face='bold'),plot.subtitle = element_text(hjust = 0.5)) +ggtitle(title)
  }else{

    p <- ggplot2::ggplot(df,aes(x=vertex, y=y))+
      geom_point(data = df, alpha=.1, color="black",shape=20,size=.5)+
      geom_point(data = df %>% filter(lim==1), alpha=.3, color="red",shape=20,size=1)+
      geom_hline(aes(yintercept=center), linetype="solid")+facet_wrap(~timepoints,nrow = 1,switch = "x", scales = "fixed") +
      geom_hline(aes(yintercept=ucl), linetype="dashed")+
      labs(x="vertex")+theme_bw()+ylab(TeX("$y_{i}^{(t)}$"))+
      theme(axis.text.x = element_text(angle=45),legend.position = "none",plot.title = element_text(hjust = 0.5,size=10, face='bold'),plot.subtitle = element_text(hjust = 0.5)) +ggtitle(title)
  }
  return(p)
  invisible()
}
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

#-----------------------------------------------------------------------------#
#                                                                             #
#                     QUALITY CONTROL CHARTS IN R                             #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Luca Scrucca                                                   #
#              Department of Statistics                                       #
#              University of Perugia, ITALY                                   #
#              luca@stat.unipg.it                                             #
#                                                                             #
#-----------------------------------------------------------------------------#

#
#'  Main function to create a 'qcc' object
#' @param plot.LCL a Boolean variable to decide whether to show the anomalies lower than lower limits (LCL \eqn{\mu^{t}-3\sigma^t}) .
#' @references see R package qcc
qcc <- function(data, type = c("xbar", "R", "S", "xbar.one", "p", "np", "c", "u", "g"), sizes, center, std.dev, limits, data.name, labels, newdata, newsizes, newdata.name, newlabels, nsigmas = 3, confidence.level, rules = shewhart.rules, plot = TRUE, plot.LCL=FALSE,...)
{
  call <- match.call()

  if (missing(data))
    stop("'data' argument is not specified")

  if(identical(type, eval(formals(qcc)$type)))
  { type <- as.character(type)[1]
  warning("chart 'type' not specified, assuming \"", type, "\"",
          immediate. = TRUE) }
  if(!exists(paste("stats.", type, sep = ""), mode="function") |
     !exists(paste("sd.", type, sep = ""), mode="function") |
     !exists(paste("limits.", type, sep = ""), mode="function"))
    stop(paste("invalid", type, "control chart. See help(qcc) "))

  if (missing(data.name))
    data.name <- deparse(substitute(data))
  data <- data.matrix(data)
  if (missing(sizes))
  { if (any(type==c("p", "np", "u")))
    stop(paste("sample 'sizes' must be given for a", type, "Chart"))
    else
      sizes <- apply(data, 1, function(x) sum(!is.na(x)))  }
  else
  { if (length(sizes)==1)
    sizes <- rep(sizes, nrow(data))
  else if (length(sizes) != nrow(data))
    stop("sizes length doesn't match with data") }

  if (missing(labels))
  { if (is.null(rownames(data))) labels <- 1:nrow(data)
  else                         labels <- rownames(data) }

  stats <- paste("stats.", type, sep = "")
  if (!exists(stats, mode="function"))
    stop(paste("function", stats, "is not defined"))
  stats <- do.call(stats, list(data, sizes))
  statistics <- stats$statistics
  if (missing(center)) center <- stats$center

  sd <- paste("sd.", type, sep = "")
  if (!exists(sd, mode="function"))
    stop(paste("function", sd, "is not defined!"))
  missing.std.dev <- missing(std.dev)
  if (missing.std.dev)
  { std.dev <- NULL
  std.dev <- switch(type,
                    "xbar" = { if(any(sizes > 25)) "RMSDF"
                      else                "UWAVE-R" },
                    "xbar.one" = "MR",
                    "R" = "UWAVE-R",
                    "S" = "UWAVE-SD",
                    NULL)
  std.dev <- do.call(sd, list(data, sizes, std.dev)) }
  else
  { if (is.character(std.dev))
  { std.dev <- do.call(sd, list(data, sizes, std.dev)) }
    else
    { if (!is.numeric(std.dev))
      stop("if provided the argument 'std.dev' must be a method available or a numerical value. See help(qcc).")  }
  }

  names(statistics) <-  rownames(data) <-  labels
  names(dimnames(data)) <- list("Group", "Samples")

  object <- list(call = call, type = type,
                 data.name = data.name, data = data,
                 statistics = statistics, sizes = sizes,
                 center = center, std.dev = std.dev)
  # check for new data provided and update object
  if (!missing(newdata))
  {   if (missing(newdata.name))
  {newdata.name <- deparse(substitute(newdata))}
    newdata <- data.matrix(newdata)
    if (missing(newsizes))
    { if (any(type==c("p", "np", "u")))
      stop(paste("sample sizes must be given for a", type, "Chart"))
      else
        newsizes <- apply(newdata, 1, function(x) sum(!is.na(x))) }
    else
    { if (length(newsizes)==1)
      newsizes <- rep(newsizes, nrow(newdata))
    else if (length(newsizes) != nrow(newdata))
      stop("newsizes length doesn't match with newdata") }
    stats <- paste("stats.", type, sep = "")
    if (!exists(stats, mode="function"))
      stop(paste("function", stats, "is not defined"))
    newstats <- do.call(stats, list(newdata, newsizes))$statistics
    if (missing(newlabels))
    { if (is.null(rownames(newdata)))
    { start <- length(statistics)
    newlabels <- seq(start+1, start+length(newstats)) }
      else
      { newlabels <- rownames(newdata) }
    }
    names(newstats) <- newlabels
    object$newstats <- newstats
    object$newdata  <- newdata
    object$newsizes <- newsizes
    object$newdata.name <- newdata.name
    statistics <- c(statistics, newstats)
    sizes <- c(sizes, newsizes)
  }

  conf <- nsigmas
  if (!missing(confidence.level))
    conf <- confidence.level
  if (conf >= 1)
  { object$nsigmas <- conf }
  else
    if (conf > 0 & conf < 1)
    { object$confidence.level <- conf }

  # get control limits
  if (missing(limits))
  { limits <- paste("limits.", type, sep = "")
  if (!exists(limits, mode="function"))
    stop(paste("function", limits, "is not defined"))
  limits <- do.call(limits, list(center = center, std.dev = std.dev,
                                 sizes = sizes, conf = conf))
  }
  else
  { if (!missing.std.dev)
    warning("'std.dev' is not used when limits is given")
    if (!is.numeric(limits))
      stop("'limits' must be a vector of length 2 or a 2-columns matrix")
    limits <- matrix(limits, ncol = 2)
    dimnames(limits) <- list(rep("",nrow(limits)), c("LCL ", "UCL"))
  }

  lcl <- limits[,1]
  ucl <- limits[,2]
  object$limits <- limits
  #number of deviation for each data point, i.e., deviation=(x-mu)/sigma
  object$deviation <- (object$data - object$center)/object$std.dev
  if (is.function(rules)) violations <- rules(object,plot.LCL=plot.LCL)
  else                    violations <- NULL
  object$violations <- violations

  class(object) <- "qcc"
  if(plot) plot(object, ...)
  return(object)
}

#' Function to perform anomaly detection for time series of graphs
#' @param glist a list of undirected simple graphs (simple graphs are graphs which do not contain self-loop and multiple edges)
#' in igraph format with same number of vertices with vertices are 1-1 matched.
#' Graphs can be weighted or binary. (Say the length of list to be tmax)
#' @param method a character variable  to be chosen among c("OMNI","MASE").
#' The code will first do OMNIbus embedding (OMNI) or Multiple Adjacency Spectrally Ebedding (MASE) with two adjcaency matrices(can be weighted or not) of all input adjacent graphs sequentially.
#' Then use latent positions to calculate test statistics \eqn{y^{t}=||X^{t}- X^{t+1}||} using operator norm.
#' Then for \eqn{t=l,...,tmax-1}, we calculate the moving means \eqn{\mu^t} and moving standard deviations \eqn{\sigma^t} at \eqn{t} by \eqn{y^{t-l+1},...,y^{t-1}}.
#' So only tmax-l time points are ploted as first l graphs have been used as estimating moving means and standard deviations.
#' @param d a fixed integer of dimension to perform OMNI and individual ASE for MASE.  If d is NULL, then dimension is chosen automatically.
#' @param l an integer of the number of graphs in time window in estimating the moving mean and moving standard deviation. l must be less than number of graphs and be greater than 3.
#' @param approx a Boolean variable to decide whether to use irlba package to find a few approximate singular values and corresponding singular vectors of a matrix. Default is TRUE.
#' @param diag.augment a Boolean variable to decide whether to do diagonal augmentation when performing adjacency spectral embedding. Default is TRUE.
#' @param par a Boolean variable to decide whether to run in parallel. Default is FALSE.
#' @param numpar an integer number to decide number of clusters for parallel implmentation. Default is 2.
#' @param dsvd An integer number of dimension only used in joint embedding for MASE. If NULL, then dimension is chosen automatically as second elbow selected in Zhu & Ghodsi method for the scree plot of the singular values of the concatenated spectral embeddings of MASE.
#' @param elbow number of elbow in Zhu & Ghodsi method for the scree plot of each individual graph singular values for MASE or of each omnibus matrix singular values for OMNI.
#' @param plot.figure a Boolean variable to decide whether plot control chart. Default is TRUE.
#' @param plot.LCL a Boolean variable to decide whether to show the anomalies lower than lower limits (LCL \eqn{\mu^{t}-3\sigma^t}) .
#' @references Zhu, Mu and Ghodsi, Ali (2006), Automatic dimensionality selection from
#' the scree plot via the use of profile likelihood, Computational
#' Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006.
#' @references Levin, K., Athreya, A., Tang, M., Lyzinski, V. and Priebe, C.E., 2017, November. A central limit theorem for an omnibus embedding of multiple random dot product graphs. In 2017 IEEE International Conference on Data Mining Workshops (ICDMW) (pp. 964-967). IEEE.
#' @references Arroyo, J., Athreya, A., Cape, J., Chen, G., Priebe, C.E. and Vogelstein, J.T., 2019. Inference for multiple heterogeneous networks with a common invariant subspace. arXiv preprint arXiv:1906.10026.
#' @return A list containing a vector GraphAD of length tmax-l which consists of control charts deviations, with the
#' a list VertexAD (length tmax-l) with vectors of anomalous vertices indices for each graph.
#' @export
#' @examples
#' glist <- list()
#' for (i in 1:5) {
#'   glist[[i]] <- sample_gnp(100,.1)
#'   }
#'   glist[[6]] <- sample_gnp(100,.9)
#'   glist[[7]] <- sample_gnp(100,.1)
#'   for (i in 8:12) {
#'     glist[[i]] <- sample_gnp(100,.1)
#'     }
#'   result<- qccAD(glist, l=4,d=1,dsvd=1,method="OMNI",
#'   diag.augment = TRUE,approx=FALSE, par=FALSE, numpar=2)
#'   print(result.OMNI$GraphAD) #print the number of deviation for GraphAD, only positive ones are meaningful
#'
#' @examples
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

#' #Do anomaly detection with OMNI in parallel
#' result.OMNI <- qccAD(glist, l=11,d=1,dsvd=NULL,method="OMNI",
#'                      diag.augment = TRUE, approx=FALSE, par=TRUE, numpar=2)
#' #print the number of deviation for GraphAD, only positive ones are meaningful
#' print(result.OMNI$GraphAD)
#'
#' # Do anomaly detection with MASE in parallel
#' result.MASE<- qccAD(glist, l=11,d=1,dsvd=2,method="MASE",
#'                     diag.augment = TRUE, approx=FALSE, par=TRUE, numpar=2)
#' #print the number of deviation for GraphAD, only positive ones are meaningful
#' print(result.MASE$GraphAD)
#'

qccAD <- function(glist, method="OMNI", diag.augment = TRUE, l=3,d=NULL, dsvd=d, approx=TRUE, par=FALSE, numpar=2, elbow=2, plot.figure=TRUE,plot.LCL=FALSE){
  tmax <- length(glist)
  n <- igraph::vcount(glist[[1]])
  tvec <- 1:tmax
  m2 <- names(gtools::running(tvec, width=2))
  #check if all the graphs are binary
  if(Reduce("+",lapply(glist, function(x) igraph::is.weighted(x)))){
    attrweight <- "weight"
  }else{
    attrweight <- NULL
  }
  if(l>(tmax-1)||l<3||tmax<4){
    print(paste0("l has to an integer between 3 and ",tmax-1, " and you need to at least have 4 graphs"))
  }
  #require(foreach)
  if (method=="OMNI"){
    #Initialize moving averages and moving standard deviation for GraphAD
    mean2omni <- rep(1, (tmax-1)-(l-1))
    std2omni <- rep(1,(tmax-1)-(l-1))
    #Initialize moving averages and moving standard deviation for VertexAD
    mean2vomni <- rep(1, (tmax-1)-(l-1))
    std2vomni <- rep(1,(tmax-1)-(l-1))
    #Do omnibus embedding
    if(par){
      #require(doParallel)
      cl <- parallel::makeCluster(numpar)
      doParallel::registerDoParallel(cl)
      parallel::clusterEvalQ(cl, library("AnomalyDetection"))#source("utilAD.R"))
      parallel::clusterExport(cl = cl, varlist = list("ase", "fast2buildOmni", "getElbows", "glist",
                                            "elbow", "d","approx", "diag.augment","attrweight"), envir = environment())
      allomni <- foreach::foreach(i =1:(tmax-1), .combine='comb', .multicombine=TRUE,.packages = "igraph",
                         .init=list(list(), list())) %dopar% {
                           O <- fast2buildOmni( glist[i:(i+1)]  , diagaug=diag.augment, attrweight=attrweight)
                           a <- ase(O, d=d, diagaug=FALSE, approx=approx,elbow=elbow)
                           list(a, a$d)
                         }
      parallel::stopCluster(cl)
    } else {
      foreach::registerDoSEQ()
      allomni <- foreach::foreach(i =1:(tmax-1), .combine='comb', .multicombine=TRUE,.packages = "igraph",
                         .init=list(list(), list())) %dopar% {
                           O <- fast2buildOmni( glist[i:(i+1)]  , diagaug=diag.augment, attrweight=attrweight)
                           Z <- ase(O, d=d, diagaug=FALSE, approx=approx,elbow=elbow)
                           list(Z, Z$d)
                         }
    }
    #Z is omnibus embedding estimates for latent positions
    Z = allomni[[1]]
    #Calculate l_2 distance between adjacent latent position estmate for graphs and vertices with OMNI
    out2omni <- fast2doOmni(n, Z=Z)
    for (w in 1:(tmax-1-(l-1))) {
      #Iterative calulating moving averages and moving deviations
      omni2v <- t(out2omni$pdist[,w:(w+l-2)])
      mean2vomni[w] <- mean(omni2v)
      std2vomni[w] <- sd.xbar(omni2v, rep(n,l-1), "UWAVE-SD")
      omni2 <- matrix(out2omni$tnorm[w:(w+l-2)], l-1, 1)
      mean2omni[w] <- mean(omni2)
      std2omni[w] <- sd.xbar.one(omni2, rep(1,l-1), "MR")
    }
    #Do vertexAD
    df <- out2omni$pdist
    omni2v <- matrix(df, tmax-1, n, byrow = TRUE)
    #c1o2v is a list of qcc (quantative control chart) object
    c1o2v <- list()
    for (i in l:(tmax-1)) {
      c1o2v[[i-l+1]] <- qcc(omni2v[i,],type = "xbar.one",center = mean2vomni[i-l+1],std.dev = std2vomni[i-l+1], nsigmas=3, plot=FALSE, plot.LCL = plot.LCL)
    }
    #Do GraphAD
    df <- out2omni$tnorm
    omni2 <- matrix(df, tmax-1, 1)
    #c1o2 is a qcc object
    c1o2 <- qcc(omni2[l:(tmax-1)],type = "xbar.one",nsigmas=3, center = mean2omni, std.dev = std2omni, plot=FALSE, plot.LCL = plot.LCL)
    if(plot.figure){
      print(plot.qcc(c1o2,title="Control Chart OMNI", l=l,plot.LCL = plot.LCL ))
      print(plot.qcc.vertex(c1o2v, l=l,title="Control Chart OMNI",plot.LCL = plot.LCL ))
    }
    #Let c1o2v be a list of anomalous vertices indices across time points
    c1o2v <- list()
    for (i in l:(tmax-1)) {
      c1o2v[[i-l+1]] <- qcc(omni2v[i,],type = "xbar.one",center = mean2vomni[i-l+1],std.dev = std2vomni[i-l+1], nsigmas=3, plot=FALSE, plot.LCL = plot.LCL )$violations$beyond.limits
    }
    result <- list(GraphAD=c1o2$deviation, VertexAD=c1o2v)
    return(result)
  }else if(method=="MASE"){
    #Initialize moving averages and moving standard deviation for GraphAD
    mean2 <- rep(1,(tmax-1)-(l-1))
    std2 <- rep(1,(tmax-1)-(l-1))
    #Initialize moving averages and moving standard deviation for VertexAD
    mean2v <- rep(1,(tmax-1)-(l-1))
    std2v <- rep(1,(tmax-1)-(l-1))
    #Do individual ASE for MASE
    if(par){
      #require(doParallel)
      cl <- parallel::makeCluster(numpar)
      doParallel::registerDoParallel(cl)
      parallel::clusterEvalQ(cl, library("AnomalyDetection"))#source("utilAD.R"))
      parallel::clusterExport(cl = cl, varlist = list("ase", "getElbows", "glist",
                                            "elbow", "d","approx", "diag.augment","attrweight"), envir = environment())
      allase <- foreach::foreach(i =1:tmax, .combine='comb', .multicombine=TRUE,.packages = "igraph",
                     .init=list(list(), list())) %dopar% {
                       latpos <- ase(igraph::get.adjacency(glist[[i]], attr=attrweight), d=d, diagaug=diag.augment, approx=approx,elbow=elbow)$Xhat
                       list(latpos, dim(latpos)[2])
                     }
      parallel::stopCluster(cl)
    } else {
      foreach::registerDoSEQ()
      allase <- foreach::foreach(i =1:tmax, .combine='comb', .multicombine=TRUE,.packages = "igraph",
                     .init=list(list(), list())) %dopar% {
                       latpos <- ase(igraph::get.adjacency(glist[[i]], attr=attrweight), d=d, diagaug=diag.augment, approx=approx,elbow=elbow)$Xhat
                       list(latpos, dim(latpos)[2])
                     }
    }
    latpos.list <- allase[[1]]
    #Do joint embedding in MASE and calculate l_2 distance between adjacent latent position estmate for graphs and vertices with MASE
    out2 <- doMase(glist,latpos.list, 2, dsvd, attrweight)
    for (w in 1:(tmax-1-(l-1))) {
      #Iterative calulating moving averages and moving deviations
      mase2v <- t(matrix(out2$pdist[,w:(w+l-2)],n,l-1))
      mean2v[w] <- mean(mase2v)
      std2v[w] <- sd.xbar(mase2v, rep(n,l-1), "UWAVE-SD")
      mase2 <- matrix(out2$tnorm[w:(w+l-2)], l-1, 1)
      mean2[w] <- mean(mase2)
      std2[w] <- sd.xbar.one(mase2, rep(1,l-1), "MR")
    }
    #Do VertexAD (vertex anomaly detection)
    df <- out2$pdist
    mase2v <- matrix(df, tmax-1, n, byrow = TRUE)
    c1m2v <- list()
    for (i in l:(tmax-1)) {
      c1m2v[[i-l+1]] <- qcc(mase2v[i,],type = "xbar.one",center = mean2v[i-l+1],std.dev = std2v[i-l+1], nsigmas=3, plot=FALSE, plot.LCL = plot.LCL )
    }
    #Do GraphAD (graph anomaly detection)
    df <- out2$tnorm
    mase2 <- matrix(df, tmax-1, 1)
    c1m2 <- qcc(mase2[l:(tmax-1)],type = "xbar.one",center = mean2 ,std.dev = std2, nsigmas=3, plot=FALSE, plot.LCL = plot.LCL )
    if(plot.figure){
      print(plot.qcc(c1m2,title="Control Chart MASE", l=l,plot.LCL = plot.LCL ))
      print(plot.qcc.vertex(c1m2v, l=l,title="Control Chart MASE",plot.LCL = plot.LCL ))
    }
    c1m2v <- list()
    for (i in l:(tmax-1)) {
      c1m2v[[i-l+1]] <- qcc(mase2v[i,],type = "xbar.one",center = mean2v[i-l+1],std.dev = std2v[i-l+1], nsigmas=3, plot=FALSE, plot.LCL = plot.LCL)$violations$beyond.limits
    }
    result <- list(GraphAD=c1m2$deviation, VertexAD=c1m2v)
    return(result)
  }
}
