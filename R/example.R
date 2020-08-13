source("utilAD.R")
source("simulation.R")
#Example 1
#Generate a time series of ER graphs of length 12,
#create graph anomaly at time point 6
n <- 100
glist <- list()
for (i in 1:5) {
  glist[[i]] <- sample_gnp(n,.1)
  }
  glist[[6]] <- sample_gnp(n,.9)
  glist[[7]] <- sample_gnp(n,.1)
  for (i in 8:12) {
    glist[[i]] <- sample_gnp(n,.1)
  }
# Do anomaly detection with OMNI, provide the quantitative control chart for GraphAD and VertexAD
result.OMNI<- qccAD(glist, l=4,d=1,dsvd=NULL,method="OMNI",
 diag.augment = TRUE, approx=FALSE, par=FALSE, numpar=2)
#print the number of deviation for GraphAD, only positive ones are meaningful
print(result.OMNI$GraphAD)

# Do anomaly detection with MASE
result.MASE<- qccAD(glist, l=4,d=2,dsvd=2,method="MASE",
               diag.augment = TRUE, approx=FALSE, par=FALSE, numpar=2)
#print the number of deviation for GraphAD, only positive ones are meaningful
print(result.MASE$GraphAD)


#Example 2
# Sample a time series of RDPG graph (length tmax > 17) with same 1-1 matched vertices unweighted
# hollow symmetric undirected graphs, the latent positions i.i.d uniform.
# Some vertices in 16-th and 17-th graphs are given perturbations so there exists anomalies at 16:17.
n <- 100 #number of vertices
nperturb <- 20 #number of perturbed vertices
cperturb <- .12 #number of perturbation, larger cperturb means more obvious anomalies.
rmin <- .2 # parameter for uniform[rmin, rmax].
rmax <- .8 # parameter for uniform[rmin, rmax].
tmax <- 22 # number of graphs must be greater than 17.
#Generate data or load the data you want
glist <- generate.tsg(n, nperturb, cperturb=NULL, rmin, rmax, tmax)$glist

#Do anomaly detection with OMNI in parallel
result.OMNI <- qccAD(glist, l=11,d=1,dsvd=NULL,method="OMNI",
                     diag.augment = TRUE, approx=FALSE, par=TRUE, numpar=2)
#print the number of deviation for GraphAD, only positive ones are meaningful
print(result.OMNI$GraphAD)

# Do anomaly detection with MASE in parallel
result.MASE<- qccAD(glist, l=11,d=1,dsvd=2,method="MASE",
                    diag.augment = TRUE, approx=FALSE, par=TRUE, numpar=2)
#print the number of deviation for GraphAD, only positive ones are meaningful
print(result.MASE$GraphAD)



#Example 3
#five of ER tsg with change point at t=6 and five at t=8.
n <- 100
dat <- matrix(0, 10, 8)
for (j in 1:5) {
  glist <- list()
  for (i in 1:5) {
    glist[[i]] <- sample_gnp(n,.1)
  }
  glist[[6]] <- sample_gnp(n,.9)
  for (i in 7:12) {
    glist[[i]] <- sample_gnp(n,.1)
  }
  # Do anomaly detection with OMNI, provide the quantitative control chart for GraphAD and VertexAD
  result.OMNI<- qccAD(glist, l=4,d=1,dsvd=NULL,method="OMNI",
                      diag.augment = TRUE, approx=FALSE, par=FALSE, numpar=2, plot.figure = FALSE)
  dat[j,] <- result.OMNI$GraphAD
}
for (j in 6:10) {
  glist <- list()
  for (i in 1:7) {
    glist[[i]] <- sample_gnp(n,.1)
  }
  glist[[8]] <- sample_gnp(n,.9)
  for (i in 9:12) {
    glist[[i]] <- sample_gnp(n,.1)
  }
  # Do anomaly detection with OMNI, provide the quantitative control chart for GraphAD and VertexAD
  result.OMNI<- qccAD(glist, l=4,d=1,dsvd=NULL,method="OMNI",
                      diag.augment = TRUE, approx=FALSE, par=FALSE, numpar=2, plot.figure = FALSE)
  dat[j,] <- result.OMNI$GraphAD
}
df <- data.frame(dat,class=factor(c(rep("t6",5),rep("t8",5))))
pca_res <- prcomp(dat, scale. = TRUE)
library(ggfortify)
autoplot(pca_res)
autoplot(pca_res, data = df,colour="class")
