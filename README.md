# AnomalyDetection R package
AnomalyDetection is a R package to detect anomalies in time series of graphs. 
Specifically, it can be applied in two related, complementary inference tasks: 
the detection of anomalous graphs within a time series, and the detection
of temporally anomalous vertices. 

## How the package works

The underlying algorithm approaches these tasks via the adaptation of statistically
principled methods for joint graph inference, specifically multiple adjacency
spectral embedding (MASE) and omnibus embedding (OMNI). It also builds upon the quantative
control chart (https://en.wikipedia.org/wiki/Control_chart) for detecting anomalies. 
Note that it can be used to detect anomalies both in overall graphs as well as individual anomalies. 
This is achieved by employing joint graph embedding for adjacent graphs sequentially, and track the
change of adjacent latent positions difference via quantative control charts. Con

The package provides some preprocessing tools to deal, for example, different number of vertices in
time series of graphs.
The user can specify the length of moving window of interest (such as seven days, twenty-four hours), 
enable/disable singular value decompostion approximation etc.

## How to get started

Install the R package using the following commands on the R console:

```
install.packages("devtools")
devtools::install_github("gdchen94/AnomalyDetection")
library(AnomalyDetection)
```
Alternatively, you can download the file AnomalyDetection_0.1.0.tar.gz
and install it using either (within R)
```
install.packages("AnomalyDetection_0.1.0.tar.gz", repos = NULL, type ="source")
library(AnomalyDetection)
```
or in terminal (say you download it in Desktop)
```
R CMD INSTALL Desktop/AnomalyDetection_0.1.0.tar.gz
```

The function qccAD is called to detect one or more statistically
significant anomalies in the input time series of graphs. The documentation of the
function qccAD, which can be seen by using the following command,
details the input arguments and the output of the function qccAD.

```
?qccAD
```

## A simple example

To get started, the user is recommended to use the example.R file to
get the results in the example.pdf


### Demo
```
require(knitr)
stitch("example.R")
```
and the result should look like example.pdf   
(This may take around 1 minutes on a typical laptop.)

![Fig 1](https://github.com/gdchen94/AnomalyDetection/blob/master/example.pdf)




