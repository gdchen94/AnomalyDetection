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
Alternatively, you can download the

The function qccAD is called to detect one or more statistically
significant anomalies in the input time series of graphs. The documentation of the
function qccAD, which can be seen by using the following command,
details the input arguments and the output of the function qccAD.

```
?qccAD
```

## A simple example

To get started, the user is recommended to use the example dataset which comes
with the packages. Execute the following commands:

```
data(raw_data)
res = AnomalyDetectionTs(raw_data, max_anoms=0.02, direction='both', plot=TRUE)
res$plot
```


All the codes are written in `R`, and here is the instruction for getting the results in the example.pdf using the R file example.R


### Demo
First, put the 
```
require(knitr)
stitch("example.R")
```
and the result should look like example.pdf   
(This may take around 1 minutes on a typical laptop.)

![Fig 1](https://github.com/twitter/AnomalyDetection/blob/master/figs/Fig1.png)

From the plot, we observe that the input time series experiences both positive 
and negative anomalies. Furthermore, many of the anomalies in the time series
are local anomalies within the bounds of the time series’ seasonality (hence,
cannot be detected using the traditional approaches). The anomalies detected
using the proposed technique are annotated on the plot. In case the timestamps 
for the plot above were not available, anomaly detection could then carried 
out using the AnomalyDetectionVec function; specifically, one can use the 
following command:

```
AnomalyDetectionVec(raw_data[,2], max_anoms=0.02, period=1440, direction='both', only_last=FALSE, plot=TRUE)
```

Often, anomaly detection is carried out on a periodic basis. For instance, at
times, one may be interested in determining whether there was any anomaly
yesterday. To this end, we support a flag only_last whereby one can subset the
anomalies that occurred during the last day or last hour. Execute the following 
command:

```
res = AnomalyDetectionTs(raw_data, max_anoms=0.02, direction='both', only_last=”day”, plot=TRUE)
res$plot
```

![Fig 2](https://github.com/twitter/AnomalyDetection/blob/master/figs/Fig2.png)

From the plot, we observe that only the anomalies that occurred during the last
day have been annotated. Further, the prior six days are included to expose the
seasonal nature of the time series but are put in the background as the window
of prime interest is the last day.

Anomaly detection for long duration time series can be carried out by setting
the longterm argument to T. 

## Copyright and License
Copyright 2015 Twitter, Inc and other contributors

Licensed under the GPLv3


