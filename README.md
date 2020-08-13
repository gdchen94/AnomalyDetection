# AnomalyDetection
## Code and Results

All the codes are written in `R`, and here is the instruction for getting the results in the example.R 

### Demo

```
if (!require(RCurl)) install.packages("RCurl")
library(RCurl)

script <- getURL("https://raw.githubusercontent.com/youngser/dhatkhat/master/R/demo.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
```
and the result (each plot on its own window) should look like example.pdf   
(This may take around 1 minutes on a typical laptop.)


