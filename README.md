# ipccwGEE
Inverse probability censoring cluster weighted generalized estimating equations for clustered longitudinal data with informative cluster size and informative drop-out

#### R installation Instructions
```
install.packages("devtools")
library(devtools)
devtools::install_github("AyaMitani/ipccwGEE")
library(ipccwGEE)
```
#### Example for using this package:
```
data(dental)
## First fit the drop model or missingness model.
outdrop <- dropmodel(toothstat ~ prevmaxcal5mm + basenumteeth, data = dental,
cluster.var = subject, unit.var = tooth, time.var = visit)
## Save the outputs.
outdat <- outdrop$outdata
Slist <- outdrop$Slist
## Restrict data set to include visits where the unit (tooth) was observed.
subdental <- subset(outdat, outdat$toothstat == 1)
## Fit marginal model with cluster-weights and IPCWs.
ipccwGEEout <- ipccwGEE(maxcal5mm ~ baseage + edu, data = subdental, cluster = subject, unit = tooth, 
ipcw = ipcw, Slist = Slist, show.iter = TRUE)
summary(ipccwGEEout)
```
