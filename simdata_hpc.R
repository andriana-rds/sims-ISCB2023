#function to simulate the data
#source: libraries.R; expit.R; cancerDataDesign.R; SimulTreat.R; SimulOut.R
#load dgm_index.csv for dgm identification 
setwd("/home/lsh2005314/New Basic/Code/")
dt <- read.csv("dgm_index.csv",header = TRUE)
source("libraries.R")
source("expit.R")
source("cancerDataDesign.R")
source("SimulTreat.R")
source("SimulOut.R")
#simdata returns a list of dataframes for each iteration; at the end we will end-up with 288 lists of data sets - or could be combined all in one list.
simdata <- function(i,dt,dgm){#i: iteration; dgm: 1-288, dgm index
  
  mydesign <- cancerDataDesign(NCluster=dt$NCluster[dt$dgm==dgm], ClusterSize=dt$ClusterSize[dt$dgm==dgm], unbal = dt$unbal[dt$dgm==dgm], 
                               correlated = dt$correlation[dt$dgm==dgm], rho = dt$rho[dt$dgm==dgm])
  simtreat <- SimulTreat(alpha0 = dt$alpha0[dt$dgm==dgm], alpha1 = 0.1,alpha2 = -0.2,alpha3 = 0.6, random = dt$trt_mod[dt$dgm==dgm],
                         true.mean0=0,true.mean1=0,true.sd0=0.5,true.sd1=0.5,true.cor01 = 0.3,
                         mydesign = mydesign)
  mydesigntr <- simtreat$mydesign
  mydesigntr$mycluster=as.factor(mydesigntr$mycluster)
  simout <- SimulOut(beta0 = dt$beta0[dt$dgm==dgm], beta1 = 0.2, beta2 = -1, beta3 = 1, beta4 = dt$beta4[dt$dgm==dgm],
                     random = dt$out_mod[dt$dgm==dgm],true.meanY0 = 0,true.meanY1 = 0,true.sdY0=0.5,true.sdY1=0.5,true.cor01= 0.3,tab2=mydesigntr)
  data <- simout$simdata
  
  data[["i"]] <- i
  data[["dgm"]] <- dgm
  return(data)
  
}
dgm <-as.numeric(Sys.getenv("SGE_TASK_ID"))
#k <- 10 #number of dgms
#dgm <- 1:k
B <- 10 #number of iterations
#set.seed(20230720)

comput.time <- system.time(dfs <- #foreach(j=dgm) %do% {
  foreach(i=1:B) %do% {
    print(i)
    set.seed(2711*dgm)
    Sys.sleep(3)
  simdata(i=i,dt=dt,dgm=dgm)
 }
#}
)
names(dfs) <- paste0("dgm", dgm)

str(dfs)
