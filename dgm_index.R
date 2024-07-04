#create a .csv description of the simulated dgms; dgm: numerical indicator corresponding to each distinct dgm.
dt <- cbind(expand.grid(alpha0 = c(-1.499,0.15),
                        beta0 = c(-2.35,-0.63),rho = c(0,0.5),clus = c("50 clusters of 20 - unbl",
                                                                      "100 clusters of 10 - unbl"),rr = c(1,0.8),
                        trt_mod=c("none","random intercept","random slope"),
                        out_mod=c("none","random intercept","random slope"),unbal = 0.20))
dt$beta4 <- ifelse(dt$beta0==-2.35,log(0.73),-0.41)
dt$beta4 <- ifelse(dt$rr==1,0,dt$beta4)
dt$true.rd <- ifelse(dt$beta0==-2.35,-0.02,-0.06) #a column for the true ate
dt$true.rd <- ifelse(dt$rr==1,0,dt$true.rd)
dt <- cbind(dt,dgm=1:288)
#the below is a string description of each dgm - set exactly in the same order as the parameter value vectors in the 
#initial line of generation of dt <- cbind(expand.grid(alpha0 = c(-1.499,0.15),
#beta0 = c(-1.35,0.50),rho = c(0,0.5),clus = c("50 clusters of 20 - unbl",
#                                              "100 clusters of 10 - unbl"),rr = c(1,0.8)...)
dt <- cbind(dt,expand.grid(prev=c("trt prev 20%","trt prev 50%"),
                           rate=c("out control 10%","out control 30%"),correlation=c("no corr","corr"),
                           design=c("50 clusters of 20 - unbl","100 clusters of 10 - unbl"),true.effect=c("no effect","moderate effect")))
dt$NCluster <- ifelse(dt$clus=="50 clusters of 20 - unbl",50,100)
dt$ClusterSize <- ifelse(dt$clus=="50 clusters of 20 - unbl",20,10)

setwd("H:/My Documents/Simulation ISCB44/New Basic/Results")
write.csv(dt,"dgm_index.csv",row.names = FALSE)