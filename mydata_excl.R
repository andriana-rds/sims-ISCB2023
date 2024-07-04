#Create a data set which excludes any clusters where units are assigned only one treatment level; required for the clustered estimator 
#depending on the seed value; we could adjust the code to avoid even this minor change in our sample;
#E.g.: set e_hk = 0.7e_hk* + 0.15 to reassure adequate numbers of treated and untreated within each cluster; see Lee et al., 2021)

exclusion <- function(data) { mydata.excl <- data

#Exclude

index <- rep(0, NROW(unique(mydata.excl$mycluster)))
for (i in 1:NROW(unique(mydata.excl$mycluster))) {
  
  index[i] <- nrow(matrix(table(mydata.excl$Z[mydata.excl$mycluster==i]))) #index= 1, if only 1 treatment level, 2 if both treatment levels
}                                                                #within hospital (cluster)

#CHECKS
#table(index)
#prop.table(table(index)) #what proportion would be alarming for us to change our methods?

index <- data.frame("index"=index, "mycluster"= 1:length(unique((mydata.excl$mycluster))))
mydata.excl <- merge(index, mydata.excl, by = "mycluster") # table(mydata$index)
mydata.excl <- mydata.excl[mydata.excl$index==2,] #keep only hospitals with both treatment levels being assigned  

#CHECKS
#length(unique(mydata$mycluster)) 
mydata.excl <- mydata.excl[order(mydata.excl$mycluster,mydata.excl$Id),]
temp <- data.frame("mycluster"= unique(mydata.excl$mycluster), "clusterid" = 1:length(unique(mydata.excl$mycluster)))
mydata.excl <- merge(mydata.excl, temp, by = "mycluster") #cluster indicator for "mydata.excl" data set: clusterid

return(mydata.excl)
}