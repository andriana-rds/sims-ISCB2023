cancerDataDesign <- function(NCluster, ClusterSize, unbal = 0.20, correlated = "corr", rho = 0.5){
  
  #  NCluster : Total number of clusters
  #  ClusterSize : Number of observations per cluster for balanced design (or number of obs per cluster on average for an unbalanced design)
  n = NCluster*ClusterSize 
  
  
  #total sample size (i.e., total number of cancer patients) for either balanced or unbalanced design
  #corresponds to an unbalanced design; if NULL, design is balanced
  if (is.null(unbal)==FALSE){cluster.sizes <- rmultinom(1,n,runif(NCluster)^(unbal)) #AK: not too sure yet if the effect of setting a specific value for unbal impacts the range of the 
  #cluster sizes or whether it is the seed that truly impacts the results - to write a detailed explanation (although, probably 'unbal' value makes the difference)
  
  if (sum(cluster.sizes==0)>0){cluster.sizes <- rmultinom(1,n,runif(NCluster)^(unbal))}#if clusters with zero obs exist, then draw again
  else {cluster.sizes <- cluster.sizes} #else, keep the cluster sizes from initial draw
  } 
  #else 
  if (is.null(unbal)==TRUE){cluster.sizes <- rep(ClusterSize,NCluster)} #corresponds to a balanced design
  #if (is.null(unbal)==TRUE)
  #discard any clusters with zero obs 
  if (sum(cluster.sizes==0)>0){cluster.sizes <- cluster.sizes[cluster.sizes!=0]}
  NCluster <- length(cluster.sizes) #total number of clusters now does not include those with zero obs, if any
  n = sum(cluster.sizes) #sanity check: total sample size redefined (anticipated to be unchanged)
  
  # --------
  # patient ID
  #Id <- c(1:n)
  Id <- foreach(i=1:NCluster, .combine = 'c') %do% rep(1:cluster.sizes[i]) #create an ID for each patient within each cluster
  
  # .combine = 'c' is to concatenate results into a vector
  # U1 = continuous individual-level covariate
  U1 <- rnorm(n, 0, 1)
  
  # U2 = binary individual-level covariate
  U2 <- rbinom(n, 1, prob=0.5)
  
  # V = Cluster-level covariate: assumed continuous within this function; AK note: could decrease number of lines by applying CL's approach to create the final data set
  V <- rnorm(NCluster, 0, 1)
  
  # mycluster = Cluster ID
  #mycluster <- sample(rep(1:NCluster,each=n/NCluster))
  mycluster <- rep(1:NCluster, times=cluster.sizes) #create an ID for each cluster 
  
  dat_cl <- data.frame(myc=sort(unique(mycluster)), V)
  
  # Create a data frame to hold the variables created
  dat_ind <- data.frame(Id, U1, U2, mycluster)
  
  # Merge the cluster data frame to the cancer data frame
  dat <- merge(dat_ind, dat_cl, by.x="mycluster", by.y="myc", all.x = TRUE)
  
  # U1 could be either correlated with V (with assumed correlation = rho) or independent of V
  dat$U1_2 <- rho*dat$V + sqrt(1 - rho**2)*dat$U1
  
  if(correlated=="corr"){dat$U1<-dat$U1_2}
  if(correlated=="no corr"){dat$U1<-dat$U1}
  
  
  # Return the data frame created
  return(dat)
  
}

#AK note: to be on the safe side, always insert "TRUE" (or "FALSE") instead of just T (or F, respectively), as the former
#are already reserved objects by R, whereas T and F are not.