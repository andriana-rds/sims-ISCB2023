# Simulate outcome values;
# general function to simulate either a random intercept or random slope model
# random: could be "none": no random effects, but V inlcuded; "random intercept": random inte, V included; 
# "random slope": random slope for U1; V included
# by default it applies a random slope for U1

SimulOut <- function(beta0,beta1,beta2,beta3,beta4,random = "random slope",
                     true.meanY0 = 0,true.meanY1 = 0,true.sdY0=0.5,true.sdY1=0.5,true.cor01= 0.3,tab2)
{ 
  # B: Outcome model
  
  beta.vec <-  c(beta0,beta1,beta2,beta3,beta4)
  n <- nrow(tab2)
  NCluster <- length(unique(tab2$mycluster))
  
  if (random == "none"){
    
    true.sigma <- NULL
    
    #linear predictor
    #X_y <- as.matrix(cbind(X0=1,tab2[,c("U1", "U2", "V","Z")])) # design matrix of the outcome Y
    
    X_y0 <- as.matrix(cbind(X0=1,tab2[,c("U1", "U2", "V")],Z=0)) # design matrix of the potential outcome under "control"
    X_y1 <- as.matrix(cbind(X0=1,tab2[,c("U1", "U2", "V")],Z=1)) # design matrix of the potential outcome under "treatment"
    
    #prob p <- 1/(1 + exp(-linpred)) where linpred = X_y %*% beta.vec + rand.effY
    #probY <- as.vector(1/(1 + exp(-(X_y %*% beta.vec + tab2$rand.effY)))) #or: plogis(X %*% alpha.vec + tab2$rand.eff)
    
    #Add two columns for the potential outcome values Y0, Y1; THESE ARE RANDOM VARIABLES
    probY0 <- as.vector(1/(1 + exp(-(X_y0 %*% beta.vec)))) 
    probY1 <- as.vector(1/(1 + exp(-(X_y1 %*% beta.vec))))
  }
  
  if (random == "random intercept"){
    true.sigma <- true.sdY0
    # simulate a random intercept (by cluster, so in the same way as we did for the cluster level covariate);
    # Normal distribution with mean = true.meanY0 and sd = true.sdY0
    temp3 <- data.frame(mycluster = rep(1:NCluster), 
                        rand.effY = rnorm(NCluster, mean = true.meanY0, sd = true.sigma))
    
    # merge random intercept column with tab2 (our simulated data set)
    tab2 <- merge(tab2, temp3, by.x = "mycluster", by.y = "mycluster", all.x=T)
    tab2 <- tab2[order(tab2$mycluster,tab2$Id),]
    
    
    #linear predictor
    #X_y <- as.matrix(cbind(X0=1,tab2[,c("U1", "U2", "V","Z")])) # @@AK design matrix of the outcome Y
    
    X_y0 <- as.matrix(cbind(X0=1,tab2[,c("U1", "U2", "V")],Z=0)) # @@AK design matrix of the potential outcome under "control"
    X_y1 <- as.matrix(cbind(X0=1,tab2[,c("U1", "U2", "V")],Z=1)) # @@AK design matrix of the potential outcome under "treatment"
    
    #prob p <- 1/(1 + exp(-linpred)) where linpred = X_y %*% beta.vec + rand.effY
    #probY <- as.vector(1/(1 + exp(-(X_y %*% beta.vec + tab2$rand.effY)))) #or: plogis(X %*% alpha.vec + tab2$rand.eff)
    
    #Add two columns for the potential outcome values Y0, Y1; THESE ARE RANDOM VARIABLES
    probY0 <- as.vector(1/(1 + exp(-(X_y0 %*% beta.vec + tab2$rand.effY)))) 
    probY1 <- as.vector(1/(1 + exp(-(X_y1 %*% beta.vec + tab2$rand.effY))))
    
  }
  
  # simulate the random effects dist (by cluster, i.e., as we did for the cluster level covariate)

  if (random == "random slope"){# Random intercept and a random slope for U1 follow: joint bivariate normal dist
    
    true.mean <- c(true.meanY0,true.meanY1) #vector of means
    true.sigma <- matrix(c(true.sdY0^2,true.cor01*true.sdY0*true.sdY1,true.cor01*true.sdY0*true.sdY1,true.sdY1^2),2,2) #sanity check: Sigma is positive definite, i.e., there is x, so that: transpose(x)*Sigma*x >0
    rand.eff <- mvrnorm(n = NCluster, mu = true.mean, Sigma = true.sigma)
    temp3 <- data.frame(mycluster = rep(1:NCluster), rand.effY0 = rand.eff[,1], rand.effY1 = rand.eff[,2])
  
  #temp3 <- data.frame(mycluster = rep(1:NCluster), 
                      #rand.effY = rnorm(NCluster, mean = true.meanY, sd = true.sdY))
  
  # merge random intercept & slope column with tab2 (our simulated data set)
  tab2 <- merge(tab2, temp3, by.x = "mycluster", by.y = "mycluster", all.x=T)
  tab2 <- tab2[order(tab2$mycluster,tab2$Id),]
  
  
  #linear predictor
  #X_y <- as.matrix(cbind(X0=1,tab2[,c("U1", "U2", "V","Z")])) # design matrix of the outcome Y
  
  X_y0 <- as.matrix(cbind(X0=1,tab2[,c("U1", "U2", "V")],Z=0)) # design matrix of the potential outcome under "control"
  X_y1 <- as.matrix(cbind(X0=1,tab2[,c("U1", "U2", "V")],Z=1)) # design matrix of the potential outcome under "treatment"
  
  #prob p <- 1/(1 + exp(-linpred)) where linpred = X_y %*% beta.vec + rand.effY
  #probY <- as.vector(1/(1 + exp(-(X_y %*% beta.vec + tab2$rand.effY)))) #or: plogis(X %*% alpha.vec + tab2$rand.eff)
  
  #Add two columns for the potential outcome values Y0, Y1; THESE ARE RANDOM VARIABLES
  probY0 <- as.vector(1/(1 + exp(-(X_y0 %*% beta.vec + tab2$rand.effY0 + tab2[,"U1"] * tab2$rand.effY1)))) 
  probY1 <- as.vector(1/(1 + exp(-(X_y1 %*% beta.vec + tab2$rand.effY0 + tab2[,"U1"] * tab2$rand.effY1)))) #AK: maybe to change names here, to avoid confusion between random intercept values and potential outcome under control level, for example, etc.
  }
  #2. use the probabilities to simulate the potential outcome values Y0, Y1 (recall: binary outcome)
  
  tab2$Y0 <- rbinom(n = n, size = 1, prob = probY0) # table(tab2$Y0)
  tab2$Y1 <- rbinom(n = n, size = 1, prob = probY1) # table(tab2$Y1)
  
  #and the observed outcome values, Y 
  tab2$Y <- ifelse(tab2$Z==0,tab2$Y0,tab2$Y1) # table(tab2$Y)
  
  return(list("simdata" = tab2, "true beta" = beta.vec, "true sigma out" = true.sigma))
  
}