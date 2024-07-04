# Simulate treatment values;
# general function to simulate either a random intercept or random slope model
# random: could be "none": neither random effects nor cluster-level covariate, V; "random intercept": random inter and V; 
# "random slope": random slope for U1; V included
# by default it applies a random slope for U1

SimulTreat <- function(alpha0,alpha1,alpha2,alpha3,random = "random slope",true.mean0=0,true.mean1=0,true.sd0=0.5,true.sd1=0.5,true.cor01 = 0.3,
                       mydesign){  
  
  #A: PS model
  #set fixed alphas; for U1: alpha1, U2: alpha2 and V: alpha3
  
  alpha.vec <-  c(alpha0,alpha1,alpha2,alpha3)
  NCluster <- length(unique(mydesign$mycluster))
  n <- nrow(mydesign)
  
  if (random == "none"){# when "none", we generate treatment to be affected only by U1, U2
    alpha.vec <-  c(alpha0,alpha1,alpha2)
    
    true.sigma <- NULL
    
    tab2 <- mydesign[order(mydesign$mycluster,mydesign$Id),]
    
    X <- as.matrix(cbind(X0=1,tab2[,c("U1", "U2")])) # Design matrix
    
    #prob p <- 1/(1 + exp(-linpred)) 
    probz <- as.vector(1/(1 + exp(-(X %*% alpha.vec)))) #or: plogis(linpred) 
    
  }
  
  if (random == "random intercept"){
    # simulate a random intercept (by cluster, so in the same way as we did for the cluster level covariate);
    # Normal dist with mean = true.mean and sd = true.sd
    true.mean <- true.mean0
    true.sigma <- true.sd0
    rand.eff <- rnorm(NCluster, mean = true.mean, sd = true.sigma)
    temp2 <- data.frame(mycluster = rep(1:NCluster), rand.eff = rand.eff)
    
    # merge random intercept column with tab2 (our simulated data set)
    tab2 <- merge(mydesign, temp2, by.x = "mycluster", by.y = "mycluster", all.x=T)
    tab2 <- tab2[order(tab2$mycluster,tab2$Id),]
    
    #1. calculate the Xalphas+rand.eff (Normal with mean=0 and sd = 2)
    #linear predictor
    X <- as.matrix(cbind(X0=1,tab2[,c("U1", "U2", "V")])) # Design matrix
    
    #prob p <- 1/(1 + exp(-linpred)) where linpred = X %*% alpha.vec + rand.eff
    probz <- as.vector(1/(1 + exp(-(X %*% alpha.vec + tab2$rand.eff)))) #or: plogis(X %*% alpha.vec + tab2$rand.eff)
  }
  
  # simulate the random effects dist (by cluster, so in the same way as we did for the cluster level covariate);
  
  if (random == "random slope"){# Random intercept and slope for U1 follow: joint bivariate normal dist
  
  true.mean <- c(true.mean0,true.mean1) #vector of means
  true.sigma <- matrix(c(true.sd0^2,true.cor01*true.sd0*true.sd1,true.cor01*true.sd0*true.sd1,true.sd1^2),2,2) #sanity check: Sigma is positive definite, i.e., there is x, so that: transpose(x)*Sigma*x >0
  rand.eff <- mvrnorm(n = NCluster, mu = true.mean, Sigma = true.sigma)
  temp2 <- data.frame(mycluster = rep(1:NCluster), rand.eff0 = rand.eff[,1], rand.eff1 = rand.eff[,2])
  
  # merge random intercept & slope columns with tab2 (simulated data set)
  tab2 <- merge(mydesign, temp2, by.x = "mycluster", by.y = "mycluster", all.x=T)
  tab2 <- tab2[order(tab2$mycluster,tab2$Id),]
  
  #1. calculate the linpred = X*alphas + rand.eff0 + rand.eff1 * U1 
  #linear predictor
  X <- as.matrix(cbind(X0=1,tab2[,c("U1", "U2", "V")])) # Design matrix
  
  #prob p <- 1/(1 + exp(-linpred)) 
  probz <- as.vector(1/(1 + exp(-(X %*% alpha.vec + tab2$rand.eff0 + tab2[,"U1"] * tab2$rand.eff1)))) #or: plogis(linpred)
}
  
  #2. use the probabilities to simulate the treatment value
  tab2$Z <- rbinom(n = n, size = 1, prob = probz) #table(Z)
  
  return(list("mydesign"=tab2, "true alpha"=alpha.vec, "true sigma trt"= true.sigma))
}

