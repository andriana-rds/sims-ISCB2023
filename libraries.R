# Install required packages
reqPcks <- c("ICCbin","boot","MASS","statmod","foreach","GLMMadaptive")#"lme4", 
for(p in reqPcks){
  if(!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)}
}

#library(lme4) #for mixed effects models
#library(ICCbin) #for unadjusted ICC calculation; AK: to revisit
#library(boot) #for bootstrap SEs
#library(MASS) #for multivariate normal distributions
#library(marginaleffects) #AK: probably will not be used
#library(statmod) #for identifying optimal quadrature points and weights when performing Gauss-Hermite quadrature
#library(foreach) #to apply the foreach function instead of a for loop