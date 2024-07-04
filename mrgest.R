# marginal estimator (Formula 7 in Li et al., 2013)
pi_marg <- function(mydata){
  
  w.1=sum(mydata$w*(mydata$Z==1))
  w.0=sum(mydata$w*(mydata$Z==0))
  Yhkwhk.1=sum(mydata$w*mydata$Y*(mydata$Z==1))
  Yhkwhk.0=sum(mydata$w*mydata$Y*(mydata$Z==0))
  pi_marg <- Yhkwhk.1/w.1-Yhkwhk.0/w.0
  return(pi_marg)
}
#(pi_marg(mydata = mydata))
#ATE