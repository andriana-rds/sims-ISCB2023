# DR parametric estimator

pi_dr <- function(mydata) {
  
  #Step i): within cluster ATE estimates
  pihk = NULL
  for (j in unique(mydata$mycluster)){
    myd = subset(mydata, mycluster==j)
    pihk.1 = myd$Z*myd$Y/myd$ps - (myd$Z - myd$ps)*myd$Y1.hat/myd$ps
    pihk.0 = (1-myd$Z)*myd$Y/(1-myd$ps) + (myd$Z - myd$ps)*myd$Y0.hat/(1 - myd$ps)
    pihk <- c(pihk, pihk.1 - pihk.0)}
  
  # Step ii): DR estimator (formula 9 in Li et al., 2013)
  pi_dr = mean(pihk)
  return(list("pi_dr"=pi_dr, "pihk"=pihk))
  
}
#(pi_dr(mydata = mydata))
#ATE
