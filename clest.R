# non-parametric clustered estimator
pi_cl <- function(mydata) {
  
  pi_h=NULL
  for (j in unique(mydata$mycluster)){
    myd=subset(mydata, mycluster==j)
    w.1=sum(myd$w*(myd$Z==1))
    w.0=sum(myd$w*(myd$Z==0))
    Yhkwhk.1=sum(myd$w*myd$Y*(myd$Z==1))
    Yhkwhk.0=sum(myd$w*myd$Y*(myd$Z==0))
    pi_h <- c(pi_h, Yhkwhk.1/w.1-Yhkwhk.0/w.0)}
  
  
  # Formula 8 (in Li et al., 2013)
  temp1=aggregate(mydata[,c("w")] , 
                  by=list("mycluster"=mydata$mycluster), FUN=sum) # temp1$x represents wh of formula 8
  sumwhpih=sum(temp1$x*pi_h)
  sumwh=sum(temp1$x)
  
  pi_cl=sumwhpih/sumwh
  return(pi_cl)
}
#(pi_cl(mydata = mydata))
#ATE