# Function to analyse the simulated data sets
# across analysis methods we omit the cluster-level covariate, V
#ate <- se <- convergence <- trt_prev <- out_rate <- NULL 
#Run in advance: expit.R; libraries.R; mydata_excl.R; mrgest.R; clest.R; drest.R
#R: number of bootstrap replicates

#dt includes the dgms description
setwd("H:/My Documents/Simulation ISCB44/New Basic/Results")
dt <- read.csv("dgm_index.csv",header = TRUE)

analysis_bin <- function(i,dt,dgm,data){#AK: to incude R version and software as inputs
  
  trt_mod <- dt$trt_mod[dt$dgm==dgm]
  out_mod <- dt$out_mod[dt$dgm==dgm]
  
  t0<-t1<-data
  t0$Z<-0
  t1$Z<-1
  
  trt_prev <- mean(data$Z) #trt prevalence
  out_rate <- mean(data$Y0) #event rate in the controls
  
  #Design matrices - omit V across analysis models
  X_y0 <- as.matrix(cbind(X0=1,data[,c("U1", "U2")],Z=0)) # design matrix of potential outcome under "control" 
  X_y1 <- as.matrix(cbind(X0=1,data[,c("U1", "U2")],Z=1)) # design matrix of potential outcome under "treatment"
  X_z <- as.matrix(cbind(X0=1,data[,c("U1", "U2")])) # design matrix of treatment
  #Design matrices - benchmark
  X_y0.bench <- as.matrix(cbind(X0=1,data[,c("U1", "U2", "V")],Z=0)) # design matrix of potential outcome under "control" 
  X_y1.bench <- as.matrix(cbind(X0=1,data[,c("U1", "U2", "V")],Z=1)) # design matrix of potential outcome under "treatment" 
  
  #boot.res <- boot(data = dat, statistic = npar_boot, R = R) #TO CHECK
  #Benchmark model(s) - AK: favorable toward G-estimation? 
  #1. we predict the cluster-specific probabilities by plugging-in the Empirical Bayes estimates (EBE) of random effects
  #2. we predict the marginal probabilities over the random effects distribution - as in Hedeker et al., 2018.
  #bench.mod <- bench(data,out_mod)
  
  
  #system.time(
  if (out_mod == "random slope"){
    
    #Empirical Bayes Estimates: benchmark  ebe
    bench.mod <- tryCatch(mixed_model(fixed = Y ~ U1 + as.factor(U2) + V + as.factor(Z), 
                                      random = ~1+U1|mycluster, data = data, family = binomial(),control = list(nAGQ=7)),
                          error = function(e) NULL)
    if(!is.null(bench.mod)){
    beta.vec <- beta.vec1 <- fixef(bench.mod)
    temp <- data.frame(mycluster = rep(1:length(unique(data$mycluster))), 
                       ebe0.bench = ranef(bench.mod)[,1],ebe1.bench = ranef(bench.mod)[,2])
    data <- merge(data, temp, by.x = "mycluster", by.y = "mycluster", all.x=T)
    y0 <- as.vector(1/(1 + exp(-(X_y0.bench %*% beta.vec + data$ebe0.bench + data$U1*data$ebe1.bench)))) 
    y1 <- as.vector(1/(1 + exp(-(X_y1.bench %*% beta.vec + data$ebe0.bench + data$U1*data$ebe1.bench))))
    #convergence check
    ate1 <- ifelse(bench.mod$converged==TRUE,mean(y1) - mean(y0),NA)
    #se1 <- 
    
    #Hedeker et al. estimation of the population average probabilities for benchmark model: benchmark  pa
    
    y0<-predict(bench.mod, newdata = t0,
                type = "marginal") #'arg': "subject_specific" does not plug-in the empirical Bayes estimates; it plugs-in 
    #the posterior means instead. This brings completely different ATE estimates. To be investigated.
    y1<-predict(bench.mod, newdata = t1,
                type = "marginal") #we apply Monte Carlo integration for the integral approximation
    ate2 <- ifelse(bench.mod$converged==TRUE,mean(y1) - mean(y0),NA)
    #se2 <- 
    convergence <- ifelse(bench.mod$converged==TRUE,1,0)
    } else {ate1 <- ate2 <- NA; se1 <- se2 <- NA; convergence <- 0}
    m1 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate1,
      #se = se1,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 1
    )
    
    m2 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate2,
      #se = se2,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 2
    )
    
    #G-computation estimators
    out.mod <- tryCatch(mixed_model(fixed = Y ~ U1 + as.factor(U2) + as.factor(Z), random = ~1+U1|mycluster, data = data, family = binomial(),
                                    control = list(nAGQ=7)),
                        error = function(e) NULL)
    if (!is.null(out.mod)){
    #Empirical Bayes Estimates: g-comp ebe
    beta.vec <- beta.vec2 <- fixef(out.mod)
    temp <- data.frame(mycluster = rep(1:length(unique(data$mycluster))), 
                       ebe0.g = ranef(out.mod)[,1],ebe1.g = ranef(out.mod)[,2])
    data <- merge(data, temp, by.x = "mycluster", by.y = "mycluster", all.x=T)
    data$Y0.hat.ebe <- as.vector(1/(1 + exp(-(X_y0 %*% beta.vec + data$ebe0.g + data$U1*data$ebe1.g)))) #AK: replace with expit() to reduce code
    data$Y1.hat.ebe <- as.vector(1/(1 + exp(-(X_y1 %*% beta.vec + data$ebe0.g + data$U1*data$ebe1.g))))
    #convergence check
    ate3 <- ifelse(out.mod$converged==TRUE,mean(data$Y1.hat.ebe) - mean(data$Y0.hat.ebe),NA)
    #se3 <-
    
    #Hedeker et al. estimation of the population average probabilities for benchmark model: g-comp pa
    #reminder: integral is approximated via Monte Carlo integration here
    data$Y0.hat.pa<-predict(out.mod, newdata = t0,
                            type = "marginal") #'arg': "subject_specific" does not plug-in the empirical Bayes estimates; it plugs-in 
    #the posterior means instead. This brings completely different ATE estimates. To be investigated.
    data$Y1.hat.pa<-predict(out.mod, newdata = t1,
                            type = "marginal") #we apply Monte Carlo integration for the integral approximation
    ate4 <- ifelse(out.mod$converged==TRUE,mean(data$Y1.hat.pa) - mean(data$Y0.hat.pa),NA)
    #se4 <- 
    convergence <- ifelse(out.mod$converged==TRUE,1,0)} else {ate3 <- ate4 <- NA; se3 <- se4 <- NA; convergence <- 0}
    m3 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate3,
      #se = se3,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 3
    )
    
    m4 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate4,
      #se = se4,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 4
    )
    
    #IPTW estimators: marginal and clustered estimators - Li et al., 2013
    #marginal ebe
    ps.mod <- tryCatch(mixed_model(fixed = Z ~ U1 + as.factor(U2), random = ~1+U1|mycluster, data = data, family = binomial(),
                                   control = list(nAGQ=7)),
                       error = function(e) NULL)
    if (!is.null(ps.mod)){
    beta.vec <- beta.vec3 <- fixef(ps.mod)
 
    temp <- data.frame(mycluster = rep(1:length(unique(data$mycluster))), 
                       ebe0.ps = ranef(ps.mod)[,1],ebe1.ps = ranef(ps.mod)[,2])
    data <- merge(data, temp, by.x = "mycluster", by.y = "mycluster", all.x=T)
    data$ps.ebe <- as.vector(1/(1 + exp(-(X_z %*% beta.vec + data$ebe0.ps + data$U1*data$ebe1.ps)))) #AK: replace with expit() to reduce code
    data$ps <- data$ps.ebe
    data$w <- ifelse(data$Z==0,1/(1-data$ps),1/data$ps) 
    ate5 <- ifelse(ps.mod$converged==TRUE,pi_marg(mydata = data),NA) #Marginal IPTW
    #se5 <- 
    
    #marginal pa
    data$ps.pa <- fitted(ps.mod, type = "marginal")
    data$ps <- data$ps.pa
    data$w <- ifelse(data$Z==0,1/(1-data$ps),1/data$ps) 
    ate6 <- ifelse(ps.mod$converged==TRUE,pi_marg(mydata = data),NA) #Marginal IPTW
    #se6 <- 
    
    convergence <- ifelse(ps.mod$converged==TRUE,1,0)} else {ate5 <- ate6 <- NA; se5 <- se6 <- NA; convergence <- 0}
    m5 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate5,
      #se = se5,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 5
    )
    
    m6 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate6,
      #se = se6,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 6
    )
    
    #Clustered IPTW
    #exclusion of clusters where all units are assigned only one treatment level
    #setwd("H:/My Documents/Simulation ISCB44/New Basic/Code")
    #source("mydata_excl.R")
    mydata.excl <- exclusion(data)
    mydata.excl$mycluster <- mydata.excl$clusterid #see function exclusion()
    
    #cl ebe
    ps.mod.ex <- tryCatch(mixed_model(fixed = Z ~ U1 + as.factor(U2), random = ~1+U1|mycluster, data = mydata.excl, family = binomial(),
                                      control = list(nAGQ=7)),
                          error = function(e) NULL)
    if (!is.null(ps.mod.ex)){
      
    beta.vec <- beta.vec4 <- fixef(ps.mod.ex)
    
    temp <- data.frame(mycluster = rep(1:length(unique(mydata.excl$mycluster))), 
                       ebe0.ps.ex = ranef(ps.mod.ex)[,1],ebe1.ps.ex = ranef(ps.mod.ex)[,2])
    mydata.excl <- merge(mydata.excl, temp, by.x = "mycluster", by.y = "mycluster", all.x=T)
    X_z.ex <- as.matrix(cbind(X0=1,mydata.excl[,c("U1", "U2")])) # design matrix of treatment for the data set of remaining clusters
    
    mydata.excl$ps <- as.vector(1/(1 + exp(-(X_z.ex %*% beta.vec + mydata.excl$ebe0.ps.ex + mydata.excl$U1*mydata.excl$ebe1.ps.ex)))) #AK: replace with expit() to reduce code
    mydata.excl$w <- ifelse(mydata.excl$Z==0,1/(1-mydata.excl$ps),1/mydata.excl$ps) 
    ate7 <- ifelse(ps.mod.ex$converged==TRUE,pi_cl(mydata = mydata.excl),NA) #clustered IPTW
    #se7 <- 
    
    #cl pa
    mydata.excl$ps <- fitted(ps.mod.ex, type = "marginal")
    mydata.excl$w <- ifelse(mydata.excl$Z==0,1/(1-mydata.excl$ps),1/mydata.excl$ps) 
    ate8 <- ifelse(ps.mod.ex$converged==TRUE,pi_cl(mydata = mydata.excl),NA) #Marginal IPTW
    #se8 <- 
    
    convergence <- ifelse(ps.mod.ex$converged==TRUE,1,0)} else {ate7 <- ate8 <- NA; se7 <- se8 <- NA; convergence <- 0}
    m7 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate7,
      #se = se7,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 7
    )
    
    m8 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate8,
      #se = se8,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 8
    )
    
  
    #DOUBLY ROBUST ESTIMATOR: AIPTW 
    if (!is.null(out.mod) & !is.null(ps.mod)){
    #aipw ebe
    data$Y0.hat <- data$Y0.hat.ebe
    data$Y1.hat <- data$Y1.hat.ebe
    data$ps <- data$ps.ebe
    ate9 <- ifelse(out.mod$converged==TRUE & ps.mod$converged==TRUE,pi_dr(mydata = data)$pi_dr,NA) #Doubly robust estimator: AIPTW
    #se9 <- 
    
    #aipw pa
    data$Y0.hat <- data$Y0.hat.pa
    data$Y1.hat <- data$Y1.hat.pa
    data$ps <- data$ps.pa
    ate10 <- ifelse(out.mod$converged==TRUE & ps.mod$converged==TRUE,pi_dr(mydata = data)$pi_dr,NA) #Doubly robust estimator: AIPTW
    #se10 <- 
    
    convergence <- ifelse(out.mod$converged==TRUE & ps.mod$converged==TRUE,1,0)} else {ate9 <- ate10 <- NA; se9 <- se10 <- NA; convergence <- 0}
    m9 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate9,
      #se = se9,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 9
    )
    
    m10 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate10,
      #se = se10,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 10
    )
    
  } else {
    
    #we analyze via random intercepts either r.i or single-level models
    #Empirical Bayes Estimates: benchmark - ebe
    bench.mod <- tryCatch(mixed_model(fixed = Y ~ U1 + as.factor(U2) + V + as.factor(Z), random = ~1|mycluster, data = data, family = binomial(),
                                      control = list(nAGQ=7)),
                          error = function(e) NULL)
    if (!is.null(bench.mod)){
    beta.vec <- beta.vec1 <- fixef(bench.mod)
    temp <- data.frame(mycluster = rep(1:length(unique(data$mycluster))), 
                       ebe0.bench = ranef(bench.mod)[,1])
    
    data <- merge(data, temp, by.x = "mycluster", by.y = "mycluster", all.x=T)
    y0 <- as.vector(1/(1 + exp(-(X_y0.bench %*% beta.vec + data$ebe0.bench)))) 
    y1 <- as.vector(1/(1 + exp(-(X_y1.bench %*% beta.vec + data$ebe0.bench))))
    #convergence check
    ate1 <- ifelse(bench.mod$converged==TRUE,mean(y1) - mean(y0),NA)
    #se1 <-
    
    #Hedeker et al. estimation of the population average probabilities for benchmark model: benchmrak - pa
    
    y0<-predict(bench.mod, newdata = t0,
                type = "marginal") #'arg': "subject_specific" does not plug-in the empirical Bayes estimates; it plugs-in 
    #the posterior means instead. This brings completely different ATE estimates. To be investigated.
    y1<-predict(bench.mod, newdata = t1,
                type = "marginal") #we apply Monte Carlo integration for the integral approximation
    ate2 <- ifelse(bench.mod$converged==TRUE,mean(y1) - mean(y0),NA)
    #se2 <- 
    
    convergence <- ifelse(bench.mod$converged==TRUE,1,0)} else {ate1 <- ate2 <- NA; se1 <- se2 <- NA; convergence <- 0}
    m1 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate1,
      #se = se1,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 1
    )
    
    m2 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate2,
      #se = se2,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 2
    )
    
    #we analyze via random intercepts either r.i or single-level models
    #Empirical Bayes Estimates: g-comp ebe
    out.mod <- tryCatch(mixed_model(fixed = Y ~ U1 + as.factor(U2) + as.factor(Z), random = ~1|mycluster, data = data, family = binomial(),
                                    control = list(nAGQ=7)),
                        error = function(e) NULL)
    if (!is.null(out.mod)){
    beta.vec <- beta.vec2 <- fixef(out.mod)
    temp <- data.frame(mycluster = rep(1:length(unique(data$mycluster))), 
                       ebe0.g = ranef(out.mod)[,1])
    
    data <- merge(data, temp, by.x = "mycluster", by.y = "mycluster", all.x=T)
    data$Y0.hat.ebe <- as.vector(1/(1 + exp(-(X_y0 %*% beta.vec + data$ebe0.g)))) 
    data$Y1.hat.ebe <- as.vector(1/(1 + exp(-(X_y1 %*% beta.vec + data$ebe0.g))))
    #convergence check
    ate3 <- ifelse(out.mod$converged==TRUE,mean(data$Y1.hat.ebe) - mean(data$Y0.hat.ebe),NA)
    #se3 <- 
    
    #Hedeker et al. estimation of the population average probabilities for benchmark model: g-comp pa
    
    data$Y0.hat.pa <-predict(out.mod, newdata = t0,
                             type = "marginal") #'arg': "subject_specific" does not plug-in the empirical Bayes estimates; it plugs-in 
    #the posterior means instead. This brings completely different ATE estimates. To be investigated.
    data$Y1.hat.pa <-predict(out.mod, newdata = t1,
                             type = "marginal") #we apply Monte Carlo integration for the integral approximation
    ate4 <- ifelse(out.mod$converged==TRUE,mean(data$Y1.hat.pa) - mean(data$Y0.hat.pa),NA)
    #se4 <- 
    
    convergence <- ifelse(out.mod$converged==TRUE,1,0)} else {ate3 <- ate4 <- NA; se3 <- se4 <- NA; convergence <- 0}
    m3 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate3,
      #se = se3,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 3
    )
  
    m4 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate4,
      #se = se4,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 4
    )
    
    
    #IPTW estimators: marginal and clustered estimators - Li et al., 2013
    #marginal ebe
    ps.mod <- tryCatch(mixed_model(fixed = Z ~ U1 + as.factor(U2), random = ~1|mycluster, data = data, family = binomial(),
                                   control = list(nAGQ=7)),
                       error = function(e) NULL)
    if (!is.null(ps.mod)){
    
    beta.vec <- beta.vec3 <- fixef(ps.mod)
    
    temp <- data.frame(mycluster = rep(1:length(unique(data$mycluster))), 
                       ebe0.ps = ranef(ps.mod)[,1])
    data <- merge(data, temp, by.x = "mycluster", by.y = "mycluster", all.x=T)
    data$ps.ebe <- as.vector(1/(1 + exp(-(X_z %*% beta.vec + data$ebe0.ps)))) #AK: replace with expit() to reduce code
    data$ps <- data$ps.ebe
    data$w <- ifelse(data$Z==0,1/(1-data$ps),1/data$ps) 
    ate5 <- ifelse(ps.mod$converged==TRUE,pi_marg(mydata = data),NA) #Marginal IPTW
    #se5 <- 
    
    #marginal pa
    data$ps.pa <- fitted(ps.mod, type = "marginal") #AK: fitted(..., type = "") should be giving the same output as predict(..., type = "")
    data$ps <- data$ps.pa
    data$w <- ifelse(data$Z==0,1/(1-data$ps),1/data$ps) 
    ate6 <- ifelse(ps.mod$converged==TRUE,pi_marg(mydata = data),NA) #Marginal IPTW
    #se6 <- 
    
    convergence <- ifelse(ps.mod$converged==TRUE,1,0)} else {ate5 <- ate6 <- NA; se5 <- se6 <- NA; convergence <- 0}
    m5 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate5,
      #se = se5,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 5
    )
    
    m6 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate6,
      #se = se6,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 6
    )
    
    #Clustered IPTW
    #exclusion of clusters where all units are assigned only one treatment level
    #setwd("H:/My Documents/Simulation ISCB44/New Basic/Code")
    #source("mydata_excl.R")
    mydata.excl <- exclusion(data)
    mydata.excl$mycluster <- mydata.excl$clusterid
    
    #cl ebe
    ps.mod.ex <- tryCatch(mixed_model(fixed = Z ~ U1 + as.factor(U2), random = ~1|mycluster, data = mydata.excl, family = binomial(),
                                      control = list(nAGQ=7)),
                          error = function(e) NULL)
    if (!is.null(ps.mod.ex)){
      
    beta.vec <- beta.vec4 <- fixef(ps.mod.ex)
    
    temp <- data.frame(mycluster = rep(1:length(unique(mydata.excl$mycluster))), 
                       ebe0.ps.ex = ranef(ps.mod.ex)[,1])
    mydata.excl <- merge(mydata.excl, temp, by.x = "mycluster", by.y = "mycluster", all.x=T)
    X_z.ex <- as.matrix(cbind(X0=1,mydata.excl[,c("U1", "U2")])) # design matrix of treatment for the data set of remaining clusters
    
    mydata.excl$ps <- as.vector(1/(1 + exp(-(X_z.ex %*% beta.vec + mydata.excl$ebe0.ps.ex)))) #AK: replace with expit() to reduce code
    mydata.excl$w <- ifelse(mydata.excl$Z==0,1/(1-mydata.excl$ps),1/mydata.excl$ps) 
    ate7 <- ifelse(ps.mod.ex$converged==TRUE,pi_cl(mydata = mydata.excl),NA) #clustered IPTW
    #se7 <- 
    
    #cl pa
    mydata.excl$ps <- fitted(ps.mod.ex, type = "marginal")
    mydata.excl$w <- ifelse(mydata.excl$Z==0,1/(1-mydata.excl$ps),1/mydata.excl$ps) 
    ate8 <- ifelse(ps.mod.ex$converged==TRUE,pi_cl(mydata = mydata.excl),NA) #clustered IPTW
    #se8 <-
    
    convergence <- ifelse(ps.mod.ex$converged==TRUE,1,0)} else {ate7 <- ate8 <- NA; se7 <- se8 <- NA; convergence <- 0}
    m7 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate7,
      #se = se7,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 7
    )
    
    m8 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate8,
      #se = se8,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 8
    )
   
    if (!is.null(out.mod) & !is.null(ps.mod)){
    #DOUBLY ROBUST ESTIMATOR: AIPTW 
    #aipw ebe
    data$Y0.hat <- data$Y0.hat.ebe
    data$Y1.hat <- data$Y1.hat.ebe
    data$ps <- data$ps.ebe
    ate9 <- ifelse(out.mod$converged==TRUE & ps.mod$converged==TRUE,pi_dr(mydata = data)$pi_dr,NA) #Doubly robust estimator: AIPTW
    #se9 <- 
    
    #aipw pa
    data$Y0.hat <- data$Y0.hat.pa
    data$Y1.hat <- data$Y1.hat.pa
    data$ps <- data$ps.pa
    ate10 <- ifelse(out.mod$converged==TRUE & ps.mod$converged==TRUE,pi_dr(mydata = data)$pi_dr,NA) #Doubly robust estimator: AIPTW
    #se10 <- 
    
    convergence <- ifelse(out.mod$converged==TRUE & ps.mod$converged==TRUE,1,0)} else {ate9 <- ate10 <- NA; se9 <- se10 <- NA; convergence <- 0}
    m9 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate9,
      #se = se9,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 9
    )
  
    m10 <- data.frame(
      i = i,
      dgm = dgm,
      true.ate = dt$true.rd[dt$dgm==dgm],
      ate = ate10,
      #se = se10,
      convergence = convergence,
      trt_prev = trt_prev,
      out_rate_0 = out_rate, #AK: maybe intra-cluster correlation calculations could be added here!
      method = 10
    )
    
  }#)
  output <- rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
  #return(list(output,"beta.vec1" = beta.vec1,"beta.vec2"=beta.vec2,"beta.vec3"=beta.vec3,"beta.vec4"=beta.vec4))
  #return(list("beta.vec1" = beta.vec1,"beta.vec2"=beta.vec2,"beta.vec3"=beta.vec3,"beta.vec4"=beta.vec4))
  
}

#quick checks
dfs <- readRDS("H:/My Documents/Simulation ISCB44/New Basic/Results/dfs.rds")
dgm=33
data <- as.data.frame(dfs$dgm33[1])
B=10
out <- NULL
set.seed(2711*33)
foreach(i=1:B) %do% {
  print(i)
   out <- rbind(out,analysis_bin(i,dt,33,data=dfs[[33]][[i]]))
   return(out)
}
system.time(t <- analysis_bin(1,dt,33,data))

#small run to estimate sets of initial values for each analysis model within each dgm
B <- 1
k=2
dgm=1:k
output <- NULL
 system.time(t <- foreach(j=dgm) %do% {
foreach(i=1:B) %do% {
  
  #print(i)
  output <- list(output,analysis_bin(i=i,dt=dt,dgm=j,data = ival.dfs[[j]][[i]]))

  return(output)
  }
})