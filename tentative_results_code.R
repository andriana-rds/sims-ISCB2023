#tentative code for results visualization
setwd("H:/My Documents/Simulation ISCB44/New Basic/Results")
dt <- read.csv("dgm_index.csv",header = TRUE)
library(foreach)

dat <- NULL
#n_dgm <- 288
foreach (i=c(1:2,5:11,13,17,18,21:30,33:35,40:43,45,46,49,53:67,71:78,81:84,86:93,95,97:108,110:112,114:118,120,122,123,126:133,
             135:138,140,142:147,149:157,
             159:163,165:181,183:185,187:192,194:196,198:200,203:204,207:208,211:213,215:216,219:220,223,225,227,228,229,
             230,231,232,233,235,236,238,243,240,242,
             244,246,247,248,251,252,255,256,257,258,260,263,272,275,265,267,268,273,276,279,283,284,287,288)) %do% {
  
  dat <- rbind(dat, read.csv(paste0("output_","dgm",i,".csv"), header = TRUE, sep = ","))
  return(dat)
             }
# nrow(dat)
# length(unique(dat$dgm))
colSums(is.na(dat))
prop.table(table(dat$convergence)) #1% of models did not converge (0: did not converge)

#merge the data set with parameters description with dat:
dat <- merge(dat,dt, by = "dgm")

#Bias calculation:
results_bias <- NULL

foreach(i=c(1:2,5:11,13,17,18,21:30,33:35,40:43,45,46,49,53:67,71:78,81:84,86:93,95,97:108,110:112,114:118,120,122,123,126:133,
            135:138,140,142:147,149:157,
            159:163,165:181,183:185,187:192,194:196,198:200,203:204,207:208,211:213,215:216,219:220,223,225,227,228,229,
            230,231,232,233,235,236,238,243,240,242,
            244,246,247,248,251,252,255,256,257,258,260,263,272,275,265,267,268,273,276,279,283,284,287,288)) %do% {
              foreach(j=1:10) %do% {
                
                bias <- mean(dat$ate[dat$dgm==i& dat$method==j & dat$convergence==1] - 
                               dat$true.ate[dat$dgm==i& dat$method==j & dat$convergence==1])
                
                results_bias <- rbind(results_bias, data.frame(bias=bias,dgm=i,method=j))
                return(results_bias)
                
              }
            }

#decide on parameters that affect bias values the most (inspired by Rubin)
results_bias <- merge(results_bias, dt, by = "dgm")
mod <- lm(bias ~ method + alpha0 + beta0 + rho + rr + trt_mod + out_mod + NCluster, data = results_bias)
summary(mod) #see which parameters are s.s. to include in our graphs
#AK note1: interestingly, the dgm of the outcome model results s.s. only at 10%...true treatment prev is not s.s. either. Number of clusters is not s.s. 
#AK note2: NCluster is significant.

library(broom)
RES <- tidy(mod) # get coefficient table as a data frame
glance(mod) # get rest of stats as a data frame
glance(mod)$p.value # get p value...

#Tables with results from bias regressed on all varying parameters of the simulation design
if (!requireNamespace("xtable")) install.packages("xtable")
library(xtable) #for the .tex format


tabtex <- xtable(cbind(RES[,1],RES[,2],RES[,3],RES[,5]), 
                 label = paste0("bias_reg"), digits = rep(3, 5), align = c("l",rep("p{1.5cm}", 4)), auto = T)
tabtex

#we repeat the same test: this time among the benchmark and g-computation methods only
results_bias_g <- results_bias[results_bias$method==1|results_bias$method==2|results_bias$method==3|results_bias$method==4,]
mod_g <- lm(bias ~ method + alpha0 + beta0 + rho + rr + trt_mod + out_mod + NCluster, data = results_bias_g)
summary(mod_g) #see which parameters are s.s. to include in our graphs
#AK note3: within g-comp methods, NCluster becomes even more significant, whereas within IPW estimators, it is no longer significant. TO DISCUSS. 
RES <- tidy(mod_g) # get coefficient table as a data frame

tabtex <- xtable(cbind(RES[,1],RES[,2],RES[,3],RES[,5]), 
                 label = paste0("bias_reg"), digits = rep(3, 5), align = c("l",rep("p{1.5cm}", 4)), auto = T)
tabtex


#include only ipw methods
results_bias_ipw <- results_bias[results_bias$method==5|results_bias$method==6|results_bias$method==7|results_bias$method==8,]
mod_ipw <- lm(bias ~ method + alpha0 + beta0 + rho + rr + trt_mod + out_mod + NCluster, data = results_bias_ipw)
summary(mod_ipw) #correlation between V and U1 is significant again - intuitively, this might be due to the fact that under corr, we have a modelling assumption of the PS analysis 
#model being violated...


RES <- tidy(mod_ipw) # get coefficient table as a data frame

tabtex <- xtable(cbind(RES[,1],RES[,2],RES[,3],RES[,5]), 
                 label = paste0("bias_reg"), digits = rep(3, 5), align = c("l",rep("p{1.5cm}", 4)), auto = T)
tabtex

#include only dr methods
results_bias_dr <- results_bias[results_bias$method==9|results_bias$method==10,]
mod_dr <- lm(bias ~ method + alpha0 + beta0 + rho + rr + trt_mod + out_mod + NCluster, data = results_bias_dr)
summary(mod_dr) #AK note4: for AIPTW estimators, NCluster is significant at .1.

RES <- tidy(mod_dr) # get coefficient table as a data frame

tabtex <- xtable(cbind(RES[,1],RES[,2],RES[,3],RES[,5]), 
                 label = paste0("bias_reg"), digits = rep(3, 5), align = c("l",rep("p{1.5cm}", 4)), auto = T)
tabtex

#Plots prep
#results_bias$prettyconfig <- ifelse(results_bias$prev=="trt prev 20%" & results_bias$rate=="out control 10%","p(Z)=20%, p(Y0)=10%","p(Z)=50%, p(Y0)=10%")
#results_bias$prettyconfig <- ifelse(results_bias$prev=="trt prev 20%" & results_bias$rate=="out control 30%",
                   #                 "p(Z)=20%, p(Y0)=30%",results_bias$prettyconfig)
#results_bias$prettyconfig <- ifelse(results_bias$prev=="trt prev 50%" & results_bias$rate=="out control 30%",
 #                                   "p(Z)=50%, p(Y0)=30%",results_bias$prettyconfig)
results_bias$prettyconfig <- ifelse(results_bias$prev=="trt prev 20%" & results_bias$rate=="out control 10%",
"list(pi[Z]==20*\"%\",pi[Y0]==10*\"%\")","list(pi[Z]==20*\"%\",pi[Y0]==30*\"%\")")
results_bias$prettyconfig <- ifelse(results_bias$prev=="trt prev 50%" & results_bias$rate=="out control 10%",
                                    "list(pi[Z]==50*\"%\",pi[Y0]==10*\"%\")",results_bias$prettyconfig)
results_bias$prettyconfig <- ifelse(results_bias$prev=="trt prev 50%" & results_bias$rate=="out control 30%",
                                    "list(pi[Z]==50*\"%\",pi[Y0]==30*\"%\")",results_bias$prettyconfig)
#results_bias$method <- factor(results_bias$method, labels = c("Bench EBE","Bench PA","G-comp EBE","G-comp PA","IPTW mar EBE","IPTW mar PA","IPTW cl EBE","IPTW cl PA","AIPTW EBE","AIPTW PA"))

#Plots
library(ggplot2)
library(ggstance)

#1. DGMS: 50 clusters of 20 patients (on average), random intercept models for trt and outcome
results_bias_riri <- results_bias[results_bias$design=="50 clusters of 20 - unbl" 
                               & results_bias$trt_mod=="random intercept" & results_bias$out_mod=="random intercept",]

p_rd_bias_riri <- ggplot(data = results_bias_riri, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 
#supremum and infimum for the bias values (so that lines on the graph don't have the same width) 
ggtitle('Random intercept models for treatment and outcome')
#pdf(paste("H:/My Documents/Simulation ISCB44/New Basic/Results/Figures/", "bias1.pdf", sep = ""), width = 12, height = 7)
#print(p_rd_bias_1)
#dev.off()
#2. DGMS: 50 clusters of 20 patients, random intercept & slope models for trt and 
#random intercept for outcome
results_bias_rirs <- results_bias[results_bias$design=="50 clusters of 20 - unbl" 
                               & results_bias$trt_mod=="random intercept" & results_bias$out_mod=="random slope",]

p_rd_bias_rirs <- ggplot(data = results_bias_rirs, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 


#3. DGMS: 50 clusters of 20 patients, random intercept and random slopes models for trt and random intercept models for outcome

results_bias_rsri <- results_bias[results_bias$design=="50 clusters of 20 - unbl" 
                               & results_bias$trt_mod=="random slope" & results_bias$out_mod=="random intercept",]
p_rd_bias_rsri <- ggplot(data = results_bias_rsri, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 


#4. DGMS: 50 clusters of 20 patients, random intercept and slope models for trt and outcome
results_bias_rsrs <- results_bias[results_bias$design=="50 clusters of 20 - unbl" 
                               & results_bias$trt_mod=="random slope" & results_bias$out_mod=="random slope",]

p_rd_bias_rsrs <- ggplot(data = results_bias_rsrs, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 

#AK note 5: consider the fact that we have a low number of dgms corresponding to this specific scenario: we end-up
# with missing dgms for cases, e.g.: no correlation; p(Z)=0.2; p(Y0)=0.3; r.s both for trt and out...Hence, we cannot 
#really check whether the trend is consistent or not for this particular case yet.

#5. DGMS: 50 clusters of 20 patients, single-level for trt and out
results_bias_ss <- results_bias[results_bias$design=="50 clusters of 20 - unbl" 
                               & results_bias$trt_mod=="none" & results_bias$out_mod=="none",]

p_rd_bias_ss <- ggplot(data = results_bias_ss, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 

#AK note6: see note 5

#6. DGMS: 50 clusters of 20 patients, single-level for trt, random intercept for out
results_bias_sri <- results_bias[results_bias$design=="50 clusters of 20 - unbl" 
                                & results_bias$trt_mod=="none" & results_bias$out_mod=="random intercept",]

p_rd_bias_sri <- ggplot(data = results_bias_sri, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define


#7. DGMS: 50 clusters of 20 patients, single-level for trt, random intercept and slope for out
results_bias_srs <- results_bias[results_bias$design=="50 clusters of 20 - unbl" 
                                 & results_bias$trt_mod=="none" & results_bias$out_mod=="random slope",]

p_rd_bias_srs <- ggplot(data = results_bias_srs, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define

#AK note7: as in 5,6, but for p(Z)=0.2; p(Y0)=0.1




#8. DGMS: 50 clusters of 20 patients, trt random intercept and single-level for out
results_bias_ris <- results_bias[results_bias$design=="50 clusters of 20 - unbl" 
                                 & results_bias$trt_mod=="random intercept" & results_bias$out_mod=="none",]

p_rd_bias_ris <- ggplot(data = results_bias_ris, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define


#9. DGMS: 50 clusters of 20 patients, trt random intercept slope and single-level for out
results_bias_rss <- results_bias[results_bias$design=="50 clusters of 20 - unbl" 
                                 & results_bias$trt_mod=="random slope" & results_bias$out_mod=="none",]

p_rd_bias_rss <- ggplot(data = results_bias_rss, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define

# Combine the above plots into one 3x3 multiplot (for 50 clusters of 20)
library(patchwork)
p_rd_bias_riri + p_rd_bias_rirs 
p_rd_bias_rsri + p_rd_bias_rsrs 
p_rd_bias_ss + p_rd_bias_sri 
p_rd_bias_srs + p_rd_bias_ris 
p_rd_bias_rss



#Do the same for 100 clusters of 10 and compare
#1. DGMS: 100 clusters of 10 patients (on average), random intercept models for trt and outcome
results_bias_riri_100 <- results_bias[results_bias$design=="100 clusters of 10 - unbl" 
                                  & results_bias$trt_mod=="random intercept" & results_bias$out_mod=="random intercept",]

p_rd_bias_riri_100 <- ggplot(data = results_bias_riri_100, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 
#supremum and infimum for the bias values (so that lines on the graph don't have the same width) 
ggtitle('Random intercept models for treatment and outcome')
#pdf(paste("H:/My Documents/Simulation ISCB44/New Basic/Results/Figures/", "bias1.pdf", sep = ""), width = 12, height = 7)
#print(p_rd_bias_1)
#dev.off()

#2. DGMS: 100 clusters of 10 patients, random intercept & slope models for trt and 
#random intercept for outcome
results_bias_rirs_100 <- results_bias[results_bias$design=="100 clusters of 10 - unbl" 
                                  & results_bias$trt_mod=="random intercept" & results_bias$out_mod=="random slope",]

p_rd_bias_rirs_100 <- ggplot(data = results_bias_rirs_100, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 


#3. DGMS: 100 clusters of 10 patients, random intercept and random slopes models for trt and random intercept models for outcome

results_bias_rsri_100 <- results_bias[results_bias$design=="100 clusters of 10 - unbl" 
                                  & results_bias$trt_mod=="random slope" & results_bias$out_mod=="random intercept",]
p_rd_bias_rsri_100 <- ggplot(data = results_bias_rsri_100, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 

#4. DGMS: 100 clusters of 10 patients, random intercept and slope models for trt and outcome
results_bias_rsrs_100 <- results_bias[results_bias$design=="100 clusters of 10 - unbl" 
                                  & results_bias$trt_mod=="random slope" & results_bias$out_mod=="random slope",]

p_rd_bias_rsrs_100 <- ggplot(data = results_bias_rsrs_100, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 

#5. DGMS: 100 clusters of 10 patients, single-level for trt and out
results_bias_ss_100 <- results_bias[results_bias$design=="100 clusters of 10 - unbl" 
                                & results_bias$trt_mod=="none" & results_bias$out_mod=="none",]

p_rd_bias_ss_100 <- ggplot(data = results_bias_ss_100, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 

#6. DGMS: 100 clusters of 10 patients, single-level for trt, random intercept for out
results_bias_sri_100 <- results_bias[results_bias$design=="100 clusters of 10 - unbl" 
                                 & results_bias$trt_mod=="none" & results_bias$out_mod=="random intercept",]

p_rd_bias_sri_100 <- ggplot(data = results_bias_sri_100, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define

#7. DGMS: 100 clusters of 10 patients, single-level for trt, random intercept and slope for out
results_bias_srs_100 <- results_bias[results_bias$design=="100 clusters of 10 - unbl" 
                                 & results_bias$trt_mod=="none" & results_bias$out_mod=="random slope",]

p_rd_bias_srs_100 <- ggplot(data = results_bias_srs_100, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define

#8. DGMS: 100 clusters of 10 patients, trt random intercept and single-level for out
results_bias_ris_100 <- results_bias[results_bias$design=="100 clusters of 10 - unbl" 
                                 & results_bias$trt_mod=="random intercept" & results_bias$out_mod=="none",]

p_rd_bias_ris_100 <- ggplot(data = results_bias_ris_100, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define

#9. DGMS: 100 clusters of 10 patients, trt random intercept slope and single-level for out
results_bias_rss_100 <- results_bias[results_bias$design=="100 clusters of 10 - unbl" 
                                 & results_bias$trt_mod=="random slope" & results_bias$out_mod=="none",]

p_rd_bias_rss_100 <- ggplot(data = results_bias_rss_100, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define

p_rd_bias_riri_100 + p_rd_bias_rirs_100 
p_rd_bias_rsri_100 + p_rd_bias_rsrs_100 
p_rd_bias_ss_100 + p_rd_bias_sri_100 
p_rd_bias_srs_100 + p_rd_bias_ris_100 
p_rd_bias_rss_100


###########################################################################################################
#Conclusion: results are similar across the two designs (i.e., 50 clusters of 20 vs. 100 clusters of 10); 
#We demonstrate results for the latter only
# 50 clusters of 20; 50% trt prev; moderate true ATE; single-level, random intercept or random slope out model 
###########################################################################################################

#trt random intercept; out single-level (correlation between V and U1 plays no role)???Correction: we don't have dgms for both categories of correlation
results_bias_ris_50_final <- results_bias[results_bias$design=="50 clusters of 20 - unbl" & results_bias$prev=="trt prev 50%"&
                                      results_bias$trt_mod=="random intercept" & results_bias$out_mod=="none"
                                      & results_bias$true.effect=="moderate effect",]

p_rd_bias_ris_50_final <- ggplot(data = results_bias_ris_50_final, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 
#supremum and infimum for the bias values (so that lines on the graph don't have the same width) 
#ggtitle('Random intercept models for treatment and outcome')
#pdf(paste("H:/My Documents/Simulation ISCB44/New Basic/Results/Figures/", "bias1.pdf", sep = ""), width = 12, height = 7)
#print(p_rd_bias_1)
#dev.off()


#trt random intercept; out random intercept 
results_bias_riri_50_final <- results_bias[results_bias$design=="50 clusters of 20 - unbl" & results_bias$prev=="trt prev 50%"&
                                            results_bias$trt_mod=="random intercept" & results_bias$out_mod=="random intercept"
                                          & results_bias$true.effect=="moderate effect",]

p_rd_bias_riri_50_final <- ggplot(data = results_bias_riri_50_final, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 


#trt random intercept; out random slope 
results_bias_rirs_50_final <- results_bias[results_bias$design=="50 clusters of 20 - unbl" & results_bias$prev=="trt prev 50%"&
                                             results_bias$trt_mod=="random intercept" & results_bias$out_mod=="random slope"
                                           & results_bias$true.effect=="moderate effect",]

p_rd_bias_rirs_50_final <- ggplot(data = results_bias_rirs_50_final, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 


#trt random slope; out single-level 
results_bias_rss_50_final <- results_bias[results_bias$design=="50 clusters of 20 - unbl" & results_bias$prev=="trt prev 50%"&
                                            results_bias$trt_mod=="random slope" & results_bias$out_mod=="none"
                                          & results_bias$true.effect=="moderate effect",]

p_rd_bias_rss_50_final <- ggplot(data = results_bias_rss_50_final, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 



#trt random slope; out random intercept 
results_bias_rsri_50_final <- results_bias[results_bias$design=="50 clusters of 20 - unbl" & results_bias$prev=="trt prev 50%"&
                                             results_bias$trt_mod=="random slope" & results_bias$out_mod=="random intercept"
                                           & results_bias$true.effect=="moderate effect",]

p_rd_bias_rsri_50_final <- ggplot(data = results_bias_rsri_50_final, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 


#trt random slope; out random slope 
results_bias_rsrs_50_final <- results_bias[results_bias$design=="50 clusters of 20 - unbl" & results_bias$prev=="trt prev 50%"&
                                             results_bias$trt_mod=="random slope" & results_bias$out_mod=="random slope"
                                           & results_bias$true.effect=="moderate effect",]

p_rd_bias_rsrs_50_final <- ggplot(data = results_bias_rsrs_50_final, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 

#dgms only for no correlation and p(Y0=1)=0.3 available here...

#trt single-level; out single-level
results_bias_ss_50_final <- results_bias[results_bias$design=="50 clusters of 20 - unbl" & results_bias$prev=="trt prev 50%"&
                                             results_bias$trt_mod=="none" & results_bias$out_mod=="none"
                                           & results_bias$true.effect=="moderate effect",]

p_rd_bias_ss_50_final <- ggplot(data = results_bias_ss_50_final, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 

#trt single-level; out random intercept
results_bias_sri_50_final <- results_bias[results_bias$design=="50 clusters of 20 - unbl" & results_bias$prev=="trt prev 50%"&
                                           results_bias$trt_mod=="none" & results_bias$out_mod=="random intercept"
                                         & results_bias$true.effect=="moderate effect",]

p_rd_bias_sri_50_final <- ggplot(data = results_bias_sri_50_final, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 

#trt single-level; out random slope
results_bias_srs_50_final <- results_bias[results_bias$design=="50 clusters of 20 - unbl" & results_bias$prev=="trt prev 50%"&
                                            results_bias$trt_mod=="none" & results_bias$out_mod=="random slope"
                                          & results_bias$true.effect=="moderate effect",]

p_rd_bias_srs_50_final <- ggplot(data = results_bias_srs_50_final, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 

p_rd_bias_srs_50_final/p_rd_bias_rirs_50_final
p_rd_bias_rsrs_50_final





#####FINAL (WE HOPE) ################################################### 
#trt random intercept; out random slope 
results_bias_rirs_50_20_50_10_final <- results_bias[results_bias$design=="50 clusters of 20 - unbl" & results_bias$prev=="trt prev 50%"&
                                             results_bias$rate=="out control 10%"&
                                             results_bias$trt_mod=="random intercept" & results_bias$out_mod=="random slope"
                                           & results_bias$true.effect=="moderate effect",]

p_rd_bias_rirs_50_20_50_10_final <- ggplot(data = results_bias_rirs_50_20_50_10_final, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrange(aes(xmin = -0.03, xmax = 0.1), size = 0.3, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 

p_rd_bias_rirs_50_20_50_10_final


results_bias_rirs_50_20_50_30_final <- results_bias[results_bias$design=="50 clusters of 20 - unbl" & results_bias$prev=="trt prev 50%"&
                                             results_bias$rate=="out control 30%"&
                                             results_bias$trt_mod=="random intercept" & results_bias$out_mod=="random slope"
                                           & results_bias$true.effect=="moderate effect",]

p_rd_bias_rirs_50_20_50_30_final <- ggplot(data = results_bias_rirs_50_20_50_30_final, aes(x = bias, y = true.effect, color = correlation, group = correlation)) +
  geom_pointrangeh(aes(xmin = -0.05, xmax = 0.2), size = 0.2, position = position_dodge(0.3)) + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_grid(method ~ prettyconfig, labeller = label_parsed) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), strip.background = element_rect(fill = "grey60", size=0),
        strip.text = element_text(colour = "white")) + 
  xlab("Bias") +
  theme(legend.position = "bottom", legend.key = element_rect(fill = "white")) + 
  scale_color_discrete("") +
  ylab("True Average Treatment effect - risk difference") #AK: to define 

p_rd_bias_rirs_50_20_50_30_final


#IDEA I: we need the descriptives of bias; we apply library rsimsum
library(rsimsum)
#data: we use dat - but only results for null ate in the first place
dat0 <- dat[dat$true.ate==0 & dat$trt_mod=="random slope" & dat$out_mod=="random intercept" &
              dat$prev=="trt prev 50%"&
              dat$rate=="out control 10%",]
dat0$method <- factor(dat0$method, levels = c(1:10), labels = c("bench ebe","bench pa",
        "g-comp ebe","g-comp pa","iptw-mar ebe","iptw-mar pa","iptw-cl ebe",
        "iptw-cl pa","aiptw ebe","aiptw pa"))
s <- simsum(data = dat0, estvarname = "ate",true = 0, 
            methodvar = "method", by = "dgm", ref = "bench ebe", x = TRUE)
print(summary(s, stats = c("bias")))

library(dplyr)
get_data(summary(s)) %>%
  filter(stat == "bias") %>%
  ggplot(aes(x = method, y = est)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 1 / 2) +
  geom_point() +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~dgm, label = label_both) +
  labs(x = "Method", y = "Bias")

autoplot(summary(s), type = "lolly", stats = "bias") + ggtitle("Bias: 50 clusters of 20")

autoplot(s, type = "forest", stats = "bias")
#autoplot(s, type = "est", stats = "bias")
#autoplot(s, type = "est_ba", stats = "bias")




################################################################################################
#TO PRINT
library(rsimsum)
#data: we use dat - but only results for null ate in the first place

# trt prev 50%; out control rate 10%
dat0_50_10 <- dat[dat$true.ate==0 & dat$trt_mod=="random slope" & dat$out_mod=="random intercept" &
              dat$prev=="trt prev 50%"&
              dat$rate=="out control 10%" & dat$design=="100 clusters of 10 - unbl",]
dat0_50_10$method <- factor(dat0_50_10$method, levels = c(1:10), labels = c("bench ebe","bench pa",
                                                                "g-comp ebe","g-comp pa","iptw-mar ebe","iptw-mar pa","iptw-cl ebe",
                                                                "iptw-cl pa","aiptw ebe","aiptw pa"))
dat0_50_10$dgm <- ifelse(dat0_50_10$dgm==170,"no correlation","correlation")
s <- simsum(data = dat0_50_10, estvarname = "ate",true = 0, 
            methodvar = "method", by = "dgm", ref = "bench ebe", x = TRUE)
print(summary(s, stats = c("bias")))

autoplot(summary(s), type = "lolly", stats = "bias") + ggtitle("100 clusters of 10 patients (on average): \n 
Pr(Z=1)=50%, Pr(Y0=1)=10%, null ATE
  \n From left to right hand side: correlation and no correlation between V and U1.")


# trt prev 50%; out control rate 10%
dat0_50_20 <- dat[dat$true.ate==0 & dat$trt_mod=="random slope" & dat$out_mod=="random intercept" &
                    dat$prev=="trt prev 50%"&
                    dat$rate=="out control 10%" & dat$design=="50 clusters of 20 - unbl",]
dat0_50_20$method <- factor(dat0_50_20$method, levels = c(1:10), labels = c("bench ebe","bench pa",
                                                                            "g-comp ebe","g-comp pa","iptw-mar ebe","iptw-mar pa","iptw-cl ebe",
                                                                            "iptw-cl pa","aiptw ebe","aiptw pa"))
dat0_50_20$dgm <- ifelse(dat0_50_20$dgm==162,"no correlation","correlation")
s <- simsum(data = dat0_50_20, estvarname = "ate",true = 0, 
            methodvar = "method", by = "dgm", ref = "bench ebe", x = TRUE)
print(summary(s, stats = c("bias")))

autoplot(summary(s), type = "lolly", stats = "bias") + ggtitle("50 clusters of 20 patients (on average): \n 
Pr(Z=1)=50%, Pr(Y0=1)=10%, null ATE
  \n From left to right hand side: correlation and no correlation between V and U1.")

# Nested loop plot:
#data("nlp", package = "rsimsum")
dat$method <- factor(dat$method, levels = c(1:10), labels = c("bench ebe","bench pa",
                                                                     "g-comp ebe","g-comp pa","iptw-mar ebe","iptw-mar pa","iptw-cl ebe",
                                                                     "iptw-cl pa","aiptw ebe","aiptw pa"))

# dat0 <- dat[dat$true.rd==0,]
# s0 <- simsum(
  # data = dat0, estvarname = "ate", true = 0,
  # methodvar = "method", ref = "bench ebe", by = c("prev", "rate", 
  # "correlation","NCluster","trt_mod","out_mod") this does not work for now, as we have missing combinations
# )

dat0_none <- dat[dat$true.rd==0 &dat$out_mod=="none",]
dat0_r_inter <- dat[dat$true.rd==0 &dat$out_mod=="random intercept",]
dat0_r_sl <- dat[dat$true.rd==0 &dat$out_mod=="random slope",]


dat6 <- dat[dat$true.rd==-0.06,]

s0_none <- simsum(
  data = dat0_none, estvarname = "ate", true = 0, #we need to specify the different true ATEs
  methodvar = "method", ref = "bench ebe", by = "dgm"
  #by = c("prev", "rate", "correlation","NCluster","true.effect","trt_mod","out_mod")
)

s0_r_inter <- simsum(
  data = dat0_r_inter, estvarname = "ate", true = 0, #we need to specify the different true ATEs
  methodvar = "method", ref = "bench ebe", by = "dgm"
  #by = c("prev", "rate", "correlation","NCluster","true.effect","trt_mod","out_mod")
)

s0_r_sl <- simsum(
  data = dat0_r_sl, estvarname = "ate", true = 0, #we need to specify the different true ATEs
  methodvar = "method", ref = "bench ebe", by = "dgm"
  #by = c("prev", "rate", "correlation","NCluster","true.effect","trt_mod","out_mod")
)
#> 'ref' method was not specified, 1 set as the reference
autoplot(s0_none, stats = "bias", type = "nlp")
p <- autoplot(s0_r_inter, stats = "bias", type = "nlp")
p
# p + scale_x_continuous(limits = c(0,90))
autoplot(s0_r_sl, stats = "bias", type = "nlp")


# dat2 <- dat[dat$true.rd==-0.02,]
dat2_none <- dat[dat$true.rd==-0.02 &dat$out_mod=="none",]
dat2_r_inter <- dat[dat$true.rd==-0.02 &dat$out_mod=="random intercept",]
dat2_r_sl <- dat[dat$true.rd==-0.02 &dat$out_mod=="random slope",]

s2_none <- simsum(
  data = dat2_none, estvarname = "ate", true =-0.02, #we need to specify the different true ATEs
  methodvar = "method", ref = "bench ebe", by = "dgm"
  #by = c("prev", "rate", "correlation","NCluster","true.effect","trt_mod","out_mod")
)

s2_r_inter <- simsum(
  data = dat2_r_inter, estvarname = "ate", true = -0.02, #we need to specify the different true ATEs
  methodvar = "method", ref = "bench ebe", by = "dgm"
  #by = c("prev", "rate", "correlation","NCluster","true.effect","trt_mod","out_mod")
)

s2_r_sl <- simsum(
  data = dat2_r_sl, estvarname = "ate", true = -0.02, #we need to specify the different true ATEs
  methodvar = "method", ref = "bench ebe", by = "dgm"
  #by = c("prev", "rate", "correlation","NCluster","true.effect","trt_mod","out_mod")
)

autoplot(s2_none, stats = "bias", type = "nlp")
autoplot(s2_r_inter, stats = "bias", type = "nlp")
autoplot(s2_r_sl, stats = "bias", type = "nlp")

dat6_none <- dat[dat$true.rd==-0.06 &dat$out_mod=="none",]
dat6_r_inter <- dat[dat$true.rd==-0.06 &dat$out_mod=="random intercept",]
dat6_r_sl <- dat[dat$true.rd==-0.06 &dat$out_mod=="random slope",]

s6_none <- simsum(
  data = dat6_none, estvarname = "ate", true =-0.06, #we need to specify the different true ATEs
  methodvar = "method", ref = "bench ebe", by = "dgm"
  #by = c("prev", "rate", "correlation","NCluster","true.effect","trt_mod","out_mod")
)

s6_r_inter <- simsum(
  data = dat6_r_inter, estvarname = "ate", true = -0.06, #we need to specify the different true ATEs
  methodvar = "method", ref = "bench ebe", by = "dgm"
  #by = c("prev", "rate", "correlation","NCluster","true.effect","trt_mod","out_mod")
)

s6_r_sl <- simsum(
  data = dat6_r_sl, estvarname = "ate", true = -0.06, #we need to specify the different true ATEs
  methodvar = "method", ref = "bench ebe", by = "dgm"
  #by = c("prev", "rate", "correlation","NCluster","true.effect","trt_mod","out_mod")
)

autoplot(s6_none, stats = "bias", type = "nlp")
autoplot(s6_r_inter, stats = "bias", type = "nlp")
autoplot(s6_r_sl, stats = "bias", type = "nlp")






s6 <- simsum(
  data = dat6, estvarname = "ate", true = -0.06, #we need to specify the different true ATEs
  methodvar = "method", ref = "bench ebe", by = "dgm"
  #by = c("prev", "rate", "correlation","NCluster","true.effect","trt_mod","out_mod")
)
#> 'ref' method was not specified, 1 set as the reference
autoplot(s6, stats = "bias", type = "nlp")



