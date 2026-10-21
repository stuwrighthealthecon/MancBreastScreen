#Calibration runs: each run draws new calibrated parameters (lc_bump and the
#genesis-age distribution) and compares predicted with observed outcomes
library(here)
library(tidyverse)
library(doParallel)
library(MASS)
library(dqrng)
library(compiler)
library(matrixStats)

for (f in list.files(here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

calibruns<-50

calibresults<-as.data.frame(matrix(nrow=calibruns,ncol=6))
colnames(calibresults)<-c("c","Mean","SD","P","LCbump","RMSE")

plotlist<-list()

for (c in 1:calibruns){

  #A different seed for each run gives different calibrated parameter draws
  results <- run_model(default_controls(seed = c))
  det_output_path <- results$det_output_path
  source(here("analysis", "validation.R"))
  calibresults[c,]<-c(c,results$params$gamma_mean,results$params$gamma_sd,
                      results$params$gamma_p,results$params$lc_bump,RMSE)

plotlist[[c]]<-incidenceplot
print(paste(c," of ",calibruns," complete"))
}

calibresults %>% arrange(RMSE)

plotlist[[20]]

load(here("output", "risksample", "possample.Rdata"))
lc_bump <- results$params$lc_bump
plot(density(risksample$liferisk),main="Distribution of Lifetime Breast Cancer Risk at Birth")
count(((risksample$liferisk-2.2-lc_bump)>40))/nrow(risksample)
