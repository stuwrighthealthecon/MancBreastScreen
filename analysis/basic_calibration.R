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

#Set number of runs
calibruns<-30

calibresults<-as.data.frame(matrix(nrow=calibruns,ncol=6))
colnames(calibresults)<-c("c","Mean","SD","P","LCbump","RMSE")

#Create plots of outcomes
plotlist<-list()

#Draw values for calibration runs
set.seed(3)   # seed for the calibration draws only
draws <- data.frame(
  lc_bump    = runif(calibruns, 0, 10),
  gamma_mean = runif(calibruns, 60, 90),
  gamma_sd   = runif(calibruns, 10, 30),
  gamma_p    = runif(calibruns, 0.1, 10)
)

#Assign value for run, run model, record results
for (c in 1:calibruns) {
  params <- get_parameters()
  params$lc_bump    <- draws$lc_bump[c]
  params$gamma_mean <- draws$gamma_mean[c]
  params$gamma_sd   <- draws$gamma_sd[c]
  params$gamma_p    <- draws$gamma_p[c]
  params$gengammaparams <- get_gengamma_params(params$gamma_mean - 18,
                                               params$gamma_sd, params$gamma_p)
  
  controls <- default_controls(desired_cases = 20000, strategies = 3)
  results  <- run_model(controls, params)

  det_output_path <- results$det_output_path
  source(here("analysis", "validation.R"))
  calibresults[c,]<-c(c,results$params$gamma_mean,results$params$gamma_sd,
                      results$params$gamma_p,results$params$lc_bump,RMSE)

plotlist[[c]]<-incidenceplot
print(paste(c," of ",calibruns," complete"))
}

#Arrange by fit
calibresults %>% arrange(RMSE)




