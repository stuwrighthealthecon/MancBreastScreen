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

calibresults<-as.data.frame(matrix(nrow=calibruns,ncol=8))
colnames(calibresults)<-c("c","Mean","SD","P","LCbump","RMSE","cdr","propscreen")

#Create plots of outcomes
plotlist<-list()

#Draw values for calibration runs
set.seed(3)
draws <- data.frame(
  lc_bump       = runif(calibruns, 1.15, 1.5),
  gamma_mean    = runif(calibruns, 65, 74),
  gamma_sd      = runif(calibruns, 15, 28),
  gamma_p       = runif(calibruns, 2.2, 8)   # variance of log growth rate
)

base_params <- get_parameters()                # calculated once
controls    <- default_controls(desired_cases = 20000, strategies = 3, n_cores = 1)

pdf(NULL)                                      # don't keep plots in RStudio during the loop
for (c in 1:calibruns) {
  params <- base_params
  params$lc_bump       <- draws$lc_bump[c]
  params$gamma_mean    <- draws$gamma_mean[c]
  params$gamma_sd      <- draws$gamma_sd[c]
  params$gamma_p       <- draws$gamma_p[c]
  
  gg <- tryCatch(get_gengamma_params(params$gamma_mean - 18, params$gamma_sd, params$gamma_p),
                 error = function(e) NULL)
  if (is.null(gg)) { message("Run ", c, ": no generalised gamma fit, skipped"); next }
  params$gengammaparams <- gg
  
  results <- run_model(controls, params)
  det_output_path <- results$det_output_path
  source(here("analysis", "validation.R"))
  
  calibresults[c, ] <- c(c, unlist(draws[c, c("gamma_mean", "gamma_sd", "gamma_p", "lc_bump")]),
                         RMSE, cdr, propscreen)
  plotlist[[c]] <- incidenceplot
  rm(results, alldata, detresults); gc()
  print(paste(c, "of", calibruns, "complete"))
}
dev.off()

#Arrange by fit
calibresults %>% arrange(RMSE)



