library("mgcv")
library("dampack")
library(gt)
library(tidyverse)
library(magrittr)
library("MASS")
library("patchwork")

#Set number of PSA runs to estimate over
mcruns <- 100000

#Set alternatives
alternative <- c(0, 1, 2, 3, 4, 9)

#Set wtp thresholds to estimate over
wtp <- seq(from = 0, to = 50000, by = 1000)

#Load GAMs
modQ <- readRDS(here::here("output", "gam_models", "QALYmodelcan.RDS"))
modQ2 <- readRDS(here::here("output", "gam_models", "QALYmodelnocan.RDS"))
modC <- readRDS(here::here("output", "gam_models", "costmodelcan.RDS"))
modC2 <- readRDS(here::here("output", "gam_models", "costmodelnocan.RDS"))

#Draw stage I to III survival parameters
survmvn <- data.frame(
  c(-5.6178, -5.2857, -6.1138),
  c(-3.8083, -3.7400, -3.8704),
  c(-2.7311, -2.6687, -2.7963)
)
survcovmat <- cov(survmvn)
survmeans <- c(survmvn[1, 1], survmvn[1, 2], survmvn[1, 3])
PSA_gamma_survival <- mvrnorm(mcruns, survmeans, survcovmat)

# Draw Metatstatic survival parameters
metmvn <- data.frame(
  c(-1.78723, -1.67922, -1.89434),
  c(-1.38762, -1.33512, -1.49956),
  c(-1.01051, -0.93338, -1.08304)
)
metmat <- cov(metmvn)
metmeans <- c(metmvn[1, 1], metmvn[1, 2], metmvn[1, 3])
PSA_meta_survival <- mvrnorm(mcruns, metmeans, metmat)

#Growth model used to generate the PSA data ("Nanaa" or "MacInnes")
growth_model <- "Nanaa"

#Draw mammography sensitivity coefficients (must match get_parameters() and
#get_psa_parameters() for the growth model used)
dens_logOR <- if (growth_model == "Nanaa") c(0.127, -0.017, -1.364, -2.167) else
  c(0.737, 0.502, -0.799, -1.639)
PSA_sens <- data.frame(
  rnorm(mcruns, -4.8344, 0.202),
  rnorm(mcruns, 0.4948, 0.0284),
  sapply(1:4, function(g) rnorm(mcruns, dens_logOR[g], c(0.289, 0.135, 0.148, 0.226)[g]))
)

#Draw supplemental Screening CDRs
PSA_MRI_cdr <- rbeta(mcruns, 99.495, 19799.5) #CDR for MRI in Mammo negative women (incremental)
PSA_US_cdr <- rbeta(mcruns, 35.89, 11927) #CDR for US in Mammo negative women (incremental)

#Draw tumour growth rate parameters
#Must match the growth model used to generate the PSA data (set above)
if (growth_model == "MacInnes") {
  PSA_log_norm_mean <- rnorm(mcruns, 0.26, 0.053) #MacInnes et al. 2020
  PSA_log_norm_sd <- pmax(rnorm(mcruns, 0.32, 0.05), 0.01)
} else {
  PSA_log_norm_mean <- rnorm(mcruns, -0.04, 0.052) #Nanaa et al. 2025
  PSA_log_norm_sd <- pmax(rnorm(mcruns, 0.48, 0.07), 0.01)
}

#Draw costs
PSA_cost_strat <- (rlnorm(mcruns, 1.8826894, 0.1015175))
PSA_costvar <- rnorm(mcruns, 0, 0.05102041)
PSA_costscreen <- rnorm(mcruns, 0, 0.05102041)
PSA_cost_follow_up <- rnorm(mcruns, 0, 0.05102041)
PSA_cost_biop <- rnorm(mcruns, 0, 0.05102041)
PSA_cost_US <- rnorm(mcruns, 0, 0.05102041)
PSA_cost_MRI <- rnorm(mcruns, 0, 0.05102041)

#Generate utility draws
utilmat <- data.frame(
  c(1 - 0.82, 1 - 0.81, 1 - 0.83),
  c(1 - 0.75, 1 - 0.73, 1 - 0.77)
)
lnutilmat <- log(utilmat)
covutil <- cov(lnutilmat)
utilmeans <- c(log(1 - 0.82), log(1 - 0.75))
PSA_util <- 1 - exp(mvrnorm(mcruns, utilmeans, covutil))

#Bind monte carlo draws
PSA_all_p <- cbind(
  PSA_gamma_survival,
  PSA_meta_survival,
  PSA_sens,
  PSA_MRI_cdr,
  PSA_US_cdr,
  PSA_log_norm_mean,
  PSA_log_norm_sd,
  PSA_cost_strat,
  PSA_costvar,
  PSA_util,
  PSA_costscreen,
  PSA_cost_follow_up,
  PSA_cost_biop,
  PSA_cost_US,
  PSA_cost_MRI
)
PSA_all_p <- as.data.frame(PSA_all_p)
colnames(PSA_all_p) <- c(
  "PSA_gamma_survival_1",
  "PSA_gamma_survival_2",
  "PSA_gamma_survival_3",
  "PSA_meta_survival_54",
  "PSA_meta_survival_74",
  "PSA_meta_survival_99",
  "PSA_sens_b0",
  "PSA_sens_bsize",
  "PSA_dens_logOR_1",
  "PSA_dens_logOR_2",
  "PSA_dens_logOR_3",
  "PSA_dens_logOR_4",
  "PSA_MRI_cdr",
  "PSA_US_cdr",
  "PSA_log_norm_mean",
  "PSA_log_norm_sd",
  "PSA_cost_strat",
  "PSA_costvar",
  "PSA_util_1to3",
  "PSA_util_4",
  "PSA_costscreen",
  "PSA_cost_follow_up",
  "PSA_cost_biop",
  "PSA_cost_US",
  "PSA_cost_MRI"
)

alt_names <- c(
  "No Screening",
  "Risk-1",
  "Risk-2",
  "3 yearly",
  "2 yearly",
  "Risk-3"
)

#Create data.frames to estimate costs and QALYs
output_costscan <- data.frame(matrix(nrow = mcruns, ncol = length(alternative)))
output_costsnocan <- data.frame(matrix(
  nrow = mcruns,
  ncol = length(alternative)
))
output_costs <- data.frame(matrix(nrow = mcruns, ncol = length(alternative)))
output_qalyscan <- data.frame(matrix(nrow = mcruns, ncol = length(alternative)))
output_qalysnocan <- data.frame(matrix(
  nrow = mcruns,
  ncol = length(alternative)
))
output_qalys <- data.frame(matrix(nrow = mcruns, ncol = length(alternative)))
colnames(output_costscan) <- alt_names
colnames(output_costsnocan) <- alt_names
colnames(output_costs) <- alt_names
colnames(output_qalyscan) <- alt_names
colnames(output_qalysnocan) <- alt_names
colnames(output_qalys) <- alt_names

#For each alternative
for (i in 1:length(alternative)) {
  PSA_all_p$alternative <- as.factor(alt_names[i])

  #Predict QALYs and Costs using GAMS
  output_qalyscan[, i] <- predict.bam(modQ, PSA_all_p)
  output_costscan[, i] <- predict.bam(modC, PSA_all_p)
  output_qalysnocan[, i] <- predict.bam(modQ2, PSA_all_p)
  output_costsnocan[, i] <- predict.bam(modC2, PSA_all_p)
  output_qalys[, i] <- (output_qalyscan[, i] * 0.0934) +
    (output_qalysnocan[, i] * (1 - 0.0934))
  output_costs[, i] <- (output_costscan[, i] * 0.0934) +
    (output_costsnocan[, i] * (1 - 0.0934))
}

#Rename columns
alt_names <- c(
  "No Screening",
  "Risk-1",
  "Risk-2",
  "Three yearly",
  "Two yearly",
  "Risk-3"
)
colnames(output_costs) <- alt_names
colnames(output_qalys) <- alt_names

#Produce PSA object
psa_obj <- make_psa_obj(
  cost = output_costs,
  effectiveness = output_qalys,
  parameters = PSA_all_p[1:25],
  strategies = alt_names,
  currency = "£"
)

#Create CEAC object
ceac_obj <- ceac(wtp, psa_obj)

#Plot CEAC
plot(
  ceac_obj,
  frontier = "FALSE",
  points = "FALSE",
  xlab = "Willingness to Pay (Thousand £ / QALY"
)

#Creat PSA results table and save
psa_sum <- summary(psa_obj, calc_sds = TRUE)
psa_sum
write.csv(psa_sum, file = here::here("output", "psa results summary.csv"))

#Plot cost-effectiveness plane for PSA results
icers <- calculate_icers(
  cost = psa_sum$meanCost,
  effect = psa_sum$meanEffect,
  strategies = psa_sum$Strategy
)
plot(icers, labels = "all")


#One-way sensitivity analysis for costs, effectiveness and net benefit
sacost <- owsa(psa_obj, outcome = "cost")
owsa_tornado(sacost, return = "plot")

saeff <- owsa(psa_obj, outcome = "eff")
owsa_tornado(saeff, return = "plot")

sanb <- owsa(psa_obj, outcome = "nmb", wtp = 20000)
owsa_tornado(sanb, return = "plot")

#Plot chart of impact of changes in parameters on optimal strategy
owsa_opt_strat(
  sanb,
  params = c(unique(a$parameter))[c(1, 2, 3, 5, 6, 7, 8, 9, 10)],
  return = "plot",
  col = "full"
)


#Plot ceac omitting infeasible strategies (too high number of scans)
feasqalys <- output_qalys[-c(3, 5)]
feascosts <- output_costs[-c(3, 5)]
feas_psa_obj <- make_psa_obj(
  cost = feascosts,
  effectiveness = feasqalys,
  parameters = PSA_all_p,
  strategies = alt_names[-c(3, 5)],
  currency = "£"
)
feas_ceac_obj <- ceac(wtp, feas_psa_obj)
plot(
  feas_ceac_obj,
  frontier = "FALSE",
  points = "FALSE",
  xlab = "Willingness to Pay (Thousand £ / QALY"
)
