# =============================================================================
# MancBreastScreen: main script
# =============================================================================
# Runs the discrete event simulation for the chosen screening strategies.
# Strategies are simulated in parallel; set n_cores = 1 to run them one at a
# time. See README.md for the strategy codes.

# ----------------------------------------------------------------------------
# Load required packages
# ----------------------------------------------------------------------------
library(doParallel)
library(MASS)
library(dqrng)
library(compiler)
library(tidyverse)
library(iterators)
library(here)
library(purrr)
library(tictoc)
library(matrixStats)

# ----------------------------------------------------------------------------
# Load model functions
# ----------------------------------------------------------------------------
here::i_am("scripts/MancBreastScreen_main_script.R")
for (f in list.files(here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

# ----------------------------------------------------------------------------
# Set model controls
# ----------------------------------------------------------------------------
controls <- default_controls(
  strategies             = c(0, 1, 2, 3, 4, 9), #Vector of strategies to simulate
  gensample              = TRUE,   #Create new sample to simulate?
  MISCLASS               = FALSE,  #Should error in risk prediction be included?
  PREVENTATIVE_DRUG      = FALSE,  #Should risk reducing medicines be used for high risk?
  supplemental_screening = FALSE,  #Should ultrasound and MRI be used as supplemental screening?
  PSA                    = FALSE,  #Run PSA?
  intervals              = FALSE,  #Run PSA with wide distributions for GAM estimation?
  desired_cases          = 300000, #Number of cancer cases required
  mcruns                 = 1,      #Number of Monte Carlo runs
  seed                   = 42,     #Random seed (NULL for an unseeded run)
  n_cores                = NULL,   #Cores to use (NULL = one per strategy, up to cores - 1)
  growth_model           = "Nanaa" #Tumour growth model: "Nanaa" or "MacInnes"
)

# ----------------------------------------------------------------------------
# Set model parameters
# ----------------------------------------------------------------------------
# Leave as NULL to use the base-case parameters (drawn after the seed is set).
# To change inputs, e.g.:
#   set.seed(controls$seed)
#   params <- get_parameters(growth_model = controls$growth_model)
#   params$cost_screen <- 45
params <- NULL

# ----------------------------------------------------------------------------
# Run the model
# ----------------------------------------------------------------------------
tic()
results <- run_model(controls, params)
toc()

# Output folders, used by the scripts in analysis/
det_output_path <- results$det_output_path
psa_output_path <- results$psa_output_path
