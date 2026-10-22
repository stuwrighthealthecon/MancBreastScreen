# =============================================================================
# MancBreastScreen model parameters
# =============================================================================
# The model is run with:
#   controls <- default_controls()          # how to run the model
#   params   <- get_parameters()            # model inputs
#   run_model(controls, params)
#
# To change an input, edit the returned list before calling run_model(), e.g.
#   params <- get_parameters(); params$cost_screen <- 45
# Note that derived values (e.g. cost look-up tables) are not recalculated when
# a single element is edited; edit the code below to change those inputs.

#' Default model run settings
#'
#' @param ... Named settings that override the defaults.
#' @return A named list of run settings.
default_controls <- function(...) {
  controls <- list(
    strategies             = c(0, 1, 2, 3, 4, 9), #Vector of strategies to simulate (see README)
    gensample              = TRUE,   #Create a new sample of women to simulate?
    MISCLASS               = FALSE,  #Include error in risk prediction?
    PREVENTATIVE_DRUG      = FALSE,  #Offer risk-reducing medicines to high-risk women?
    supplemental_screening = FALSE,  #Use ultrasound and MRI as supplemental screening?
    PSA                    = FALSE,  #Run a probabilistic sensitivity analysis?
    intervals              = FALSE,  #Use wide PSA distributions (for GAM estimation)?
    desired_cases          = 300000, #Approximate number of cancer cases to simulate
    mcruns                 = 1,      #Number of Monte Carlo runs (PSA only)
    seed                   = 42,     #Random seed; NULL for an unseeded run
    n_cores                = NULL,   #Cores to use; NULL = one per strategy, up to cores - 1
    output_dir             = here::here("output"), #Folder for samples and results
    growth_model           = "Nanaa" #Tumour growth model: "Nanaa" or "MacInnes" (see get_parameters())
  )
  overrides <- list(...)
  unknown <- setdiff(names(overrides), names(controls))
  if (length(unknown) > 0)
    stop("Unknown control(s): ", paste(unknown, collapse = ", "))
  utils::modifyList(controls, overrides, keep.null = TRUE)
}

#' Base-case model parameters
#'
#' Some calibrated parameters (lc_bump, gamma_mean, gamma_sd, gamma_p) are drawn
#' from uniform distributions each time this function is called, so set a seed
#' first (or let run_model() call it after seeding) for reproducible values.
#'
#' @param screen_strategy Strategy used to initialise the drug matrices.
#' @param growth_model Tumour growth model, "Nanaa" (default) or "MacInnes". Sets the
#'   growth-rate parameters and the natural history parameters calibrated with
#'   them (lc_bump and the tumour genesis age distribution).
#' @return A named list of model parameters.
get_parameters <- function(screen_strategy = 0L, growth_model = "Nanaa") {
  #Tumour growth model and the parameters calibrated with it (22/10/2026:
  #cumulative hazard scaling of lifetime risk, Strandberg sensitivity with
  #age-related density decline; calibrated to incidence by age, 3-yearly screening)
  #  MacInnes et al. 2020: incidence RMSE 18.5, 7.3% undiagnosed at death,
  #    CDR 7.0 per 1,000 screens, 49% of cancers aged 50-70 screen-detected
  #  Nanaa et al. 2025: incidence RMSE 22.4, 10.0% undiagnosed at death,
  #    CDR 8.4 per 1,000 screens, 54% of cancers aged 50-70 screen-detected
  growth_model <- match.arg(growth_model, c("MacInnes", "Nanaa"))
  calibrated <- list(
    MacInnes = list(log_norm_mean = 0.26, log_norm_sd = 0.32, lc_bump = 0.52,
                    gamma_mean = 82.25, gamma_sd = 24.38, gamma_p = 5.53),
    Nanaa    = list(log_norm_mean = -0.04, log_norm_sd = 0.48, lc_bump = 0.81,
                    gamma_mean = 69.67, gamma_sd = 20.19, gamma_p = 6.06)
  )[[growth_model]]

  #Latent cancer bump up (cumulative hazard scaling of lifetime risk)
  lc_bump <- calibrated$lc_bump
  # Previous (Cambridge) 0.67
  # Leeds 6.00
  # WF 4.34

  # Global option for whether to use corrected age-at-detection distribution
  CORRECT_BC_AGE <- TRUE

  #Screening uptake
  uptakefirstscreen <- 0.625 #Uptake for the first screen
  uptakeotherscreen <- 0.887 #Uptake if woman has attended >=1 screen
  uptakenoscreen <- 0.208 #Uptake if woman has not previously attended screening

  #Uptake for risk stratification
  risk_uptake <- 1 #Uptake for risk prediction
  risk_feedback <- 1 #Uptake for attending risk feedback consultation
  screen_change <- 1 #Uptake for changing screening intervals based on risk

  #Age of an individual at start of simulation
  start_age <- 18 #Default is 38 to give time for tumours to develop pre-screening

  #Set time horizon
  time_horizon <- 100

  #Set health and cost discount rates
  discount_health <- 0.035
  discount_cost <- 0.035

  ##############Set clinical parameters for screening########################

  #Set parameters of a Weibull survival curve to represent all cause mortality
  acmmortality_wb_a <- 7.937

  acmmortality_wb_b <- 86.788

  #Cancer genesis age
  gamma_mean <- calibrated$gamma_mean #MacInnes 82.25, Nanaa 69.67
  #Cambridge 77.18
  #Leeds 81
  #WF 85.04


  gamma_sd <- calibrated$gamma_sd #MacInnes 24.38, Nanaa 20.19
  #Cambridge 20.73
  #Leeds 21
  #WF 23.07

  gamma_p <- calibrated$gamma_p #MacInnes 5.53, Nanaa 6.06
  #Cambridge 6.03
  #Leeds 6
  #WF 5.45


  gengammaparams<-get_gengamma_params(gamma_mean-18,gamma_sd,gamma_p)

  #Set parameters for all cause mortality following breast cancer
  gamma_survival_1 <- exp(-5.618) #Exponential distribution scale parameter stage 1
  gamma_survival_2 <- exp(-3.808) #Exponential distribution scale parameter stage 2
  gamma_survival_3 <- exp(-2.731) #Exponential distribution scale parameter stage 3
  gamma_stage <- c(gamma_survival_1, gamma_survival_2, gamma_survival_3)

  # Expected population prevalence for doing Bayes correction
  expected_prev <- .12

  #Read in distribution of cancer incidence by age
  Incidence_Mortality <- if (CORRECT_BC_AGE) {
    read.csv(here::here("inst", "extdata", "Incidence_Mortality_ONS2.csv")) %>%
      mutate(
        BC_age = expected_prev *
          Cond.on.getting.BC..prob.of.getting.cancer.at.age.t /
          (surv / surv[1])
      ) %>%
      mutate(BC_age = replace(BC_age, age == 100, 1 - sum(BC_age[1:100])))
  } else {
    read.csv(here::here("inst", "extdata", "Incidence_Mortality_ONS2.csv")) %>%
      mutate(BC_age = Cond.on.getting.BC..prob.of.getting.cancer.at.age.t)
  }

  #Set metastatic cancer probabilities by cancer size
  metastatic_prob <- data.frame(
    c(25, 35, 45, 55, 65, 75, 85),
    c(
      0.046218154,
      0.086659039,
      0.109768116,
      0.127099924,
      0.142505975,
      0.159837783,
      1.73E-01
    )
  )


  #Create matrix of probability of cancer stage by cancer size
  stage_by_size_mat <- data.frame(
    "v1" = c(0.383, 0.567, 0.611, 0.557, 0, 0),
    "v2" = c(0.033, 0.111, 0.180, 0.208, 0.723, 0.563),
    "v3" = c(0.058, 0.057, 0.089, 0.147, 0.206, 0.351),
    "v5" = c(0.526, 0.265, 0.120, 0.088, 0.071, 0.086)
  )

  #Set mean and sd of tumour doublings at clinical detection
  clin_detection_m <- 6.5
  clin_detection_sd <- 0.535

  #Mammography sensitivity by tumour diameter and breast density (Strandberg &
  #Humphreys 2022): sensitivity = logistic(sens_b0 + sens_bsize * diameter (mm)
  #+ sens_bdens * percent density (as a proportion, 0-1))
  sens_b0    <- -4.99 #Intercept
  sens_bsize <-  0.49 #Per mm of tumour diameter
  sens_bdens <- -2.09 #Per unit percent density (proportion)
    #Previously Weedon-Fekjaer size-only curve: beta1 1.47, beta2 6.51

  #Volpara density grade cut-points (density grades are used for false-positive
  #rates and supplemental screening eligibility)
  VDG_interval <- c(4.5, 7.5, 15.5)
  #Mammography sensitivity by volpara density grade (no longer used for
  #mammography sensitivity, which now uses sens_bdens; kept for reference)
  Sen_VDG <- c(0.75, 0.735, 0.598, 0.513)
  Sen_VDG_av <- 0.666

  #Age-related decline in Volpara volumetric percent density. Each woman's
  #density is scaled from the age it was measured by the ratio of mean density
  #at the two ages (BreastScreen Norway, CEBP 2018; cross-sectional,
  #interpolated; held flat outside 49-69)
  density_age     <- c(49,   50,   55,   60,   62.5, 65,   69)
  density_rel     <- c(1.04, 1.00, 0.93, 0.86, 0.82, 0.80, 0.77)
  density_ref_age <- 50 #Measurement age assumed if a woman's age is missing

  #Supplemental screening sensitivity parameters
  Mammo_cdr <- 4.2 #Cancer detection rate per 1000 high dense screens
  MRI_cdr <- 5 #CDR for MRI in Mammogram negative women (incremental)
  US_cdr <- 3 #CDR for US in Mammogram negative women (incremental)

  density_cutoff <- 3

  #Set tumour growth rate parameters
  #Mean of log growth rate: MacInnes et al. 2020 0.26 (median doubling time
  #196 days); Nanaa et al. 2025 -0.04 (264 days); WF 1.07
  log_norm_mean <- calibrated$log_norm_mean
  #Variance of log growth rate: MacInnes 0.32 (approx.); Nanaa 0.48; WF 1.31
  log_norm_sd <- calibrated$log_norm_sd
  max_size <- 128 #Maximum size of tumours, mm diameter
  start_size <- 0.25 #Starting size of tumours, diameter in mm
  Vc = (4 / 3) * pi * (start_size / 2)^3 #Volume at start
  Vm = (4 / 3) * pi * (max_size / 2)^3 #Max volume

  #Metatstatic survival parameters
  meta_survival_54 <- exp(-1.787) #Age <= 54
  meta_survival_74 <- exp(-1.388) #Age 55-74
  meta_survival_99 <- exp(-1.011) #Age 75+
  metastatic_survival <- c(meta_survival_54, meta_survival_74, meta_survival_99)

  #Set screening ages
  screen_startage <- 50 #Starting age
  screen_endage <- 70 #Ending age
  high_risk_screentimes <- seq(screen_startage, screen_endage, 1) #Yearly screening
  med_risk_screentimes <- seq(screen_startage, screen_endage, 2) #Two yearly screening
  low_risk_screentimes <- seq(screen_startage, screen_endage, 3) #Three yearly screening

  #Set maximum screening sensitivity
  sensitivity_max <- 0.95

  #Risk cut-offs for different screening approaches
  risk_cutoffs_procas <- c(1.5, 3.5, 5, 8, 100) #PROCAS and PROCAS Full
  risk_cutoffs_tert <- c(1.946527, 2.942792) #Tertiles of risk
  low_risk_cut <- 1.5 #Cut off in low risk only strategies

  #Cancer size cut-points
  ca_size_cut <- c(0.025, 5, 10, 15, 20, 30, 128) #Size cut-points for deciding stage

  ##########False Positive and Overdiagnosis parameters################
  fp_rate <- c(0.0239,0.0355,0.0459,0.0439) #UK false-positive rate by breast density grade 1-4
  biopsy_rate <- 0.506 #Proportion of referrals without cancer that have biopsy

  #### Drug data ####

  # First bring in log hazard ratios from networked analysis
  loghaz_ests <- readRDS(here::here("inst", "extdata", "PreventionOutputs.RDS"))
  efficacy_ests <- loghaz_ests[1]
  dropout_ests <- loghaz_ests[4]

  # Extract parameters for multivariate normal draws
  efficacy_mu <- efficacy_ests$AnyBC$means %>% as.numeric()
  efficacy_sigma <- efficacy_ests$AnyBC$vcov %>% as.matrix()
  dropout_mu <- dropout_ests$Adherence$means %>% as.numeric()
  dropout_sigma <- dropout_ests$Adherence$vcov %>% as.matrix()

  # Fit a Weibull distribution to data in Incidence_Mortality. Drug acts to change scale parameter.
  weibull_fit <- sample(
    Incidence_Mortality$age,
    size = 100000,
    prob = Incidence_Mortality$BC_age,
    replace = TRUE
  ) %>%
    fitdistr("weibull", lower = c(0, 0))
  inc_scale <- weibull_fit$estimate["scale"]
  inc_shape <- weibull_fit$estimate["shape"]

  ana_eff <- exp(efficacy_ests$AnyBC$means$Anastrozole)
  tam_eff <- exp(efficacy_ests$AnyBC$means$Tamoxifen)

  full_course_len <- 5

  # Assume constant drop out rate with 77% of individuals reaching 5yr mark
  tam_completion_prob <- .77

  # Estimate Anastrozole completion probability from Tamoxifen estimate and hazard ratios:
  ana_completion_prob <- exp(
    dropout_ests$Adherence$means$Anastrozole -
      dropout_ests$Adherence$means$Tamoxifen
  ) *
    tam_completion_prob

  completion_prob <- c(ana_completion_prob, tam_completion_prob)

  # Now estimate per-unit-time dropout rates based on exponential time to dropout:
  tam_dropout_rate <- (1. / full_course_len) * log(1 / (tam_completion_prob))
  ana_dropout_rate <- (1. / full_course_len) * log(1 / (ana_completion_prob))

  # Work out mean time taking each drug
  mean_tam_length <- full_course_len *
    (1 - tam_completion_prob) /
    log(1 / tam_completion_prob)
  mean_ana_length <- full_course_len *
    (1 - ana_completion_prob) /
    log(1 / ana_completion_prob)

  # Quick check: quantities below give efficacy of taking full five-year course
  # assuming linear relationship between time taking and hazard ratio. If both of
  # these are less than one then the linear model is safe to use because no one
  # goes past this point.
  tam_full_course_eff <- 1 -
    (1 - tam_eff) * log(1 / tam_completion_prob) / (1 - tam_completion_prob)
  ana_full_course_eff <- 1 -
    (1 - ana_eff) * log(1 / ana_completion_prob) / (1 - ana_completion_prob)
  if ((tam_full_course_eff > 1) | (ana_full_course_eff > 1)) {
    print(
      "Assumption of linear change to hazard ratio with time taking drug will
          not work for one or both drugs being simulated."
    )
  }

  course_length <- c(5., 5.)

  #Drug matrices for the strategy given in screen_strategy (recomputed for each
  #strategy by get_drug_matrices() during a model run)
  if (screen_strategy == 1 | screen_strategy == 9) {
    risk_red <- matrix(
      c(
        ana_eff,
        tam_eff,
        ana_eff,
        tam_eff,
        ana_eff,
        tam_eff,
        ana_eff,
        tam_eff,
        ana_eff,
        tam_eff
      ),
      nrow = 5,
      ncol = 2
    )

    course_length <- c(5., 5.)

    uptake <- rbind(c(0., 0.), c(0., 0.), c(0., 0.), c(.71, .71), c(.71, .71))

    persistence <- matrix(
      c(
        ana_dropout_rate,
        tam_dropout_rate,
        ana_dropout_rate,
        tam_dropout_rate,
        ana_dropout_rate,
        tam_dropout_rate,
        ana_dropout_rate,
        tam_dropout_rate,
        ana_dropout_rate,
        tam_dropout_rate
      ),
      nrow = 5,
      ncol = 2
    )
  } else if (screen_strategy == 2) {
    risk_red <- matrix(
      c(ana_eff, tam_eff, ana_eff, tam_eff, ana_eff, tam_eff),
      nrow = 3,
      ncol = 2
    )

    uptake <- rbind(c(0., 0.), c(0., 0.), c(.71, .71))

    persistence <- matrix(
      c(
        ana_dropout_rate,
        tam_dropout_rate,
        ana_dropout_rate,
        tam_dropout_rate,
        ana_dropout_rate,
        tam_dropout_rate
      ),
      nrow = 3,
      ncol = 2
    )
  } else if (screen_strategy == 7 | screen_strategy == 8) {
    risk_red <- matrix(c(ana_eff, tam_eff, ana_eff, tam_eff), nrow = 2, ncol = 2)

    uptake <- rbind(c(0., 0.), c(.71, .71))

    persistence <- matrix(
      c(ana_dropout_rate, tam_dropout_rate, ana_dropout_rate, tam_dropout_rate),
      nrow = 2,
      ncol = 2
    )
  } else {
    risk_red <- matrix(c(ana_eff, tam_eff), nrow = 1, ncol = 2)

    uptake <- matrix(c(0., 0.), nrow = 1, ncol = 2)

    persistence <- matrix(
      c(ana_dropout_rate, tam_dropout_rate),
      nrow = 1,
      ncol = 2
    )
  }

  age_prescribed <- 50

  median_age_at_men <- 51
  prob_premen <- exp(-age_prescribed * (log(2) / median_age_at_men)) # Assuming exponential distribution (not correct!)

  cost_in_full_courses <- TRUE # If TRUE then each person to be prescribed drug incurs cost of full course, otherwise cost is proportional to time taking

  #######################Cost Data#########################################

  cost_strat <- 7.40 #Cost of risk prediction
  cost_screen_base <- 39.93 #Cost of Mammography
  cost_follow_up_base <- 130.97 #Cost of follow-up
  cost_biop_base <- 558 #Cost of biopsyy
  cost_DCIS_base <- 12346.83 #Cost of treating DCIS
  cost_US_base <- 102.00 #Cost of ultrasound
  cost_MRI_base <- 162.00 #Cost of MRI
  cost_drug_base <- c(100., 100.) # Cost of full course of drug

  #If deterministic analysis then set costs as base costs
    cost_DCIS <- cost_DCIS_base
    cost_screen <- cost_screen_base
    cost_follow_up <- cost_follow_up_base
    cost_biop <- cost_biop_base
    cost_US <- cost_US_base
    cost_MRI <- cost_MRI_base
    cost_drug <- cost_drug_base


  #Set up look-up table for treatment costs
  tbl <- tribble(
    ~Yr , ~Early_18.64 , ~Late_18.64 , ~Diff1 , ~Early_65plus , ~Late_65plus , ~Diff2 ,
      0 ,          464 ,         607 ,    143 ,          1086 ,         1324 ,    238 ,
      1 ,        10746 ,       13315 ,   2569 ,          7597 ,         8804 ,   1207 ,
      2 ,         3357 ,        5785 ,   2429 ,          2529 ,         3650 ,   1121 ,
      3 ,         1953 ,        3782 ,   1829 ,          2156 ,         3170 ,   1014 ,
      4 ,         1627 ,        2932 ,   1305 ,          2230 ,         2924 ,    693 ,
      5 ,         1617 ,        2841 ,   1225 ,          2077 ,         2957 ,    880 ,
      6 ,         1547 ,        2645 ,   1099 ,          2174 ,         2783 ,    609 ,
      7 ,         1394 ,        2618 ,   1225 ,          2063 ,         2903 ,    840 ,
      8 ,         1376 ,        2559 ,   1183 ,          2134 ,         2454 ,    320 ,
      9 ,         1279 ,        1848 ,    569 ,          2204 ,         2932 ,    728
  ) %>%
    dplyr::select(-Diff1, -Diff2) %>%
    pivot_longer(
      cols = contains("6"),
      names_to = c("Stage", "Age"),
      names_sep = "_",
      values_to = "Cost"
    ) %>%
    group_by(Stage, Age) %>%
    mutate(
      DCost = Cost - first(Cost),
      DCost.i = DCost * 1.402093, # NHSCII inflator for 2010/11-->2023/2024
      disc = 1 / 1.035^(Yr - 0.5),
      DCost.i.d = DCost.i * disc,
      CDCost.i.d = cumsum(DCost.i.d),
      Yr1 = as.factor(Yr == 1),
      Yr2 = as.factor(Yr == 2),
      Yr3 = as.factor(Yr == 3)
    ) %>%
    filter(Yr > 0) %>%
    arrange(Stage, Age, Yr)

  #Create predictive model of cost by age of diagnosis, stage, life expectancy
  mod <- lm(
    data = tbl,
    formula = log(DCost) ~ (Yr1 + Yr2 + Yr3 + Yr) * Stage * Age
  )

  #Prediction matrix
  tblNewDat <- crossing(
    Yr = 1:50,
    Stage = c("Early", "Late"),
    Age = c("18.64", "65plus")
  ) %>%
    mutate(
      Yr1 = as.factor(Yr == 1),
      Yr2 = as.factor(Yr == 2),
      Yr3 = as.factor(Yr == 3)
    )

  #Generate cost predictions
  tblNewDat %>%
    bind_cols(pred = mod %>% predict(newdata = tblNewDat)) %>%
    mutate(DCost.p = exp(pred)) -> tblPred

  #Make lookup table for costs
  tblLookup <- tblPred %>%
    filter(Yr == 1) %>%
    mutate(across(c(Yr, pred, DCost.p), ~0)) %>%
    bind_rows(tblPred) %>%
    group_by(Stage, Age) %>%
    mutate(
      DCost.p.i = DCost.p * 1.402093, # NHSCII inflator for 2010/11-->2023/2024
      disc = 1 / 1.035^(Yr - 0.5),
      DCost.p.i.d = DCost.p.i * disc,
      CDCost.p.i.d = cumsum(DCost.p.i.d),
      StageEarly = Stage == "Early",
      AgeYoung = Age == "18.64"
    ) %>%
    arrange(Stage, Age, Yr) %>%
    ungroup()

  tblLookup$stage_flag <- ifelse(tblLookup$Stage == "Early", 0L, 1L)
  tblLookup$age_flag   <- ifelse(tblLookup$Age   == "18.64", 0L, 1L)
  tblLookup$lookup_key <- paste(tblLookup$stage_flag,
                                    tblLookup$age_flag,
                                    tblLookup$Yr, sep = "|")

  #######################Utility Weights#########################################

  #Set age adjusted utility values
  utility_ages <- read.csv(here::here("inst", "extdata", "age_related_utilities.csv"))

  #Set time independent utility decrements
  #Utility of DCIS
  utility_DCIS <- 1 #assumes no effect

  #Set first year cancer utilities:
  utility_stage_cat_y1 <- c(
    "stage1" = 0.82 / 0.822,
    "stage2" = 0.82 / 0.822,
    "stage3" = 0.75 / 0.822,
    "Metastatic" = 0.75 / 0.822,
    "DCIS" = utility_DCIS
  )

  #Set following year cancer utilities:
  utility_stage_cat_follow <- c(
    "stage1" = 0.82 / 0.822,
    "stage2" = 0.82 / 0.822,
    "stage3" = 0.75 / 0.822,
    "Metastatic" = 0.75 / 0.822,
    "DCIS" = utility_DCIS
  )


  params <- as.list(environment())
  params$screen_strategy <- NULL
  params$calibrated <- NULL
  params
}

#' Parameters for the probabilistic sensitivity analysis
#'
#' @param intervals If TRUE, return the wide distributions used to generate
#'   data for the GAM metamodels.
#' @return A named list of PSA distribution parameters.
get_psa_parameters <- function(intervals = FALSE, growth_model = "Nanaa") {
  growth_model <- match.arg(growth_model, c("MacInnes", "Nanaa"))
  if (!intervals) {
    #Stage I to III survival parameters
    survmvn <- data.frame(
      c(-5.6178, -5.2857, -6.1138),
      c(-3.8083, -3.7400, -3.8704),
      c(-2.7311, -2.6687, -2.7963)
    )
    survcovmat <- cov(survmvn)
    survmeans <- c(survmvn[1, 1], survmvn[1, 2], survmvn[1, 3])

    #Metastatic survival parameters
    metmvn <- data.frame(
      c(-1.78723, -1.67922, -1.89434),
      c(-1.38762, -1.33512, -1.49956),
      c(-1.01051, -0.93338, -1.08304)
    )
    metmat <- cov(metmvn)
    metmeans <- c(metmvn[1, 1], metmvn[1, 2], metmvn[1, 3])

    #Mammography sensitivity coefficients (normal: mean, SE). Means from
    #Strandberg & Humphreys 2022; SEs approximate, from the 95% CIs reported for
    #the KARMA cohort model (Sci Rep 2023). Coefficients are drawn independently
    #(their correlation is not available), which overstates uncertainty
    PSA_sens_b0    <- c(-4.99, 0.20)
    PSA_sens_bsize <- c(0.49, 0.028)
    PSA_sens_bdens <- c(-2.09, 0.46)

    #Draw supplemental Screening CDRs
    PSA_MRI_cdr <- c(99.495, 19799.5)
    PSA_US_cdr <- c(35.89, 11927)

    #Tumour growth rate parameters (normal: mean, SE); log_norm_sd is a variance
    if (growth_model == "MacInnes") {
      #MacInnes et al. 2020: SE from the 95% CI of the mean doubling time;
      #variance approximate, SE assumed
      PSA_log_norm_mean <- c(0.26, 0.053)
      PSA_log_norm_sd <- c(0.32, 0.05)
    } else {
      #Nanaa et al. 2025: SE of the log median approximated from the IQR
      #(n = 281); SE of the variance assumed
      PSA_log_norm_mean <- c(-0.04, 0.052)
      PSA_log_norm_sd <- c(0.48, 0.07)
    }

    # Chemoprevention Drug parameters

    # First bring in log hazard ratios from networked analysis
    loghaz_ests <- readRDS(here::here("inst", "extdata", "PreventionOutputs.RDS"))
    efficacy_ests <- loghaz_ests[1]
    dropout_ests <- loghaz_ests[4]

    # Extract parameters for multivariate normal draws
    efficacy_mu <- efficacy_ests$AnyBC$means %>% as.numeric()
    efficacy_sigma <- efficacy_ests$AnyBC$vcov %>% as.matrix()
    dropout_mu <- dropout_ests$Adherence$means %>% as.numeric()
    dropout_sigma <- dropout_ests$Adherence$vcov %>% as.matrix()

    #Drug uptake parameters
    PSA_uptake_1 <- c(.71, .1)
    PSA_uptake_2 <- c(.71, .1)

    #Draw costs
    PSA_cost_strat <- c(1.8826894, 0.1015175)
    cost_inflator <- 0.1
    cost_dist_sd <- 0.05102041
    PSA_costvar <- c(0, cost_dist_sd)
    PSA_costscreen <- c(0, cost_dist_sd)
    PSA_cost_follow_up <- c(0, cost_dist_sd)
    PSA_cost_biop <- c(0, cost_dist_sd)
    PSA_cost_US <- c(0, cost_dist_sd)
    PSA_cost_MRI <- c(0, cost_dist_sd)
    PSA_cost_drug <- c(0, cost_dist_sd)

    #Generate utility draws
    utilmat <- data.frame(
      c(1 - 0.82, 1 - 0.81, 1 - 0.83),
      c(1 - 0.75, 1 - 0.73, 1 - 0.77)
    )
    lnutilmat <- log(utilmat)
    covutil <- cov(lnutilmat)
    utilmeans <- c(log(1 - 0.82), log(1 - 0.75))

  } else {
    #Stage I to III survival parameters
    PSA_gamma_survival1 <- c(-6, -3)
    PSA_gamma_survival2 <- c(-5, -1)
    PSA_gamma_survival3 <- c(-3, -0.01)

    #Metastatic survival parameters
    PSA_meta_survival1 <- c(-2.5, -0.5)
    PSA_meta_survival2 <- c(-2, -0.3)
    PSA_meta_survival3 <- c(-1.5, -0.1)

    #Mammography sensitivity coefficients: twice the standard PSA SEs
    PSA_sens_b0    <- c(-4.99, 0.40)
    PSA_sens_bsize <- c(0.49, 0.056)
    PSA_sens_bdens <- c(-2.09, 0.92)

    #Draw supplemental Screening CDRs
    PSA_MRI_cdr <- c(99.495, 19799.5)
    PSA_US_cdr <- c(35.89, 11927)

    #Tumour growth rate parameters: wide distributions centred on the chosen model
    if (growth_model == "MacInnes") {
      PSA_log_norm_mean <- c(0.26, 0.15)
      PSA_log_norm_sd <- c(0.32, 0.1)
    } else {
      PSA_log_norm_mean <- c(-0.04, 0.15)
      PSA_log_norm_sd <- c(0.48, 0.1)
    }

    # Chemoprevention Drug parameters

    # First bring in log hazard ratios from networked analysis
    loghaz_ests <- readRDS(here::here("inst", "extdata", "PreventionOutputs.RDS"))
    efficacy_ests <- loghaz_ests[1]
    dropout_ests <- loghaz_ests[4]

    # Extract parameters for multivariate normal draws
    efficacy_mu <- efficacy_ests$AnyBC$means %>% as.numeric()
    efficacy_sigma <- efficacy_ests$AnyBC$vcov %>% as.matrix()
    dropout_mu <- dropout_ests$Adherence$means %>% as.numeric()
    dropout_sigma <- dropout_ests$Adherence$vcov %>% as.matrix()

    #Drug uptake parameters
    PSA_uptake_1 <- c(.71, .1)
    PSA_uptake_2 <- c(.71, .1)

    #Draw costs
    PSA_cost_strat <- c(1.8826894, 0.1015175)
    cost_inflator <- 0.2
    cost_dist_sd <- 0.1020408
    PSA_costvar <- c(0, cost_dist_sd)
    PSA_costscreen <- c(0, cost_dist_sd)
    PSA_cost_follow_up <- c(0, cost_dist_sd)
    PSA_cost_biop <- c(0, cost_dist_sd)
    PSA_cost_US <- c(0, cost_dist_sd)
    PSA_cost_MRI <- c(0, cost_dist_sd)
    PSA_cost_drug <- c(0, cost_dist_sd)

    #Generate utility draws
    PSA_util_1 <- c(0.6, 0.9)
    PSA_util_4 <- c(0.5, 0.8)

  }
  params <- as.list(environment())
  params$intervals <- NULL
  params$growth_model <- NULL
  params
}
