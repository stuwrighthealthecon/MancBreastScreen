# =============================================================================
# run_model(): run the MancBreastScreen discrete event simulation
# =============================================================================
#
# How it works
# ------------
# The model code was written to read its parameters from the R workspace.
# run_model() keeps that code unchanged but runs it inside a private "model
# environment" instead of the global workspace:
#   1. controls and params are copied into a new environment;
#   2. the model functions are copied into it, so they look up parameters there;
#   3. the cohort is created, then each strategy is simulated in a fresh copy
#      of the environment (in parallel), so strategies cannot affect each other.
# Nothing is read from or written to the global workspace.

# Model functions that look up parameters in the model environment
.model_functions <- c(
  "get_drug_matrices", "get_screen_times", "assign_risk_groups",
  "vec_fnLookupBase", "vec_stage_by_size", "vec_screening_result",
  "vec_ca_survival_time", "vec_QALY_counter_core", "vec_QALY_counter",
  "vec_QALY_counter_fast", "run_des_vectorised",
  "create_sample", "create_sample_with_misclass", "get_drug_adj_IM",
  "redraw_drug_pars", "get_psa_parameters", "vbd_at_age",
  "rgengamma_custom", "get_gengamma_params", "pgengamma_custom",
  "qgengamma_custom"
)

#' Run the MancBreastScreen model
#'
#' @param controls A list of run settings; see default_controls(). Settings not
#'   supplied take their default values.
#' @param params A list of model parameters; see get_parameters(). If NULL,
#'   get_parameters(growth_model = controls$growth_model) is called after the
#'   random seed is set. If supplied, it must have been created with the same
#'   growth model as controls$growth_model.
#' @return Invisibly, a list with the controls and parameters used, the output
#'   folders and the result files written. Results are saved as .rds files,
#'   one per strategy and risk group.
run_model <- function(controls = default_controls(), params = NULL) {
  controls <- utils::modifyList(default_controls(), controls, keep.null = TRUE)

  # ---- Random seed ----------------------------------------------------------
  # Seeds both base R and dqrng generators. Each strategy is then reseeded
  # with seed + its position in controls$strategies.
  if (!is.null(controls$seed)) {
    set.seed(controls$seed)
    dqrng::dqset.seed(controls$seed)
  }
  if (is.null(params)) params <- get_parameters(growth_model = controls$growth_model)
  if (!identical(params$growth_model, controls$growth_model)) {
    stop("params were created for growth model '", params$growth_model,
         "' but controls$growth_model is '", controls$growth_model,
         "'. Use get_parameters(growth_model = \"", controls$growth_model, "\").")
  }

  env <- .new_model_env(params)

  # ---- Values derived from controls ----------------------------------------
  env$controls          <- controls
  env$MISCLASS          <- controls$MISCLASS
  env$PREVENTATIVE_DRUG <- controls$PREVENTATIVE_DRUG
  env$PSA               <- as.integer(controls$PSA)
  env$intervals         <- as.integer(controls$intervals)
  env$desired_cases     <- controls$desired_cases
  env$mcruns            <- controls$mcruns
  env$inum              <- ceiling(controls$desired_cases / params$expected_prev) #Women to simulate
  env$output_dir        <- controls$output_dir
  env$sample_fname      <- "possample"
  env$screen_strategies <- unlist(controls$strategies)

  paths <- .output_paths(controls)
  dir.create(paths$det, showWarnings = FALSE, recursive = TRUE)
  dir.create(paths$psa, showWarnings = FALSE, recursive = TRUE)
  env$det_output_path <- paths$det
  env$psa_output_path <- paths$psa

  # ---- Create (or load) the cohort of women ---------------------------------
  eval(.cohort_code, envir = env)

  # ---- Simulate each strategy in parallel ----------------------------------
  n_strat <- length(env$screen_strategies)
  n_cores <- if (is.null(controls$n_cores)) {
    max(1L, parallel::detectCores() - 1L)
  } else {
    controls$n_cores
  }
  n_cores <- max(1L, min(n_strat, n_cores)) #Limit cores for memory management

  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores, outfile = "")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
  } else {
    # One core: run in this R session, avoiding the time and memory cost of
    # starting a separate R process (results are identical)
    foreach::registerDoSEQ()
  }

  `%dopar%` <- foreach::`%dopar%`
  foreach::foreach(
    r         = seq_len(n_strat),
    .packages = c("dqrng", "matrixStats", "tidyverse", "here"),
    .export   = c(".simulate_strategy", ".clone_model_env", ".strategy_code"),
    .noexport = setdiff(ls(), c("env"))
  ) %dopar% {
    .simulate_strategy(r, env, .strategy_code)
  }

  files <- list.files(
    if (controls$PSA) paths$psa else paths$det,
    pattern = "\\.rds$", full.names = TRUE
  )
  invisible(list(
    controls        = controls,
    params          = params,
    det_output_path = paths$det,
    psa_output_path = paths$psa,
    files           = files
  ))
}

# Build the model environment from a parameter list
.new_model_env <- function(params) {
  env <- list2env(params, parent = globalenv())
  for (fn in .model_functions) {
    f <- get(fn, mode = "function")
    environment(f) <- env
    assign(fn, f, envir = env)
  }
  env$cmp_create_sample <- compiler::cmpfun(env$create_sample)
  env
}

# Copy a model environment so a strategy run cannot alter the original
.clone_model_env <- function(env) {
  new_env <- list2env(as.list(env, all.names = TRUE), parent = parent.env(env))
  for (nm in ls(new_env, all.names = TRUE)) {
    f <- get(nm, envir = new_env)
    if (is.function(f) && identical(environment(f), env)) {
      environment(f) <- new_env
      assign(nm, f, envir = new_env)
    }
  }
  new_env
}

# Simulate strategy number r (position in controls$strategies)
.simulate_strategy <- function(r, env, code) {
  e <- .clone_model_env(env)
  e$r <- r
  if (!is.null(e$controls$seed)) {
    set.seed(e$controls$seed + r)
    dqrng::dqset.seed(e$controls$seed + r)
  }
  eval(code, envir = e)
  invisible(NULL)
}

# Output folders for a given set of controls
.output_paths <- function(controls) {
  sub <- if (controls$MISCLASS && controls$PREVENTATIVE_DRUG) {
    c("misclassification_and_preventative_drug", "misclassification_and_preventative_drug")
  } else if (controls$MISCLASS) {
    c("misclassification", "misclassification")
  } else if (controls$PREVENTATIVE_DRUG) {
    c("preventative_drug", "misclassification_and_preventative_drug")
  } else {
    c("nomisclassification", "nomisclassification")
  }
  list(
    det = paste0(file.path(controls$output_dir, "deterministic_results", sub[1]), "/"),
    psa = paste0(file.path(controls$output_dir, "psa_results", sub[2]), "/")
  )
}

# -----------------------------------------------------------------------------
# Cohort creation (evaluated inside the model environment)
# -----------------------------------------------------------------------------
.cohort_code <- quote({
  screen_strategies <- unlist(controls$strategies)

  if (controls$gensample) {
    # Use first strategy for sample generation — sample is shared across strategies
    screen_strategy <- screen_strategies[1]

    if (MISCLASS) {
      dir.create(file.path(output_dir, "risksample_misclass"), showWarnings = FALSE, recursive = TRUE)
      create_sample_with_misclass(PSA, intervals, controls$seed, screen_strategy)
    } else {
      dir.create(file.path(output_dir, "risksample"), showWarnings = FALSE, recursive = TRUE)
      cmp_create_sample(PSA, intervals, controls$seed, screen_strategy)
    }
  }

  # Load sample
  if (MISCLASS) {
    load(file.path(output_dir, "risksample_misclass", paste0(sample_fname, ".Rdata")))
  } else {
    load(file.path(output_dir, "risksample", paste0(sample_fname, ".Rdata")))
  }

  #Fix column names for risksample
  prefix <- paste0("^", "X", 1, ".")
  names(risksample) <- sub(prefix, "", names(risksample))

  # Risk group assignment
  if (MISCLASS) {
    risk_col <- "tenyrrisk_est"
  } else {
    risk_col <- "tenyrrisk"
  }

  # Initialise risk_group — will be overwritten per-strategy inside loop
  # for strategies that use risk stratification; set to 0 (average) as default
  risksample$risk_group <- 0L

  #Assign supplemental screening if used
  if (controls$supplemental_screening) {
    supp_col        <- ifelse(MISCLASS, "tenyrrisk_true", "tenyrrisk_est")
    high_risk_dense <- risksample$VDG >= density_cutoff &
                       risksample[[supp_col]] >= 8
    low_risk_dense  <- risksample$VDG >= density_cutoff &
                       risksample[[supp_col]] < 8
    risksample$MRI_screen[high_risk_dense] <- 1L
    risksample$US_screen[low_risk_dense]   <- 1L
  }

  # Drug setup for PREVENTATIVE_DRUG runs
  if (PREVENTATIVE_DRUG) {
    nsample <- nrow(risksample)
    risksample$starting_menses_status <- ifelse(
      dqrunif(nsample, 0, 1) < prob_premen, 1L, 2L
    )
    risksample$takes_drug        <- logical(nsample)
    risksample$time_taking_drug  <- numeric(nsample)
  }

  # Keep a clean master copy for reuse across strategies
  risksample_master <- risksample

  #Drop un-needed columns for memory management, some now some later
  if (MISCLASS) {
    drop_cols_now  <- c("cancer", "feedback", "liferisk_est", "liferisk_true")
    drop_cols_late <- c("tenyrrisk_est", "tenyrrisk_true")
  } else {
    drop_cols_now  <- c("cancer", "feedback", "liferisk")
    drop_cols_late <- c("tenyrrisk")
  }
  splitmaster_base <- risksample_master[, !names(risksample_master) %in% drop_cols_now]

  keep <- c("risk_group", "VDG", "VBD", "vbd_age", "MRI_screen", "US_screen", "risk_predicted", "interval_change",
            "life_expectancy", "ca_incidence", "clin_detect_size_g", "genage", "growth_rate",
            "invite_jitter_years", "tenyrrisk", "tenyrrisk_est", "tenyrrisk_true",
            "starting_menses_status", "time_taking_drug", grep("^PSA_|^mcid$", names(splitmaster_base), value = TRUE))
  splitmaster_base <- splitmaster_base[, intersect(keep, names(splitmaster_base))]

  for (nm in c("risk_group", "VDG", "MRI_screen", "US_screen", "risk_predicted", "interval_change"))
    if (nm %in% names(splitmaster_base)) splitmaster_base[[nm]] <- as.integer(splitmaster_base[[nm]])



  if (PSA == 0L) {
    psa_cols     <- grep("^PSA_", names(splitmaster_base), value = TRUE)
    splitmaster_base <- splitmaster_base[, !names(splitmaster_base) %in% psa_cols]
  }
  rm(risksample_master, risksample); gc()
})

# -----------------------------------------------------------------------------
# Simulation of one strategy (evaluated inside a copy of the model environment)
# -----------------------------------------------------------------------------
.strategy_code <- quote({
  screen_strategy        <- screen_strategies[r]
  supplemental_screening <- as.integer(controls$supplemental_screening)

  # ---- Recompute strategy-dependent drug matrices --------------------------
  dm          <- get_drug_matrices(screen_strategy)
  risk_red    <- dm$risk_red
  uptake      <- dm$uptake
  persistence <- dm$persistence
  course_length <- dm$course_length

  # ---- Assign risk groups for this strategy --------------------------------
  splitmaster <- splitmaster_base
  splitmaster$risk_group <- assign_risk_groups(splitmaster, screen_strategy,
                                               MISCLASS)
  risk_groups <- sort(unique(splitmaster$risk_group))

  # Drop the risk score columns now that risk groups are assigned
  splitmaster <- splitmaster[, !names(splitmaster) %in% drop_cols_late]

  # Pre-split data.frame before ii loop (for memory management)
  risk_group_list <- split(splitmaster, splitmaster$risk_group)
  rm(splitmaster); gc()

  # ---- ii loop: one iteration per risk group --------------------------------
  for (ii in seq_along(risk_groups)) {
    tryCatch({
      risksample <- risk_group_list[[as.character(risk_groups[ii])]]

    # Catch NA risk_group before it causes errors
    if (is.na(risksample$risk_group[1]))
      stop(sprintf("strategy %d ii %d: risk_group[1] is NA",
                   screen_strategy, ii))

    # ---- Drug incidence adjustment (PREVENTATIVE_DRUG only) ----------------
    if (PREVENTATIVE_DRUG && risksample$risk_group[1] != 0) {

      if (PSA == 1L) {
        drug_matrix_list <- redraw_drug_pars(risksample[1, ])
        risk_red    <- drug_matrix_list[[1]]
        uptake      <- drug_matrix_list[[2]]
        persistence <- drug_matrix_list[[3]]
      }

      idx         <- cbind(risksample$risk_group, risksample$starting_menses_status)
      uptake_probs <- uptake[idx]
      risksample$uptake <- dqrunif(nrow(risksample), 0, 1) < uptake_probs

      #Create column of time taking drugs for those who take them
      risksample$time_taking_drug <- risksample$uptake * pmin(
        rexp(nrow(risksample), rate = persistence[idx]),
        course_length
      )

      risksample$weibullrisk <- (1 - risk_red[idx]) *
        risksample$time_taking_drug *
        log(1 / completion_prob[risksample$starting_menses_status])
      risksample$weibullrisk <- inc_scale / risksample$weibullrisk

      ages <- Incidence_Mortality$age[start_age:101]
      n    <- nrow(risksample)

      prob_matrix <- outer(
        risksample$weibullrisk, ages,
        FUN = function(scale, age) dweibull(age, shape = inc_shape, scale = scale)
      )

      age_mask    <- outer(risksample$life_expectancy, ages, FUN = ">")
      prob_matrix <- prob_matrix * age_mask

      row_sums    <- rowSums(prob_matrix)
      prob_matrix <- prob_matrix / row_sums

      cum_probs   <- matrixStats::rowCumsums(prob_matrix)
      u           <- dqrunif(n, 0, 1)
      col_indices <- rowSums(cum_probs < u) + 1L
      col_indices <- pmin(col_indices, length(ages))

      risksample$ca_incidence <- ifelse(
        risksample$time_taking_drug > 0,
        ages[col_indices] + dqrunif(n, 0, 1),
        risksample$ca_incidence
      )
    }

    # ---- Cap incidence and life expectancy ----------------------------------
    risksample$life_expectancy <- ifelse(
      risksample$life_expectancy >= time_horizon,
      99.99,
      risksample$life_expectancy
    )

    # ---- Screen times for this strategy + risk group -----------------------
    screen_times <- get_screen_times(screen_strategy, risksample$risk_group[1])

    # ---- Screen attendance columns -----------------------------------------
    if (length(screen_times) > 1) {
      n_screens <- length(screen_times)

      # Initialise columns one at a time
      for (k in seq_len(n_screens)) {
        col_k <- paste0("screen_", k)
        risksample[[col_k]] <- rep(0L, nrow(risksample))
      }

      # First screen attendance
      risksample[["screen_1"]] <- rbinom(nrow(risksample), 1L, uptakefirstscreen)

      ever_attended <- risksample[["screen_1"]]

      for (i in seq_len(n_screens - 1L)) {
        next_col <- paste0("screen_", i + 1L)
        new_vals <- ifelse(
          ever_attended >= 1L,
          rbinom(nrow(risksample), 1L, uptakeotherscreen),
          rbinom(nrow(risksample), 1L, uptakenoscreen)
        )
        risksample[[next_col]] <- new_vals
        ever_attended <- pmax(ever_attended, risksample[[next_col]])
      }

    } else {
      risksample[["screen_1"]] <- 999L
    }

    # ---- PSA parameter overrides -------------------------------------------
    if (PSA == 1L) {
      risk_data <- risksample[1, ]

      sens_b0    <- risk_data$PSA_sens_b0
      sens_bsize <- risk_data$PSA_sens_bsize
      dens_logOR <- c(risk_data$PSA_dens_logOR_1, risk_data$PSA_dens_logOR_2,
                      risk_data$PSA_dens_logOR_3, risk_data$PSA_dens_logOR_4)

      # Growth-rate PSA draws are applied to each woman when the cohort is
      # created (create_sample), so they are not overridden here

      gamma_stage <- c(
        exp(risk_data$PSA_gamma_survival_1),
        exp(risk_data$PSA_gamma_survival_2),
        exp(risk_data$PSA_gamma_survival_3)
      )
      metastatic_survival <- c(
        exp(risk_data$PSA_meta_survival_54),
        exp(risk_data$PSA_meta_survival_74),
        exp(risk_data$PSA_meta_survival_99)
      )

      MRI_cdr <- risk_data$PSA_MRI_cdr
      US_cdr  <- risk_data$PSA_US_cdr

      utility_stage_cat_y1 <- c(
        "stage1"     = risk_data$PSA_util_1to3 / 0.822,
        "stage2"     = risk_data$PSA_util_1to3 / 0.822,
        "stage3"     = risk_data$PSA_util_1to3 / 0.822,
        "Metastatic" = risk_data$PSA_util_4    / 0.822,
        "DCIS"       = utility_DCIS
      )
      utility_stage_cat_follow <- utility_stage_cat_y1

      cost_strat    <- risk_data$PSA_cost_strat
      cost_DCIS     <- cost_DCIS_base     * (1 + risk_data$PSA_costvar)
      cost_screen   <- cost_screen_base   * (1 + risk_data$PSA_costscreen)
      cost_follow_up <- cost_follow_up_base * (1 + risk_data$PSA_cost_follow_up)
      cost_biop     <- cost_biop_base     * (1 + risk_data$PSA_cost_biop)
      cost_US       <- cost_US_base       * (1 + risk_data$PSA_cost_US)
      cost_MRI      <- cost_MRI_base      * (1 + risk_data$PSA_cost_MRI)
      cost_drug     <- cost_drug_base     * (1 + risk_data$PSA_cost_drug)
    }


    # ---- Run DES ------------------------------------------------------------
    run_des_vectorised(
      risk_df            = risksample,
      screen_strategy    = screen_strategy,
      start_age          = start_age,
      screen_startage    = screen_startage,
      discount_cost      = discount_cost,
      fp_rate        = fp_rate,
      Vm                 = Vm,
      Vc                 = Vc,
      cost_screen        = cost_screen,
      cost_strat         = cost_strat,
      cost_US            = cost_US,
      cost_MRI           = cost_MRI,
      cost_follow_up     = cost_follow_up,
      cost_biop          = cost_biop,
      cost_DCIS          = cost_DCIS,
      cost_drug          = cost_drug,
      biopsy_rate        = biopsy_rate,
      course_length      = course_length,
      PREVENTATIVE_DRUG  = PREVENTATIVE_DRUG,
      PSA                = PSA,
      cost_in_full_courses = cost_in_full_courses,
      uptake             = uptake,
      persistence        = persistence
    )

    }, error = function(e) {
      message(sprintf(
        "FAILED: strategy=%s ii=%d risk_group=%s nrow=%d PSA=%d mcruns=%d — %s",
        screen_strategy, ii, risk_groups[ii],
        nrow(risk_group_list[[as.character(risk_groups[ii])]]),
        PSA, mcruns, conditionMessage(e)
      ))
      stop(e)
    })
  } # end ii loop

  message(paste("Strategy", r, "complete"))
})
