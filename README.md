# MancBreastScreen

*(Previously called MANC-RISK-SCREEN.)*

This model was developed by Ewan Gray, Anna Donten, and Katherine Payne and validated by Stuart Wright, Gabriel Rogers, Katherine Payne, and Rob Hainsworth at the Manchester Centre for Health Economics. The model is currently maintained by Stuart Wright (stuart.j.wright@manchester.ac.uk)

This discrete event simulation model estimates the QALYs, costs, and clinical outcomes of women attending breast cancer screening under a range of screening strategies, or no screening. Strategies are selected by their code in `controls$strategies`:

| Code | Strategy | Screening schedule (ages 50–70 unless stated) |
|---|---|---|
| 0 | No screening | None |
| 1 | Risk-stratified (PROCAS risk cut-offs), using the Tyrer-Cuzick V8 questionnaire and Volpara breast density | 10-year risk < 5%: 3-yearly; 5–8%: 2-yearly; ≥ 8%: annual |
| 2 | Risk-stratified (risk tertiles) | Lowest third: 3-yearly; middle third: 2-yearly; highest third: annual |
| 3 | 3-yearly screening | Every 3 years |
| 4 | 2-yearly screening | Every 2 years |
| 5 | 5-yearly screening | Every 5 years |
| 6 | 10-yearly screening | At ages 50 and 60 |
| 7 | Low risk 5-yearly | 10-year risk < 1.5%: 5-yearly; others: 3-yearly |
| 8 | Low risk 6-yearly | 10-year risk < 1.5%: 6-yearly; others: 3-yearly |
| 9 | Fully stratified (PROCAS risk cut-offs) | 10-year risk < 1.5%: 5-yearly; 1.5–5%: 3-yearly; 5–8%: 2-yearly; ≥ 8%: annual |

Any other code is treated as no screening. The analysis scripts in `analysis/` are set up for strategies 0, 1, 2, 3, 4 and 9, as used in the most recent analysis.

The model was originally published as an early economic evaluation in the journal Value in Health:
Gray E, Donten A, Karssemeijer N, et al. Evaluation of a Stratified National Breast Screening Program in the United Kingdom: An Early Model-Based Cost-Effectiveness Analysis. Value Heal 2017;20:1100–9.

The most recent reference for this model is:
Wright SJ, Hilton J, Rogers G, Gray E, Donten A, McWilliams L, et al. Implementing a Risk-Stratified National Breast Screening Programme in the United Kingdom: A Model-Based Cost-effectiveness Analysis. Cancer Epidemiol Biomarkers Prev. 2026 Sep 2;OF1–11. doi:10.1158/1055-9965.EPI-26-0103


The `docs/` folder contains information about the model. An overall summary is provided in the Model Summary word document. The research team recommend reading this file and the published economic evaluations to gain a full understanding of the model. The documentation folder also contains diagrams outlining the structure of the model, tables of model parameters and their sources, and details about the process of updating the model parameters between version 0 and version 1.

Below are some key notes about the model for new users:
1) The model simulates each strategy listed in `controls$strategies` separately, in parallel. The results for each strategy are saved separately and are not incremental; strategies are compared using the scripts in `analysis/`.
2) Like most discrete event simulation models this model is computationally expensive. Following vectorisation of the code, the model now runs very quickly.
3) DES models create natural variation in the outcomes and so require a large population to be sampled in the model to achieve stable results. Currently we recommend that users simulate results for at least 3 million women, particularly for the more variable risk-stratified strategies. We are investigating variance reduction techniques to reduce this number.
4) Results are saved separately for each strategy and risk group (`Determ_<strategy>_<risk group>.rds`, or `PSA_...` for a PSA). If the model is interrupted mid-run, the results already saved can be retrieved to save running the whole sample again.
5) Key data in this model include a sample of ~15,000 women's Volpara breast density estimates, estimated 10 year risk of breast cancer, and estimated lifetime risk of breast cancer. This data was based on real data from the PROCAS 2 study. In order to create a model version which could widely be shared while being sensitive to data protection concerns, the researchers created a synthetic dataset of these variables basedon the original data. This used the synthpop package in R which maintains the original structure of the data. 

## Folder structure

```
MancBreastScreen/
├── MancBreastScreen.Rproj   Open this first: all paths are found from this folder with here()
├── DESCRIPTION              Package metadata (authors, licence, dependencies)
├── R/                       Model functions
│   ├── run_model.R                    run_model(): runs the whole model
│   ├── parameters.R                   default_controls(), get_parameters(), get_psa_parameters()
│   ├── MancBreastScreen_functions.R   Core model functions (screening, stage, survival, QALYs)
│   ├── des_vectorised_df.R            Vectorised discrete event simulation (run_des_vectorised)
│   └── risksample_functions.R         Creation of the simulated cohort of women
├── scripts/
│   └── MancBreastScreen_main_script.R Sets the controls and runs the model
├── analysis/                Post-processing of model output
│   ├── base_case_analysis.R, psa_analysis.R, validation.R
│   ├── basic_calibration.R
│   └── gam_estimation.R               Fits GAM metamodels to PSA output
├── data-raw/                Scripts used to create the input datasets
├── inst/extdata/            Input data (ONS incidence/mortality, utilities,
│                            synthetic risk data, prevention and posterior outputs)
├── docs/                    Model documentation (summary, diagrams, parameters, user guide)
├── output/                  Model results (not tracked by git)
│   ├── deterministic_results/, psa_results/
│   ├── risksample/, risksample_misclass/
│   └── gam_models/
└── archive/                 Superseded scripts and duplicate data kept for reference
```

## Running the model

1. Open `MancBreastScreen.Rproj` in RStudio.
2. Set the options in the `controls` list in `scripts/MancBreastScreen_main_script.R`.
3. Run the script. Results are written to `output/`.
4. Summarise the results with the scripts in `analysis/`.

The model can also be run directly as a function:

```r
for (f in list.files(here::here("R"), pattern = "\\.R$", full.names = TRUE)) source(f)

controls <- default_controls(strategies = c(0, 3), desired_cases = 30000)
results  <- run_model(controls)          # base-case parameters

# Run with modified parameters
set.seed(controls$seed)
params <- get_parameters()
params$cost_screen <- 45
results <- run_model(controls, params)
```

`run_model()` returns the controls and parameters used, the output folders and the list of result files.

**Reproducibility.** With `seed` set (default 42), runs are reproducible: the seed is set before the parameters and the cohort are drawn, and each strategy is reseeded with `seed` + its position in `strategies`. Results do not depend on the number of cores used. Set `seed = NULL` for an unseeded run.

## AI usage statement

Generative AI tools were used in the development of this repository:

- **Code vectorisation:** Claude Sonnet 4.8 (Anthropic) was used to help vectorise the
  discrete event simulation code (e.g. `R/des_vectorised_df.R` and the vectorised
  functions in `R/MancBreastScreen_functions.R`) to reduce model run-time.
- **Repository reorganisation:** In October 2026, Claude (Anthropic; model
  `claude-opus-5-5`) was used with the HEPackageR skill
  (<https://github.com/HealthyUncertainty/HEPackageR>) to reorganise the repository
  into an R package-style structure, rename files from the previous model name
  (MANC-RISK-SCREEN) to MancBreastScreen, update file paths, and remove references
  to the retired `negsample` function. No model logic or parameter values were
  changed in this step.
- **Refactoring into functions:** Claude was also used to restructure the model so
  it runs through `run_model(controls, params)`, with parameters returned by
  `get_parameters()`; to use `here()` for file paths; to merge the sequential and
  parallel run scripts; and to add random seeding so runs are reproducible. Model
  logic was not changed: with seeding added to the original code, outputs from
  the original and refactored versions were identical for deterministic,
  misclassification, PSA and preventative drug runs on a small test cohort.

All AI-assisted changes were reviewed by the research team, who take full
responsibility for the code and its outputs.
