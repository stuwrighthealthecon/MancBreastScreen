# MancBreastScreen

*(Previously called MANC-RISK-SCREEN.)*

This model was developed by Ewan Gray, Anna Donten, and Katherine Payne and validated by Stuart Wright, Gabriel Rogers, Katherine Payne, and Rob Hainsworth at the Manchester Centre for Health Economics. The model is currently maintained by Stuart Wright (stuart.j.wright@manchester.ac.uk)

This discrete event simulation model estimates the QALYs, costs, and clinical outcomes of women attending breast cancer screening delivered using one of 6 strategies or for#no screening. The available screening programmes are 
   1) risk stratified screening using the Tyrer-Cuzick V8 questionnaire and Volpara breast density scan using risk cut-offs defined in the PROCAS project.
   2) risk stratified screening using the Tyrer-Cuzick V8 questionnaire and Volpara breast density scan using risk cut-offs dividing the population into tertiles
   3) 3 yearly breast cancer screening
   4) 2 yearly breast cancer screening
   5) 5 yearly breast cancer screening
   6) 10 yearly breast cancer screening (at age 50 and 60)
   7) No breast cancer screening
   8) 5 yearly breast cancer screening for women at low risk predict using Tyrer-Cuzick V8 questionnaire and Volpara breast density scan using risk cut-offs defined in the PROCAS project
   9) 6 yearly breast cancer screening for women at low risk predict using Tyrer-Cuzick V8 questionnaire and Volpara breast density scan using risk cut-offs defined in the PROCAS project
   10) A fully stratified breast cancer screening programme with women at high or moderate risk receiving annual and bi-annual screening respectively and lower risk women receiving 5 yearly screening.

The model was originally published as an early economic evaluation in the journal Value in Health:
Gray E, Donten A, Karssemeijer N, et al. Evaluation of a Stratified National Breast Screening Program in the United Kingdom: An Early Model-Based Cost-Effectiveness Analysis. Value Heal 2017;20:1100–9.

The `docs/` folder contains information about the model. An overall summary is provided in the Model Summary word document. The accompanying file titled "Model Text Algorithm" was created by Stuart Wright and checked by Ewan Gray in the model validation process. This file contains a step by step text based explanation of the processes taking place in the model. The research team recommend reading this file and the published economic evaluation to gain a full understanding of the model. The documentation folder also contains diagrams outlining the structure of the model, tables of model parameters and their sources, and details about the process of updating the model parameters between version 0 and version 1.

Below are some key notes about the model for new users:
1) The model evalutes one screening strategy at a time (as defined by the user in the screen_strategy parameter). The results are therefore not incremental and need to be compared to those from another strategy which must also be run seperately
2) Like most discrete event simulation models this model is computationally expensive. We have taken steps to decrease run-time and will be looking into this further. 
3) DES models create natural variation in the outcomes and so require a large population to be sampled in the model to achieve stable results. Currently we recommend that users simulate results for at least 7 million women, particularly for the more variable risk-stratified strategies. We are investigating variance reduction techniques to reduce this number.
4) Given the large number of required model runs, the model splits the total sample into 10 sub-samples. If the model is interrupted mid-run then the results from some of these sub-samples may be retrieved to save running the whole sample again
5) Key data in this model include a sample of ~15,000 women's Volpara breast density estimates, estimated 10 year risk of breast cancer, and estimated lifetime risk of breast cancer. This data was based on real data from the PROCAS 2 study. In order to create a model version which could widely be shared while being sensitive to data protection concerns, the researchers created a synthetic dataset of these variables basedon the original data. This used the synthpop package in R which maintains the original structure of the data. 

## Folder structure

```
MancBreastScreen/
├── MancBreastScreen.Rproj   Open this first: all paths are relative to this folder
├── DESCRIPTION              Package metadata (authors, licence, dependencies)
├── R/                       Function definitions
│   ├── MancBreastScreen_functions.R   Core model functions (screening, stage, survival, QALYs)
│   ├── des_vectorised_df.R            Vectorised discrete event simulation (run_des_vectorised)
│   └── risksample_functions.R         Creation of the simulated cohort of women
├── parameters/              Model inputs
│   ├── params.R                       Base-case parameters
│   ├── psa_params.R                   PSA distributions
│   └── psa_params_intervals.R         Wide PSA distributions for GAM estimation
├── scripts/                 Scripts that run the model
│   ├── MancBreastScreen_main_script.R Runs strategies sequentially
│   └── MancBreastScreen_parallel.R    Runs strategies in parallel (faster)
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

1. Open `MancBreastScreen.Rproj` in RStudio so the working directory is the project root.
2. Set the options in the `controls` list at the top of `scripts/MancBreastScreen_parallel.R`
   (or `scripts/MancBreastScreen_main_script.R`).
3. Run the script. Results are written to `output/`.
4. Summarise the results with the scripts in `analysis/`.

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

All AI-assisted changes were reviewed by the research team, who take full
responsibility for the code and its outputs.
