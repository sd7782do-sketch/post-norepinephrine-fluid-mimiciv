Post-norepinephrine fluid exposure in MIMIC-IV

Reproducible cohort-construction, statistical-analysis, and figure-generation code accompanying the manuscript:

Post-norepinephrine fluid volume and clinical outcomes in septic shock: a landmark retrospective cohort analysis of MIMIC-IV

Study overview

This repository contains the code used to reconstruct an adult MIMIC-IV cohort with coded sepsis and norepinephrine initiated during the first 24 hours of an ICU stay. Cumulative crystalloid and albumin exposure during the first 6 hours after norepinephrine initiation was compared using a 6-hour landmark design.

The primary contrast was lower-volume exposure (<=1,000 mL) versus higher-volume exposure (>1,000 mL). The primary outcome was 28-day mortality after the landmark. Analyses used stabilized inverse-probability weights truncated at the 1st and 99th percentiles, covariate-adjusted survey logistic regression, and patient-clustered robust standard errors.

Repository contents

MIMIC_IV_Fluid_04_MASTER_REBUILD.R: source-level cohort reconstruction, covariate and outcome construction, propensity-score weighting, balance assessment, primary models, sensitivity analyses, threshold analyses, quartile summaries, and spline analysis.

HEART_LUNG_NEJM_PLOTS.R: generation of the four publication figures from the definitive analysis outputs.

Data source and access

The analysis was developed for MIMIC-IV version 3.0. MIMIC-IV is available through PhysioNet to users who have completed the required credentialing, training, and data-use agreement.

No patient-level MIMIC-IV data or derived individual-level datasets are included in this repository. Users must obtain MIMIC-IV independently through PhysioNet and comply with the applicable credentialed health-data license.

Software requirements

The analysis requires R and the following packages:

install.packages(c(
  "DBI", "duckdb", "data.table", "survey",
  "ggplot2", "scales"
))

The scripts were developed for execution in R/RStudio on Windows. Paths are configurable and are not hard-coded.

Required MIMIC-IV files

The following compressed source tables must be placed directly in the directory assigned to MIMIC_ROOT:

patients.csv.gz
admissions.csv.gz
icustays.csv.gz
diagnoses_icd.csv.gz
d_items.csv.gz
inputevents.csv.gz
chartevents.csv.gz
labevents.csv.gz
procedureevents.csv.gz

The reconstruction script retrieves the official MIT-LCP Charlson Comorbidity Index SQL definition from the public mimic-code repository. If it cannot be downloaded, a local copy named MIT_LCP_charlson_source.sql may be supplied in the report parent directory.

Running the analysis

Open MIMIC_IV_Fluid_04_MASTER_REBUILD.R.

Set the two paths at the beginning of the script:

MIMIC_ROOT <- "path/to/mimiciv"
REPORT_PARENT <- "path/to/output_parent"

Run the entire script.

The script creates a time-stamped results directory under:

MIMIC_IV_Fluid_MASTER_REBUILD/

The expected definitive reconstruction yields 2,696 eligible ICU stays from 2,493 patients, with 1,117 primary outcome events. These values are reproducibility checks, not distributed patient-level data.

Generating the figures

After the master reconstruction has completed:

Open HEART_LUNG_NEJM_PLOTS.R.

Set RESULTS_DIR to the definitive time-stamped output directory.

Optionally set FIGURE_DIR; if left empty, a NEJM_Figures subdirectory is created.

Run the entire script.

RESULTS_DIR <- "path/to/definitive/results"
FIGURE_DIR <- ""

The script produces four LZW-compressed TIFF figures at 600 dpi.

Reproducibility and data protection

Generated .duckdb, .csv, .csv.gz, .rds, .RData, and ZIP files may contain restricted individual-level information and must not be committed to a public repository. Only code and non-disclosive aggregate documentation should be shared publicly.

License

The original code in this repository is released under the MIT License. The externally retrieved Charlson SQL is maintained by the MIT-LCP MIMIC Code Repository and remains subject to its original attribution and license.

Citation

If using this code, please cite the accompanying Heart & Lung manuscript and MIMIC-IV according to the citation instructions provided on PhysioNet. Publication details will be added after acceptance.
