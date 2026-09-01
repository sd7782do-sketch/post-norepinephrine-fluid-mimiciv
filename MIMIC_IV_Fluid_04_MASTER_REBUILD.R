# Heart & Lung major revision: complete MIMIC-IV cohort rebuild and analysis
# Version 1.1 - 2026-09-01
#
# Primary design
#   Population: adults with a sepsis ICD code and norepinephrine initiated
#               during ICU hours 0-24.
#   Landmark:   6 hours after norepinephrine initiation; patients must be alive
#               and still observable in the ICU at the landmark.
#   Exposure:   crystalloid/albumin volume during [NE, NE+6h), using the
#               explicit item map below; primary contrast <=1000 vs >1000 mL.
#   Baseline:   information strictly before NE initiation.
#   Primary:    death from the landmark through day 28 after NE initiation.
#   Inference:  stabilized IPTW truncated at the 1st/99th
#               percentiles, doubly adjusted survey logistic regression with
#               subject_id as the clustering unit.
#
# Required packages: DBI, duckdb, data.table, survey
# Raw files expected directly under MIMIC_ROOT as *.csv.gz.

MIMIC_ROOT <- ""
REPORT_PARENT <- ""

run_fluid_master_rebuild <- function(root = "", report_parent = "") {
  pkgs <- c("DBI", "duckdb", "data.table", "survey")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Installare prima: ", paste(missing_pkgs, collapse = ", "))

  profile <- Sys.getenv("USERPROFILE", unset = path.expand("~"))
  cloud <- Sys.getenv(c("OneDrive", "OneDriveConsumer", "OneDriveCommercial"), unset = "")
  desktops <- unique(c(file.path(profile, "Desktop"), file.path(profile, "Scrivania"),
                       file.path(cloud[nzchar(cloud)], "Desktop"),
                       file.path(cloud[nzchar(cloud)], "Scrivania")))
  if (!nzchar(root)) {
    candidates <- unlist(lapply(desktops[dir.exists(desktops)], function(z) {
      x <- list.dirs(z, recursive = FALSE, full.names = TRUE)
      x[tolower(basename(x)) == "mimiciv"]
    }), use.names = FALSE)
    candidates <- unique(candidates)
    if (length(candidates) != 1) stop("Impostare MIMIC_ROOT all'inizio dello script.")
    root <- candidates
  }
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!nzchar(report_parent)) {
    report_parent <- if (dir.exists(file.path(profile, "Documents")))
      file.path(profile, "Documents") else dirname(root)
  }
  parent <- file.path(report_parent, "MIMIC_IV_Fluid_MASTER_REBUILD")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(parent, format(Sys.time(), "%Y%m%d_%H%M%S"))
  if (!dir.create(out)) stop("Impossibile creare la cartella di output: ", out)
  db_path <- file.path(out, "master_rebuild.duckdb")

  log_lines <- character()
  say <- function(...) {
    z <- paste0(...)
    log_lines <<- c(log_lines, z)
    message(z)
  }
  report_file <- file.path(out, "00_run_report.txt")
  on.exit(writeLines(log_lines, report_file, useBytes = TRUE), add = TRUE)
  say("MIMIC-IV FLUID MASTER REBUILD v1.1")
  say("Root: ", root)
  say("Output: ", out)
  say("All baseline covariates use [NE-24h, NE); exposure uses [NE, NE+6h).")

  raw_names <- c("patients", "admissions", "icustays", "diagnoses_icd",
                 "d_items", "inputevents", "chartevents", "labevents",
                 "procedureevents")
  raw <- setNames(file.path(root, paste0(raw_names, ".csv.gz")), raw_names)
  if (any(!file.exists(raw))) {
    stop("File MIMIC mancanti: ", paste(basename(raw[!file.exists(raw)]), collapse = ", "))
  }
  qpath <- function(x) gsub("'", "''", normalizePath(x, winslash = "/", mustWork = TRUE), fixed = TRUE)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
  for (nm in raw_names) {
    DBI::dbExecute(con, sprintf(
      "CREATE OR REPLACE VIEW %s AS SELECT * FROM read_csv_auto('%s', header=true, sample_size=20000)",
      nm, qpath(raw[[nm]])))
  }
  get <- function(sql) DBI::dbGetQuery(con, sql)
  put <- function(x, filename) data.table::fwrite(data.table::as.data.table(x),
                                                   file.path(out, filename), na = "")

  # Fixed, auditable exposure map. No dextrose carriers, hypotonic/hypertonic
  # saline, irrigants, nutrition, blood products or medication carrier volumes.
  fluid_map <- data.frame(
    itemid = c(225158, 225828, 220861, 220863, 220862, 220864),
    fluid_class = c("0.9% sodium chloride", "lactated Ringer's",
                    "albumin 20%", "albumin 4%", "albumin 25%", "albumin 5%"),
    include_primary = 1L,
    stringsAsFactors = FALSE
  )
  DBI::dbWriteTable(con, "fluid_map", fluid_map, overwrite = TRUE)
  map_check <- get("SELECT f.*,d.label,d.category,d.unitname FROM fluid_map f LEFT JOIN d_items d USING(itemid) ORDER BY itemid")
  put(map_check, "01_fluid_item_mapping.csv")
  if (any(is.na(map_check$label))) stop("Uno o più itemid della mappa fluidi non sono presenti in d_items.")

  say("1/7 - Costruzione del denominatore, indice NE e landmark...")
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE source_index AS
    WITH sepsis_adm AS (
      SELECT DISTINCT hadm_id FROM diagnoses_icd
      WHERE (icd_version=9 AND upper(replace(icd_code,'.','')) IN ('99591','99592','78552'))
         OR (icd_version=10 AND (upper(replace(icd_code,'.','')) LIKE 'A40%'
             OR upper(replace(icd_code,'.','')) LIKE 'A41%'
             OR upper(replace(icd_code,'.','')) IN ('R6520','R6521')))
    ), ne AS (
      SELECT stay_id,min(starttime) AS ne_starttime
      FROM inputevents
      WHERE itemid=221906
        AND lower(coalesce(CAST(statusdescription AS VARCHAR),'')) NOT LIKE '%rewritten%'
      GROUP BY stay_id
    )
    SELECT i.subject_id,i.hadm_id,i.stay_id,i.intime,i.outtime,
      p.anchor_age+date_part('year',i.intime)-p.anchor_year AS age,
      CASE WHEN lower(CAST(p.gender AS VARCHAR)) IN ('m','male') THEN 1 ELSE 0 END AS male,
      a.admittime,a.dischtime,a.deathtime,p.dod,ne.ne_starttime,
      ne.ne_starttime+INTERVAL 6 HOUR AS landmark_time,
      epoch(ne.ne_starttime-i.intime)/3600.0 AS hours_icu_to_ne
    FROM icustays i JOIN patients p USING(subject_id)
    JOIN admissions a USING(hadm_id) JOIN sepsis_adm s USING(hadm_id)
    JOIN ne USING(stay_id)")
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE timed_adults AS
    SELECT * FROM source_index WHERE age>=18 AND ne_starttime>=intime
      AND ne_starttime<intime+INTERVAL 25 HOUR")

  say("2/7 - Estrazione dei valori basali pre-NE...")
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE baseline_vitals AS
    SELECT t.stay_id,
      median(CASE WHEN ce.itemid IN (220052,220181,225312) AND ce.valuenum BETWEEN 20 AND 250 THEN ce.valuenum END) AS base_map,
      median(CASE WHEN ce.itemid=220045 AND ce.valuenum BETWEEN 20 AND 300 THEN ce.valuenum END) AS base_hr,
      min(CASE WHEN ce.itemid IN (220052,220181,225312) AND ce.valuenum BETWEEN 20 AND 250 THEN ce.valuenum END) AS map_min,
      min(CASE WHEN ce.itemid=220739 THEN ce.valuenum END)
       +min(CASE WHEN ce.itemid=223900 THEN ce.valuenum END)
       +min(CASE WHEN ce.itemid=223901 THEN ce.valuenum END) AS gcs_min
    FROM timed_adults t LEFT JOIN chartevents ce ON ce.stay_id=t.stay_id
      AND ce.charttime>=t.ne_starttime-INTERVAL 24 HOUR AND ce.charttime<t.ne_starttime
      AND ce.itemid IN (220052,220181,225312,220045,220739,223900,223901)
    GROUP BY t.stay_id")
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE baseline_labs AS
    WITH x AS (
      SELECT t.stay_id,le.itemid,le.valuenum,le.charttime,
        row_number() OVER(PARTITION BY t.stay_id,le.itemid ORDER BY le.charttime DESC) AS rn
      FROM timed_adults t JOIN labevents le ON le.hadm_id=t.hadm_id
      AND le.charttime>=t.ne_starttime-INTERVAL 24 HOUR AND le.charttime<t.ne_starttime
      AND ((le.itemid IN (50813,52442) AND le.valuenum BETWEEN 0.1 AND 40)
        OR (le.itemid=50912 AND le.valuenum BETWEEN 0.1 AND 30)
        OR (le.itemid=51265 AND le.valuenum BETWEEN 1 AND 2000)
        OR (le.itemid=50885 AND le.valuenum BETWEEN 0.1 AND 80))
    )
    SELECT stay_id,
      max(CASE WHEN itemid IN (50813,52442) AND rn=1 THEN valuenum END) AS base_lactate,
      max(CASE WHEN itemid=50912 AND rn=1 THEN valuenum END) AS base_creatinine,
      min(CASE WHEN itemid=51265 THEN valuenum END) AS platelet_min,
      max(CASE WHEN itemid=50885 THEN valuenum END) AS bilirubin_max,
      max(CASE WHEN itemid=50912 THEN valuenum END) AS creatinine_max
    FROM x GROUP BY stay_id")
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE baseline_ne AS
    SELECT t.stay_id,
      max(CASE WHEN i.itemid=221906 AND lower(coalesce(CAST(i.rateuom AS VARCHAR),''))='mcg/kg/min'
          THEN i.rate END) AS ne_rate_preindex_max,
      max(CASE WHEN i.itemid=221906 AND i.starttime=t.ne_starttime
          AND lower(coalesce(CAST(i.rateuom AS VARCHAR),''))='mcg/kg/min' THEN i.rate END) AS initial_ne_rate
    FROM timed_adults t LEFT JOIN inputevents i ON i.stay_id=t.stay_id
      AND i.itemid=221906 AND i.starttime<=t.ne_starttime
      AND coalesce(i.endtime,i.starttime)>=t.ne_starttime-INTERVAL 24 HOUR
      AND lower(coalesce(CAST(i.statusdescription AS VARCHAR),'')) NOT LIKE '%rewritten%'
    GROUP BY t.stay_id")
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE baseline_vent AS
    SELECT t.stay_id,max(CASE WHEN lower(d.label)='invasive ventilation'
      AND p.starttime<t.ne_starttime
      AND coalesce(p.endtime,p.starttime)>=t.ne_starttime THEN 1 ELSE 0 END) AS baseline_invasive_vent
    FROM timed_adults t LEFT JOIN procedureevents p ON p.stay_id=t.stay_id
    LEFT JOIN d_items d ON d.itemid=p.itemid AND lower(d.label)='invasive ventilation'
    GROUP BY t.stay_id")

  # Pro-rate infusions crossing a window boundary; boluses are assigned by start time.
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE fluid_windows AS
    WITH e AS (
      SELECT t.stay_id,t.ne_starttime,i.starttime,coalesce(i.endtime,i.starttime) AS endtime,
        i.amount,greatest(0,epoch(coalesce(i.endtime,i.starttime)-i.starttime))/3600.0 AS duration_h
      FROM timed_adults t JOIN inputevents i USING(stay_id) JOIN fluid_map f USING(itemid)
      WHERE i.starttime<t.ne_starttime+INTERVAL 6 HOUR
        AND coalesce(i.endtime,i.starttime)>=t.ne_starttime-INTERVAL 6 HOUR
        AND i.amount>=0
        AND lower(coalesce(CAST(i.statusdescription AS VARCHAR),'')) NOT LIKE '%rewritten%'
    ), z AS (
      SELECT *,
        greatest(0,epoch(least(endtime,ne_starttime)-greatest(starttime,ne_starttime-INTERVAL 6 HOUR)))/3600.0 AS pre_h,
        greatest(0,epoch(least(endtime,ne_starttime+INTERVAL 6 HOUR)-greatest(starttime,ne_starttime)))/3600.0 AS post_h
      FROM e)
    SELECT t.stay_id,
      coalesce(sum(CASE WHEN z.duration_h>0 THEN z.amount*z.pre_h/z.duration_h
        WHEN z.starttime>=z.ne_starttime-INTERVAL 6 HOUR AND z.starttime<z.ne_starttime THEN z.amount ELSE 0 END),0) AS pre_ne_fluid_6h_ml,
      coalesce(sum(CASE WHEN z.duration_h>0 THEN z.amount*z.post_h/z.duration_h
        WHEN z.starttime>=z.ne_starttime AND z.starttime<z.ne_starttime+INTERVAL 6 HOUR THEN z.amount ELSE 0 END),0) AS post_ne_fluid_6h_ml
    FROM timed_adults t LEFT JOIN z USING(stay_id) GROUP BY t.stay_id")

  say("3/7 - Charlson, sede di infezione e modified SOFA...")
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE infection_flags AS
    SELECT t.stay_id,
      max(CASE WHEN regexp_matches(upper(replace(d.icd_code,'.','')),
        '^(480|481|482|483|484|485|486|J09|J10|J11|J12|J13|J14|J15|J16|J17|J18)') THEN 1 ELSE 0 END) AS pneumonia,
      max(CASE WHEN regexp_matches(upper(replace(d.icd_code,'.','')),
        '^(590|595|5990|N10|N11|N12|N30|N39)') THEN 1 ELSE 0 END) AS urinary_source,
      max(CASE WHEN regexp_matches(upper(replace(d.icd_code,'.','')),
        '^(540|541|542|K35|K36|K37|K55|K57|K63|K65|K80|K81|K83)') THEN 1 ELSE 0 END) AS abdominal_source,
      max(CASE WHEN regexp_matches(upper(replace(d.icd_code,'.','')),
        '^(038|A40|A41|R652)') THEN 1 ELSE 0 END) AS bloodstream_sepsis_code
    FROM timed_adults t LEFT JOIN diagnoses_icd d USING(hadm_id) GROUP BY t.stay_id")

  char_file <- file.path(out, "MIT_LCP_charlson_source.sql")
  char_url <- "https://raw.githubusercontent.com/MIT-LCP/mimic-code/refs/heads/main/mimic-iv/concepts_postgres/comorbidity/charlson.sql"
  dl <- try(utils::download.file(char_url, char_file, mode = "wb", quiet = TRUE), silent = TRUE)
  if (inherits(dl, "try-error") || !file.exists(char_file)) {
    local_candidates <- list.files(report_parent, pattern = "MIT_LCP_charlson_source[.]sql$",
                                   recursive = TRUE, full.names = TRUE)
    if (!length(local_candidates)) stop("Charlson SQL non scaricabile e nessuna copia locale trovata.")
    file.copy(local_candidates[1], char_file, overwrite = TRUE)
  }
  char_lines <- readLines(char_file, warn = FALSE)
  char_start <- grep("CREATE TABLE mimiciv_derived\\.charlson AS", char_lines)[1]
  if (is.na(char_start)) stop("Formato inatteso del codice Charlson MIT-LCP.")
  char_lines[char_start] <- sub("^.*CREATE TABLE mimiciv_derived\\.charlson AS", "", char_lines[char_start])
  char_sql <- paste(char_lines[char_start:length(char_lines)], collapse = "\n")
  char_sql <- gsub("mimiciv_hosp\\.diagnoses_icd", "diagnoses_icd", char_sql)
  char_sql <- gsub("mimiciv_hosp\\.admissions", "admissions", char_sql)
  char_sql <- gsub("mimiciv_derived\\.age", "age_derived", char_sql)
  DBI::dbExecute(con, "CREATE OR REPLACE VIEW age_derived AS SELECT a.hadm_id,
    p.anchor_age+date_part('year',a.admittime)-p.anchor_year AS age
    FROM admissions a JOIN patients p USING(subject_id)")
  DBI::dbExecute(con, paste0("CREATE OR REPLACE TABLE charlson_all AS ", char_sql))

  # Five-component pre-index score: coagulation, liver, cardiovascular, CNS,
  # renal-creatinine. Respiration and urine output are intentionally omitted.
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE preindex_score AS
    SELECT t.stay_id,
      CASE WHEN l.platelet_min IS NULL THEN NULL WHEN l.platelet_min<20 THEN 4 WHEN l.platelet_min<50 THEN 3 WHEN l.platelet_min<100 THEN 2 WHEN l.platelet_min<150 THEN 1 ELSE 0 END AS score_coag,
      CASE WHEN l.bilirubin_max IS NULL THEN NULL WHEN l.bilirubin_max>=12 THEN 4 WHEN l.bilirubin_max>=6 THEN 3 WHEN l.bilirubin_max>=2 THEN 2 WHEN l.bilirubin_max>=1.2 THEN 1 ELSE 0 END AS score_liver,
      CASE WHEN n.ne_rate_preindex_max>0.1 THEN 4 WHEN n.ne_rate_preindex_max>0 THEN 3 WHEN v.map_min<70 THEN 1 WHEN v.map_min IS NOT NULL THEN 0 ELSE NULL END AS score_cv,
      CASE WHEN v.gcs_min IS NULL THEN NULL WHEN v.gcs_min<6 THEN 4 WHEN v.gcs_min<10 THEN 3 WHEN v.gcs_min<13 THEN 2 WHEN v.gcs_min<15 THEN 1 ELSE 0 END AS score_cns,
      CASE WHEN l.creatinine_max IS NULL THEN NULL WHEN l.creatinine_max>=5 THEN 4 WHEN l.creatinine_max>=3.5 THEN 3 WHEN l.creatinine_max>=2 THEN 2 WHEN l.creatinine_max>=1.2 THEN 1 ELSE 0 END AS score_renal
    FROM timed_adults t LEFT JOIN baseline_vitals v USING(stay_id)
    LEFT JOIN baseline_labs l USING(stay_id) LEFT JOIN baseline_ne n USING(stay_id)")

  say("4/7 - Applicazione dei criteri di eleggibilità e costruzione outcome...")
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE analytic_sql AS
    WITH joined AS (
      SELECT t.*,v.base_map,v.base_hr,l.base_lactate,l.base_creatinine,
        n.initial_ne_rate,coalesce(b.baseline_invasive_vent,0) AS baseline_invasive_vent,
        f.pre_ne_fluid_6h_ml,f.post_ne_fluid_6h_ml,
        ch.charlson_comorbidity_index AS charlson,
        inf.pneumonia,inf.urinary_source,inf.abdominal_source,inf.bloodstream_sepsis_code,
        s.score_coag,s.score_liver,s.score_cv,s.score_cns,s.score_renal,
        coalesce(s.score_coag,0)+coalesce(s.score_liver,0)+coalesce(s.score_cv,0)+coalesce(s.score_cns,0)+coalesce(s.score_renal,0) AS modified_sofa_5c,
        (s.score_coag IS NOT NULL)::INTEGER+(s.score_liver IS NOT NULL)::INTEGER+(s.score_cv IS NOT NULL)::INTEGER+(s.score_cns IS NOT NULL)::INTEGER+(s.score_renal IS NOT NULL)::INTEGER AS modified_sofa_components,
        CASE WHEN t.deathtime IS NOT NULL THEN t.deathtime
             WHEN t.dod IS NOT NULL THEN CAST(t.dod AS TIMESTAMP)+INTERVAL 1 DAY-INTERVAL 1 SECOND
             ELSE NULL END AS death_time
      FROM timed_adults t LEFT JOIN baseline_vitals v USING(stay_id)
      LEFT JOIN baseline_labs l USING(stay_id) LEFT JOIN baseline_ne n USING(stay_id)
      LEFT JOIN baseline_vent b USING(stay_id) LEFT JOIN fluid_windows f USING(stay_id)
      LEFT JOIN charlson_all ch USING(hadm_id) LEFT JOIN infection_flags inf USING(stay_id)
      LEFT JOIN preindex_score s USING(stay_id)
    )
    SELECT *,
      CASE WHEN post_ne_fluid_6h_ml<=1000 THEN 1 ELSE 0 END AS restrictive_1000,
      CASE WHEN death_time>=landmark_time AND death_time<ne_starttime+INTERVAL 28 DAY THEN 1 ELSE 0 END AS death28_landmark,
      CASE WHEN death_time>=landmark_time AND death_time<ne_starttime+INTERVAL 30 DAY THEN 1 ELSE 0 END AS death30_landmark,
      CASE WHEN deathtime>=landmark_time AND deathtime<=dischtime THEN 1 ELSE 0 END AS hospital_death_landmark
    FROM joined
    WHERE base_map IS NOT NULL AND base_hr IS NOT NULL AND base_lactate IS NOT NULL
      AND base_creatinine IS NOT NULL AND initial_ne_rate IS NOT NULL
      AND outtime>=landmark_time
      AND (death_time IS NULL OR death_time>=landmark_time)")

  flow <- get("SELECT * FROM (VALUES
    (1,'All ICU stays',(SELECT count(*) FROM icustays)),
    (2,'Adults age >=18',(SELECT count(*) FROM icustays i JOIN patients p USING(subject_id) WHERE p.anchor_age+date_part('year',i.intime)-p.anchor_year>=18)),
    (3,'Adults with sepsis ICD and norepinephrine',(SELECT count(*) FROM source_index WHERE age>=18)),
    (4,'NE initiated during ICU hours 0-24',(SELECT count(*) FROM timed_adults)),
    (5,'Complete pre-NE MAP, HR, lactate, creatinine and initial NE rate',(SELECT count(*) FROM timed_adults t JOIN baseline_vitals v USING(stay_id) JOIN baseline_labs l USING(stay_id) JOIN baseline_ne n USING(stay_id) WHERE v.base_map IS NOT NULL AND v.base_hr IS NOT NULL AND l.base_lactate IS NOT NULL AND l.base_creatinine IS NOT NULL AND n.initial_ne_rate IS NOT NULL)),
    (6,'Alive and in ICU at 6-hour landmark',(SELECT count(*) FROM analytic_sql))) x(step_order,step,n) ORDER BY step_order")
  put(flow, "02_screening_flow.csv")

  # Strict incident invasive ventilation after the landmark and through hour 72.
  DBI::dbExecute(con, "CREATE OR REPLACE TABLE incident_vent AS
    SELECT a.stay_id,min(p.starttime) AS incident_vent_time
    FROM analytic_sql a JOIN procedureevents p USING(stay_id) JOIN d_items d USING(itemid)
    WHERE lower(d.label)='invasive ventilation' AND p.starttime>=a.landmark_time
      AND p.starttime<a.ne_starttime+INTERVAL 72 HOUR GROUP BY a.stay_id")
  analytic <- data.table::as.data.table(get("SELECT a.*,
    CASE WHEN a.baseline_invasive_vent=0 AND v.incident_vent_time IS NOT NULL THEN 1 ELSE 0 END AS incident_vent_72h
    FROM analytic_sql a LEFT JOIN incident_vent v USING(stay_id) ORDER BY a.stay_id"))
  if (nrow(analytic) != data.table::uniqueN(analytic$stay_id)) stop("Errore: dataset non univoco per stay_id.")
  put(analytic, "03_analytic_dataset_unweighted.csv")

  say("5/7 - Propensity score, pesi, bilanciamento e modelli...")
  # Missing modified-score components are represented by a zero contribution plus
  # the observed-component count. Charlson and infection flags are complete by design.
  analytic[, `:=`(
    charlson = data.table::fifelse(is.na(charlson), 0, charlson),
    pneumonia = data.table::fifelse(is.na(pneumonia), 0, pneumonia),
    urinary_source = data.table::fifelse(is.na(urinary_source), 0, urinary_source),
    abdominal_source = data.table::fifelse(is.na(abdominal_source), 0, abdominal_source),
    bloodstream_sepsis_code = data.table::fifelse(is.na(bloodstream_sepsis_code), 0, bloodstream_sepsis_code)
  )]
  analytic[, log1p_pre_ne_fluid_6h_ml := log1p(pre_ne_fluid_6h_ml)]
  covars <- c("age", "male", "hours_icu_to_ne", "base_map", "base_hr",
              "base_lactate", "base_creatinine", "baseline_invasive_vent",
              "charlson", "pre_ne_fluid_6h_ml", "log1p_pre_ne_fluid_6h_ml", "initial_ne_rate",
              "modified_sofa_5c", "modified_sofa_components", "pneumonia",
              "urinary_source", "abdominal_source", "bloodstream_sepsis_code")
  rhs <- paste(c("scale(age)", "male", "scale(hours_icu_to_ne)", "scale(base_map)",
                 "scale(base_hr)", "scale(log1p(base_lactate))", "scale(log1p(base_creatinine))",
                 "baseline_invasive_vent", "scale(charlson)", "scale(pre_ne_fluid_6h_ml)",
                 "scale(log1p_pre_ne_fluid_6h_ml)",
                 "scale(log1p(initial_ne_rate))", "scale(modified_sofa_5c)",
                 "modified_sofa_components", "pneumonia", "urinary_source",
                 "abdominal_source", "bloodstream_sepsis_code"), collapse = "+")
  ps_fit <- stats::glm(stats::as.formula(paste("restrictive_1000~", rhs)),
                       data = analytic, family = stats::binomial())
  analytic[, ps := pmin(pmax(stats::predict(ps_fit, type = "response"), 0.001), 0.999)]
  p_treat <- mean(analytic$restrictive_1000)
  analytic[, sw := data.table::fifelse(restrictive_1000 == 1, p_treat/ps, (1-p_treat)/(1-ps))]
  limits <- stats::quantile(analytic$sw, c(.01, .99), na.rm = TRUE)
  analytic[, sw_trunc := pmin(pmax(sw, limits[1]), limits[2])]
  put(analytic, "04_analytic_dataset_weighted.csv")

  weight_diag <- data.table::data.table(
    n = nrow(analytic), patients = data.table::uniqueN(analytic$subject_id),
    restrictive_n = sum(analytic$restrictive_1000 == 1),
    liberal_n = sum(analytic$restrictive_1000 == 0),
    ps_min = min(analytic$ps), ps_max = max(analytic$ps),
    weight_p01 = limits[1], weight_median = median(analytic$sw_trunc),
    weight_p99 = limits[2], weight_max = max(analytic$sw_trunc),
    effective_sample_size = sum(analytic$sw_trunc)^2/sum(analytic$sw_trunc^2)
  )
  put(weight_diag, "05_weight_diagnostics.csv")

  weighted_mean <- function(x, w) sum(x*w, na.rm=TRUE)/sum(w[!is.na(x)])
  balance_one <- function(v) {
    x <- analytic[[v]]; g <- analytic$restrictive_1000; w <- analytic$sw_trunc
    m1 <- mean(x[g==1],na.rm=TRUE); m0 <- mean(x[g==0],na.rm=TRUE)
    s <- sqrt((stats::var(x[g==1],na.rm=TRUE)+stats::var(x[g==0],na.rm=TRUE))/2)
    wm1 <- weighted_mean(x[g==1],w[g==1]); wm0 <- weighted_mean(x[g==0],w[g==0])
    wv <- function(z,ww) sum(ww*(z-weighted_mean(z,ww))^2,na.rm=TRUE)/sum(ww[!is.na(z)])
    ws <- sqrt((wv(x[g==1],w[g==1])+wv(x[g==0],w[g==0]))/2)
    data.frame(variable=v,mean_restrictive=m1,mean_liberal=m0,smd_unweighted=(m1-m0)/s,
               weighted_mean_restrictive=wm1,weighted_mean_liberal=wm0,smd_weighted=(wm1-wm0)/ws)
  }
  balance <- data.table::rbindlist(lapply(covars, balance_one), fill=TRUE)
  put(balance, "06_covariate_balance.csv")

  fit_outcome <- function(data, outcome, label, subset_expr = NULL) {
    d <- data
    if (!is.null(subset_expr)) d <- d[eval(substitute(subset_expr), d, parent.frame())]
    design <- survey::svydesign(ids=~subject_id, weights=~sw_trunc, data=d)
    form <- stats::as.formula(paste(outcome, "~restrictive_1000+", rhs))
    fit <- survey::svyglm(form, design=design, family=stats::quasibinomial())
    b <- stats::coef(fit)["restrictive_1000"]
    se <- sqrt(stats::vcov(fit)["restrictive_1000","restrictive_1000"])
    p <- summary(fit)$coefficients["restrictive_1000","Pr(>|t|)"]
    pred1 <- d; pred1$restrictive_1000 <- 1
    pred0 <- d; pred0$restrictive_1000 <- 0
    r1 <- weighted_mean(stats::predict(fit,newdata=pred1,type="response"),d$sw_trunc)
    r0 <- weighted_mean(stats::predict(fit,newdata=pred0,type="response"),d$sw_trunc)
    data.frame(analysis=label,outcome=outcome,n=nrow(d),events=sum(d[[outcome]]==1,na.rm=TRUE),
      OR=exp(b),CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se),p_value=p,
      adjusted_risk_restrictive=r1,adjusted_risk_liberal=r0,risk_difference=r1-r0)
  }
  results <- data.table::rbindlist(list(
    fit_outcome(analytic,"death28_landmark","Primary: 28-day mortality"),
    fit_outcome(analytic,"death30_landmark","Sensitivity: 30-day mortality"),
    fit_outcome(analytic,"hospital_death_landmark","Sensitivity: hospital mortality"),
    fit_outcome(analytic[baseline_invasive_vent==0],"incident_vent_72h","Incident invasive ventilation")
  ), fill=TRUE)
  put(results, "07_primary_and_secondary_models.csv")

  # First eligible ICU stay per patient.
  first_stay <- analytic[order(subject_id, ne_starttime), .SD[1], by=subject_id]
  first_result <- fit_outcome(first_stay,"death28_landmark","Sensitivity: first eligible stay")
  put(first_result, "08_first_stay_sensitivity.csv")

  say("6/7 - Soglie, quartili e dose-response...")
  threshold_results <- list()
  for (cut in c(500,1000,1500)) {
    d <- data.table::copy(analytic)
    d[, restrictive_1000 := as.integer(post_ne_fluid_6h_ml<=cut)]
    pf <- stats::glm(stats::as.formula(paste("restrictive_1000~",rhs)),data=d,family=stats::binomial())
    d[, ps:=pmin(pmax(stats::predict(pf,type="response"),.001),.999)]
    pt <- mean(d$restrictive_1000); d[,sw:=data.table::fifelse(restrictive_1000==1,pt/ps,(1-pt)/(1-ps))]
    lim <- stats::quantile(d$sw,c(.01,.99));d[,sw_trunc:=pmin(pmax(sw,lim[1]),lim[2])]
    z <- fit_outcome(d,"death28_landmark",paste0("Threshold <=",cut," mL"))
    z$threshold_ml <- cut; threshold_results[[as.character(cut)]] <- z
  }
  put(data.table::rbindlist(threshold_results,fill=TRUE), "09_threshold_sensitivity.csv")

  analytic[, fluid_quartile := cut(post_ne_fluid_6h_ml,
    breaks=unique(stats::quantile(post_ne_fluid_6h_ml,probs=seq(0,1,.25),na.rm=TRUE)),
    include.lowest=TRUE,ordered_result=TRUE)]
  quartiles <- analytic[,.(n=.N,events=sum(death28_landmark),
    fluid_min=min(post_ne_fluid_6h_ml),fluid_median=median(post_ne_fluid_6h_ml),
    fluid_max=max(post_ne_fluid_6h_ml)),by=fluid_quartile]
  put(quartiles,"10_fluid_quartiles.csv")
  spline_design <- survey::svydesign(ids=~subject_id,weights=~sw_trunc,data=analytic)
  spline_fit <- survey::svyglm(death28_landmark ~ splines::ns(log1p(post_ne_fluid_6h_ml),df=3)
    + scale(age)+male+scale(hours_icu_to_ne)+scale(base_map)+scale(base_hr)
    + scale(log1p(base_lactate))+scale(log1p(base_creatinine))+baseline_invasive_vent
    + scale(charlson)+scale(pre_ne_fluid_6h_ml)+scale(log1p_pre_ne_fluid_6h_ml)
    + scale(log1p(initial_ne_rate))
    + scale(modified_sofa_5c)+modified_sofa_components+pneumonia+urinary_source
    + abdominal_source+bloodstream_sepsis_code,
    design=spline_design,family=stats::quasibinomial())
  spline_summary <- data.frame(term=names(stats::coef(spline_fit)),estimate=stats::coef(spline_fit),
    std_error=sqrt(diag(stats::vcov(spline_fit))),row.names=NULL)
  spline_summary$OR <- exp(spline_summary$estimate)
  spline_summary$CI_low <- exp(spline_summary$estimate-1.96*spline_summary$std_error)
  spline_summary$CI_high <- exp(spline_summary$estimate+1.96*spline_summary$std_error)
  spline_terms <- grep("^splines::ns",names(stats::coef(spline_fit)))
  spline_beta <- stats::coef(spline_fit)[spline_terms]
  spline_vcov <- stats::vcov(spline_fit)[spline_terms,spline_terms,drop=FALSE]
  spline_wald <- as.numeric(t(spline_beta)%*%solve(spline_vcov,spline_beta))
  spline_summary$overall_spline_p <- NA_real_
  spline_summary$overall_spline_p[1] <- stats::pchisq(spline_wald,df=length(spline_terms),lower.tail=FALSE)
  put(spline_summary,"11_spline_model_cluster_robust.csv")

  # Marginal standardized risk curve. Confidence limits use the delta method
  # with the patient-clustered covariance matrix from svyglm.
  grid <- unique(as.numeric(stats::quantile(analytic$post_ne_fluid_6h_ml,
    probs=seq(.01,.99,length.out=60),na.rm=TRUE)))
  curve_rows <- lapply(grid,function(volume) {
    nd <- data.table::copy(analytic)
    nd[,post_ne_fluid_6h_ml:=volume]
    X <- stats::model.matrix(stats::delete.response(stats::terms(spline_fit)),data=nd)
    X <- X[,names(stats::coef(spline_fit)),drop=FALSE]
    eta <- as.numeric(X%*%stats::coef(spline_fit)); pr <- stats::plogis(eta)
    ww <- nd$sw_trunc; risk <- sum(ww*pr)/sum(ww)
    gradient <- colSums(X*(ww*pr*(1-pr)))/sum(ww)
    se_risk <- sqrt(as.numeric(t(gradient)%*%stats::vcov(spline_fit)%*%gradient))
    se_logit <- se_risk/(risk*(1-risk))
    data.frame(fluid_ml=volume,adjusted_risk=risk,
      CI_low=stats::plogis(stats::qlogis(risk)-1.96*se_logit),
      CI_high=stats::plogis(stats::qlogis(risk)+1.96*se_logit))
  })
  put(data.table::rbindlist(curve_rows),"12_spline_adjusted_risk_curve.csv")

  say("7/7 - Controlli finali e archivio...")
  outcome_audit <- analytic[,.(
    n=.N,patients=data.table::uniqueN(subject_id),death28=sum(death28_landmark),death30=sum(death30_landmark),
    hospital_death=sum(hospital_death_landmark),baseline_invasive_vent=sum(baseline_invasive_vent),
    baseline_nonventilated=sum(baseline_invasive_vent==0),incident_vent_72h=sum(incident_vent_72h),
    fluid_median=median(post_ne_fluid_6h_ml),fluid_q1=quantile(post_ne_fluid_6h_ml,.25),
    fluid_q3=quantile(post_ne_fluid_6h_ml,.75))]
  put(outcome_audit,"13_outcome_and_exposure_audit.csv")
  missingness <- data.table::rbindlist(lapply(names(analytic), function(v)
    data.frame(variable=v,n_missing=sum(is.na(analytic[[v]])),pct_missing=mean(is.na(analytic[[v]]))*100)))
  put(missingness,"14_missingness.csv")
  if (any(analytic$death28_landmark > analytic$death30_landmark)) stop("Controllo mortalità 28/30 giorni fallito.")
  if (any(analytic$post_ne_fluid_6h_ml < 0)) stop("Volume fluido negativo rilevato.")
  max_smd <- max(abs(balance$smd_weighted),na.rm=TRUE)
  say("Maximum absolute weighted SMD: ",sprintf("%.3f",max_smd),".")
  if (max_smd>.10) say("ATTENZIONE: almeno uno SMD pesato supera 0.10; non usare i risultati come definitivi.")
  say("Final cohort: ", nrow(analytic), " stays; ", data.table::uniqueN(analytic$subject_id), " patients.")
  say("Primary deaths: ", sum(analytic$death28_landmark), ".")
  say("Finished: ", format(Sys.time(), tz="UTC", usetz=TRUE))
  writeLines(log_lines, report_file, useBytes=TRUE)

  archive_files <- list.files(out, full.names=FALSE)
  archive <- file.path(parent,paste0(basename(out),".zip"))
  oldwd <- getwd(); on.exit(setwd(oldwd),add=TRUE); setwd(out)
  zip_result <- try(utils::zip(archive,files=archive_files),silent=TRUE)
  if (!inherits(zip_result,"try-error") && file.exists(archive)) {
    message("\nAnalisi completata. Invia questo ZIP:\n",archive)
  } else message("\nAnalisi completata. Comprimi e invia la cartella:\n",out)
  invisible(out)
}

run_fluid_master_rebuild(MIMIC_ROOT, REPORT_PARENT)
