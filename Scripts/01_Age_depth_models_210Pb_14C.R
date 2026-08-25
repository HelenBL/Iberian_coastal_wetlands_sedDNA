# ======================================================================
# COMBINED 210Pb + AMS 14C AGE-DEPTH MODELLING WITH rplum
# Compatible with: RPLUM_INPUT_AUDITED_WITH_SENSITIVITY.xlsx
# ======================================================================
#
# IMPORTANT
# ---------
# This script is written specifically for the audited workbook created
# from ALL_Datings(1).xlsx.
#
# It DOES NOT guess missing metadata.
# It will stop cleanly until the following fields are completed:
#
# CONFIG sheet:
#   - Confirmed_rPlum_date_sample_AD
#   - Confirmed_activity_unit
#
# Each Pb core sheet:
#   - rPlum_bottom_depth_cm_REQUIRED
#   - dry_bulk_density_g_cm3_REQUIRED
#   - sample_thickness_cm_REQUIRED
#
# The raw AMS 14C determinations are supplied directly to rplum via
# `otherdates` with cc = 1 (IntCal20).
#
# No marine reservoir correction or empirical old-carbon correction is
# applied.
#
# ======================================================================


# ----------------------------------------------------------------------
# 0. INSTALL PACKAGES ONCE
# ----------------------------------------------------------------------
# install.packages(c("readxl", "rplum", "rbacon"))


# ----------------------------------------------------------------------
# 1. LOAD PACKAGES
# ----------------------------------------------------------------------

library(readxl)
library(rplum)
library(rbacon)


# ----------------------------------------------------------------------
# 2. USER PATHS
# ----------------------------------------------------------------------

setwd("C:/Users/Elena Baños/Desktop/Doctorat/Doctorat/ALL_DATA/DATINGS/rPLUM")

INPUT_XLSX <- "RPLUM_INPUT_AUDITED_WITH_SENSITIVITY.xlsx"
RUN_DIR <- file.path(getwd(), "Plum_runs_combined")

if (!file.exists(INPUT_XLSX)) {
  stop(
    paste0(
      "Cannot find '", INPUT_XLSX, "'. ",
      "Put the audited Excel file in the working directory or ",
      "change INPUT_XLSX to the correct path."
    )
  )
}

dir.create(RUN_DIR, showWarnings = FALSE, recursive = TRUE)


# ----------------------------------------------------------------------
# 3. GLOBAL ANALYTICAL SETTINGS
# ----------------------------------------------------------------------

SSIZE <- 8000
BURNIN <- 1000
PROB <- 0.95
SEED_BASE <- 20260812

# Execution switches:
#
# First run:
#   RUN_BASELINE    <- TRUE
#   RUN_SENSITIVITY <- FALSE
#
# After the baseline models have been checked, sensitivity only:
#   RUN_BASELINE    <- FALSE
#   RUN_SENSITIVITY <- TRUE
#
# If both are TRUE, baseline models are run first and sensitivity runs follow.

RUN_BASELINE <- TRUE
RUN_SENSITIVITY <- TRUE


# ----------------------------------------------------------------------
# 4. READ AND VALIDATE CONFIG
# ----------------------------------------------------------------------

cfg <- read_excel(INPUT_XLSX, sheet = "CONFIG")
names(cfg) <- trimws(names(cfg))

required_cfg <- c(
  "Core",
  "Site",
  "Core_base_cm",
  "C14_sample",
  "C14_depth_cm",
  "C14_age_BP",
  "C14_error_1sigma",
  "C14_cc_IntCal20",
  "CFCS_formula_date_reference_AD",
  "Confirmed_rPlum_date_sample_AD",
  "Confirmed_activity_unit",
  "n_supp",
  "ra_case",
  "Model_thick_cm",
  "acc_mean_yr_per_cm",
  "acc_shape",
  "mem_mean",
  "mem_strength",
  "Source_date_note",
  "Run_status"
)

missing_cfg <- setdiff(required_cfg, names(cfg))

if (length(missing_cfg) > 0) {
  stop(
    paste0(
      "CONFIG sheet structure does not match the audited workbook.\n",
      "Missing column(s): ",
      paste(missing_cfg, collapse = ", ")
    )
  )
}

if (anyDuplicated(cfg$Core)) {
  stop("Duplicated Core values found in CONFIG.")
}

# Explicit numeric conversion.
numeric_cfg <- c(
  "Core_base_cm",
  "C14_depth_cm",
  "C14_age_BP",
  "C14_error_1sigma",
  "C14_cc_IntCal20",
  "CFCS_formula_date_reference_AD",
  "Confirmed_rPlum_date_sample_AD",
  "n_supp",
  "ra_case",
  "Model_thick_cm",
  "acc_mean_yr_per_cm",
  "acc_shape",
  "mem_mean",
  "mem_strength"
)

for (nm in numeric_cfg) {
  cfg[[nm]] <- suppressWarnings(as.numeric(cfg[[nm]]))
}


# ----------------------------------------------------------------------
# 5. VALIDATE THE 14C CONFIGURATION
# ----------------------------------------------------------------------

if (anyNA(cfg$C14_sample) ||
    anyNA(cfg$C14_depth_cm) ||
    anyNA(cfg$C14_age_BP) ||
    anyNA(cfg$C14_error_1sigma)) {
  stop("Missing raw AMS 14C information in CONFIG.")
}

if (any(cfg$C14_age_BP <= 0)) {
  stop("All conventional 14C ages must be > 0 yr BP.")
}

if (any(cfg$C14_error_1sigma <= 0)) {
  stop("All 14C analytical errors must be > 0.")
}

# cc = 1 means IntCal20 in Bacon/rplum.
if (anyNA(cfg$C14_cc_IntCal20) ||
    any(cfg$C14_cc_IntCal20 != 1)) {
  stop(
    "All C14_cc_IntCal20 values must equal 1 for this IntCal20 analysis."
  )
}


# ----------------------------------------------------------------------
# 6. VALIDATE THE REQUIRED rPLUM METADATA
# ----------------------------------------------------------------------
#
# IMPORTANT:
# CFCS_formula_date_reference_AD is retained only for provenance.
# It is NOT automatically substituted for the confirmed rplum date.sample.

missing_date <- is.na(cfg$Confirmed_rPlum_date_sample_AD)

if (any(missing_date)) {
  stop(
    paste0(
      "Confirmed_rPlum_date_sample_AD is still missing for: ",
      paste(cfg$Core[missing_date], collapse = ", "),
      ".\nDo not use CFCS_formula_date_reference_AD automatically. ",
      "Enter the confirmed calendar date at which the Pb-210 samples ",
      "were measured."
    )
  )
}

missing_unit <- is.na(cfg$Confirmed_activity_unit) |
  trimws(as.character(cfg$Confirmed_activity_unit)) == ""

if (any(missing_unit)) {
  stop(
    paste0(
      "Confirmed_activity_unit is still missing for: ",
      paste(cfg$Core[missing_unit], collapse = ", "),
      ".\nAllowed values: Bq/kg, Bq/g, dpm/g."
    )
  )
}

allowed_units <- c("Bq/kg", "Bq/g", "dpm/g")

if (any(!cfg$Confirmed_activity_unit %in% allowed_units)) {
  bad <- unique(
    cfg$Confirmed_activity_unit[
      !cfg$Confirmed_activity_unit %in% allowed_units
    ]
  )
  
  stop(
    paste0(
      "Unknown Confirmed_activity_unit value(s): ",
      paste(bad, collapse = ", "),
      ". Allowed values: Bq/kg, Bq/g, dpm/g."
    )
  )
}


# ----------------------------------------------------------------------
# 7. HELPER FUNCTIONS
# ----------------------------------------------------------------------

as_logical_include <- function(x) {
  
  if (is.logical(x)) {
    if (anyNA(x)) {
      stop("'include_in_rPlum' contains missing TRUE/FALSE values.")
    }
    return(x)
  }
  
  y <- tolower(trimws(as.character(x)))
  
  out <- rep(NA, length(y))
  
  out[y %in% c("true", "t", "yes", "y", "1")] <- TRUE
  out[y %in% c("false", "f", "no", "n", "0")] <- FALSE
  
  if (anyNA(out)) {
    stop(
      "'include_in_rPlum' must contain only TRUE/FALSE ",
      "(or equivalent yes/no values)."
    )
  }
  
  out
}


prepare_pb_data <- function(core, unit) {
  
  dat <- read_excel(INPUT_XLSX, sheet = core)
  names(dat) <- trimws(names(dat))
  
  # EXACT column names in RPLUM_INPUT_AUDITED_FROM_ALL_Datings.xlsx
  required <- c(
    "labID",
    "clientID",
    "reported_depth_cm",
    "total_Po210_activity_source",
    "lab_error_1sigma_source",
    "unsupported_Pb_source",
    "rPlum_bottom_depth_cm_REQUIRED",
    "dry_bulk_density_g_cm3_REQUIRED",
    "sample_thickness_cm_REQUIRED",
    "include_in_rPlum"
  )
  
  miss <- setdiff(required, names(dat))
  
  if (length(miss) > 0) {
    stop(
      paste0(
        core,
        " sheet does not match the audited workbook structure.\n",
        "Missing column(s): ",
        paste(miss, collapse = ", ")
      )
    )
  }
  
  dat$include_in_rPlum <- as_logical_include(dat$include_in_rPlum)
  dat <- dat[dat$include_in_rPlum, , drop = FALSE]
  
  if (nrow(dat) == 0) {
    stop(core, ": no Pb observations have include_in_rPlum = TRUE.")
  }
  
  numeric_cols <- c(
    "reported_depth_cm",
    "total_Po210_activity_source",
    "lab_error_1sigma_source",
    "unsupported_Pb_source",
    "rPlum_bottom_depth_cm_REQUIRED",
    "dry_bulk_density_g_cm3_REQUIRED",
    "sample_thickness_cm_REQUIRED"
  )
  
  for (nm in numeric_cols) {
    dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
  }
  
  # rPlum uses TOTAL measured Pb/Po proxy activity, not the unsupported
  # Pb values calculated previously by the CF:CS spreadsheet.
  #
  # Therefore unsupported_Pb_source is retained only for provenance and
  # is NOT passed to Plum.
  
  required_values <- c(
    "labID",
    "rPlum_bottom_depth_cm_REQUIRED",
    "dry_bulk_density_g_cm3_REQUIRED",
    "total_Po210_activity_source",
    "lab_error_1sigma_source",
    "sample_thickness_cm_REQUIRED"
  )
  
  incomplete <- which(
    !complete.cases(dat[, required_values, drop = FALSE])
  )
  
  if (length(incomplete) > 0) {
    
    cat(
      "\n============================================================\n",
      core,
      ": REQUIRED rPlum INPUT IS INCOMPLETE\n",
      "============================================================\n",
      sep = ""
    )
    
    print(
      dat[
        incomplete,
        c(
          "labID",
          "clientID",
          "reported_depth_cm",
          "rPlum_bottom_depth_cm_REQUIRED",
          "dry_bulk_density_g_cm3_REQUIRED",
          "sample_thickness_cm_REQUIRED"
        ),
        drop = FALSE
      ],
      row.names = FALSE
    )
    
    stop(
      paste0(
        core,
        ": fill all yellow REQUIRED fields in the audited workbook ",
        "before running rPlum. No value will be guessed."
      )
    )
  }
  
  # Physical validity checks.
  if (any(dat$rPlum_bottom_depth_cm_REQUIRED <= 0)) {
    stop(core, ": rPlum bottom depths must be > 0 cm.")
  }
  
  if (any(diff(dat$rPlum_bottom_depth_cm_REQUIRED) <= 0)) {
    stop(
      core,
      ": rPlum bottom depths must increase monotonically down-core."
    )
  }
  
  if (any(dat$dry_bulk_density_g_cm3_REQUIRED <= 0)) {
    stop(core, ": dry bulk density must be > 0 g cm^-3.")
  }
  
  if (any(dat$sample_thickness_cm_REQUIRED <= 0)) {
    stop(core, ": sample thickness must be > 0 cm.")
  }
  
  if (any(dat$total_Po210_activity_source <= 0)) {
    stop(core, ": total Po-210/Pb-210 proxy activity must be > 0.")
  }
  
  if (any(dat$lab_error_1sigma_source <= 0)) {
    stop(core, ": all activity 1-sigma errors must be > 0.")
  }
  
  
  # --------------------------------------------------------------------
  # Unit conversion
  # --------------------------------------------------------------------
  #
  # rplum accepts:
  # Bqkg = TRUE  -> Bq/kg
  # Bqkg = FALSE -> dpm/g
  #
  # If the confirmed source unit is Bq/g, convert to Bq/kg.
  
  if (unit == "Bq/kg") {
    
    activity <- dat$total_Po210_activity_source
    activity_error <- dat$lab_error_1sigma_source
    Bqkg_flag <- TRUE
    
  } else if (unit == "Bq/g") {
    
    activity <- dat$total_Po210_activity_source * 1000
    activity_error <- dat$lab_error_1sigma_source * 1000
    Bqkg_flag <- TRUE
    
  } else if (unit == "dpm/g") {
    
    activity <- dat$total_Po210_activity_source
    activity_error <- dat$lab_error_1sigma_source
    Bqkg_flag <- FALSE
    
  } else {
    
    stop(core, ": unsupported activity unit.")
  }
  
  
  # --------------------------------------------------------------------
  # Construct EXACT rplum 210Pb input
  # --------------------------------------------------------------------
  
  out <- data.frame(
    labID = as.character(dat$labID),
    depth.cm. = dat$rPlum_bottom_depth_cm_REQUIRED,
    density.g.cm.3. = dat$dry_bulk_density_g_cm3_REQUIRED,
    X210Pb.Bq.kg. = activity,
    sd.210Pb. = activity_error,
    thickness.cm. = dat$sample_thickness_cm_REQUIRED,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Ensure rows are ordered by bottom depth.
  out <- out[order(out$depth.cm.), , drop = FALSE]
  row.names(out) <- NULL
  
  list(
    data = out,
    Bqkg = Bqkg_flag
  )
}


prepare_c14_data <- function(config_row) {
  
  # Raw conventional AMS 14C age.
  # DO NOT pass the previously calibrated BCE/CE HPD ranges here.
  #
  # cc = 1 -> IntCal20.
  
  data.frame(
    labID = as.character(config_row$C14_sample),
    Age = as.numeric(config_row$C14_age_BP),
    Error = as.numeric(config_row$C14_error_1sigma),
    Depth = as.numeric(config_row$C14_depth_cm),
    cc = as.integer(config_row$C14_cc_IntCal20),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


# ----------------------------------------------------------------------
# 8. WRITE rPLUM INPUT FILES
# ----------------------------------------------------------------------

write_core_files <- function(config_row) {
  
  core <- as.character(config_row$Core)
  unit <- as.character(config_row$Confirmed_activity_unit)
  
  prepared <- prepare_pb_data(core, unit)
  c14 <- prepare_c14_data(config_row)
  
  core_dir <- file.path(RUN_DIR, core)
  dir.create(core_dir, showWarnings = FALSE, recursive = TRUE)
  
  pb_file <- file.path(core_dir, paste0(core, ".csv"))
  c14_file <- file.path(core_dir, paste0(core, "_C14.csv"))
  
  write.csv(
    prepared$data,
    file = pb_file,
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
  
  write.csv(
    c14,
    file = c14_file,
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
  
  list(
    core_dir = core_dir,
    pb_file = pb_file,
    c14_file = c14_file,
    Bqkg = prepared$Bqkg
  )
}


# ----------------------------------------------------------------------
# 9. PRE-FLIGHT VALIDATION
# ----------------------------------------------------------------------
#
# No MCMC is launched unless EVERY included Pb row has valid:
# - bottom depth
# - dry bulk density
# - sample thickness
# and each core has a confirmed measurement date and activity unit.

cat("\n============================================================\n")
cat("PRE-FLIGHT VALIDATION\n")
cat("============================================================\n")

prepared_files <- vector("list", nrow(cfg))
names(prepared_files) <- cfg$Core

for (i in seq_len(nrow(cfg))) {
  
  core <- as.character(cfg$Core[i])
  
  prepared_files[[core]] <- write_core_files(cfg[i, ])
  
  cat(
    "PASS:",
    core,
    "- rPlum Pb file and raw 14C file created.\n"
  )
}

cat("\nAll cores passed pre-flight validation.\n")


# ----------------------------------------------------------------------
# 10. RUN ONE COMBINED rPLUM MODEL
# ----------------------------------------------------------------------

run_combined_model <- function(config_row, seed_offset = 0) {
  
  core <- as.character(config_row$Core)
  files <- prepared_files[[core]]
  
  cat("\n============================================================\n")
  cat("RUNNING COMBINED MODEL:", core, "\n")
  cat("============================================================\n")
  
  cat("Site:", config_row$Site, "\n")
  cat(
    "Confirmed Pb measurement date:",
    config_row$Confirmed_rPlum_date_sample_AD,
    "AD\n"
  )
  cat(
    "Activity unit:",
    config_row$Confirmed_activity_unit,
    "\n"
  )
  cat(
    "AMS 14C:",
    config_row$C14_age_BP,
    "+/-",
    config_row$C14_error_1sigma,
    "14C yr BP at",
    config_row$C14_depth_cm,
    "cm\n"
  )
  cat("14C calibration curve: IntCal20 (cc = 1)\n")
  cat("Reservoir correction: 0 yr\n")
  cat("Model base:", config_row$Core_base_cm, "cm\n")
  
  result <- Plum(
    core = core,
    otherdates = paste0(core, "_C14.csv"),
    coredir = RUN_DIR,
    
    # Pb-210 metadata
    date.sample = as.numeric(
      config_row$Confirmed_rPlum_date_sample_AD
    ),
    n.supp = as.integer(config_row$n_supp),
    ra.case = as.integer(config_row$ra_case),
    Bqkg = files$Bqkg,
    
    # Depth domain
    d.min = 0,
    d.max = as.numeric(config_row$Core_base_cm),
    d.by = 1,
    thick = as.numeric(config_row$Model_thick_cm),
    
    # Accumulation-rate / memory priors
    acc.mean = as.numeric(config_row$acc_mean_yr_per_cm),
    acc.shape = as.numeric(config_row$acc_shape),
    mem.mean = as.numeric(config_row$mem_mean),
    mem.strength = as.numeric(config_row$mem_strength),
    
    # Radiocarbon calibration
    cc1 = "IntCal20",
    delta.R = 0,
    delta.STD = 0,
    
    # MCMC
    seed = SEED_BASE + seed_offset,
    prob = PROB,
    ssize = SSIZE,
    burnin = BURNIN,
    
    # Plot/output
    BCAD = TRUE,
    plot.pdf = TRUE,
    save.info = TRUE,
    save.elbowages = TRUE,
    
    # Because the 14C date lies below the recent 210Pb horizon.
    remove.tail = FALSE,
    
    # Non-interactive/reproducible execution
    ask = FALSE,
    suggest = FALSE,
    verbose = TRUE
  )
  
  saveRDS(
    result,
    file.path(
      files$core_dir,
      paste0(core, "_Plum_result_object.rds")
    )
  )
  
  invisible(result)
}


# ----------------------------------------------------------------------
# 11. RUN BASELINE MODELS
# ----------------------------------------------------------------------
#
# Baseline priors come ONLY from the CONFIG sheet.

model_results <- vector("list", nrow(cfg))
names(model_results) <- cfg$Core

if (RUN_BASELINE) {
  
  for (i in seq_len(nrow(cfg))) {
    
    core <- as.character(cfg$Core[i])
    
    model_results[[core]] <- run_combined_model(
      cfg[i, ],
      seed_offset = i
    )
  }
  
} else {
  
  cat("
Baseline runs skipped because RUN_BASELINE = FALSE.
")
}


# ----------------------------------------------------------------------
# 12. OPTIONAL SENSITIVITY ANALYSIS
# ----------------------------------------------------------------------
#
# Sensitivity scenarios are read DIRECTLY from the Excel sheet:
#   SENSITIVITY
#
# CONFIG remains the baseline model.
#
# The SENSITIVITY sheet contains one-at-a-time prior perturbations:
#   - acc.mean lower / higher
#   - mem.mean lower / higher
#
# The rows labelled "baseline" in SENSITIVITY are retained for
# documentation only and are NOT re-run here, because the baseline model
# already comes from CONFIG.
#
# This is a robustness analysis. Do NOT select a sensitivity run simply
# because it gives a preferred chronology.

if (RUN_SENSITIVITY) {
  
  # ------------------------------------------------------------
  # 12.1 READ SENSITIVITY SHEET
  # ------------------------------------------------------------
  
  if (!"SENSITIVITY" %in% excel_sheets(INPUT_XLSX)) {
    stop(
      "RUN_SENSITIVITY = TRUE, but the workbook has no sheet called ",
      "'SENSITIVITY'."
    )
  }
  
  sens <- read_excel(INPUT_XLSX, sheet = "SENSITIVITY")
  names(sens) <- trimws(names(sens))
  
  required_sens <- c(
    "Core",
    "run_id",
    "parameter_varied",
    "acc_mean_yr_per_cm",
    "acc_shape",
    "mem_mean",
    "mem_strength",
    "Use_in_sensitivity",
    "Notes"
  )
  
  missing_sens <- setdiff(required_sens, names(sens))
  
  if (length(missing_sens) > 0) {
    stop(
      paste0(
        "SENSITIVITY sheet is missing column(s): ",
        paste(missing_sens, collapse = ", ")
      )
    )
  }
  
  # Convert inclusion flag robustly.
  sens$Use_in_sensitivity <-
    as_logical_include(sens$Use_in_sensitivity)
  
  # Explicit numeric conversion.
  sens_numeric <- c(
    "acc_mean_yr_per_cm",
    "acc_shape",
    "mem_mean",
    "mem_strength"
  )
  
  for (nm in sens_numeric) {
    sens[[nm]] <- suppressWarnings(as.numeric(sens[[nm]]))
  }
  
  # ------------------------------------------------------------
  # 12.2 VALIDATE SENSITIVITY SETTINGS
  # ------------------------------------------------------------
  
  if (anyNA(sens$Core) || anyNA(sens$run_id)) {
    stop("Missing Core or run_id in SENSITIVITY.")
  }
  
  unknown_cores <- setdiff(unique(sens$Core), cfg$Core)
  
  if (length(unknown_cores) > 0) {
    stop(
      "SENSITIVITY contains unknown core(s): ",
      paste(unknown_cores, collapse = ", ")
    )
  }
  
  if (anyDuplicated(paste(sens$Core, sens$run_id, sep = "__"))) {
    stop(
      "Duplicated Core + run_id combinations found in SENSITIVITY."
    )
  }
  
  if (anyNA(sens$acc_mean_yr_per_cm) ||
      any(sens$acc_mean_yr_per_cm <= 0)) {
    stop("All sensitivity acc_mean_yr_per_cm values must be > 0.")
  }
  
  if (anyNA(sens$acc_shape) ||
      any(sens$acc_shape <= 0)) {
    stop("All sensitivity acc_shape values must be > 0.")
  }
  
  if (anyNA(sens$mem_mean) ||
      any(sens$mem_mean <= 0) ||
      any(sens$mem_mean >= 1)) {
    stop("All sensitivity mem_mean values must lie strictly between 0 and 1.")
  }
  
  if (anyNA(sens$mem_strength) ||
      any(sens$mem_strength <= 0)) {
    stop("All sensitivity mem_strength values must be > 0.")
  }
  
  # Keep only requested sensitivity scenarios.
  #
  # Baseline rows are deliberately excluded here because the true
  # baseline comes from CONFIG.
  sens_to_run <- sens[
    sens$Use_in_sensitivity &
      tolower(trimws(sens$parameter_varied)) != "baseline",
    ,
    drop = FALSE
  ]
  
  if (nrow(sens_to_run) == 0) {
    stop(
      "RUN_SENSITIVITY = TRUE but no non-baseline sensitivity rows ",
      "have Use_in_sensitivity = TRUE."
    )
  }
  
  cat("
============================================================
")
  cat("SENSITIVITY ANALYSIS
")
  cat("============================================================
")
  cat("Runs requested:", nrow(sens_to_run), "
")
  
  print(
    sens_to_run[
      ,
      c(
        "Core",
        "run_id",
        "parameter_varied",
        "acc_mean_yr_per_cm",
        "acc_shape",
        "mem_mean",
        "mem_strength"
      )
    ],
    row.names = FALSE
  )
  
  
  # ------------------------------------------------------------
  # 12.3 RUN EXACTLY THE SCENARIOS LISTED IN EXCEL
  # ------------------------------------------------------------
  
  sensitivity_manifest <- data.frame()
  
  for (j in seq_len(nrow(sens_to_run))) {
    
    core <- as.character(sens_to_run$Core[j])
    run_id <- as.character(sens_to_run$run_id[j])
    
    cfg_index <- match(core, cfg$Core)
    
    if (is.na(cfg_index)) {
      stop("Internal error: cannot match sensitivity core ", core)
    }
    
    files <- prepared_files[[core]]
    
    # Prefix with underscore so rplum output names remain readable,
    # e.g. DEE_acc_low_74.pdf rather than DEEacc_low_74.pdf.
    plum_runname <- paste0("_", run_id)
    
    cat("
------------------------------------------------------------
")
    cat("Sensitivity run:", core, "/", run_id, "
")
    cat("Varied parameter:", sens_to_run$parameter_varied[j], "
")
    cat(
      "acc.mean =", sens_to_run$acc_mean_yr_per_cm[j], "yr/cm | ",
      "acc.shape =", sens_to_run$acc_shape[j], " | ",
      "mem.mean =", sens_to_run$mem_mean[j], " | ",
      "mem.strength =", sens_to_run$mem_strength[j], "
"
    )
    cat("------------------------------------------------------------
")
    
    sens_result <- Plum(
      core = core,
      otherdates = paste0(core, "_C14.csv"),
      coredir = RUN_DIR,
      
      # Same raw chronological data and metadata as baseline
      date.sample = as.numeric(
        cfg$Confirmed_rPlum_date_sample_AD[cfg_index]
      ),
      n.supp = as.integer(cfg$n_supp[cfg_index]),
      ra.case = as.integer(cfg$ra_case[cfg_index]),
      Bqkg = files$Bqkg,
      
      # Same model depth domain as baseline
      d.min = 0,
      d.max = as.numeric(cfg$Core_base_cm[cfg_index]),
      d.by = 1,
      thick = as.numeric(cfg$Model_thick_cm[cfg_index]),
      
      # Priors read from SENSITIVITY sheet
      acc.mean = as.numeric(
        sens_to_run$acc_mean_yr_per_cm[j]
      ),
      acc.shape = as.numeric(
        sens_to_run$acc_shape[j]
      ),
      mem.mean = as.numeric(
        sens_to_run$mem_mean[j]
      ),
      mem.strength = as.numeric(
        sens_to_run$mem_strength[j]
      ),
      
      # Same 14C calibration as baseline
      cc1 = "IntCal20",
      delta.R = 0,
      delta.STD = 0,
      
      # Same MCMC settings as baseline
      seed = SEED_BASE + 10000 + j,
      prob = PROB,
      ssize = SSIZE,
      burnin = BURNIN,
      
      # Save separate output for every scenario
      BCAD = TRUE,
      plot.pdf = TRUE,
      save.info = TRUE,
      save.elbowages = TRUE,
      remove.tail = FALSE,
      
      runname = plum_runname,
      
      ask = FALSE,
      suggest = FALSE,
      verbose = TRUE
    )
    
    # Save the returned R object as an additional reproducibility record.
    saveRDS(
      sens_result,
      file.path(
        files$core_dir,
        paste0(
          core,
          plum_runname,
          "_Plum_result_object.rds"
        )
      )
    )
    
    sensitivity_manifest <- rbind(
      sensitivity_manifest,
      data.frame(
        Core = core,
        run_id = run_id,
        parameter_varied =
          as.character(sens_to_run$parameter_varied[j]),
        acc_mean_yr_per_cm =
          as.numeric(sens_to_run$acc_mean_yr_per_cm[j]),
        acc_shape =
          as.numeric(sens_to_run$acc_shape[j]),
        mem_mean =
          as.numeric(sens_to_run$mem_mean[j]),
        mem_strength =
          as.numeric(sens_to_run$mem_strength[j]),
        seed = SEED_BASE + 10000 + j,
        Notes = as.character(sens_to_run$Notes[j]),
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Exact record of all sensitivity runs that were actually executed.
  write.csv(
    sensitivity_manifest,
    file = file.path(
      RUN_DIR,
      "SENSITIVITY_runs_executed.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  cat("
Sensitivity runs completed.
")
  cat(
    "Manifest saved as: ",
    file.path(RUN_DIR, "SENSITIVITY_runs_executed.csv"),
    "
",
sep = ""
  )
  
} else {
  
  cat("
Sensitivity analysis skipped because RUN_SENSITIVITY = FALSE.
")
}


# ----------------------------------------------------------------------
# 13. REPRODUCIBILITY RECORD
# ----------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    RUN_DIR,
    "sessionInfo_combined_rplum_models.txt"
  )
)

capture.output(
  citation("rplum"),
  file = file.path(RUN_DIR, "citation_rplum.txt")
)

capture.output(
  citation("rbacon"),
  file = file.path(RUN_DIR, "citation_rbacon.txt")
)


# ----------------------------------------------------------------------
# 14. FINAL MESSAGE
# ----------------------------------------------------------------------

cat("\n============================================================\n")
cat("REQUESTED rPLUM RUNS COMPLETED\n")
cat("============================================================\n")

cat("Output directory:\n", RUN_DIR, "\n")
cat("RUN_BASELINE =", RUN_BASELINE, "\n")
cat("RUN_SENSITIVITY =", RUN_SENSITIVITY, "\n")

cat("\nFor every core inspect:\n")
cat("  1. MCMC mixing / stability\n")
cat("  2. prior versus posterior accumulation rate\n")
cat("  3. prior versus posterior memory\n")
cat("  4. measured versus modelled 210Pb\n")
cat("  5. 14C constraint versus posterior chronology\n")
cat("  6. 95% age-depth uncertainty envelope\n")
cat("  7. sensitivity to priors and model thickness\n")

cat("\nIMPORTANT:\n")
cat(
  "The CF:CS-derived unsupported Pb values and extrapolated basal ages\n",
  "are NOT used as chronological constraints in this combined model.\n",
  "The model uses the raw total Pb/Po proxy activities plus the raw AMS\n",
  "14C determination.\n",
  sep = ""
)

cat("============================================================\n")


# ======================================================================
# END OF SCRIPT
# ======================================================================