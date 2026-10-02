#!/usr/bin/env Rscript
# =============================================================================
# test_results.R
#
# Validates the output CSVs produced by the XevaDB seeding scripts for every
# subfolder inside a given results directory.
#
# Usage:
#   Rscript test_results.R [results_dir]
#
# If no argument is supplied, the script defaults to:
#   ../results_2026_sep_23   (relative to the script location)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

results_log <- list()   # collects all test outcomes

log_test <- function(subfolder, file, test_name, passed, message = "") {
  cat(sprintf("  [%s] %s | %s%s\n",
              if (isTRUE(passed)) "PASS" else "FAIL",
              file, test_name,
              if (nzchar(message)) paste0(" -- ", message) else ""))
  results_log[[length(results_log) + 1]] <<- list(
    subfolder  = subfolder,
    file       = file,
    test       = test_name,
    passed     = isTRUE(passed),
    message    = message
  )
}

log_warn <- function(subfolder, file, test_name, message = "") {
  cat(sprintf("  [WARN] %s | %s%s\n",
              file, test_name,
              if (nzchar(message)) paste0(" -- ", message) else ""))
  results_log[[length(results_log) + 1]] <<- list(
    subfolder  = subfolder,
    file       = file,
    test       = test_name,
    passed     = NA,   # NA = warning, not a hard failure
    message    = message
  )
}

# Safe CSV reader -- returns NULL on error.
# R's write.csv always writes a leading unnamed row-index column.
# We detect and drop it so downstream tests can reference columns by name.
safe_read <- function(path) {
  dt <- tryCatch(fread(path, showProgress = FALSE, header = TRUE), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  # Drop leading row-index column: named "" or "V1", containing integers 1..nrow
  first_col <- names(dt)[1]
  if (first_col %in% c("", "V1")) {
    vals <- dt[[1]]
    if (is.numeric(vals) && !anyNA(vals) && all(vals == seq_len(nrow(dt)))) {
      dt[[1]] <- NULL
    }
  }
  dt
}

# ---------------------------------------------------------------------------
# Per-file test definitions
# ---------------------------------------------------------------------------

test_model_response <- function(dt, subfolder) {
  f <- "model_response.csv"
  required_cols <- c("drug", "model.id", "response_type", "value")
  log_test(subfolder, f, "required columns present",
           all(required_cols %in% names(dt)),
           paste(setdiff(required_cols, names(dt)), collapse = ", "))

  log_test(subfolder, f, "no NA model.id",
           !any(is.na(dt$model.id)))

  log_test(subfolder, f, "no NA drug",
           !any(is.na(dt$drug)))

  valid_rtypes <- c("mRECIST", "best.average.response", "slope", "AUC", "survival")
  bad_rtypes   <- setdiff(unique(dt$response_type), valid_rtypes)
  log_test(subfolder, f, "only known response_type values",
           length(bad_rtypes) == 0,
           if (length(bad_rtypes)) paste(bad_rtypes, collapse = ", ") else "")

  mrecist_vals <- dt[response_type == "mRECIST", value]
  # Treat blank strings the same as NA (written that way by write.csv for R NA)
  mrecist_vals[mrecist_vals == ""] <- NA
  valid_mr     <- c("CR", "PR", "SD", "PD")
  bad_mr       <- setdiff(unique(na.omit(mrecist_vals)), valid_mr)
  log_test(subfolder, f, "mRECIST values are CR/PR/SD/PD (or NA)",
           length(bad_mr) == 0,
           if (length(bad_mr)) paste(bad_mr, collapse = ", ") else "")

  n_na_mrecist <- sum(is.na(mrecist_vals))
  if (n_na_mrecist > 0)
    log_warn(subfolder, f, "NA mRECIST values",
             sprintf("%d model(s) could not be classified", n_na_mrecist))

  numeric_types <- c("best.average.response", "slope", "AUC", "survival")
  for (rt in numeric_types) {
    raw  <- dt[response_type == rt, value]
    # write.csv writes R's NA as the string "NA"; treat those as missing
    vals <- suppressWarnings(as.numeric(ifelse(raw == "NA", NA, raw)))
    if (length(vals) > 0) {
      real_bad <- !is.na(raw) & raw != "NA" & is.na(vals)
      log_test(subfolder, f, paste0(rt, " values are numeric"),
               !any(real_bad),
               if (any(real_bad)) "non-numeric, non-NA entries found" else "")
    }
  }
}

test_model_information <- function(dt, subfolder) {
  f <- "model_information.csv"
  required_cols <- c("model.id", "tissue", "patient.id", "drug", "dataset")
  log_test(subfolder, f, "required columns present",
           all(required_cols %in% names(dt)),
           paste(setdiff(required_cols, names(dt)), collapse = ", "))

  log_test(subfolder, f, "no NA model.id",     !any(is.na(dt$model.id)))
  log_test(subfolder, f, "no NA patient.id",   !any(is.na(dt$patient.id)))
  log_test(subfolder, f, "no NA tissue",       !any(is.na(dt$tissue)))
  log_test(subfolder, f, "no duplicate model.id",
           !any(duplicated(dt$model.id)),
           sprintf("%d duplicate(s)", sum(duplicated(dt$model.id))))
}

test_batch_information <- function(dt, subfolder) {
  f <- "batch_information.csv"
  required_cols <- c("batch.id", "model.id", "type")
  log_test(subfolder, f, "required columns present",
           all(required_cols %in% names(dt)),
           paste(setdiff(required_cols, names(dt)), collapse = ", "))

  log_test(subfolder, f, "no NA batch.id",  !any(is.na(dt$batch.id)))
  log_test(subfolder, f, "no NA model.id",  !any(is.na(dt$model.id)))

  valid_types <- c("control", "treatment")
  bad_types   <- setdiff(unique(dt$type), valid_types)
  log_test(subfolder, f, "type is 'control' or 'treatment'",
           length(bad_types) == 0,
           if (length(bad_types)) paste(bad_types, collapse = ", ") else "")

  # Every batch should have at least one control
  batches_no_ctrl <- dt[, .(has_ctrl = any(type == "control")), by = batch.id][has_ctrl == FALSE, batch.id]
  log_test(subfolder, f, "each batch has at least one control",
           length(batches_no_ctrl) == 0,
           sprintf("%d batch(es) lack a control row", length(batches_no_ctrl)))
}

test_batch_response <- function(dt, subfolder) {
  f <- "batch_response.csv"
  required_cols <- c("batch.id", "response_type", "value")
  log_test(subfolder, f, "required columns present",
           all(required_cols %in% names(dt)),
           paste(setdiff(required_cols, names(dt)), collapse = ", "))

  log_test(subfolder, f, "no NA batch.id", !any(is.na(dt$batch.id)))

  raw  <- dt$value
  vals <- suppressWarnings(as.numeric(ifelse(raw %in% c("NA", "NaN"), NA, raw)))
  real_bad <- !is.na(raw) & !raw %in% c("NA", "NaN") & is.na(vals)
  log_test(subfolder, f, "value column is numeric",
           !any(real_bad),
           if (any(real_bad)) "non-numeric, non-NA entries found" else "")

  n_na <- sum(is.na(vals))
  if (n_na > 0)
    log_warn(subfolder, f, "NA/NaN batch response values",
             sprintf("%d entries could not be computed (NA/NaN)", n_na))
}

test_drug_screening <- function(dt, subfolder) {
  f <- "drug_screening.csv"
  required_cols <- c("model.id", "drug", "time", "volume", "volume.normal")
  log_test(subfolder, f, "required columns present",
           all(required_cols %in% names(dt)),
           paste(setdiff(required_cols, names(dt)), collapse = ", "))

  log_test(subfolder, f, "no NA model.id", !any(is.na(dt$model.id)))
  log_test(subfolder, f, "no NA time",     !any(is.na(dt$time)))
  log_test(subfolder, f, "no NA volume",   !any(is.na(dt$volume)))

  log_test(subfolder, f, "time values are non-negative",
           all(dt$time >= 0, na.rm = TRUE))

  log_test(subfolder, f, "volume values are non-negative",
           all(dt$volume >= 0, na.rm = TRUE))

  # Each model should have more than one time-point (otherwise no growth curve)
  n_points <- dt[, .N, by = model.id]
  single_pt <- n_points[N == 1, model.id]
  if (length(single_pt) > 0)
    log_warn(subfolder, f, "models with only one time-point",
             sprintf("%d model(s) -- cannot fit growth curve", length(single_pt)))
}

test_modelid_mapping <- function(dt, subfolder) {
  f <- "modelid_moleculardata_mapping.csv"
  required_cols <- c("model.id", "biobase.id", "mDataType")
  log_test(subfolder, f, "required columns present",
           all(required_cols %in% names(dt)),
           paste(setdiff(required_cols, names(dt)), collapse = ", "))

  log_test(subfolder, f, "no NA model.id",   !any(is.na(dt$model.id)))
  log_test(subfolder, f, "no NA biobase.id", !any(is.na(dt$biobase.id)))
  log_test(subfolder, f, "no NA mDataType",  !any(is.na(dt$mDataType)))

  valid_types <- c("mutation", "RNASeq", "CNV", "cnv", "rnaseq", "RNAseq")
  bad_types   <- setdiff(unique(dt$mDataType), valid_types)
  if (length(bad_types) > 0)
    log_warn(subfolder, f, "unexpected mDataType values",
             paste(bad_types, collapse = ", "))
}

test_mutation <- function(dt, subfolder) {
  f <- "mutation.csv"
  required_cols <- c("gene.id", "sequencing.uid", "value")
  log_test(subfolder, f, "required columns present",
           all(required_cols %in% names(dt)),
           paste(setdiff(required_cols, names(dt)), collapse = ", "))

  log_test(subfolder, f, "file contains data rows",
           nrow(dt) > 0,
           if (nrow(dt) == 0) "0 data rows (empty file)" else "")

  log_test(subfolder, f, "no NA gene.id",        !any(is.na(dt$gene.id)))
  log_test(subfolder, f, "no NA sequencing.uid", !any(is.na(dt$sequencing.uid)))

  # Mutation values may be:
  #   - 0/1 binary flags (PDXE-style)
  #   - Consequence labels e.g. Missense_Mutation, Silent, Frame_Shift_Ins (Tsao-style)
  # Just ensure no completely empty / all-NA column.
  all_empty <- all(is.na(dt$value) | dt$value == "")
  log_test(subfolder, f, "value column is not entirely empty",
           !all_empty)

  # Report the style detected
  uniq_vals <- unique(na.omit(dt$value))
  is_binary <- all(uniq_vals %in% c("0", "1", 0L, 1L, 0, 1))
  if (is_binary) {
    pct_nonzero <- 100 * mean(dt$value != 0 & !is.na(dt$value))
    log_warn(subfolder, f, "mutation style: binary 0/1",
             sprintf("%.2f%% of entries are non-zero", pct_nonzero))
  } else {
    n_nonzero <- sum(!is.na(dt$value) & dt$value != "0" & dt$value != 0)
    log_warn(subfolder, f, "mutation style: consequence labels",
             sprintf("%d non-zero/non-silent entries found", n_nonzero))
  }
}

test_cnv <- function(dt, subfolder) {
  f <- "copy_number_variation.csv"
  required_cols <- c("gene.id", "sequencing.uid", "value")
  log_test(subfolder, f, "required columns present",
           all(required_cols %in% names(dt)),
           paste(setdiff(required_cols, names(dt)), collapse = ", "))

  log_test(subfolder, f, "file contains data rows",
           nrow(dt) > 0,
           if (nrow(dt) == 0) "0 data rows (empty file)" else "")

  log_test(subfolder, f, "no NA gene.id",        !any(is.na(dt$gene.id)))
  log_test(subfolder, f, "no NA sequencing.uid", !any(is.na(dt$sequencing.uid)))

  # Expected: discrete CNV category strings (e.g. Deletion, Shallow Deletion, 0, Gain, Amplification) or numeric calls
  valid_cnv <- c("Deletion", "Deep Deletion", "deletion", "Shallow Deletion", "0", "Gain", "Amplification", "-2", "-1", "1", "2")
  vals <- as.character(dt$value)
  bad_vals <- setdiff(unique(na.omit(vals)), valid_cnv)
  log_test(subfolder, f, "valid CNV categories",
           length(bad_vals) == 0,
           if (length(bad_vals)) paste(head(bad_vals, 5), collapse = ", ") else "")
}

test_rna_sequencing <- function(dt, subfolder) {
  f <- "rna_sequencing.csv"
  required_cols <- c("gene.id", "sequencing.uid", "value")
  log_test(subfolder, f, "required columns present",
           all(required_cols %in% names(dt)),
           paste(setdiff(required_cols, names(dt)), collapse = ", "))

  log_test(subfolder, f, "file contains data rows",
           nrow(dt) > 0,
           if (nrow(dt) == 0) "0 data rows (empty file)" else "")

  log_test(subfolder, f, "no NA gene.id",        !any(is.na(dt$gene.id)))
  log_test(subfolder, f, "no NA sequencing.uid", !any(is.na(dt$sequencing.uid)))

  raw  <- dt$value
  vals <- suppressWarnings(as.numeric(ifelse(raw %in% c("NA", "NaN"), NA, raw)))
  real_bad <- !is.na(raw) & !raw %in% c("NA", "NaN") & is.na(vals)
  log_test(subfolder, f, "value column is numeric",
           !any(real_bad),
           if (any(real_bad)) "non-numeric, non-NA entries found" else "")

  n_na <- sum(is.na(vals))
  if (n_na > 0)
    log_warn(subfolder, f, "NA/NaN expression values",
             sprintf("%d entries are NA or NaN", n_na))
}

# ---------------------------------------------------------------------------
# Cross-file consistency checks
# ---------------------------------------------------------------------------

test_cross_file <- function(dts, subfolder) {
  f <- "cross-file"

  # model.id consistency: every model in model_response should be in model_information
  if (!is.null(dts$model_response) && !is.null(dts$model_information)) {
    mr_ids  <- unique(dts$model_response$model.id)
    mi_ids  <- unique(dts$model_information$model.id)
    missing <- setdiff(mr_ids, mi_ids)
    log_test(subfolder, f,
             "all model_response model.ids in model_information",
             length(missing) == 0,
             sprintf("%d id(s) missing: %s", length(missing),
                     paste(head(missing, 5), collapse = ", ")))
    extra <- setdiff(mi_ids, mr_ids)
    if (length(extra) > 0)
      log_warn(subfolder, f,
               "model_information has model.ids not in model_response",
               sprintf("%d extra id(s)", length(extra)))
  }

  # batch_information model.ids should all appear in drug_screening
  if (!is.null(dts$batch_information) && !is.null(dts$drug_screening)) {
    bi_ids  <- unique(dts$batch_information$model.id)
    ds_ids  <- unique(dts$drug_screening$model.id)
    missing <- setdiff(bi_ids, ds_ids)
    log_test(subfolder, f,
             "all batch_information model.ids in drug_screening",
             length(missing) == 0,
             sprintf("%d id(s) missing", length(missing)))
  }

  # batch_response batch.ids should match batch_information batch.ids
  if (!is.null(dts$batch_response) && !is.null(dts$batch_information)) {
    br_batches <- unique(dts$batch_response$batch.id)
    bi_batches <- unique(dts$batch_information$batch.id)
    missing    <- setdiff(br_batches, bi_batches)
    log_test(subfolder, f,
             "all batch_response batch.ids in batch_information",
             length(missing) == 0,
             sprintf("%d id(s) missing", length(missing)))
  }

  # modelid_moleculardata_mapping: biobase.ids should appear in molecular tables
  if (!is.null(dts$modelid_mapping)) {
    if (!is.null(dts$mutation)) {
      map_mut  <- unique(dts$modelid_mapping[mDataType == "mutation", biobase.id])
      mut_uids <- unique(dts$mutation$sequencing.uid)
      missing  <- setdiff(map_mut, mut_uids)
      log_test(subfolder, f,
               "mutation mapping biobase.ids in mutation.csv sequencing.uid",
               length(missing) == 0,
               sprintf("%d id(s) missing", length(missing)))
    }
    if (!is.null(dts$rna_sequencing)) {
      map_rna  <- unique(dts$modelid_mapping[tolower(mDataType) == "rnaseq", biobase.id])
      rna_uids <- unique(dts$rna_sequencing$sequencing.uid)
      missing  <- setdiff(map_rna, rna_uids)
      log_test(subfolder, f,
               "RNASeq mapping biobase.ids in rna_sequencing.csv sequencing.uid",
               length(missing) == 0,
               sprintf("%d id(s) missing", length(missing)))
    }
    if (!is.null(dts$cnv)) {
      map_cnv  <- unique(dts$modelid_mapping[tolower(mDataType) == "cnv", biobase.id])
      cnv_uids <- unique(dts$cnv$sequencing.uid)
      missing  <- setdiff(map_cnv, cnv_uids)
      log_test(subfolder, f,
               "CNV mapping biobase.ids in copy_number_variation.csv sequencing.uid",
               length(missing) == 0,
               sprintf("%d id(s) missing", length(missing)))
    }
  }
}

# ---------------------------------------------------------------------------
# Main runner -- iterate over subfolders
# ---------------------------------------------------------------------------

EXPECTED_FILES <- c(
  "model_response.csv",
  "model_information.csv",
  "batch_information.csv",
  "batch_response.csv",
  "drug_screening.csv",
  "modelid_moleculardata_mapping.csv",
  "mutation.csv",
  "copy_number_variation.csv",
  "rna_sequencing.csv"
)

run_tests_for_subfolder <- function(subfolder_path) {
  subfolder <- basename(subfolder_path)
  cat(sprintf("\n%s\n## Testing: %s\n%s\n",
              strrep("=", 60), subfolder, strrep("=", 60)))

  # 1. Check which files are present
  present <- list.files(subfolder_path, pattern = "\\.csv$")
  missing_files <- setdiff(EXPECTED_FILES, present)
  if (length(missing_files) > 0) {
    cat(sprintf("  [WARN] Missing expected files: %s\n",
                paste(missing_files, collapse = ", ")))
    for (mf in missing_files) {
      results_log[[length(results_log) + 1]] <<- list(
        subfolder = subfolder, file = mf,
        test = "file exists", passed = FALSE,
        message = "file not found in output folder"
      )
    }
  }

  # 2. Load each file that is present
  load_if_present <- function(fname) {
    path <- file.path(subfolder_path, fname)
    if (!file.exists(path)) return(NULL)
    cat(sprintf("  Loading %s ...\n", fname))
    dt <- safe_read(path)
    if (is.null(dt)) {
      cat(sprintf("  [FAIL] %s -- could not be read\n", fname))
      results_log[[length(results_log) + 1]] <<- list(
        subfolder = subfolder, file = fname,
        test = "file readable", passed = FALSE,
        message = "fread() returned an error"
      )
    }
    dt
  }

  dts <- list(
    model_response    = load_if_present("model_response.csv"),
    model_information = load_if_present("model_information.csv"),
    batch_information = load_if_present("batch_information.csv"),
    batch_response    = load_if_present("batch_response.csv"),
    drug_screening    = load_if_present("drug_screening.csv"),
    modelid_mapping   = load_if_present("modelid_moleculardata_mapping.csv"),
    mutation          = load_if_present("mutation.csv"),
    cnv               = load_if_present("copy_number_variation.csv"),
    rna_sequencing    = load_if_present("rna_sequencing.csv")
  )

  # 3. Per-file tests
  cat("\n--- Per-file tests ---\n")
  if (!is.null(dts$model_response))    test_model_response(dts$model_response, subfolder)
  if (!is.null(dts$model_information)) test_model_information(dts$model_information, subfolder)
  if (!is.null(dts$batch_information)) test_batch_information(dts$batch_information, subfolder)
  if (!is.null(dts$batch_response))    test_batch_response(dts$batch_response, subfolder)
  if (!is.null(dts$drug_screening))    test_drug_screening(dts$drug_screening, subfolder)
  if (!is.null(dts$modelid_mapping))   test_modelid_mapping(dts$modelid_mapping, subfolder)
  if (!is.null(dts$mutation))          test_mutation(dts$mutation, subfolder)
  if (!is.null(dts$cnv))               test_cnv(dts$cnv, subfolder)
  if (!is.null(dts$rna_sequencing))    test_rna_sequencing(dts$rna_sequencing, subfolder)

  # 4. Cross-file consistency
  cat("\n--- Cross-file consistency tests ---\n")
  test_cross_file(dts, subfolder)
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else {
  # Default: sibling results_2026_sep_23 relative to this script
  script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                         error = function(e) getwd())
  file.path(script_dir, "..", "results_2026_sep_23")
}

results_dir <- normalizePath(results_dir, mustWork = FALSE)

if (!dir.exists(results_dir)) {
  stop("Results directory not found: ", results_dir)
}

cat(sprintf("Results directory: %s\n", results_dir))

subfolders <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)
if (length(subfolders) == 0) {
  stop("No subfolders found in: ", results_dir)
}

for (sf in subfolders) {
  run_tests_for_subfolder(sf)
}

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------

cat(sprintf("\n%s\n## SUMMARY\n%s\n", strrep("=", 60), strrep("=", 60)))

log_dt <- rbindlist(lapply(results_log, as.data.table), fill = TRUE)

failures <- log_dt[passed == FALSE]
warnings <- log_dt[is.na(passed)]
passes   <- log_dt[passed == TRUE]

cat(sprintf("  Total tests : %d\n", nrow(log_dt)))
cat(sprintf("  PASS        : %d\n", nrow(passes)))
cat(sprintf("  WARN        : %d\n", nrow(warnings)))
cat(sprintf("  FAIL        : %d\n", nrow(failures)))

if (nrow(failures) > 0) {
  cat("\nFailed tests:\n")
  for (i in seq_len(nrow(failures))) {
    r <- failures[i]
    cat(sprintf("  [%s] %s > %s > %s\n",
                r$subfolder, r$file, r$test,
                ifelse(nzchar(r$message), r$message, "")))
  }
}

# Exit with non-zero status if any hard failures (CI-friendly)
if (nrow(failures) > 0) quit(status = 1L)
