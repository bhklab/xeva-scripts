#!/usr/bin/env Rscript
# Export *all* models to model_response.csv, guaranteeing one row per model.id
# even when Xeva cannot compute response metrics. Where metrics cannot be
# computed, they are left as NA. The output is intentionally unconstrained:
# - does not subset away experiments
# - does not trim by max.time (uses max.time=NULL)
#
# Usage (CLI):
#   Rscript export_model_response_all.R \
#     --rds /path/to/your.xeva.rds \
#     --outdir ./results
#
# If args are omitted, defaults below are used.

suppressMessages({
    library(optparse)
})

# --------- CLI args ---------
opt <- list(
    rds = "/Users/mattbocc/uhn/xeva-scripts/xevasets_obj_2026/UHN_Tsao_Lung_DrugResponse_2022_v1.rds",
    outdir = "../results_2026_sep_23"
)

if (interactive() == FALSE) {
    option_list <- list(
        make_option(c("--rds"),
            type = "character", default = opt$rds,
            help = "Path to Xeva RDS file", metavar = "file"
        ),
        make_option(c("--outdir"),
            type = "character", default = opt$outdir,
            help = "Output directory", metavar = "dir"
        )
    )
    parser <- OptionParser(option_list = option_list)
    args <- parse_args(parser)
    opt$rds <- args$rds
    opt$outdir <- args$outdir
}

# If user accidentally passed a CSV path as outdir, use its directory
if (grepl("\\.csv$", opt$outdir, ignore.case = TRUE)) {
    message("NOTE: --outdir looks like a file path; using its directory instead.")
    opt$outdir <- dirname(opt$outdir)
}

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# --------- Helpers ---------
flatten_val <- function(v) {
    # Turn Xeva's @drug which may be list(join.name=...) into a scalar
    if (is.null(v)) {
        return(NA_character_)
    }
    if (is.list(v)) {
        if (!is.null(v$join.name)) {
            return(as.character(v$join.name))
        }
        return(paste0(unlist(v, use.names = FALSE), collapse = "|"))
    }
    as.character(v)
}

ensure_drug_slots <- function(obj) {
    for (id in names(obj@experiment)) {
        ex <- obj@experiment[[id]]
        if (is.null(ex)) next
        if (!("drug" %in% slotNames(ex))) next
        val <- try(ex@drug, silent = TRUE)
        if (!inherits(val, "try-error") && is.list(val) && !is.null(val[["join.name"]]) &&
            nzchar(as.character(val[["join.name"]]))) {
            next
        }
        if (!inherits(val, "try-error") && !is.null(val) && !is.list(val) && nzchar(as.character(val))) {
            ex@drug <- list(join.name = as.character(val))
            obj@experiment[[id]] <- ex
            next
        }
        # last-resort inference
        fill <- if (grepl("control", id, ignore.case = TRUE)) "control" else NA_character_
        if (!nzchar(fill) || is.na(fill)) {
            toks <- unlist(strsplit(id, "[.-]"))
            if (length(toks)) fill <- toks[length(toks)]
        }
        if (!nzchar(fill) || is.na(fill)) fill <- "UNKNOWN"
        ex@drug <- list(join.name = fill)
        obj@experiment[[id]] <- ex
    }
    obj
}

safe_setResponse <- function(x) {
    # Use max.time=NULL to avoid trimming; swallow errors
    try(
        {
            x <- setResponse(x,
                res.measure = c("mRECIST", "slope", "AUC", "angle", "abc", "TGI"),
                max.time = NULL, verbose = FALSE
            )
        },
        silent = TRUE
    )
    x
}

safe_model_response <- function(x) {
    out <- try(model_response(x), silent = TRUE)
    if (!inherits(out, "try-error") && is.data.frame(out)) {
        return(out)
    }
    sm <- try(x@sensitivity$model, silent = TRUE)
    if (!inherits(sm, "try-error") && is.data.frame(sm)) {
        sm <- as.data.frame(sm, stringsAsFactors = FALSE, check.names = FALSE)
        if (!"model.id" %in% names(sm)) sm$model.id <- rownames(sm)
        keep <- intersect(
            c(
                "model.id", "drug", "mRECIST", "best.response", "best.response.time",
                "best.average.response", "best.average.response.time", "slope", "AUC", "angle", "abc", "TGI"
            ),
            names(sm)
        )
        return(sm[, keep, drop = FALSE])
    }
    # fall back to just model ids
    data.frame(model.id = names(x@experiment), stringsAsFactors = FALSE)
}

safe_write_csv <- function(df, path) {
    if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
    rn <- rownames(df)
    if (is.null(rn)) rn <- seq_len(nrow(df))
    rownames(df) <- rn
    utils::write.csv(df, file = path, row.names = TRUE, na = "")
}

# --------- Main ---------
suppressMessages({
    library(Xeva)
})

message("Reading: ", normalizePath(opt$rds))
x <- readRDS(opt$rds)

# Guarantee we can surface drug names even if Xeva can't compute metrics
x <- ensure_drug_slots(x)

# Attempt to compute responses without constraints (won't abort if it fails)
x <- safe_setResponse(x)

# Candidate response table (may be empty)
mres <- safe_model_response(x)

# Build the complete model list
all_models <- data.frame(model.id = names(x@experiment), stringsAsFactors = FALSE)

# Attach drug from experiments if missing/partial
ex_drug <- data.frame(
    model.id = names(x@experiment),
    drug = vapply(x@experiment, function(ex) {
        if (!is.null(ex) && "drug" %in% slotNames(ex)) flatten_val(ex@drug) else NA_character_
    }, ""),
    stringsAsFactors = FALSE
)

# Left-join to include every model exactly once
out <- merge(all_models, mres, by = "model.id", all.x = TRUE)
if (!"drug" %in% names(out) || all(is.na(out$drug))) {
    out <- merge(out, ex_drug, by = "model.id", all.x = TRUE)
} else {
    # fill missing drug only
    if ("drug" %in% names(out)) {
        idx <- which(is.na(out$drug) | out$drug == "")
        if (length(idx)) {
            map <- setNames(ex_drug$drug, ex_drug$model.id)
            out$drug[idx] <- map[out$model.id[idx]]
        }
    }
}

# Ensure canonical column order (present ones only)
col_order <- c(
    "model.id", "drug", "mRECIST", "best.response", "best.response.time",
    "best.average.response", "best.average.response.time", "slope", "AUC", "angle", "abc", "TGI"
)
ord <- intersect(col_order, names(out))
out <- out[, c(ord, setdiff(names(out), ord)), drop = FALSE]

# Write outputs (WIDE -> LONG, matching required format)
# Desired columns: drug, model.id, response_type, value
# response_types we emit (in this order): mRECIST, best.average.response, slope, AUC, survival
# where survival := best.response.time if available
make_long <- function(df) {
    pieces <- list()
    add <- function(type, col) {
        if (col %in% names(df)) {
            pieces[[length(pieces) + 1]] <<- data.frame(
                drug = df$drug,
                model.id = df$model.id,
                response_type = type,
                value = df[[col]],
                stringsAsFactors = FALSE
            )
        }
    }
    add("mRECIST", "mRECIST")
    add("best.average.response", "best.average.response")
    add("slope", "slope")
    add("AUC", "AUC")
    if ("best.response.time" %in% names(df)) {
        pieces[[length(pieces) + 1]] <- data.frame(
            drug = df$drug,
            model.id = df$model.id,
            response_type = "survival",
            value = df[["best.response.time"]],
            stringsAsFactors = FALSE
        )
    }
    if (!length(pieces)) {
        return(data.frame(drug = character(0), model.id = character(0), response_type = character(0), value = character(0)))
    }
    long <- do.call(rbind, pieces)
    # Order consistently: by drug, model.id, response_type
    long <- long[order(as.character(long$drug), as.character(long$model.id), match(long$response_type, c("mRECIST", "best.average.response", "slope", "AUC", "survival"), nomatch = 999L)), , drop = FALSE]
    rownames(long) <- seq_len(nrow(long))
    long
}

long_out <- make_long(out)

csv_path <- file.path(opt$outdir, "model_response.csv")
safe_write_csv(long_out, csv_path)

# Diagnostics: which models lack *all* of the emitted metrics (value all NA across types)
if (nrow(long_out)) {
    # pivot-like check: for each model.id, are all values NA?
    miss <- aggregate(is.na(long_out$value), by = list(model.id = long_out$model.id), FUN = all)
    miss <- miss[miss$x, , drop = FALSE]
    diag_path <- file.path(opt$outdir, "models_without_computable_response.csv")
    if (nrow(miss)) {
        # attach drug for readability (first non-NA drug per model)
        first_drug <- tapply(as.character(long_out$drug), long_out$model.id, function(v) v[which.max(!is.na(v))[1]])
        miss$drug <- unname(first_drug[miss$model.id])
        miss$x <- NULL
        safe_write_csv(miss, diag_path)
    } else {
        # write empty diag file
        safe_write_csv(data.frame(model.id = character(0), drug = character(0)), diag_path)
    }
}

message("Wrote: ", normalizePath(csv_path))
