#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressMessages({
    library(Xeva)
    library(readxl)
    library(data.table)
})

# legacy helper functions
source("xevaDB_fun.R")

# rds_path <- "/Users/mattbocc/uhn/xeva-scripts/xevasets-obj/Updated_Cescon_TNBC_Xeva_Obj_2.rds"  # Cescon TNBC
rds_path <- "/Users/mattbocc/uhn/xeva-scripts/xevasets_obj_2026/UHN_Tsao_Lung_DrugResponse_2022_v1.rds" # Tsao Lung
x <- readRDS(rds_path)

is_tnbc <- grepl("TNBC|Cescon", basename(rds_path), ignore.case = TRUE)
OUTDIR <- if (is_tnbc) "../results_2026_sep_23/TNBC_v2" else "../results_2026_sep_23/KRAS_LUNG_v2"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("Xeva version:", as.character(packageVersion("Xeva")), "\n")

# assay fetcher (list or MultiAssayExperiment)
.get_assay <- function(xeva, keys) {
    mp <- xeva@molecularProfiles
    keys <- tolower(keys)
    if (inherits(mp, "MultiAssayExperiment")) {
        if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
            stop("MultiAssayExperiment not available")
        }
        exps <- MultiAssayExperiment::experiments(mp)
        nm <- tolower(names(exps))
        hit <- which(nm %in% keys)
        if (!length(hit)) stop("Assay not found: ", paste(keys, collapse = "/"))
        exps[[hit[1]]]
    } else {
        nm <- tolower(names(mp))
        hit <- which(nm %in% keys)
        if (!length(hit)) stop("Assay not found: ", paste(keys, collapse = "/"))
        mp[[names(mp)[hit[1]]]]
    }
}

# experiments that have at least 1 data row
nonempty_ids <- function(obj) {
    if (is.null(obj@experiment) || !length(obj@experiment)) {
        return(character(0))
    }
    names(obj@experiment)[vapply(obj@experiment, function(e) {
        if (is.null(e)) {
            return(FALSE)
        }
        ok <- try(nrow(e@data), silent = TRUE)
        if (inherits(ok, "try-error")) {
            return(FALSE)
        }
        ok > 0
    }, logical(1))]
}

# make sure every experiment has a usable @drug of the *list(join.name=...)* form
ensure_drug_slots <- function(obj) {
    for (id in names(obj@experiment)) {
        ex <- obj@experiment[[id]]
        if (is.null(ex)) next
        if (!("drug" %in% slotNames(ex))) next

        val <- try(ex@drug, silent = TRUE)

        # already list(join.name=...)
        if (!inherits(val, "try-error") && is.list(val) && !is.null(val[["join.name"]]) &&
            nzchar(as.character(val[["join.name"]]))) {
            next
        }

        # scalar/string → wrap into list(join.name=...)
        if (!inherits(val, "try-error") && !is.null(val) && !is.list(val) && nzchar(as.character(val))) {
            ex@drug <- list(join.name = as.character(val))
            obj@experiment[[id]] <- ex
            next
        }

        # infer from id (control if present, else last token)
        fill <- NA_character_
        if (grepl("control", id, ignore.case = TRUE)) fill <- "control"
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

# safe subset that never aborts the pipeline
safe_subset <- function(obj, ids, keep.batch = FALSE) {
    out <- try(subsetXeva(obj, ids = ids, id.name = "model.id", keep.batch = keep.batch), silent = TRUE)
    if (inherits(out, "try-error")) {
        message("WARN: subsetXeva -> ", conditionMessage(out))
        obj
    } else {
        out
    }
}

# log which experiments are NULL
log_null_experiments <- function(obj, tag = "initial") {
    null_ids <- names(obj@experiment)[vapply(obj@experiment, is.null, logical(1))]
    cat("NULL experiments:", length(null_ids), "\n")
    if (length(null_ids)) {
        write.csv(data.frame(model.id = null_ids),
            file = file.path(OUTDIR, sprintf("%s_null_experiments.csv", tag)),
            row.names = FALSE
        )
    }
}

# batches that actually have control+treatment with enough time/volume points
# NOTE: This is the *only* place we require ≥2 time/volume points.
collect_valid_batch_ids <- function(obj) {
    keep_batches <- list()
    for (b in obj@expDesign) {
        ctrl <- b$control
        trt <- b$treatment
        if (length(ctrl) == 0 || length(trt) == 0) next
        ok <- function(id) {
            ex <- obj@experiment[[id]]
            if (is.null(ex)) {
                return(FALSE)
            }
            d <- ex@data
            is.data.frame(d) && nrow(d) >= 2 && all(!is.na(d$time)) && all(!is.na(d$volume))
        }
        ctrl_ok <- unique(Filter(ok, ctrl))
        trt_ok <- unique(Filter(ok, trt))
        if (length(ctrl_ok) > 0 && length(trt_ok) > 0) {
            keep_batches[[b$batch.name]] <- c(ctrl_ok, trt_ok)
        }
    }
    keep_batches
}

# compute angle/abc/TGI robustly (per-batch with isolated batch design)
compute_batch_metrics_robust <- function(obj) {
    kb <- collect_valid_batch_ids(obj)
    if (!length(kb)) {
        message("No valid batches with control+treatment data; batch metrics will be empty.")
        return(data.frame())
    }

    out_list <- list()
    for (bn in names(kb)) {
        b <- obj@expDesign[[bn]]
        ids <- kb[[bn]]
        ctrl_ids <- intersect(b$control, ids)
        trt_ids <- intersect(b$treatment, ids)
        if (!length(ctrl_ids) || !length(trt_ids)) next

        # Clean batch definition restricted to valid models
        b_clean <- b
        b_clean$control <- ctrl_ids
        b_clean$treatment <- trt_ids

        # Construct single-batch sub-object directly (bypasses subsetXeva multi-arm leakage and + name bugs)
        xb <- obj
        xb@experiment <- obj@experiment[ids]
        xb@model <- obj@model[obj@model$model.id %in% ids, , drop = FALSE]
        xb@expDesign <- setNames(list(b_clean), bn)

        br <- try(
            {
                xb <- setResponse(xb, res.measure = c("angle", "abc", "TGI"), max.time = NULL, verbose = FALSE)
                res <- batch_response(xb)
                res[res$batch.id == bn, ]
            },
            silent = TRUE
        )

        if (inherits(br, "try-error") || is.null(br) || nrow(br) == 0) {
            message("Skipping batch (cannot compute): ", bn)
        } else {
            out_list[[length(out_list) + 1]] <- br
        }
    }

    if (!length(out_list)) {
        return(data.frame())
    }
    data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
}

# flatten any weird column shape to a simple atomic vector (for CSV safety)
flatten_col <- function(v) {
    if (is.null(v)) {
        return(v)
    }
    if (is.matrix(v)) {
        return(as.vector(v))
    }
    if (is.data.frame(v)) {
        return(as.vector(as.matrix(v)))
    }
    if (is.list(v)) {
        return(vapply(v, function(u) {
            if (is.null(u)) {
                return(NA_character_)
            }
            if (length(u) == 0) {
                return(NA_character_)
            }
            if (is.list(u)) {
                if (!is.null(u$join.name)) {
                    return(as.character(u$join.name))
                }
                return(paste(unlist(u, use.names = FALSE), collapse = "|"))
            }
            if (length(u) > 1) {
                return(paste(as.character(u), collapse = "|"))
            }
            as.character(u)
        }, ""))
    }
    v
}

# write CSV safely (coerce ragged/list columns; keep numerics numeric where obvious)
safe_write_csv <- function(df, path, add_index = TRUE) {
    if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
    for (nm in names(df)) {
        df[[nm]] <- flatten_col(df[[nm]])
        if (nm %in% c("time", "volume", "volume.normal", "slope", "AUC", "angle", "abc", "TGI")) {
            suppressWarnings({
                df[[nm]] <- as.numeric(df[[nm]])
            })
        }
    }
    if (isTRUE(add_index)) {
        rownames(df) <- seq_len(nrow(df)) # force numeric row names 1..N
    } else {
        rownames(df) <- NULL
    }
    write.csv(df, file = path, row.names = isTRUE(add_index), na = "")
}


# -- safe model_response: try Xeva's, else build from x@sensitivity$model
safe_model_response <- function(obj) {
    out <- try(model_response(obj), silent = TRUE)
    if (!inherits(out, "try-error") && is.data.frame(out)) {
        return(out)
    }

    sm <- try(obj@sensitivity$model, silent = TRUE)
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
        sm <- sm[, keep, drop = FALSE]
        return(sm)
    }

    data.frame(model.id = names(obj@experiment), stringsAsFactors = FALSE)
}

## ---------- helpers specific to assays ----------
# Convert SummarizedExperiment -> ExpressionSet (when needed)
as_eset <- function(se) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
        !requireNamespace("Biobase", quietly = TRUE)) {
        return(NULL)
    }
    mat <- SummarizedExperiment::assay(se)
    ph <- as.data.frame(SummarizedExperiment::colData(se))
    fe <- as.data.frame(SummarizedExperiment::rowData(se))
    Biobase::ExpressionSet(
        assayData = mat,
        phenoData = Biobase::AnnotatedDataFrame(ph),
        featureData = Biobase::AnnotatedDataFrame(fe)
    )
}

# Safely get a numeric matrix from any assay-like object
assay_to_matrix <- function(obj) {
    if (inherits(obj, "ExpressionSet")) {
        Biobase::exprs(obj)
    } else if (inherits(obj, "SummarizedExperiment")) {
        SummarizedExperiment::assay(obj)
    } else if (is.matrix(obj)) {
        obj
    } else {
        as.matrix(obj)
    }
}

## ========== PREP ==========
log_null_experiments(x, "initial")

x <- ensure_drug_slots(x)

keep <- nonempty_ids(x)
cat("Non-empty experiments (found):", length(keep), "of", length(x@experiment), "\n")
if (!length(keep)) stop("All experiments are NULL/empty.")
if (length(keep) < length(x@experiment)) x <- safe_subset(x, keep, keep.batch = FALSE)
log_null_experiments(x, "post_nonempty_subset")

## ========== Normalize time (start at 0) ==========
max.time <- 60
cutAtMaxTime <- FALSE
for (id in names(x@experiment)) {
    ex <- x@experiment[[id]]
    if (is.null(ex)) next
    d <- ex@data
    if (!is.null(d$time) && length(d$time)) {
        if (!is.na(d$time[1]) && d$time[1] > 0) d$time <- d$time - d$time[1]
        if (cutAtMaxTime) d <- d[d$time <= max.time, , drop = FALSE]
    }
    x@experiment[[id]]@data <- d
}

## ========== Responses ==========
max.time <- 30
x <- setResponse(x, res.measure = "mRECIST", max.time = NULL, verbose = FALSE)
x <- setResponse(x, res.measure = c("mRECIST", "slope", "AUC"), max.time = NULL, verbose = FALSE)

## ========== Batch metrics (robust) ==========
brf <- compute_batch_metrics_robust(x)

## ========== Exports ==========
# model_information.csv
m <- try(get_model_info(x), silent = TRUE)
if (inherits(m, "try-error")) m <- try(modelInfo(x), silent = TRUE)
if (inherits(m, "try-error") || is.null(m)) m <- data.frame(model.id = names(x@experiment))

# fix drug column if it's list(join.name=...)
if ("drug" %in% names(m) && is.list(m$drug)) {
    m$drug <- vapply(m$drug, function(v) if (is.list(v) && !is.null(v[["join.name"]])) v[["join.name"]] else as.character(v), "")
} else if (!"drug" %in% names(m)) {
    m$drug <- vapply(x@experiment, function(ex) {
        if (!is.null(ex) && "drug" %in% slotNames(ex)) {
            v <- ex@drug
            if (is.list(v) && !is.null(v[["join.name"]])) as.character(v[["join.name"]]) else as.character(v)
        } else {
            NA_character_
        }
    }, "")
}
if (!"tissue" %in% names(m)) m$tissue <- if (is_tnbc) "Breast Cancer" else "Lung Cancer"
if (!"patient.id" %in% names(m)) m$patient.id <- NA
m$dataset <- if (is_tnbc) "TNBC_v2" else paste0("Tsao_v2 (", m$tissue, ")")

mi_cols <- intersect(c("model.id", "tissue", "patient.id", "drug", "dataset"), names(m))
safe_write_csv(m[, mi_cols, drop = FALSE], file.path(OUTDIR, "model_information.csv"))

# batch_information.csv
b <- get_batch_info(x)
mi <- try(modelInfo(x), silent = TRUE)
if (!inherits(mi, "try-error") && "model.id" %in% colnames(b)) b <- b[b$model.id %in% mi$model.id, ]
safe_write_csv(b, file.path(OUTDIR, "batch_information.csv"))

## ---------- drug_screening.csv (robust) ----------
cat("Before drug screening\n")
mdf <- try(drug_screening(x), silent = TRUE)
if (inherits(mdf, "try-error")) {
    # fallback: build from experiments
    mdfl <- lapply(names(x@experiment), function(id) {
        d <- x@experiment[[id]]@data
        d$model.id <- id
        d
    })
    mdf <- data.table::rbindlist(mdfl, fill = TRUE)
}
mdf <- as.data.frame(mdf, stringsAsFactors = FALSE, check.names = FALSE)

# normalize column name variants
if (!"model.id" %in% names(mdf) && "modelid" %in% names(mdf)) {
    names(mdf)[names(mdf) == "modelid"] <- "model.id"
}

# add/repair drug column from experiments if missing
if (!"drug" %in% names(mdf)) {
    ex_drug <- data.frame(
        model.id = names(x@experiment),
        drug = vapply(x@experiment, function(ex) {
            if (!is.null(ex) && "drug" %in% slotNames(ex)) {
                v <- ex@drug
                if (is.list(v) && !is.null(v[["join.name"]])) as.character(v[["join.name"]]) else as.character(v)
            } else {
                NA_character_
            }
        }, ""),
        stringsAsFactors = FALSE
    )
    if ("model.id" %in% names(mdf)) {
        mdf <- merge(mdf, ex_drug, by = "model.id", all.x = TRUE)
    }
}

# Build the output explicitly (avoid ragged subsetting)
mk <- function(nm, numeric = FALSE) {
    if (!nm %in% names(mdf)) {
        return(NULL)
    }
    v <- flatten_col(mdf[[nm]])
    if (numeric) suppressWarnings(v <- as.numeric(v))
    v
}
ds_out <- data.frame(
    model.id = mk("model.id"),
    drug = mk("drug"),
    time = mk("time", numeric = TRUE),
    volume = mk("volume", numeric = TRUE),
    volume.normal = mk("volume.normal", numeric = TRUE),
    stringsAsFactors = FALSE, check.names = FALSE
)

# compute volume.normal if missing
if (!"volume.normal" %in% names(mdf) || all(is.na(ds_out$volume.normal))) {
    if (!is.null(ds_out$time) && !is.null(ds_out$volume) && !is.null(ds_out$model.id)) {
        dt <- data.table::as.data.table(ds_out)
        dt[, volume.normal := {
            base <- volume[which.min(ifelse(is.na(time), Inf, time))]
            if (is.na(base) || base == 0) NA_real_ else volume / base
        }, by = model.id]
        ds_out <- as.data.frame(dt, stringsAsFactors = FALSE, check.names = FALSE)
    }
}

# prune obviously unusable rows (missing required fields)
req <- intersect(c("model.id", "time", "volume"), names(ds_out))
if (length(req)) {
    keep_rows <- rowSums(is.na(ds_out[req])) == 0
    ds_out <- ds_out[keep_rows, , drop = FALSE]
}

# debug dump of shapes we saw
classes <- vapply(mdf, function(col) class(col)[1], "")
utils::write.table(data.frame(name = names(mdf), class = classes),
    file = file.path(OUTDIR, "drug_screening_column_classes.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

cat("After drug screening\n")
safe_write_csv(ds_out, file.path(OUTDIR, "drug_screening.csv"))

# model_response.csv (safe)
mres <- safe_model_response(x)
safe_write_csv(mres, file.path(OUTDIR, "model_response.csv"))

# batch_response.csv  (from robust computation)
safe_write_csv(brf, file.path(OUTDIR, "batch_response.csv"))

# modelid_moleculardata_mapping.csv
mmap <- modelid_moleculardata_mapping(x)
safe_write_csv(mmap, file.path(OUTDIR, "modelid_moleculardata_mapping.csv"))

## ========== Expression (robust to SE/Eset/shapes) ==========
write_empty_expr <- function(path) {
    safe_write_csv(data.frame(feature = character(0), sample = character(0), value = double(0)), path)
}

expr_path <- file.path(OUTDIR, "rna_sequencing.csv")
expr_ok <- FALSE
try(
    {
        rna_assay <- .get_assay(x, c("RNASeq", "RNAseq", "rna"))
        # Coerce to ExpressionSet if needed
        if (inherits(rna_assay, "SummarizedExperiment")) {
            eset <- as_eset(rna_assay)
        } else if (inherits(rna_assay, "ExpressionSet")) {
            eset <- rna_assay
        } else {
            # as last resort, try to wrap a matrix into an ExpressionSet
            mat <- assay_to_matrix(rna_assay)
            if (!requireNamespace("Biobase", quietly = TRUE)) stop("Biobase not available")
            eset <- Biobase::ExpressionSet(assayData = mat)
        }

        # Filter to protein_coding if that annotation exists
        fd <- Biobase::fData(eset)
        keep_rows <- rep(TRUE, nrow(eset))
        if (nrow(fd)) {
            if ("gene_type" %in% colnames(fd)) {
                keep_rows <- fd$gene_type == "protein_coding"
            } else if ("gene_biotype" %in% colnames(fd)) {
                keep_rows <- fd$gene_biotype == "protein_coding"
            }
        }
        if (length(keep_rows) && any(keep_rows %in% TRUE)) {
            eset <- eset[which(keep_rows), ]
        }

        # Tidy feature names if present
        fd2 <- Biobase::fData(eset)
        if (nrow(fd2)) {
            if ("gene_name" %in% colnames(fd2)) {
                Biobase::featureNames(eset) <- make.names(fd2$gene_name, unique = TRUE)
            } else if ("hgnc_symbol" %in% colnames(fd2)) {
                Biobase::featureNames(eset) <- make.names(fd2$hgnc_symbol, unique = TRUE)
            }
        }

        # Map to modToBiobase (rnaseq) if possible; otherwise keep as-is
        idmap <- x@modToBiobaseMap
        idmap <- try(idmap[tolower(idmap$mDataType) == "rnaseq", , drop = FALSE], silent = TRUE)
        emat <- Biobase::exprs(eset)

        if (!inherits(idmap, "try-error") && nrow(idmap) && "biobase.id" %in% colnames(idmap)) {
            comm <- intersect(idmap$biobase.id, colnames(emat))
            if (length(comm)) {
                emat <- emat[, comm, drop = FALSE]
            } else {
                message("RNASeq: no overlap with modToBiobaseMap; keeping original columns.")
            }
        }

        # Scale per-gene if we still have data
        if (nrow(emat) && ncol(emat)) {
            emat_scaled <- try(t(scale(t(emat))[, ]), silent = TRUE)
            if (!inherits(emat_scaled, "try-error")) {
                emat <- emat_scaled
            }
            dfexp <- getFlatDF(emat)
            safe_write_csv(dfexp, expr_path)
            expr_ok <- TRUE
        }
    },
    silent = TRUE
)

if (!expr_ok) {
    message("RNASeq assay missing or could not be processed; writing empty rna_sequencing.csv")
    write_empty_expr(expr_path)
}

## ========== Mutation ==========
mut_path <- file.path(OUTDIR, "mutation.csv")
mut_ok <- FALSE
try(
    {
        mut_assay <- .get_assay(x, c("mutation", "mut"))
        mut_mat <- assay_to_matrix(mut_assay) # handles list/MAE
        dfm <- getFlatDF(mut_mat) # usually: feature, sample, value

        # Standardize column names to what your older datasets show
        nm <- names(dfm)
        nm[nm == "feature"] <- "gene.id"
        nm[nm == "sample"] <- "sequencing.uid"
        names(dfm) <- nm

        # Force character, then map empties to "0"
        dfm$value <- as.character(dfm$value)
        dfm$value[is.na(dfm$value) | dfm$value == ""] <- "0"

        # (Optional) quick sanity logs
        message("Mutation rows total: ", nrow(dfm))
        message("Non-zero mutations: ", sum(dfm$value != "0", na.rm = TRUE))

        safe_write_csv(dfm, mut_path)
        mut_ok <- TRUE
    },
    silent = TRUE
)

if (!mut_ok) {
    message("Mutation assay missing or could not be processed; writing empty mutation.csv")
    safe_write_csv(
        data.frame(gene.id = character(0), sequencing.uid = character(0), value = character(0)),
        mut_path
    )
}


## ========== CNV ==========
cnv_path <- file.path(OUTDIR, "copy_number_variation.csv")
cnv_ok <- FALSE
try(
    {
        cnvCat <- c(
            "-2" = "Deletion",
            "-1" = "Shallow Deletion",
            "0"  = "0",
            "1"  = "Gain",
            "2"  = "Amplification"
        )
        cnv_assay <- .get_assay(x, c("CNV", "cnv", "copy"))
        cnv_mat <- assay_to_matrix(cnv_assay)
        dfc <- getFlatDF(cnv_mat)
        dfc$value <- as.character(dfc$value)
        dfc$value <- cnvCat[dfc$value]
        dfc$value[is.na(dfc$value) | dfc$value == "" | dfc$value == "Diploid"] <- "0"
        safe_write_csv(dfc, cnv_path)
        cnv_ok <- TRUE
    },
    silent = TRUE
)
if (!cnv_ok) {
    message("CNV assay missing or could not be processed; writing empty copy_number_variation.csv")
    safe_write_csv(data.frame(feature = character(0), sample = character(0), value = character(0)), cnv_path)
}

message("Export complete -> ", normalizePath(OUTDIR))
