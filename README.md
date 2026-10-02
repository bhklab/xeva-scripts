# xeva-scripts — PDX Data Curation & Extraction Pipeline

This repository provides reproducible R scripts to extract, standardize, and export Patient-Derived Xenograft (PDX) pharmacogenomic data from raw **Xeva** R objects (`.rds`) into flat, query-ready CSV files for ingestion into **XevaDB**.

---
### Prerequisites
- R ≥ 4.0
- Required packages:
  ```r
  install.packages(c("data.table", "reshape2"))
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("Xeva", "Biobase"))
  ```

### Step 1: Place Input RDS Objects
The scripts expect the raw RDS objects in `xevasets_obj_2026/`:
- `xevasets_obj_2026/UHN_Tsao_Lung_DrugResponse_2022_v1.rds` (Tsao KRAS Lung)
- `xevasets_obj_2026/Xeva_McGill_v2.rds` (McGill TNBC)
- `xevasets_obj_2026/Xeva_PDXE_v2.rds` (Novartis PDXE)
- `xevasets_obj_2026/UHN_Cescon_Breast_DrugResponse_2025_v1.rds` (Cescon Breast Cancer)

### Step 2: Run Extraction Scripts
Run each script from the `src/` directory:

```bash
cd src/

# 1. Tsao UHN KRAS Lung Cancer
Rscript Tsao_to_XevaDB.R

# 2. McGill Triple-Negative Breast Cancer (TNBC)
Rscript MP_to_XevaDB.R

# 3. Novartis PDXE (Multi-tissue collection)
Rscript pdxe_to_XevaDB.R

# 4. Cescon UHN Breast Cancer
Rscript tnbc_to_XevaDB.R
```

Outputs will be saved into separate subdirectories in `results_2026_sep_23/`.

### Step 3: Run the Automated Validation Suite
Verify data integrity across all generated tables:

```bash
Rscript test_results.R
```
This tests for missing files, required columns, valid categorical ranges, numeric validity, and cross-file foreign key consistency.

---

## Curated Datasets Overview

| Dataset | Tissue | Source Script | Output Folder | Models | Batches | Drugs |
|---|---|---|---|---|---|---|
| **Tsao (KRAS Lung)** | Non-Small Cell Lung Cancer | `src/Tsao_to_XevaDB.R` | `results_2026_sep_23/KRAS_LUNG_v2/` | 1,415 | 225 | 36 |
| **McGill (TNBC)** | Triple-Negative Breast Cancer | `src/MP_to_XevaDB.R` | `results_2026_sep_23/McGill_TNBC_v2/` | 241 | 135 | 6 |
| **PDXE (Novartis)** | Multi-tissue (Colorectal, Breast, Lung, Pancreas, Stomach, Skin) | `src/pdxe_to_XevaDB.R` | `results_2026_sep_23/pdxe_v2/` | 4,706 | 3,968 | 62 |
| **Cescon (Breast)** | Breast Cancer | `src/tnbc_to_XevaDB.R` | `results_2026_sep_23/TNBC_v2/` | 907 | 105 | 12 |

---

## Standardized Output Files

For each dataset, the pipeline generates 9 core CSV files and 1 schema metadata file:

```
results_2026_sep_23/<dataset>/
├── model_information.csv             # Model metadata (patient ID, tissue, drug)
├── model_response.csv                # Model-level efficacy metrics (mRECIST, AUC, slope, survival)
├── batch_information.csv             # Batch-to-model mapping (control vs. treatment)
├── batch_response.csv                # Batch-level response metrics (angle, abc, TGI)
├── drug_screening.csv                # Full raw longitudinal tumor volume time series
├── drug_screening_column_classes.tsv # Data types for time series columns
├── modelid_moleculardata_mapping.csv # Links model IDs to molecular sample IDs
├── mutation.csv                      # Somatic mutation calls (gene × sample)
├── copy_number_variation.csv         # Gene copy number calls (gene × sample)
└── rna_sequencing.csv                # Normalized gene expression (gene × sample)
```

### Table Definitions

| File | Primary Key / Index | Description |
|---|---|---|
| `model_information.csv` | `model.id` | One row per mouse/model. Contains `model.id`, `patient.id`, `tissue`, `drug`, and `dataset`. |
| `model_response.csv` | `model.id` × `response_type` | Long-format response metrics: `mRECIST`, `best.average.response`, `slope`, `AUC`, and `survival`. |
| `batch_information.csv` | `batch.id` × `model.id` | Groups control (untreated) and treatment arms into experimental batches. |
| `batch_response.csv` | `batch.id` × `response_type` | Batch metrics comparing treated vs. control growth curves: `angle`, `abc` (area between curves), and `TGI` (tumor growth inhibition). |
| `drug_screening.csv` | `model.id` × `time` | **Untruncated** raw measurement history: day of measurement (`time`), tumor size in mm³ (`volume`), and volume normalized to baseline (`volume.normal`). |
| `modelid_moleculardata_mapping.csv` | `model.id` × `mDataType` | Maps PDX models to molecular assay sample IDs (`biobase.id`). |
| `mutation.csv` | `gene.id` × `sequencing.uid` | Long-format mutation matrix (`0` = wild-type/silent, `1` or mutation consequence label = mutated). |
| `copy_number_variation.csv` | `gene.id` × `sequencing.uid` | Categorical copy number calls: `'Deletion'`, `'Shallow Deletion'`, `'0'` (diploid), `'Gain'`, `'Amplification'`. |
| `rna_sequencing.csv` | `gene.id` × `sequencing.uid` | Z-score standardized gene expression values. |

---

## Core Design Choices & Key Parameters

> [!IMPORTANT]
> **`max.time = 30` is applied across ALL FOUR datasets (`KRAS_LUNG`, `McGill_TNBC`, `pdxe`, and `Cescon_Breast`).**

#### What does `max.time` do?
In PDX studies, cancer tissue from a patient is implanted into mice. Some mice are monitored for only 14–21 days because their tumors grow very aggressively and reach the humane endpoint quickly. Other mice are treated with effective drugs and their tumors shrink or stay dormant for 100 to 400+ days.

If we calculate curve metrics like **AUC (Area Under the Curve)** or **Slope** across the entire duration of every experiment:
- A mouse tracked for 14 days will have an AUC calculated over 14 days.
- A mouse tracked for 150 days will have an AUC calculated over 150 days.

The 150-day mouse would naturally have a much larger AUC simply because the experiment lasted 10 times longer, not because the drug was less effective, this is not our desired outcome.

#### Why 30 days?
By truncating growth curve calculations at **30 days** (`max.time = 30`), every model's response is assessed over the **exact same standardized 30-day window**. This makes response metrics directly comparable across drugs, models, and datasets (across datasets isn't available yet in the platform but will be eventually) in XevaDB.

#### What does `max.time = 30` affect vs. not affect?
- **Affected:** `mRECIST`, `slope`, `AUC`, `best.average.response`, `angle`, `abc`, and `TGI`.
- **NOT Affected:**
  - `survival`: Always reports the animal's true maximum observation day (up to 253 days in Tsao, 404 days in PDXE).
  - `drug_screening.csv`: Retains 100% of all raw timepoints with zero truncation.

---

### 2. `min.time = 10` and Unclassified (`NA`) mRECIST Calls

#### What is mRECIST?
mRECIST is an adaptation of clinical RECIST criteria for mice. It classifies response into:
- **CR** (Complete Response): Tumor shrinks completely.
- **PR** (Partial Response): Tumor shrinks significantly.
- **SD** (Stable Disease): Tumor does not grow or shrink substantially.
- **PD** (Progressive Disease): Tumor grows rapidly.

#### Why do some models have `NA` for mRECIST?
To distinguish genuine stable disease from an experiment that simply ended before the tumor had time to grow, Xeva requires at least **10 days** of data (`min.time = 10`). Models that finished or were sacrificed in fewer than 10 days receive `NA`. These `NA` values are statistically expected and correct.

---

### 3. Safety Patch: Handling Single-Point Models in `Xeva::TGI`

> [!NOTE]
> **Dataset affected:** Handled inside `src/xevaDB_fun.R` (used by **PDXE** and **McGill_TNBC**).

#### The Problem
In the Novartis PDXE collection (4,706 models and 3,968 batches), some treatment arms have only 1 data point measured before day 30 (or only 1 data point in their entire experimental lifetime, such as model `X.2088.CG97` at day 15).

`Xeva` requires at least 2 points to compute a treatment trajectory. When fewer than 2 points exist:
1. `Xeva` sets `treat.volume = NULL`.
2. The built-in `TGI()` function attempts:
   ```R
   tgi = (contr.volume[length(contr.volume)] - treat.volume[length(treat.volume)]) / ...
   ```
   Because `treat.volume` is `NULL`, this collapses into `numeric(0)` (a zero-length vector).
3. Assigning `numeric(0)` into a data frame cell crashes R:
   ```text
   Error in x[[jj]][iseq] <- vjj : replacement has length zero
   ```

#### The Solution (Safety Patch)
In `src/xevaDB_fun.R`, we patch `Xeva::TGI` so that whenever an arm has $< 1$ valid timepoint or produces `numeric(0)`, it safely returns `NA_real_`:

```R
try({
  ns <- asNamespace("Xeva")
  if (bindingIsLocked("TGI", ns)) unlockBinding("TGI", ns)
  orig_TGI <- get("TGI", ns)
  assign("TGI", function(contr.volume, treat.volume) {
    if (is.null(contr.volume) || is.null(treat.volume) || 
        length(contr.volume) < 1 || length(treat.volume) < 1) {
      return(ns$batch_response_class(name = "TGI", value = NA_real_))
    }
    res <- orig_TGI(contr.volume, treat.volume)
    if (length(res$value) == 0) res$value <- NA_real_
    return(res)
  }, ns)
  lockBinding("TGI", ns)
}, silent = TRUE)
```
This allows `setResponse()` to finish processing all 3,968 PDXE batches cleanly.

---

### 4. Zero-Variance Standardization in RNA-seq

> [!NOTE]
> **Dataset affected:** Handled in `src/xevaDB_fun.R` for **PDXE** (and applicable to any unexpressed genes).

#### Plain English Explanation
RNA sequencing measures the activity level of thousands of genes. To make expression levels comparable across different samples, each gene's expression is normalized using a **z-score**:

$$\text{Z-score} = \frac{\text{Value} - \text{Mean}}{\text{Standard Deviation}}$$

If a gene is completely inactive in every single tumor sample (expression is 0 everywhere), its standard deviation across samples is **0**. Dividing by 0 in R produces `NaN` ("Not a Number").

#### How it is handled
In `src/xevaDB_fun.R` (`expression()`), any gene with zero variance across all samples has its z-score set to **`0`** (baseline/average) instead of `NaN`:
```R
scaled <- t(scale(t(df))[, ])
scaled[is.nan(scaled)] <- 0
```
In PDXE, this cleanly resolves 1,395 unexpressed genes across 375 samples into valid numeric values.

---

### 5. Molecular Annotation & PDX Filtering (McGill TNBC)

> [!NOTE]
> **Dataset affected:** Handled in `src/MP_to_XevaDB.R` for **McGill TNBC**.

#### PDX vs. Patient Filtering
The raw `Xeva_McGill_v2.rds` object contains molecular profiles from both the original human patient tumors and the mouse PDX models. The script filters for `Source == "PDX"` across RNA-seq, CNV, and mutation data so only xenograft molecular profiles are mapped.

#### Modern GENCODE Column Mapping
In `Xeva_McGill_v2.rds`, feature annotations follow modern GENCODE conventions:
- Biotype is stored under `gencode.biotype` (filtered for `"protein_coding"`).
- Gene symbols are stored under `gencode.symbol` (e.g. `TSPAN6`, `DPM1`).

`src/MP_to_XevaDB.R` dynamically resolves these column names with fallback to Ensembl IDs for any unnamed loci.

---

### 6. Copy Number Variation (CNV) Categorization

In **XevaDB**, all copy number calls are standardized into discrete categorical strings:
- **Deletion**: Complete / homozygous deletion (loss of both alleles).
- **Shallow Deletion**: Heterozygous loss (loss of one allele).
- **0**: Neutral / diploid state (normal copy number).
- **Gain**: Low-level copy number gain (broad / low copy increase).
- **Amplification**: High-level focal amplification (multiple additional copies).

Depending on whether the source dataset provides continuous total copy numbers or pre-categorized GISTIC calls, categorization is handled as follows:

| Dataset | Input Format | Categorization Rule | Output States |
|---|---|---|---|
| **PDXE (Novartis)** | Continuous total copy number | Continuous cutoffs applied:<br>• $\text{Copy Number} > 5 \rightarrow$ `"Amplification"`<br>• $\text{Copy Number} \le 0.8 \rightarrow$ `"Deletion"`<br>• $\text{Otherwise} \rightarrow$ `"0"` (Diploid / neutral) | `"Deletion"`, `"0"`, `"Amplification"` |
| **McGill (TNBC)** | Pre-categorized GISTIC calls | Integer states mapped directly:<br>• `-2` $\rightarrow$ `"Deletion"`<br>• `-1` $\rightarrow$ `"Shallow Deletion"`<br>• `0` (or `"Diploid"`) $\rightarrow$ `"0"`<br>• `1` $\rightarrow$ `"Gain"`<br>• `2` $\rightarrow$ `"Amplification"` | `"Deletion"`, `"Shallow Deletion"`, `"0"`, `"Gain"`, `"Amplification"` |
| **Cescon (Breast / TNBC)** | Pre-categorized GISTIC calls | Integer states mapped directly:<br>• `-2` $\rightarrow$ `"Deletion"`<br>• `-1` $\rightarrow$ `"Shallow Deletion"`<br>• `0` (or `"Diploid"`) $\rightarrow$ `"0"`<br>• `1` $\rightarrow$ `"Gain"`<br>• `2` $\rightarrow$ `"Amplification"` | `"Deletion"`, `"Shallow Deletion"`, `"0"`, `"Gain"`, `"Amplification"` |
| **Tsao (KRAS Lung)** | Pre-categorized copy number calls | Discretized states mapped directly:<br>• `-2` $\rightarrow$ `"Deletion"`<br>• `-1` $\rightarrow$ `"Shallow Deletion"`<br>• `0` (or `"Diploid"`) $\rightarrow$ `"0"`<br>• `1` $\rightarrow$ `"Gain"`<br>• `2` $\rightarrow$ `"Amplification"` | `"Deletion"`, `"Shallow Deletion"`, `"0"`, `"Gain"`, `"Amplification"` |

---

## Known Source Data Exceptions

### KRAS Lung: 6 Batches Without Matched Control Arms
In `results_2026_sep_23/KRAS_LUNG_v2/batch_information.csv`, 6 treatment batches do not have a corresponding vehicle control arm in the source publication data:
- `f7.PHLC137.TPX.erlotinib`
- `f18.PHLC164.TPX.afatinib`
- `f19.PHLC164.TPX.dacomitinib`
- `f67.PHLC77.TPX.erlotinib`
- `f70.PHLC77.TPX.afatinib`
- `f71.PHLC77.TPX.afatinib`

These batches are exported as-is. Batch-level comparison metrics (`angle`, `abc`, `TGI`) that require a control curve cannot be computed for them and are left blank/`NA`.

---

## File Schema Summary

```
model_information.csv:
  - model.id    (character): Unique PDX model identifier
  - patient.id  (character): De-identified patient identifier
  - tissue      (character): Tissue or cancer type
  - drug        (character): Administered drug or drug combination
  - dataset     (character): Originating dataset name

model_response.csv:
  - drug          (character): Drug name
  - model.id      (character): Unique PDX model identifier
  - response_type (character): Metric name (mRECIST, AUC, slope, best.average.response, survival)
  - value         (character/numeric): Calculated value or category (e.g., CR, PR, SD, PD)

batch_information.csv:
  - batch.id  (character): Unique batch identifier (patient.drug)
  - model.id  (character): PDX model assigned to batch
  - type      (character): Arm type ("control" or "treatment")

batch_response.csv:
  - batch.id      (character): Unique batch identifier
  - response_type (character): Metric name (angle, abc, TGI)
  - value         (numeric): Calculated batch response value

drug_screening.csv:
  - model.id      (character): Unique PDX model identifier
  - drug          (character): Drug name
  - time          (numeric): Day of measurement since treatment start
  - volume        (numeric): Absolute tumor volume in mm³
  - volume.normal (numeric): Volume relative to baseline (day 0 = 1.0)
```
