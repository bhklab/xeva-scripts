options(stringsAsFactors = FALSE)
library(Xeva)
source("xevaDB_fun.R")
#mpx = readRDS("~/CXP/XG/MP/Data/MP_XevaSet.rds")
#mpx = readRDS("/Users/mattbocc/uhn/xeva-scripts/xevasets-obj/McGill_XevaSet_2024.rds")
mpx = readRDS("/Users/mattbocc/uhn/xeva-scripts/xevasets_obj_2026/Xeva_McGill_v2.rds")

out_dir <- "../results_2026_sep_23/McGill_TNBC_v2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max.time = 30

mpx  <- setResponse(mpx, res.measure = c("mRECIST", "slope", "AUC"), 
                     max.time = max.time, verbose=FALSE)
mpx  <- setResponse(mpx, res.measure = c("angle", "abc", "TGI"), 
                     max.time = max.time, verbose=FALSE)

##=======subset only PDX mol data ==============================================
mpx@molecularProfiles$RNAseq = 
  mpx@molecularProfiles$RNAseq[,mpx@molecularProfiles$RNAseq$Source=="PDX"]

mpx@molecularProfiles$CNV = 
  mpx@molecularProfiles$CNV[,mpx@molecularProfiles$CNV$Source=="PDX"]

mpx@molecularProfiles$mutation = 
  mpx@molecularProfiles$mutation[,mpx@molecularProfiles$mutation$Source=="PDX"]


mbmap = mpx@modToBiobaseMap
mbmap = mbmap[!is.na(mbmap$biobase.id),]
#mbmap$mDataType[mbmap$mDataType=="CNV"] <- "cnv"

updateMBMAP <- function(mpx, mbmap, type)
{
  mx = mbmap[mbmap$mDataType==type, ]
  mx = mx[mx$biobase.id %in% sampleNames(mpx@molecularProfiles[[type]]),]
  mx = mx[!duplicated(mx$model.id), ]
  mol= mpx@molecularProfiles[[type]][,unique(mx$biobase.id)]
  return(list(mol=mol, map=mx))
}

rna = updateMBMAP(mpx, mbmap, type="RNAseq")
mut = updateMBMAP(mpx, mbmap, type="mutation")
cnv = updateMBMAP(mpx, mbmap, type="CNV")

mpx@molecularProfiles$RNAseq <- rna$mol
mpx@molecularProfiles$mutation <- mut$mol
mpx@molecularProfiles$CNV <- cnv$mol

mpx@modToBiobaseMap <- rbind(rna$map, mut$map, cnv$map)

##==============================================================================
##------------------------------------------------------------------------------
###response is already done
if(1==2){
max.time = 60
mpx  <- setResponse(mpx, res.measure = c("mRECIST", "slope", "AUC"), 
                         max.time = max.time, verbose=FALSE)

mpx  <- setResponse(mpx, res.measure = c("angle", "abc", "TGI"), 
                         max.time = max.time, verbose=F)
}
##---------------model_info-----------------------------------------
m = get_model_info(mpx)
m$tissue <- "Breast Cancer"
m$dataset <- paste0("McGill_v2 (", m$tissue, ")")
m <- m[,c("model.id","tissue","patient.id","drug","dataset")]
write.csv(m, file = file.path(out_dir, "model_information.csv"))

##---------------batch_information----------------------------------
b = get_batch_info(mpx)
write.csv(b, file = file.path(out_dir, "batch_information.csv"))

##-----drug_screening-----------------------------------
mdf = drug_screening(mpx)
write.csv(mdf, file = file.path(out_dir, "drug_screening.csv"))

classes <- vapply(mdf, function(col) class(col)[1], "")
utils::write.table(data.frame(name = names(mdf), class = classes),
    file = file.path(out_dir, "drug_screening_column_classes.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

##----------model_response---------------------------------
mres <- model_response(mpx)
write.csv(mres, file = file.path(out_dir, "model_response.csv"))

##----------batch_response---------------------------------
brf <- batch_response(mpx)
write.csv(brf, file = file.path(out_dir, "batch_response.csv"))

###---------modelid_moleculardata_mapping ------------------
mmap = modelid_moleculardata_mapping(mpx, dt = c("RNAseq","mutation","CNV"))
write.csv(mmap, file = file.path(out_dir, "modelid_moleculardata_mapping.csv"))

##-----------expression ----------
fd <- fData(mpx@molecularProfiles$RNAseq)

# Identify biotype column (new RDS uses 'gencode.biotype', old uses 'gene_type')
biotype_col <- if ("gene_type" %in% colnames(fd)) {
  "gene_type"
} else if ("gencode.biotype" %in% colnames(fd)) {
  "gencode.biotype"
} else if ("gene_biotype" %in% colnames(fd)) {
  "gene_biotype"
} else {
  NULL
}

if (!is.null(biotype_col)) {
  pc <- fd[!is.na(fd[[biotype_col]]) & fd[[biotype_col]] == "protein_coding", ]
  mpx@molecularProfiles$RNAseq <- mpx@molecularProfiles$RNAseq[rownames(pc), ]
}

# Identify gene symbol column (new RDS uses 'gencode.symbol', old uses 'gene_name')
symbol_col <- if ("gene_name" %in% colnames(fd)) {
  "gene_name"
} else if ("gencode.symbol" %in% colnames(fd)) {
  "gencode.symbol"
} else if ("hgnc_symbol" %in% colnames(fd)) {
  "hgnc_symbol"
} else {
  NULL
}

if (!is.null(symbol_col)) {
  syms <- fData(mpx@molecularProfiles$RNAseq)[[symbol_col]]
  missing_syms <- is.na(syms) | syms == ""
  syms[missing_syms] <- rownames(fData(mpx@molecularProfiles$RNAseq))[missing_syms]
  featureNames(mpx@molecularProfiles$RNAseq) <- make.names(syms, unique = TRUE)
}

df <- expression(mpx, dt = "RNAseq")
write.csv(df, file = file.path(out_dir, "rna_sequencing.csv"))

##----------- mutation ----------
##mut=exprs(mpx@molecularProfiles$mutation)
df <- getFlatDF(exprs(mpx@molecularProfiles$mutation))
#mutMap=rep("mutation", length(unique(df$value)))
#names(mutMap)=unique(df$value)
#mutMap[c("0", "Silent")]="0"

#df$value = mutMap[df$value]
write.csv(df, file = file.path(out_dir, "mutation.csv"))

##----------- cnv ----------
cnvCat <- c("-2"= "Deletion", # indicates a deep loss, possibly a homozygous deletion
            "-1"= "Shallow Deletion",#indicates a shallow loss, possibley a heterozygous deletion
            "0" = "0",
            "1" = "Gain", #indicates a low-level gain (a few additional copies, often broad)
            "2" = "Amplification" # indicate a high-level amplification (more copies, often focal)
)

df <- getFlatDF(exprs(mpx@molecularProfiles$CNV))
df$value <- as.character(df$value)
df$value <- cnvCat[df$value]
df$value[is.na(df$value) | df$value == "" | df$value == "Diploid"] <- "0"
write.csv(df, file = file.path(out_dir, "copy_number_variation.csv"))

