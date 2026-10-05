options(stringsAsFactors = FALSE)
library(Xeva)
library(readxl)
source("xevaDB_fun.R")

tnbc <- readRDS("/Users/mattbocc/uhn/xeva-scripts/xevasets_obj_2026/UHN_Cescon_Breast_DrugResponse_2025_v1.rds")
out_dir <- "../results_2026_sep_23/TNBC_v2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# In-memory sanitization: make.names ensures model IDs (e.g. leading numbers or "-"/"+")
# match the rownames that Xeva::modelInfo generates, preventing subsetXeva from failing.
# Every slot keyed by model.id must be sanitized together. subsetXeva() indexes
# @model BY ROWNAME, so if rownames(tnbc@model) keep the raw IDs, every model whose ID
# changes under make.names (leading digit, "-", "+") becomes an NA row and is silently
# dropped by the is.na(patient.id) filter below.
tnbc@model$model.id <- make.names(tnbc@model$model.id)
rownames(tnbc@model) <- tnbc@model$model.id
names(tnbc@experiment) <- make.names(names(tnbc@experiment))
for (i in seq_along(tnbc@expDesign)) {
  tnbc@expDesign[[i]]$treatment <- make.names(tnbc@expDesign[[i]]$treatment)
  tnbc@expDesign[[i]]$control <- make.names(tnbc@expDesign[[i]]$control)
}
tnbc@sensitivity$model$model.id <- make.names(tnbc@sensitivity$model$model.id)
rownames(tnbc@sensitivity$model) <- tnbc@sensitivity$model$model.id
tnbc@modToBiobaseMap$model.id <- make.names(tnbc@modToBiobaseMap$model.id)
for (i in seq_along(tnbc@experiment)) {
  tnbc@experiment[[i]]@model.id <- names(tnbc@experiment)[i]
}

##-----for non RES and remove DMSO -------
mi <- modelInfo(tnbc)
mid <- mi[!grepl("_RES", mi$model.id, ignore.case=TRUE) & mi$drug != "DMSO", ]

tnbc.nor <- subsetXeva(tnbc, ids=mid$model.id, id.name="model.id", keep.batch = FALSE)


##----------------------------------
mi <- modelInfo(tnbc.nor)
mi <- mi[!is.na(mi$patient.id),]

# mi$patient.id <- paste0("P.", mi$patient.id)
# pid2remove=c("P.1", "P.100534", "P.2018-12-28", "P.2019-01-16", "P.44721", 
#              "P.5/18", "P.85227", "P.REF020", "P.2018-12-13", "P.5")
# mi <- mi[! mi$patient.id %in% pid2remove,]
tnbc.nor <- subsetXeva(tnbc.nor, ids=mi$model.id, id.name="model.id", 
                       keep.batch = F)
###-------------------------------
max.time = 30; cutAtMaxTime = FALSE
##----make sure curve start at 0 & cut everything at max.time ---
for(i in names(tnbc.nor@experiment))
{
  d = tnbc.nor@experiment[[i]]@data
  if(d$time[1]>0) { d$time <- d$time-d$time[1] }
  if(cutAtMaxTime==TRUE) {  d = d[d$time <= max.time, ] }
  
  tnbc.nor@experiment[[i]]@data = d
}

##-------------------------------------

max.time = 30
tnbc.nor  <- setResponse(tnbc.nor, res.measure = "mRECIST",max.time = max.time, verbose=FALSE)

mr <- summarizeResponse(tnbc.nor, response.measure="mRECIST", group.by="patient.id")

nonNAcol = sort(apply(mr, 2, function(i)sum(!is.na(i))))
mr2 <- mr[, names(nonNAcol)[nonNAcol>0]]

nonNArow = sort(apply(mr2, 1, function(i)sum(!is.na(i))))
mr3 <- mr2[names(nonNArow)[nonNArow>0], ]

#plotmRECIST(mr3, control.name = "H2O")

mi=modelInfo(tnbc.nor)
mi=mi[mi$patient.id%in%colnames(mr3), ]
mi=mi[mi$drug%in%c(rownames(mr3), "H2O"), ]

tnbc.nor <- subsetXeva(tnbc.nor, ids=mi$model.id, id.name="model.id", 
                       keep.batch = F)

###-----------------------------------------------------------------
max.time = 30
tnbc.nor  <- setResponse(tnbc.nor, res.measure = c("mRECIST", "slope", "AUC"), 
                         max.time = max.time, verbose=FALSE)
tnbc.nor  <- setResponse(tnbc.nor, res.measure = c("angle", "abc", "TGI"), 
                         max.time = max.time, verbose=FALSE)

##---------------model_info-----------------------------------------
m = get_model_info(tnbc.nor)
m$dataset <- "TNBC_v2"; m$tissue <- "Breast Cancer"
mo <- m[,c("model.id","tissue","patient.id","drug","dataset")]
write.csv(mo, file = file.path(out_dir, "model_information.csv"))

##---------------batch_information----------------------------------
b = get_batch_info(tnbc.nor)
mi = modelInfo(tnbc.nor); b = b[b$model.id%in%mi$model.id,]
valid_batches <- names(which(tapply(b$type == "control", b$batch.id, any) & 
                             tapply(b$type == "treatment", b$batch.id, any)))
b <- b[b$batch.id %in% valid_batches, ]
write.csv(b, file = file.path(out_dir, "batch_information.csv"))

##-----drug_screening-----------------------------------

# Removed by Matthew because "dose1.value", "dose1.unit", "dose2.value", "dose2.unit" do not seem to exist
# mdf = drug_screening(tnbc.nor, col2take=c("model.id", "drug.join.name", "time",
#                                 "volume", "volume.normal",
#                                 "dose1.value", "dose1.unit",
#                                 "dose2.value", "dose2.unit"))
# mdf[mdf$drug=="H2O", c("dose1.value", "dose1.unit","dose2.value", "dose2.unit")] <- NA

print("Before drug screening")
mdf = drug_screening(tnbc.nor)
print("After drug screening")

write.csv(mdf, file = file.path(out_dir, "drug_screening.csv"))

classes <- vapply(mdf, function(col) class(col)[1], "")
utils::write.table(data.frame(name = names(mdf), class = classes),
    file = file.path(out_dir, "drug_screening_column_classes.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

##----------model_response---------------------------------
mres <- model_response(tnbc.nor)
write.csv(mres, file = file.path(out_dir, "model_response.csv"))

##----------batch_response---------------------------------
brf <- batch_response(tnbc.nor)
brf <- brf[brf$batch.id %in% b$batch.id, ]
write.csv(brf, file = file.path(out_dir, "batch_response.csv"))

###---------modelid_moleculardata_mapping ------------------
mmap = modelid_moleculardata_mapping(tnbc.nor, dt = c("RNASeq", "mutation", "CNV"))
mmap = unique(mmap)
write.csv(mmap, file = file.path(out_dir, "modelid_moleculardata_mapping.csv"))

##-----------expression ----------
fd = fData(tnbc.nor@molecularProfiles$RNASeq)
pc = fd[!is.na(fd$description) & fd$description == "protein_coding",]
tnbc.nor@molecularProfiles$RNASeq = tnbc.nor@molecularProfiles$RNASeq[rownames(pc),]
featureNames(tnbc.nor@molecularProfiles$RNASeq) = 
  make.names(fData(tnbc.nor@molecularProfiles$RNASeq)$hugo.id, unique = T)

df=expression(tnbc.nor)
write.csv(df, file = file.path(out_dir, "rna_sequencing.csv"))

##-----------mutation ----------
mut <- tnbc.nor@molecularProfiles$mutation
df = getFlatDF(exprs(mut))
df$value <- as.character(df$value)
df$value[is.na(df$value)] <- "0"
df$value[df$value==""] <- "0"
write.csv(df, file = file.path(out_dir, "mutation.csv"))

##-----------copy_number_variation ----------

# cnvCat <- c("-2"= "Deep Deletion", # indicates a deep loss, possibly a homozygous deletion
#             "-1"= "Shallow Deletion",#indicates a shallow loss, possibley a heterozygous deletion
#             "0" = "Diploid",
#             "1" = "Gain", #indicates a low-level gain (a few additional copies, often broad)
#             "2" = "Amplification" # indicate a high-level amplification (more copies, often focal)
# )

cnvCat <- c("-2" = "Deletion", # indicates a deep loss, possibly a homozygous deletion
            "-1" = "Shallow Deletion",#indicates a shallow loss, possibley a heterozygous deletion
            "0"  = "0",
            "1"  = "Gain", #indicates a low-level gain (a few additional copies, often broad)
            "2"  = "Amplification" # indicate a high-level amplification (more copies, often focal)
)

cnv <- tnbc.nor@molecularProfiles$CNV
df = getFlatDF(exprs(cnv))
df$value <- as.character(df$value)
df$value <- cnvCat[df$value]
df$value[is.na(df$value) | df$value == "" | df$value == "Diploid"] <- "0"
write.csv(df, file = file.path(out_dir, "copy_number_variation.csv"))

####-----------------------