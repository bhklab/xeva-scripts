library(Xeva)
source("xevaDB_fun.R")
#pdxe=readRDS("~/CXP/Xeva_dataset/data/XevaObjects/XevaSets/data/PDXE.rds")
pdxe=readRDS("/Users/mattbocc/uhn/xeva-scripts/xevasets_obj_2026/Xeva_PDXE_v2.rds")

out_dir <- "../results_2026_sep_23/pdxe_v2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max.time = 30

pdxe  <- setResponse(pdxe, res.measure = c("mRECIST", "slope", "AUC"), 
                     max.time = max.time, verbose=FALSE)

pdxe  <- setResponse(pdxe, res.measure = c("angle", "abc", "TGI"), 
                     max.time = max.time, verbose=FALSE)

##---------------model_info-----------------------------------------
m = get_model_info(pdxe)
write.csv(m, file = file.path(out_dir, "model_information.csv"))

##---------------batch_information----------------------------------
b = get_batch_info(pdxe)
write.csv(b, file = file.path(out_dir, "batch_information.csv"))

##-----drug_screening-----------------------------------
mdf = drug_screening(pdxe)
write.csv(mdf, file = file.path(out_dir, "drug_screening.csv"))

classes <- vapply(mdf, function(col) class(col)[1], "")
utils::write.table(data.frame(name = names(mdf), class = classes),
    file = file.path(out_dir, "drug_screening_column_classes.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

##----------model_response---------------------------------
mres <- model_response(pdxe)
write.csv(mres, file = file.path(out_dir, "model_response.csv"))

##----------batch_response---------------------------------
brf <- batch_response(pdxe)
write.csv(brf, file = file.path(out_dir, "batch_response.csv"))

###---------modelid_moleculardata_mapping ------------------
mmap = modelid_moleculardata_mapping(pdxe)
write.csv(mmap, file = file.path(out_dir, "modelid_moleculardata_mapping.csv"))

##-----------expression ----------
df=expression(pdxe)
write.csv(df, file = file.path(out_dir, "rna_sequencing.csv"))

##-----------mutation ----------
df=mutation(pdxe)
write.csv(df, file = file.path(out_dir, "mutation.csv"))

##-----------copy_number_variation ----------
df=cnv(pdxe)
write.csv(df, file = file.path(out_dir, "copy_number_variation.csv"))

####-----------------------

