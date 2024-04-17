library(Xeva)
source("xevaDB_fun.R")
#pdxe=readRDS("~/CXP/Xeva_dataset/data/XevaObjects/XevaSets/data/PDXE.rds")
pdxe=readRDS("~/CXP/XG/PDXE_xevaSet/data/PDXE_xevaSet/PDXE_XevaSet_All.rds")

max.time = 30

pdxe  <- setResponse(pdxe, res.measure = c("mRECIST", "slope", "AUC"), 
                     max.time = max.time, verbose=FALSE)

pdxe  <- setResponse(pdxe, res.measure = c("angle", "abc", "TGI"), 
                     max.time = max.time, verbose=FALSE)

##---------------model_info-----------------------------------------
m = get_model_info(pdxe)
write.csv(m, file = "results/pdxe/model_information.csv")

##---------------batch_information----------------------------------
b = get_batch_info(pdxe)
write.csv(b, file = "results/pdxe/batch_information.csv")

##-----drug_screening-----------------------------------
mdf = drug_screening(pdxe)
write.csv(mdf, file = "results/pdxe/drug_screening.csv")

##----------model_response---------------------------------
mres <- model_response(pdxe)
write.csv(mres, file = "results/pdxe/model_response.csv")

##----------batch_response---------------------------------
brf <- batch_response(pdxe)
write.csv(brf, file = "results/pdxe/batch_response.csv")

###---------modelid_moleculardata_mapping ------------------
mmap = modelid_moleculardata_mapping(pdxe)
write.csv(mmap, file = "results/pdxe/modelid_moleculardata_mapping.csv")

##-----------expression ----------
df=expression(pdxe)
write.csv(df, file = "results/pdxe/rna_sequencing.csv")

##-----------mutation ----------
df=mutation(pdxe)
write.csv(df, file = "results/pdxe/mutation.csv")

##-----------copy_number_variation ----------
df=cnv(pdxe)
write.csv(df, file = "results/pdxe/copy_number_variation.csv")

####-----------------------

