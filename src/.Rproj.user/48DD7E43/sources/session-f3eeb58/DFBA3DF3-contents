options(stringsAsFactors = FALSE)
library(Xeva)
source("xevaDB_fun.R")

###lpdx <- readRDS("~/CXP/KRAS_P53/data/xeva_data/KRAS_P53_XevaData.Rda")

getFixed_LPDX_Data <- function()
{
  lpdx <- readRDS("~/CXP/KRAS_P53/data/xeva_data/KRAS_P53_XevaData.Rda")
  
  ##----fix error --------------
  for(i in names(lpdx@experiment))
  {
    dx <- lpdx@experiment[[i]]@data
    lpdx@experiment[[i]]@data <- dx[!is.na(dx$width), ]
  }
  
  ###-----------subset models --------------
  goodMod <- read.csv("~/CXP/KRAS_P53/data/Table cutoff for Arvind.csv", 
                      stringsAsFactors = F)
  goodMod$Name.model <- paste0("PHLC", goodMod$Name.model)
  mi <- modelInfo(lpdx)
  m2take <- mi[mi$patient.id %in% goodMod$Name.model,]
  
  ##-----add PHLC191 (P5) ------
  phlc191 <- mi[mi$patient.id=="PHLC191",]
  m2take <- rbind(m2take, phlc191[grep("_P5", phlc191$model.id), ])
  goodMod$Name.model[goodMod$Name.model=="PHLC191 (P5)"] <- "PHLC191"
  
  ###------subset the object -----
  #unique(modelInfo(lpdx)$patient.id)
  lpdx2 <- Xeva::subsetXeva(lpdx, ids = m2take$model.id, id.name = "model.id",
                            keep.batch = T)
  mi2 <- modelInfo(lpdx2)
  #goodMod$Name.model[!goodMod$Name.model%in%mi2$patient.id]
  #unique(mi2$patient.id)
  
  ##----subset by treatment start date ------
  for(i in names(lpdx2@experiment))
  {
    pid=mi2$patient.id[mi2$model.id==i]
    v = goodMod[goodMod$Name.model==pid, ]
    
    mintm <- v$Day.1.treatment
    if(mi2[i, "drug"]=="Control")
    { maxtm <- v$Control.cutoff[1] } else {maxtm <- v$Treatment.cutoff[1]}
    
    exp <- lpdx2@experiment[[i]]
    #exp@data <- exp@data[exp@data$time>=mintm & exp@data$time <=maxtm, ]
    ##---get the closet min time ---------------------------
    mintm2I = which.min(abs(mintm - exp@data$time))[1]
    mintm2 = exp@data$time[mintm2I]
    exp@data$time <- exp@data$time - mintm2
    lpdx2@experiment[[i]] <- exp
  }
  
  ##----------remove models with low data --------
  #expLen <- sapply(lpdx2@experiment, function(exp) nrow(exp@data))
  
  ###======= fix batch ================
  #allbat <- Xeva::batchInfo(lpdx2)
  for(i in 1:length(lpdx2@expDesign))
  {
    ed <- lpdx2@expDesign[[i]]
    ed$treatment<- ed$treatment[ed$treatment%in%modelInfo(lpdx2)$model.id]
    ed$control  <- ed$control[ed$control%in%modelInfo(lpdx2)$model.id]
    lpdx2@expDesign[[i]] <- ed
  }
  
  #mi <- modelInfo(lpdx2)
  #mi$model.id[!mi$model.id %in% lpdx2@modToBiobaseMap$model.id]
  
  ###-----normalize at t0
  #for(m in names(lpdx2@experiment))
  #{
  #  vn <- getExperiment(lpdx2, m, vol.normal = T)
  #  lpdx2@experiment[[m]]@data$volume.normal <- vn$volume
  #}
  #m="PHLC110_P5.501.A1"
  #vn <- getExperiment(lpdx2, m, vol.normal = T)
  return(lpdx2)
}


lpdx <- getFixed_LPDX_Data()
##-----------------
mi <- modelInfo(lpdx)
max.time = 60

lpdx  <- setResponse(lpdx, res.measure = c("mRECIST", "slope", "AUC"), 
                         max.time = max.time, verbose=FALSE)

lpdx  <- setResponse(lpdx, res.measure = c("angle", "abc", "TGI"), 
                         max.time = max.time, verbose=FALSE)



##---------------model_info-----------------------------------------
m = get_model_info(lpdx)
m$dataset <- "UHN Lung Cancer"; m$tissue <- "Lung Cancer"
m <- m[,c("model.id","tissue","patient.id","drug","dataset")]
write.csv(m, file = "results/UHN_Lung/model_information.csv")

##---------------batch_information----------------------------------
b = get_batch_info(lpdx)
write.csv(b, file = "results/UHN_Lung/batch_information.csv")

##-----drug_screening-----------------------------------
mdf = drug_screening(lpdx)
write.csv(mdf, file = "results/UHN_Lung/drug_screening.csv")

##----------model_response---------------------------------
mres <- model_response(lpdx)
write.csv(mres, file = "results/UHN_Lung/model_response.csv")

##----------batch_response---------------------------------
brf <- batch_response(lpdx)
write.csv(brf, file = "results/UHN_Lung/batch_response.csv")

###---------modelid_moleculardata_mapping ------------------
mmap = modelid_moleculardata_mapping(lpdx, dt=c("microarray", "target_mut"))
mmap$mDataType[mmap$mDataType=="target_mut"] = "mutation"
###------this need to be changed 
mmap$mDataType[mmap$mDataType=="microarray"] = "RNASeq"

write.csv(mmap, file = "results/UHN_Lung/modelid_moleculardata_mapping.csv")

##-----------expression ----------
df = exprs(lpdx@molecularProfiles$microarray)
df = t(scale(t(df))[,])
df = getFlatDF(df)
write.csv(df, file = "results/UHN_Lung/rna_sequencing.csv")

##-----------mutation ----------

df = exprs(lpdx@molecularProfiles$target_mut)
df = getFlatDF(df)
df$value <- as.character(df$value)
df$value[df$value=="1"] = "mutation"
write.csv(df, file = "results/UHN_Lung/mutation.csv")

