library(Xeva)
suppressMessages(library(Biobase))
library(reshape2)
library(data.table)
options(stringsAsFactors = FALSE) 

# Safety patch for Xeva::TGI:
# When control or treatment volume has < 1 timepoints (e.g. truncated by max.time,
# or experiments with single observations), orig_TGI calculates numeric(0), which crashes
# setResponse data.frame assignment with 'replacement has length zero'.
# This patch ensures TGI returns NA_real_ instead of numeric(0).
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
    if (length(res$value) == 0) {
      res$value <- NA_real_
    }
    return(res)
  }, ns)
  lockBinding("TGI", ns)
}, silent = TRUE) 
##---------------model_info-----------------------------------------
get_model_info <- function(pdxe)
{
  m = modelInfo(pdxe)
  m$dataset = paste0("PDXE_v2 (", m$tissue.name, ")")
  m$tissue = m$tissue.name
  m$tissue.name=NULL
  return(m)
}

##---------------batch_info-----------------------------------------
get_batch_info <- function(pdxe)
{
  bdf=data.frame()
  for(b in pdxe@expDesign)
  {
    if(length(b$control)>0)
    {
      bdf = rbind(bdf, data.frame(batch.id=b$batch.name, model.id=b$control, 
                                  type="control"))
    }
    
    if(length(b$treatment)>0)
    {
      bdf = rbind(bdf, data.frame(batch.id=b$batch.name, model.id=b$treatment, 
                                  type="treatment"))
    }
  }
  return(bdf)
}

##-----drug_screening-----------------------------------
drug_screening <- function(pdxe, col2take=c("model.id", "drug.join.name", "time",
                                            "volume", "volume.normal"))
{
  mid <- modelInfo(pdxe)
  mdfl<- list()
  #col2take=c("model.id", "drug.join.name", "time","volume", "volume.normal")
  for(i in 1:nrow(mid))
  {
    mx <- getExperiment(pdxe, mid$model.id[i])
    mdfl[[length(mdfl)+1]] <- mx
    #mdf <- rbind(mdf, mx[, col2take])
  }
  
  #colnames(mdf)[colnames(mdf)=="drug.join.name"] <- "drug"
  allcl <- unique(unlist(sapply(mdfl, colnames)))
  
  for(i in 1:length(mdfl))
  {
    ncl <- allcl[!allcl%in%colnames(mdfl[[i]])]
    if(length(ncl)>0){ mdfl[[i]][, ncl] <- NA}
    mdfl[[i]] <- mdfl[[i]][, allcl]
  }

  # Was previously used but has now been changed by Matthew for efficiency sake, looking for same functionality rbindlist
  # mdf <- do.call(rbind.data.frame, mdfl)
  mdf <- rbindlist(mdfl, fill = TRUE)
  mdf <- mdf[, ..col2take]
  
  colnames(mdf)[colnames(mdf)=="drug.join.name"] <- "drug"

  return(mdf)
}

getNumToString <- function(v)
{ sapply(v, function(i)sprintf("%1.4f",i)) }


##----------model_response---------------------------------
model_response <- function(pdxe)
{
  mr = pdxe@sensitivity$model[, c("mRECIST", "best.average.response", "slope", "AUC")]
  mr$survival <- sapply(rownames(pdxe@sensitivity$model), function(i){max(pdxe@experiment[[i]]@data$time)})
  
  mr$best.average.response <- getNumToString(mr$best.average.response)
  mr$slope <- getNumToString(mr$slope)
  mr$AUC <- getNumToString(mr$AUC)
  mr$survival <- getNumToString(mr$survival)
  
  mid <- modelInfo(pdxe)
  drug_map <- setNames(mid$drug, mid$model.id)
  mr$model.id <- rownames(mr)
  mr$drug <- drug_map[mr$model.id]
  
  mrf <- reshape2::melt(mr, id.vars = c("drug", "model.id"),
                        variable.name = "response_type", value.name = "value")
  mrf$response_type <- as.character(mrf$response_type)
  mrf <- mrf[, c("drug", "model.id", "response_type", "value")]
  return(mrf)
}


##----------batch_response---------------------------------
batch_response <- function(pdxe)
{
  br = pdxe@sensitivity$batch[, c("angle", "abc", "TGI")]
  br$angle <- getNumToString(br$angle)
  br$abc <- getNumToString(br$abc)
  br$TGI <- getNumToString(br$TGI)
  br$batch.id <- rownames(br)
  
  brf <- reshape2::melt(br, id.vars = "batch.id",
                        variable.name = "response_type", value.name = "value")
  brf$response_type <- as.character(brf$response_type)
  brf <- brf[, c("batch.id", "response_type", "value")]
  return(brf)
}

###------modelid_moleculardata_mapping-----------------------
modelid_moleculardata_mapping <- function(pdxe, dt=c("RNASeq","mutation", "cnv"))
{
  df = pdxe@modToBiobaseMap
  df = df[df$mDataType%in%dt, ]
  return(df)
}
##----------- expression data ----------
getFlatDF <- function(df)
{
  rtx = reshape2::melt(df, varnames = c("gene.id", "sequencing.uid"),
                       value.name="value")
  rtx$gene.id = as.character(rtx$gene.id)
  rtx$sequencing.uid = as.character(rtx$sequencing.uid)
  return(rtx)
}

expression <- function(pdxe, dt="RNASeq")
{
  idmap = pdxe@modToBiobaseMap
  idmap = idmap[idmap$mDataType==dt, ]
  
  df = exprs(pdxe@molecularProfiles[[dt]])
  commanSample <- intersect(idmap$biobase.id, colnames(df))
  df = df[,commanSample]
  scaled <- t(scale(t(df))[,])
  # Genes with zero variance across all samples produce NaN when dividing by sd=0; set to 0
  scaled[is.nan(scaled)] <- 0
  rtx = getFlatDF(scaled)
  return(rtx)
}

##----------- mutation data ----------
mutation <- function(pdxe)
{
  idmap = pdxe@modToBiobaseMap
  idmap = idmap[idmap$mDataType=="mutation", ]
  
  df = exprs(pdxe@molecularProfiles$mutation)
  commanSample <- intersect(idmap$biobase.id, colnames(df))
  if (length(commanSample) > 0) {
    df = df[, commanSample, drop = FALSE]
  }
  df <- getFlatDF(df)
  
  df$value <- as.character(df$value)
  df$value[is.na(df$value) | df$value == ""] <- "0"
  df$value[df$value %in% c("1", 1) | grepl("MUT", toupper(df$value))] <- "mutation"
  return(df)
}

##----------- cnv data ----------
cnv <- function(pdxe)
{
  idmap = pdxe@modToBiobaseMap
  idmap = idmap[idmap$mDataType=="cnv", ]
  
  df = exprs(pdxe@molecularProfiles$cnv)
  df = df[!rownames(df) %in% c("ArmLevelCNScore", "FocalCNScore"), , drop = FALSE]
  commanSample <- intersect(idmap$biobase.id, colnames(df))
  if (length(commanSample) > 0) {
    df = df[, commanSample, drop = FALSE]
  }
  df <- getFlatDF(df)

  # PDXE CNV Classification:
  # Amp5 / Amp8 (copy number > 5)  -> "Amplification"
  # Del0.8 (copy number <= 0.8)     -> "Deletion"
  # Normal / diploid range          -> "0"
  val <- as.numeric(df$value)
  cat_val <- rep("0", length(val))
  cat_val[!is.na(val) & val > 5] <- "Amplification"
  cat_val[!is.na(val) & val <= 0.8] <- "Deletion"
  df$value <- cat_val

  return(df)
}




#####++++++++++++++++++++++++++++++++++
cutDataAtMaxTime <- function(pdxe, max.time)
{
  for(i in names(pdxe@experiment))
  {
    d = pdxe@experiment[[i]]@data
    pdxe@experiment[[i]]@data = d[d$time <= max.time, ]
  }
  return(pdxe)
}




