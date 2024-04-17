library(Xeva)
suppressMessages(library(Biobase))
library(reshape2)
library(data.table)
options(stringsAsFactors = FALSE) 
##---------------model_info-----------------------------------------
get_model_info <- function(pdxe)
{
  m = modelInfo(pdxe)
  m$dataset = paste0("PDXE (", m$tissue.name, ")")
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
  mrf <- data.frame()
  for(m in rownames(mr))
  {
    r=data.frame(drug=mid$drug[mid$model.id==m], model.id=m, response_type=colnames(mr),
                 stringsAsFactors = F)
    r$value = sapply(r$response_type, function(i)mr[m,i])
    mrf <- rbind(mrf, r)
  }
  return(mrf)
}


##----------batch_response---------------------------------
batch_response <- function(pdxe)
{
  br = pdxe@sensitivity$batch[, c("angle", "abc", "TGI")]
  br$angle <- getNumToString(br$angle)
  br$abc <- getNumToString(br$abc)
  br$TGI <- getNumToString(br$TGI)
  brf <- data.frame()
  for(b in rownames(br))
  {
    r=data.frame(batch.id=b, response_type=colnames(br), stringsAsFactors = F)
    r$value = sapply(r$response_type, function(i)br[b,i])
    brf <- rbind(brf, r)
  }
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
  df = t(scale(t(df))[,])
  rtx = getFlatDF(df)
  return(rtx)
}

##----------- mutation data ----------
mutation <- function(pdxe)
{
  idmap = pdxe@modToBiobaseMap
  idmap = idmap[idmap$mDataType=="mutation", ]
  
  df = exprs(pdxe@molecularProfiles$mutation)
  commanSample <- intersect(idmap$biobase.id, colnames(df))
  df = df[,commanSample]
  df <- getFlatDF(df)
  
  unqVal <- unique(df$value)
  vl <- as.list(rep("0", length(unqVal))); names(vl)=unqVal
  mutNames = names(vl)[grepl("MUT", toupper(names(vl)))]
  vl[mutNames] = "mutation"
  df$value <- unlist(vl[df$value])
  return(df)
}

##----------- mutation data ----------
cnv <- function(pdxe)
{
  # idmap = pdxe@modToBiobaseMap
  # idmap = idmap[idmap$mDataType=="cnv", ]
  # 
  # df = exprs(pdxe@molecularProfiles$cnv)
  # commanSample <- intersect(idmap$biobase.id, colnames(df))
  # df = df[,commanSample]
  # df <- getFlatDF(df)
  
  idmap = pdxe@modToBiobaseMap
  idmap = idmap[idmap$mDataType=="mutation", ]
  
  df = exprs(pdxe@molecularProfiles$mutation)
  commanSample <- intersect(idmap$biobase.id, colnames(df))
  df = df[,commanSample]
  df <- getFlatDF(df)
  
  unqVal <- unique(df$value)
  vl <- as.list(rep("0", length(unqVal))); names(vl)=unqVal
  #cnvTxt = c("Amp5", "Amp8", "Del0.8")
  for(i in names(vl))
  {
    if(grepl("Amp5|Amp8", i)==TRUE){vl[[i]]="Amplification"}
    if(grepl("Del0.8", i)==TRUE){vl[[i]]="Deletion"}
  }
  df$value <- unlist(vl[df$value])
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




