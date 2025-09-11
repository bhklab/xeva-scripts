options(stringsAsFactors = FALSE)
library(Xeva)
source("xevaDB_fun.R")
#mpx = readRDS("~/CXP/XG/MP/Data/MP_XevaSet.rds")
mpx = readRDS("/Users/mattbocc/uhn/xeva-scripts/xevasets-obj/McGill_XevaSet_2024.rds")

max.time = 30

mpx  <- setResponse(mpx, res.measure = c("mRECIST", "slope", "AUC"), 
                     max.time = max.time, verbose=FALSE)
mpx  <- setResponse(mpx, res.measure = c("angle", "abc", "TGI"), 
                     max.time = max.time, verbose=FALSE)

##=======subset only PDX mol data ==============================================
mpx@molecularProfiles$RNAseq = 
  mpx@molecularProfiles$RNAseq[,mpx@molecularProfiles$RNAseq$Source=="PDX"]

mpx@molecularProfiles$CNV = 
  mpx@molecularProfiles$CNV[,mpx@molecularProfiles$CNV$type=="PDX"]

mpx@molecularProfiles$mutation = 
  mpx@molecularProfiles$mutation[,mpx@molecularProfiles$mutation$type=="PDX"]


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
m$dataset <- "McGill TNBC"; m$tissue <- "Breast Cancer"
m <- m[,c("model.id","tissue","patient.id","drug","dataset")]
write.csv(m, file = "../results/McGill_TNBC/model_information.csv")

##---------------batch_information----------------------------------
b = get_batch_info(mpx)
write.csv(b, file = "../results/McGill_TNBC/batch_information.csv")

##-----drug_screening-----------------------------------
mdf = drug_screening(mpx)
write.csv(mdf, file = "../results/McGill_TNBC/drug_screening.csv")

##----------model_response---------------------------------
mres <- model_response(mpx)
write.csv(mres, file = "../results/McGill_TNBC/model_response.csv")

##----------batch_response---------------------------------
brf <- batch_response(mpx)
write.csv(brf, file = "../results/McGill_TNBC/batch_response.csv")

###---------modelid_moleculardata_mapping ------------------
mmap = modelid_moleculardata_mapping(mpx, dt = c("RNAseq","mutation","CNV"))
write.csv(mmap, file = "../results/McGill_TNBC/modelid_moleculardata_mapping.csv")

##-----------expression ----------
fd = fData(mpx@molecularProfiles$RNAseq)
pc = fd[fd$gene_type == "protein_coding",]
mpx@molecularProfiles$RNAseq = mpx@molecularProfiles$RNAseq[rownames(pc),]
featureNames(mpx@molecularProfiles$RNAseq) = 
  make.names(fData(mpx@molecularProfiles$RNAseq)$gene_name, unique = T)
df=expression(mpx, dt = "RNAseq")
write.csv(df, file = "../results/McGill_TNBC/rna_sequencing.csv")

##----------- mutation ----------
##mut=exprs(mpx@molecularProfiles$mutation)
df <- getFlatDF(exprs(mpx@molecularProfiles$mutation))
#mutMap=rep("mutation", length(unique(df$value)))
#names(mutMap)=unique(df$value)
#mutMap[c("0", "Silent")]="0"

#df$value = mutMap[df$value]
write.csv(df, file = "../results/McGill_TNBC/mutation.csv")

##----------- cnv ----------
cnvCat <- c("-2"= "Deep Deletion", # indicates a deep loss, possibly a homozygous deletion
            "-1"= "Shallow Deletion",#indicates a shallow loss, possibley a heterozygous deletion
            "0" = "Diploid",
            "1" = "Gain", #indicates a low-level gain (a few additional copies, often broad)
            "2" = "Amplification" # indicate a high-level amplification (more copies, often focal)
)

df <- getFlatDF(exprs(mpx@molecularProfiles$CNV))
#df$value <- cnvCat[as.character(df$value)]
df$value[df$value=="Diploid"] <- "0"
write.csv(df, file = "../results/McGill_TNBC/copy_number_variation.csv")

