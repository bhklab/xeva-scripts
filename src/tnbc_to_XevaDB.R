options(stringsAsFactors = FALSE)
library(Xeva)
library(readxl)
source("xevaDB_fun.R")

#tnbc <- readRDS("~/CXP/XG/TNBC_Xeva_obj/Data/XevaData/XEVA_OFFICIAL_SHARING/TNBC_Xeva_Obj.Rda")
tnbc <- readRDS("/Users/mattbocc/uhn/xeva-scripts/xevasets-obj/Updated_Cescon_TNBC_Xeva_Obj_2.rds")
# tnbcOld <- readRDS("/Users/mattbocc/uhn/xeva-scripts/xevasets-obj/Xeva_TNBC_Obj.rds")

##-----for non RES -------
mi <- modelInfo(tnbc)
tnbc.nor <- subsetXeva(tnbc, ids=mi$model.id[mi$resistant=="NO"], id.name="model.id")

##------subset by drug --------------------
mi <- modelInfo(tnbc.nor)
ndr <- sort(table(mi$drug))
drug2take <- names(ndr)[ndr>3]
mid <- mi[mi$drug%in%drug2take, ]

print(colnames(tnbc.nor@sensitivity$model))

tnbc.nor <- subsetXeva(tnbc.nor, ids=mid$model.id, id.name="model.id", 
                       keep.batch = F)


##----------------------------------
mi <- modelInfo(tnbc.nor)
mi <- mi[!is.na(mi$patient.id),]

mi$patient.id <- paste0("P.", mi$patient.id)
# pid2remove=c("P.1", "P.100534", "P.2018-12-28", "P.2019-01-16", "P.44721", 
#              "P.5/18", "P.85227", "P.REF020", "P.2018-12-13", "P.5")
# mi <- mi[! mi$patient.id %in% pid2remove,]
tnbc.nor <- subsetXeva(tnbc.nor, ids=mi$model.id, id.name="model.id", 
                       keep.batch = F)
###-------------------------------
max.time = 60; cutAtMaxTime = FALSE
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
mi=mi[mi$drug%in%rownames(mr3), ]

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
m$dataset <- "TNBC"; m$tissue <- "Breast Cancer"
mo <- m[,c("model.id","tissue","patient.id","drug","dataset")]
write.csv(mo, file = "../results/TNBC/model_information.csv")

fileLoc <- m[, c("model.id", "file.url")]
fileLoc$row <- sapply(fileLoc$model.id, function(i)gsub("m", "", strsplit(i, "\\.")[[1]][2]))
write.csv(fileLoc, file = "../results/TNBC/model_information_FileLink.csv")

##---------------batch_information----------------------------------
b = get_batch_info(tnbc.nor)
mi =modelInfo(tnbc.nor); b = b[b$model.id%in%mi$model.id,]
write.csv(b, file = "../results/TNBC/batch_information.csv")

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


write.csv(mdf, file = "../results/TNBC/drug_screening.csv")

##----------model_response---------------------------------
mres <- model_response(tnbc.nor)
write.csv(mres, file = "../results/TNBC/model_response.csv")

##----------batch_response---------------------------------
brf <- batch_response(tnbc.nor)
write.csv(brf, file = "../results/TNBC/batch_response.csv")

###---------modelid_moleculardata_mapping ------------------
mmap = modelid_moleculardata_mapping(tnbc.nor)
write.csv(mmap, file = "../results/TNBC/modelid_moleculardata_mapping.csv")

##-----------expression ----------
fd = fData(tnbc.nor@molecularProfiles$RNASeq)
pc = fd[fd$gene_type == "protein_coding",]
tnbc.nor@molecularProfiles$RNASeq = tnbc.nor@molecularProfiles$RNASeq[rownames(pc),]
featureNames(tnbc.nor@molecularProfiles$RNASeq) = 
  make.names(fData(tnbc.nor@molecularProfiles$RNASeq)$gene_name, unique = T)

df=expression(tnbc.nor)
write.csv(df, file = "../results/TNBC/rna_sequencing.csv")

##-----------mutation ----------

mutInfo <- read_excel("/Users/mattbocc/uhn/xeva-scripts/data/mutation.xlsx", 
                           sheet = 1)
mutCat <- mutInfo$IMPACT; names(mutCat) <- mutInfo$SO.term
mutCat[mutCat %in% c("HIGH", "MODERATE")] <- "Mutation"
mutCat[mutCat!="Mutation"] <- ""

mut <- tnbc.nor@molecularProfiles$mutation
df = getFlatDF(exprs(mut))
df$value[is.na(df$value)] <- ""
df$value[df$value==""] <- "0"
mutCat2 <- rep("", length(unique(df$value)))
names(mutCat2)<- sort(unique(df$value))

for(i in names(mutCat2))
{
  v <- strsplit(i, ",")[[1]]
  if(length(v)>0)
  {
    if("Mutation" %in% mutCat[v]){mutCat2[i] <-"Mutation"} else
    {mutCat2[i] <-"0"}
  }
}

df$value <- mutCat2[df$value]
df$value[df$value=="0"] <- ""
###df=mutation(tnbc.nor)
write.csv(df, file = "../results/TNBC/mutation.csv")

##-----------copy_number_variation ----------

cnvCat <- c("-2"= "Deep Deletion", # indicates a deep loss, possibly a homozygous deletion
            "-1"= "Shallow Deletion",#indicates a shallow loss, possibley a heterozygous deletion
            "0" = "Diploid",
            "1" = "Gain", #indicates a low-level gain (a few additional copies, often broad)
            "2" = "Amplification" # indicate a high-level amplification (more copies, often focal)
)

cnv <- tnbc.nor@molecularProfiles$CNV
df = getFlatDF(exprs(cnv))
df$value <- as.character(df$value)
df$value <- cnvCat[df$value]
df$value[df$value=="Diploid"] <- "0"
write.csv(df, file = "../results/TNBC/copy_number_variation.csv")

####-----------------------