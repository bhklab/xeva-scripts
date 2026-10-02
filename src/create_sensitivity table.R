library(Xeva)
pdxe <- readRDS("~/CXP/Xeva_dataset/data/XevaObjects/XevaSets/data/PDXE.rds")
mi <- modelInfo(pdxe)
res <- c(
    "best.response_published", "time.best.response_published",
    "best.avg.response_published", "time.best.avg.response_published",
    "timeToDouble_published", "time.last_published", "mRECIST_published",
    "mRECIST", "best.response", "best.response.time", "best.average.response",
    "best.average.response.time", "slope", "AUC"
)
res <- c("mRECIST", "best.average.response", "slope", "AUC")
rt <- data.frame()
for (i in pdxe@sensitivity$model$model.id)
{
    for (type in res)
    {
        value <- pdxe@sensitivity$model[i, type]
        if (class(value) == "numeric") {
            value <- sprintf("%1.4f", value)
        }
        rt <- rbind(rt, data.frame(
            drug = mi[i, "drug"], model.id = i, type = type,
            value = value
        ))
    }
}
rownames(rt) <- NULL
write.csv(rt, file = "../data/model_sensitivity.csv")

#### --------for batch -----------

# bres=c("slope.control", "slope.treatment", "angle", "auc.control", "auc.treatment", "abc")
bres <- c("angle", "abc")
bt <- data.frame()
for (i in pdxe@sensitivity$batch$batch.name)
{
    for (type in bres)
    {
        value <- pdxe@sensitivity$batch[i, type]
        if (class(value) == "numeric") {
            value <- sprintf("%1.4f", value)
        }
        bt <- rbind(bt, data.frame(batch.name = i, type = type, value = value))
    }
}

rownames(bt) <- NULL
write.csv(bt, file = "../data/batch_sensitivity.csv")
