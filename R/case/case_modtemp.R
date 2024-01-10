#------------------------------------------------------------------------------#
#                      Study cases: modelling temperature                      #
#------------------------------------------------------------------------------#

# Preparation ----

library("sf")
library("caret")
library("randomForest")
library("terra")
library("CAST")
library("twosamples")

# HPC
pathroot <- "****"

# Data pre-processing
source(paste0(pathroot, "R/case/case_functions.R"))
temp_sf <- st_read(paste0(pathroot, "data/temp/temp_train.gpkg"))
temp_df <- st_drop_geometry(temp_sf)
temp_df$outcome <- temp_df$temp
aoi <- st_read(paste0(pathroot, "data/boundaries/spain.gpkg"))
preds <- rast(paste0(pathroot, "data/predictors.tif"))
RFsp <- rast(paste0(pathroot, "data/proxies/RFsp_temp.tif"))
tempstack <- c(preds, RFsp)
rm("preds", "RFsp")

# Prepare random and spatial CV indices
set.seed(123)
sfolds <- knndm(temp_sf, aoi, k = 10, maxp = 0.5)
sindx <- CreateSpacetimeFolds(data.frame(ID = sfolds$clusters), spacevar = "ID", k = 10)
sctrl <- trainControl(method="cv", index = sindx$index, savePredictions='final')
rfolds <- sample(1:nrow(temp_df)%%10, nrow(temp_df))
rindx <- CreateSpacetimeFolds(data.frame(ID = rfolds), spacevar = "ID", k = 10)
rctrl <- trainControl(method="cv", index = rindx$index, savePredictions='final')
saveRDS(sfolds, paste0(pathroot, "results/case/temp_sfolds.rds"))
saveRDS(rfolds, paste0(pathroot, "results/case/temp_rfolds.rds"))

# Naive model ----

# Baseline model
set.seed(123)
basecovs <- "dem"
proxycovs <- NULL
basemod <- fitvalpred_rf(basecovs, proxycovs, rctrl, sctrl, temp_df, tempstack[[c(basecovs, proxycovs)]])
write.csv(basemod$tab, paste0(pathroot, "results/case/temp_naive_baseline.csv"))
writeRaster(basemod$preds, paste0(pathroot, "results/case/temp_naive_baseline.tif"), overwrite=TRUE)
saveRDS(basemod$mod, paste0(pathroot, "results/case/temp_naive_baseline.rds"))
rm("basemod", "proxycovs")

# Coordinate model
set.seed(123)
proxycovs <- c("X", "Y")
coordmod <- fitvalpred_rf(basecovs, proxycovs, rctrl, sctrl, temp_df, tempstack[[c(basecovs, proxycovs)]])
write.csv(coordmod$tab, paste0(pathroot, "results/case/temp_naive_coord.csv"))
writeRaster(coordmod$preds, paste0(pathroot, "results/case/temp_naive_coord.tif"), overwrite=TRUE)
saveRDS(coordmod$mod, paste0(pathroot, "results/case/temp_naive_coord.rds"))
rm("coordmod", "proxycovs")

# EDF model
set.seed(123)
proxycovs <- c("X", "Y", paste0("EDF", 1:5))
EDFmod <- fitvalpred_rf(basecovs, proxycovs, rctrl, sctrl, temp_df, tempstack[[c(basecovs, proxycovs)]])
write.csv(EDFmod$tab, paste0(pathroot, "results/case/temp_naive_EDF.csv"))
writeRaster(EDFmod$preds, paste0(pathroot, "results/case/temp_naive_EDF.tif"), overwrite=TRUE)
saveRDS(EDFmod$mod, paste0(pathroot, "results/case/temp_naive_EDF.rds"))
rm("EDFmod", "proxycovs")

# RFsp model
set.seed(123)
proxycovs <- paste0("RFsp_temp_", 1:nrow(temp_df))
RFspmod <- fitvalpred_rf(basecovs, proxycovs, rctrl, sctrl, temp_df, tempstack[[c(basecovs, proxycovs)]])
write.csv(RFspmod$tab, paste0(pathroot, "results/case/temp_naive_RFsp.csv"))
writeRaster(RFspmod$preds, paste0(pathroot, "results/case/temp_naive_RFsp.tif"), overwrite=TRUE)
saveRDS(RFspmod$mod, paste0(pathroot, "results/case/temp_naive_RFsp.rds"))
rm("RFspmod", "proxycovs", "basecovs")

# Complete model ----

# Baseline model
set.seed(123)
basecovs <- c("dem", "coast", "imd", "ndvi", "lst_day", "lst_night")
proxycovs <- NULL
basemod <- fitvalpred_rf(basecovs, proxycovs, rctrl, sctrl, temp_df, tempstack[[c(basecovs, proxycovs)]])
write.csv(basemod$tab, paste0(pathroot, "results/case/temp_complete_baseline.csv"))
writeRaster(basemod$preds, paste0(pathroot, "results/case/temp_complete_baseline.tif"), overwrite=TRUE)
saveRDS(basemod$mod, paste0(pathroot, "results/case/temp_complete_baseline.rds"))
rm("basemod", "proxycovs")

# Coordinate model
set.seed(123)
proxycovs <- c("X", "Y")
coordmod <- fitvalpred_rf(basecovs, proxycovs, rctrl, sctrl, temp_df, tempstack[[c(basecovs, proxycovs)]])
write.csv(coordmod$tab, paste0(pathroot, "results/case/temp_complete_coord.csv"))
writeRaster(coordmod$preds, paste0(pathroot, "results/case/temp_complete_coord.tif"), overwrite=TRUE)
saveRDS(coordmod$mod, paste0(pathroot, "results/case/temp_complete_coord.rds"))
rm("coordmod", "proxycovs")

# EDF model
set.seed(123)
proxycovs <- c("X", "Y", paste0("EDF", 1:5))
EDFmod <- fitvalpred_rf(basecovs, proxycovs, rctrl, sctrl, temp_df, tempstack[[c(basecovs, proxycovs)]])
write.csv(EDFmod$tab, paste0(pathroot, "results/case/temp_complete_EDF.csv"))
writeRaster(EDFmod$preds, paste0(pathroot, "results/case/temp_complete_EDF.tif"), overwrite=TRUE)
saveRDS(EDFmod$mod, paste0(pathroot, "results/case/temp_complete_EDF.rds"))
rm("EDFmod", "proxycovs")

# RFsp model
set.seed(123)
proxycovs <- paste0("RFsp_temp_", 1:nrow(temp_df))
RFspmod <- fitvalpred_rf(basecovs, proxycovs, rctrl, sctrl, temp_df, tempstack[[c(basecovs, proxycovs)]])
write.csv(RFspmod$tab, paste0(pathroot, "results/case/temp_complete_RFsp.csv"))
writeRaster(RFspmod$preds, paste0(pathroot, "results/case/temp_complete_RFsp.tif"), overwrite=TRUE)
saveRDS(RFspmod$mod, paste0(pathroot, "results/case/temp_complete_RFsp.rds"))
rm("RFspmod", "proxycovs")
