#------------------------------------------------------------------------------#
#                          Study cases: extract data                           #
#------------------------------------------------------------------------------#

library("terra")

pathroot <- "****"

# Common predictors ----
spain_preds <- rast(paste0(pathroot, "predictors.tif"))

# Temperature ----
temp_stations <- vect(paste0(pathroot, "temp/tempdata.gpkg"))
temp_RFsp <- rast(paste0(pathroot, "proxies/RFsp_temp.tif"))
temp_stack <- c(spain_preds, temp_RFsp)
temp_extr <- subset(extract(temp_stack, temp_stations), select = -c(ID))
temp_stations <- cbind(temp_stations, temp_extr)
writeVector(temp_stations, paste0(pathroot, "temp/temp_train.gpkg"), filetype = "GPKG")
rm("temp_stations", "temp_RFsp", "temp_stack", "temp_extr")


# PM25 ----
pm25_stations <- vect(paste0(pathroot, "AP/PM25_stations.gpkg"))
pm25_RFsp <- rast(paste0(pathroot, "proxies/RFsp_pm25.tif"))
pm25_stack <- c(spain_preds, pm25_RFsp)
pm25_extr <- subset(extract(pm25_stack, pm25_stations), select = -c(ID))
pm25_stations <- cbind(pm25_stations, pm25_extr)
writeVector(pm25_stations, paste0(pathroot, "AP/PM25_train.gpkg"), filetype = "GPKG")
rm("pm25_stations", "pm25_RFsp", "pm25_stack", "pm25_extr")
