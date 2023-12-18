#------------------------------------------------------------------------------#
#                     Study cases: RFsp predictors (HPC)                       #
#------------------------------------------------------------------------------#

library("sf")
library("stars")

pathroot <- "****"

# boundary and grid objects
grid_points <- st_read(paste0(pathroot, "boundaries/grid_points.gpkg"))
grid_raster <- read_stars(paste0(pathroot, "boundaries/grid_raster.tif"))
spain <- st_read(paste0(pathroot, "boundaries/spain.gpkg"))
spain500 <- st_buffer(spain, 500)

# temp
tempdata <- st_read(paste0(pathroot, "temp/tempdata.gpkg"))
RFsp_temp <- st_distance(grid_points, tempdata)/1000
RFsp_temp <- units::drop_units(RFsp_temp)
RFsp_temp <- as.data.frame(RFsp_temp)
names(RFsp_temp) <- paste0("RFsp_temp_", 1:nrow(tempdata))
RFsp_temp <- st_sf(cbind(RFsp_temp, grid_points))
RFsp_temp_stack <- c( # split to avoid stars error
  st_rasterize(RFsp_temp[,1:100], grid_raster),
  st_rasterize(RFsp_temp[,101:nrow(tempdata)], grid_raster)
)
RFsp_temp_stack <- st_crop(RFsp_temp_stack, spain500)
write_stars(merge(RFsp_temp_stack), paste0(pathroot, "proxies/RFsp_temp.tif"))
rm("RFsp_temp", "tempdata", "RFsp_temp_stack")

# pm25
pm25data <- st_read(paste0(pathroot, "AP/PM25_stations.gpkg"))
RFsp_pm25 <- st_distance(grid_points, pm25data)/1000
RFsp_pm25 <- units::drop_units(RFsp_pm25)
RFsp_pm25 <- as.data.frame(RFsp_pm25)
names(RFsp_pm25) <- paste0("RFsp_pm25_", 1:nrow(pm25data))
RFsp_pm25 <- st_sf(cbind(RFsp_pm25, grid_points))
RFsp_pm25_stack <- c( # split to avoid stars error
  st_rasterize(RFsp_pm25[,1:100], grid_raster),
  st_rasterize(RFsp_pm25[,101:nrow(pm25data)], grid_raster)
)
RFsp_pm25_stack <- st_crop(RFsp_pm25_stack, spain500)
write_stars(merge(RFsp_pm25_stack), paste0(pathroot, "proxies/RFsp_pm25.tif"))
rm("RFsp_pm25", "pm25data", "RFsp_pm25_stack")
