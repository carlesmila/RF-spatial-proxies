#------------------------------------------------------------------------------#
#                          Study cases: prepare data                           #
#------------------------------------------------------------------------------#

library("sf")
library("terra")
library("raster")
library("stars")
library("starsExtra")
library("tidyverse")
library("climaemet")
library("rnaturalearth")
library("readxl")

# Boundaries ----
spain <- ne_states(country = 'spain',returnclass = 'sf') |>
  st_transform(crs = 23030) # EPSG:23030 - ED50 / UTM zone 30N
spain <- spain[!spain$name %in% c("Ceuta", "Melilla", "Santa Cruz de Tenerife",
                                  "Las Palmas", "Baleares"),] # mainland
spain <- st_union(spain)
spain500 <- st_buffer(spain, 500)
# st_write(spain, "data/boundaries/spain.gpkg")

# Prediction grid ----
grid_raster <- st_as_stars(spain500, dx = 1000, dy = 1000)  # 1km grid
# write_stars(grid_raster, "data/boundaries/grid_raster.tif")
grid_polys <- st_as_sf(grid_raster)
grid_polys$values <- NULL
grid_polys <- grid_polys[c(st_intersects(grid_polys, spain, sparse = F)),]
# st_write(grid_polys, "data/boundaries/grid_polys.gpkg")
grid_points <- st_centroid(grid_polys)
# st_write(grid_points, "data/boundaries/grid_points.gpkg")

# Temp stations ----
aemet_api_key("****")
tempstations <- aemet_daily_clim(start = as.Date("2019-01-01"),
                                 end = as.Date("2019-12-31"), return_sf = T)
tempstations <- rename(tempstations, ID=indicativo)
templocs <- tempstations["ID"]
templocs <- templocs[!duplicated(templocs$ID),] |>
  st_transform(crs = st_crs(spain)) |>
  st_intersection(spain)
tempreads <- tempstations |>
  st_drop_geometry() |>
  filter(!is.na(tmred)) |>
  group_by(ID) |>
  summarise(temp = mean(tmed),
            Ndays = n()) |>
  filter(Ndays >= 365*0.75) |> # Min 75% temporal coverage
  dplyr::select(-Ndays)
tempdata <- inner_join(tempreads, templocs,  by = "ID") |>
  st_as_sf(crs = st_crs(spain))
tempdata_dup <- duplicated(st_coordinates(tempdata))
tempdata <- tempdata[!tempdata_dup,]
# st_write(tempdata, "data/temp/tempdata.gpkg")
rm("tempdata", "tempstations", "templocs", "tempreads")

# PM25 stations ----
pm25stations <- read_excel("data/AP/raw/metainformacion_2019_tcm30-513561.xlsx", 2) |>
  dplyr::select(COD_LOCAL, PROVINCIA, MUNICIPIO, ESTACION, LONGITUD_G, LATITUD_G) |>
  rename(ID = COD_LOCAL) |>
  st_as_sf(coords = c("LONGITUD_G", "LATITUD_G"), crs = 4326) |>
  st_transform(crs = st_crs(spain))
pm25stations <- st_intersection(pm25stations, spain)
pm25data <- read.csv("data/AP/raw/PM25_HH_2019.csv", sep = ";", dec = ".") |>
  pivot_longer(cols = starts_with("H"), names_to = "hour", values_to = "PM25") |>
  mutate(PM25 = as.numeric(PM25)) |>
  filter(!is.na(PM25))
pm25data_daily1 <- pm25data |> # daily averages
  group_by(PROVINCIA, MUNICIPIO, ESTACION, ANNO, MES, DIA) |>
  summarise(PM25 = mean(PM25),
            Nday = n()) |>
  filter(Nday >= 24*0.75) |> # 75% availability
  ungroup()
pm25data_yearly1 <- pm25data_daily1 |> # yearly averages
  group_by(PROVINCIA, MUNICIPIO, ESTACION) |>
  summarise(PM25 = mean(PM25),
            Nyear = n()) |>
  filter(Nyear >= 365*0.75) |> # 75% availability
  ungroup() |>
  inner_join(pm25stations, by = join_by(PROVINCIA, MUNICIPIO, ESTACION))
pm25data_daily2 <- read.csv("data/AP/raw/PM25_DD_2019.csv", sep = ";", dec = ".") |>
  pivot_longer(cols = starts_with("D"), names_to = "hour", values_to = "PM25") |>
  mutate(PM25 = as.numeric(PM25)) |>
  filter(!is.na(PM25))
pm25data_yearly2 <- pm25data_daily2 |> # yearly averages
  group_by(PROVINCIA, MUNICIPIO, ESTACION) |>
  summarise(PM25 = mean(PM25),
            Nyear = n()) |>
  filter(Nyear >= 365*0.75) |> # 75% availability
  ungroup() |>
  inner_join(pm25stations, by = join_by(PROVINCIA, MUNICIPIO, ESTACION))
pm25data <- rbind(pm25data_yearly1, pm25data_yearly2) |>
  st_as_sf() |>
  arrange(ID, -Nyear)
pm25data <- pm25data[!duplicated(pm25data$ID),] |>
  dplyr::select(ID, PM25)
# st_write(pm25data, "data/AP/PM25_stations.gpkg")
rm("pm25stations", "pm25data_daily1", "pm25data_yearly1",
   "pm25data_daily2", "pm25data_yearly2", "pm25data")

# Population density ----
popdens <- read_stars("data/population/JRC_1K_POP_2018.tif")
names(popdens) <- "popdens"
popdens <- mutate(popdens, popdens = ifelse(is.na(popdens), 0, popdens))
popdens <- st_warp(popdens, grid_raster, method = "bilinear", use_gdal = T,
                   no_data_value = -9999)
popdens <- st_crop(popdens, spain500)
names(popdens) <- "popdens"
# write_stars(popdens, "data/population/popdens.tif")
rm("popdens")

# Distance to coast ----
grid_points <- st_read("data/boundaries/grid_points.gpkg")
hydro <- st_read("data/coast/EUHYDRO_Coastline_EEA39_v013.shp")
hydro <- st_crop(hydro, st_transform(st_buffer(spain, 200000), crs = st_crs(hydro)))
hydro <- st_union(hydro) |>
  st_transform(crs = st_crs(grid_points))
hydrodist <- st_distance(grid_points, hydro)
hydrodist <- units::drop_units(c(hydrodist))/1000
grid_points$coast <- hydrodist
coast <- st_rasterize(grid_points["coast"], grid_raster)
coast <- st_crop(coast, spain500)
# write_stars(coast, "data/coast/coast.tif")
rm("grid_points", "hydro", "hydrodist", "coast")

# DEM ----
rasterio <- list(nBufXSize = 4000, nBufYSize = 4000, resample = "average")
dem1 <- read_stars("data/DEM/eu_dem_v11_E20N10.TIF", RasterIO = rasterio, proxy = F)
dem2 <- read_stars("data/DEM/eu_dem_v11_E20N20.TIF", RasterIO = rasterio, proxy = F)
dem3 <- read_stars("data/DEM/eu_dem_v11_E30N10.TIF", RasterIO = rasterio, proxy = F)
dem4 <- read_stars("data/DEM/eu_dem_v11_E30N20.TIF", RasterIO = rasterio, proxy = F)
dem <- st_mosaic(dem1, dem2, dem3, dem4)
rm("dem1", "dem2", "dem3", "dem4")
dem <- st_warp(dem, grid_raster, method = "average", use_gdal = T, no_data_value = -9999)
dem <- st_crop(dem, spain500)
# write_stars(dem, "data/DEM/dem.tif")
rm("dem")

# IMD ----
imd <- raster("data/IMD/IMD_2018_100m_eu_03035_V2_0.tif")
imd <- raster::crop(imd, st_sf(st_transform(st_buffer(spain, 5000), crs = st_crs(imd))))
imd <- reclassify(imd, cbind(255,255,NA), include.lowest=T, right=T)
imd <- st_as_stars(imd)
imd <- st_warp(imd, grid_raster, method = "average", use_gdal = T, no_data_value = -9999)
names(imd) <- "imd"
imd <- st_crop(imd, spain500)
# write_stars(imd, "data/IMD/imd.tif")
rm("imd")

# NTL ----
rasterio <- list(nXOff = 40000, nYOff = 7000, nXSize = 5000, nYSize = 3000)
ntl <- read_stars("data/NTL/VNL_v21_npp_2019_global_vcmslcfg_c202205302300.median_masked.dat.tif",
                   RasterIO = rasterio, proxy = F)
ntl <- st_warp(ntl, grid_raster, method = "average", use_gdal = T,
               no_data_value = -9999)
names(ntl) <- "ntl"
ntl <- st_crop(ntl, spain500)
# write_stars(ntl, "data/NTL/ntl.tif")
rm("ntl")

# NDVI ----
ndvi <- read_stars("data/NDVI/MYD13A1.006_500m_aid0001.nc", proxy = F)
ndvi <- ndvi[1]
ndvi <- st_apply(ndvi, 1:2, median, na.rm = T)
ndvi[is.na(ndvi)] <- min(ndvi$median, na.rm=T)  # Assign minimum ndvi to sea
ndvi <- st_warp(ndvi, grid_raster, method = "average", use_gdal = T, no_data_value = -9999)
names(ndvi) <-"ndvi"
ndvi <- st_crop(ndvi, spain500)
# write_stars(ndvi, "data/NDVI/ndvi.tif")
rm("ndvi")

# Roads ----

# Filter raw files per region
allfiles <- list.files("data/roads", pattern = "*.shp", all.files = TRUE,
                       recursive = TRUE, full.names = TRUE)
allfiles <- allfiles[!grepl(".zip", allfiles, fixed= TRUE)]
allroads <- NULL
for(f in allfiles){
  rfile <- read_sf(f) |>
    dplyr::select(fclass)
  rfile <- rfile[rfile$fclass %in%
                   c("motorway", "motorway_link", "primary",
                     "primary_link", "trunk", "trunk_link",
                     "secondary", "secondary_link", "tertiary",
                     "tertiary_link", "living_street", "unclassified",
                     "service", "residential", "unclassified"),]
  allroads <- rbind(allroads, rfile)
  rm("rfile")
}
# st_write(allroads, "data/roads/allroads.gpkg")
rm("f", "allroads")

# Compute road density
roads <- read_sf("data/roads/allroads.gpkg") |>
  st_transform(crs = st_crs(spain))
primary <- roads[roads$fclass %in% c("motorway", "motorway_link", "primary",
                                     "primary_link", "trunk", "trunk_link"),]
others <- roads[roads$fclass %in% c("secondary", "secondary_link", "tertiary",
                                    "tertiary_link", "living_street", "unclassified",
                                    "service", "residential", "unclassified"),]
rm("roads")
road_dens <- function(road_lines, gridpoly){

  # Compute intersection and road length per cell
  road_lines_inter <- st_intersection(road_lines, gridpoly)
  road_lines_inter$rdens <- as.numeric(st_length(road_lines_inter))
  road_lines_inter <- st_drop_geometry(road_lines_inter) %>%
    group_by(cellID) %>%
    summarise(rdens=sum(rdens))

  # Drop geometries and merge
  road_gridpoly <- st_drop_geometry(gridpoly)
  road_gridpoly <- left_join(road_gridpoly, road_lines_inter, by="cellID")

  # Fill with 0s if not covered by any road
  road_gridpoly$rdens <- ifelse(is.na(road_gridpoly$rdens), 0, road_gridpoly$rdens)
  road_gridpoly$rdens
}
grid_roads <- st_read("data/boundaries/grid_polys.gpkg")
grid_raster <- read_stars("data/boundaries/grid_raster.tif")
grid_roads$cellID <- 1:nrow(grid_roads)
grid_roads$primaryroads <- road_dens(primary, grid_roads)/1000
grid_roads$otherroads <- road_dens(others, grid_roads)/1000
primaryroads <- st_rasterize(grid_roads[c("primaryroads")], grid_raster)
primaryroads <- st_crop(primaryroads, spain500)
otherroads <- st_rasterize(grid_roads[c("otherroads")], grid_raster)
otherroads <- st_crop(otherroads, spain500)
# write_stars(primaryroads, "data/roads/primaryroads.tif")
# write_stars(otherroads, "data/roads/otherroads.tif")
rm("grid_roads", "primaryroads", "otherroads", "primary", "others")

# LULC ----
lulc_aux <- function(lulcr, gridr, codes_lulc, bound){
  lulc_bin <- lulcr
  lulc_bin <- lulc_bin%in%codes_lulc
  lulc_bin <- st_as_stars(lulc_bin, point=FALSE)
  lulc_bin <- st_warp(lulc_bin, gridr, method="average", use_gdal=TRUE,
                      no_data_value = -9999)
  lulc_bin <- st_crop(lulc_bin, bound)
  lulc_bin
}
lulc <- raster("data/LULC/U2018_CLC2018_V2020_20u1.tif")
lulc <- crop(lulc, st_sf(st_transform(st_buffer(spain, 1000), crs = st_crs(lulc))))
lulc <- as.integer(lulc)
lulc[lulc==128] <- 44 # Sea
urban <- lulc_aux(lulc, grid_raster, c(1:2, 9, 11), spain500)
# write_stars(urban, "data/LULC/urban.tif"); rm("urban")
industry <- lulc_aux(lulc, grid_raster, 3:8, spain500)
# write_stars(industry, "data/LULC/industry.tif"); rm("industry")
agriculture <- lulc_aux(lulc, grid_raster, 12:21, spain500)
# write_stars(agriculture, "data/LULC/agriculture.tif"); rm("agriculture")
natural <- lulc_aux(lulc, grid_raster, c(10, 22:44), spain500)
# write_stars(natural, "data/LULC/natural.tif"); rm("natural")
rm("lulc", "urban", "industry", "agriculture", "natural")


# CAMS PM2.5 ----
pm25_list <- list.files("data/CAMS", pattern = ".nc", recursive = TRUE,
                       full.names = TRUE, all.files = TRUE)
pm25_list <- pm25_list[grepl("pm2p5", pm25_list)]
pm25 <- read_stars(pm25_list, proxy = FALSE,
                  RasterIO = list(nXOff = 140, nYOff = 270, nXSize = 150, nYSize = 100))
st_crs(pm25) <- 4326
pm25 <- st_apply(pm25, 1:2, mean)
pm25 <- st_warp(pm25, grid_raster, method = "bilinear", use_gdal = T, no_data_value = -9999)
names(pm25) <- "CAMSrean_pm25"
# write_stars(pm25, "data/CAMS/CAMSpm25.tif")
rm("pm25_list", "pm25")

# LST ----
lst_day <- st_warp(read_stars("data/LST/LST_day.tif"), grid_raster,
                   method = "bilinear", use_gdal = T, no_data_value = -9999)
focalr <- focal2(lst_day, matrix(1, 5, 5), "mean", na.rm = T) # Fill coastal
valuesr <- ifelse(is.na(lst_day[[1]]), focalr[[1]], lst_day[[1]])
lst_day[[1]] <- valuesr

lst_night <- st_warp(read_stars("data/LST/LST_night.tif"), grid_raster,
                   method = "bilinear", use_gdal = T, no_data_value = -9999)
focalr <- focal2(lst_night, matrix(1, 5, 5), "mean", na.rm = T) # Fill coastal
valuesr <- ifelse(is.na(lst_night[[1]]), focalr[[1]], lst_night[[1]])
lst_night[[1]] <- valuesr

lst <- c(lst_day, lst_night)
names(lst) <- c("lst_day", "lst_night")
lst <- st_crop(lst, spain500)
# write_stars(merge(lst), "data/LST/lst.tif")
rm("lst_day", "lst_night", "lst", "focalr", "valuesr")

# Coordinates ----
coords <- as.data.frame(st_coordinates(grid_points))
coords <- st_sf(cbind(coords, grid_points))
coords_stack <- st_rasterize(coords, grid_raster)
coords_stack <- st_crop(coords_stack, spain500)
# write_stars(merge(coords_stack), "data/proxies/coords.tif")
rm("coords_stack")

# EDF ----
EDF <- rbind(st_sf(geom = st_sfc(st_point(c(min(coords$X),min(coords$Y))))),
             st_sf(geom = st_sfc(st_point(c(min(coords$X),max(coords$Y))))),
             st_sf(geom = st_sfc(st_point(c(max(coords$X),min(coords$Y))))),
             st_sf(geom = st_sfc(st_point(c(max(coords$X),max(coords$Y))))),
             st_sf(geom = st_sfc(st_point(c(median(coords$X),median(coords$Y))))))
EDF <- st_set_crs(EDF, st_crs(grid_points))
EDF <- st_distance(grid_points, EDF)/1000
EDF <- units::drop_units(EDF)
EDF <- as.data.frame(EDF)
names(EDF) <- paste0("EDF", 1:5)
EDF <- st_sf(cbind(EDF, grid_points))
EDF_stack <- st_rasterize(EDF, grid_raster)
EDF_stack <- st_crop(EDF_stack, spain500)
# write_stars(merge(EDF_stack), "data/proxies/EDF.tif")
rm("EDF_stack")

# Stack ----
rstack <- c(
  rast("data/population/popdens.tif"),
  rast("data/coast/coast.tif"),
  rast("data/DEM/dem.tif"),
  rast("data/IMD/imd.tif"),
  rast("data/NTL/ntl.tif"),
  rast("data/NDVI/ndvi.tif"),
  rast("data/roads/primaryroads.tif"),
  rast("data/roads/otherroads.tif"),
  rast("data/LULC/urban.tif"),
  rast("data/LULC/industry.tif"),
  rast("data/LULC/agriculture.tif"),
  rast("data/LULC/natural.tif"),
  rast("data/CAMS/CAMSpm25.tif"),
  rast("data/LST/lst.tif"),
  rast("data/proxies/coords.tif"),
  rast("data/proxies/EDF.tif")
)

# Write ----
# writeRaster(rstack, "data/predictors.tif", overwrite=TRUE)
