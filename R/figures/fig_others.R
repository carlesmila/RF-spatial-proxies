#------------------------------------------------------------------------------#
#                               Extra figures                                  #
#------------------------------------------------------------------------------#

library("tidyverse")
library("gstat")
library("raster")
library("stars")
library("sf")
library("tmap")
library("colorspace")
source("R/simulations/sim_utils.R")

# Random field simulation ----

# Create an empty results object
res <- data.frame()

# Create grids (raster and point format) and sampling area
rast_grid <- raster(ncols=100, nrows=100, xmn=0, xmx=100, ymn=0, ymx=100)
point_grid <- st_as_sf(rasterToPoints(rast_grid, spatial = TRUE))
study_area <- matrix(c(0,0,100,0,100,100,0,100,0,0), ncol=2, byrow=TRUE)
study_area <- st_sfc(st_polygon(list(study_area)))
set.seed(1234)

# Simulate random field: range 10
cov_mod <- vgm(model="Sph", psill=1, range=10, nugget=0)
cov_mod <- gstat(formula=z~1, dummy=TRUE, beta=0, model=cov_mod, nmax=100)
cov_points <- predict(cov_mod, point_grid, nsim=1)
cov10 <- rasterise_stack(cov_points, 1, "Range 10")

# Simulate random field: range 40
cov_mod <- vgm(model="Sph", psill=1, range=40, nugget=0)
cov_mod <- gstat(formula=z~1, dummy=TRUE, beta=0, model=cov_mod, nmax=100)
cov_points <- predict(cov_mod, point_grid, nsim=1)
cov40 <- rasterise_stack(cov_points, 1, "Range 40")

# Simulate random field: independent
rnoise <- raster(ncols=100, nrows=100, xmn=0, xmx=100, ymn=0, ymx=100)
vals <- rnorm(100*100, sd=1)
rnoise <- setValues(rnoise, vals)
names(rnoise) <- "Random noise"

# Plot stack
rstack <- st_as_stars(stack(cov10, cov40, rnoise))
rstack <- st_set_dimensions(rstack, "band", values = c(
  "Range = 10", "Range = 40", "Random noise"))
p1 <- ggplot() +
  geom_stars(data = rstack) +
  facet_wrap(~band, ncol = 3) +
  scale_fill_continuous_diverging(palette = "Blue-Red 3", mid = 0,
                                  breaks = seq(-3, 3, 1)) +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "bottom", legend.key.width = unit(2, 'cm'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank()) +
  labs(fill = "") +
  coord_equal()
ggsave("figures/sims_examplesrange.png", p1, width = 8, height = 6)


# Stations measurements ----

# Read data
boundaries <- read_sf("data/boundaries/spain.gpkg")
tempdata <- read_sf("data/temp/temp_train.gpkg")
pm25data <- read_sf("data/AP/PM25_train.gpkg")

# Generate maps
tempmap <- tm_shape(boundaries) +
  tm_borders() +
  tm_shape(tempdata) +
  tm_dots(col = "temp", size = 0.1,
          palette = hcl.colors(5, "Blue-Red"),
          style = "equal",
          title = "Temp. (ºC)", midpoint = NA,
          legend.format = list(digits=1)) +
  tm_graticules(alpha = 0.2)
pm25map <- tm_shape(boundaries) +
  tm_borders() +
  tm_shape(pm25data) +
  tm_dots(col = "PM25", size = 0.1,
          palette = hcl.colors(5, "Cividis"),
          style = "equal",
          title = "PM\U02082.\U02085 (\u03bcg/m\U000b3)", midpoint = NA,
          legend.format = list(digits=1)) +
  tm_graticules(alpha = 0.2)

# Join and export
allmaps <- tmap_arrange(tempmap, pm25map, nrow = 1)
tmap_save(allmaps, "figures/case_maps.png", width = 9, height = 4)
