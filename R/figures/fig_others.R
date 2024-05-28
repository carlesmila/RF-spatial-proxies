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

set.seed(1)

# Create an empty results object
res <- data.frame()

# Create grids (raster and point format) and sampling area
rast_grid <- raster(ncols=300, nrows=100, xmn=0, xmx=300, ymn=0, ymx=100)
point_grid <- st_as_sf(rasterToPoints(rast_grid, spatial = TRUE))

# Simulate covariate field: range 10
cov_mod <- vgm(model="Sph", psill=1, range=10, nugget=0)
cov_mod <- gstat(formula=z~1, dummy=TRUE, beta=0, model=cov_mod, nmax=100)
cov_points <- predict(cov_mod, point_grid, nsim=1)
cov10 <- rasterise_stack(cov_points, 1, "Range 10")

# Simulate covariate field: range 40
cov_mod <- vgm(model="Sph", psill=1, range=40, nugget=0)
cov_mod <- gstat(formula=z~1, dummy=TRUE, beta=0, model=cov_mod, nmax=100)
cov_points <- predict(cov_mod, point_grid, nsim=1)
cov40 <- rasterise_stack(cov_points, 1, "Range 40")

# Simulate error field: random
rnoise <- raster(ncols=300, nrows=100, xmn=0, xmx=300, ymn=0, ymx=100)
vals <- rnorm(300*100, sd=1)
rnoise <- setValues(rnoise, vals)
names(rnoise) <- "Random error"

# Simulate error field: autocorrelated
error_mod <- vgm(model="Sph", psill=1, range=25, nugget=0)
error_mod <- gstat(formula=z~1, dummy=TRUE, beta=0, model=error_mod, nmax=100)
error_points <- quiet(predict(error_mod, point_grid, nsim=1))
snoise <- rasterise_stack(error_points, 1, "Spatial error")

# Plot stack
rstack <- st_as_stars(stack(cov10, cov40, snoise, rnoise))
rstack <- st_set_dimensions(rstack, "band", values = c(
  "Predictor with range = 10", "Predictor with range = 40",
  "Autocorrelated error with range = 25", "Random error"))
p1 <- ggplot() +
  geom_stars(data = rstack) +
  facet_wrap(~band, ncol = 1) +
  colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3", mid = 0,
                                  breaks = seq(-3, 3, 1)) +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.key.height = unit(2, 'cm'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank()) +
  labs(fill = "") +
  coord_equal()
# ggsave("figures/sims_examplesrange.pdf", p1, width = 6, height = 8)


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
          title = expression(PM[2.5]~(mu*g/m^3)), midpoint = NA,
          legend.format = list(digits=1)) +
  tm_graticules(alpha = 0.2) +
  tm_layout(legend.position = c(0.73, 0.03))

# Join and export
allmaps <- tmap_arrange(tempmap, pm25map, nrow = 1)
tmap_save(allmaps, "figures/case_maps.pdf", width = 8.5, height = 4)
