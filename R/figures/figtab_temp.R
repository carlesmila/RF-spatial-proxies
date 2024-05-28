#------------------------------------------------------------------------------#
#                           Figures: Temperature                               #
#------------------------------------------------------------------------------#

library("tidyverse")
library("knitr")
library("CAST")
library("caret")
library("cowplot")
library("gstat")
library("gridExtra")
library("sf")
library("terra")
library("tmap")
library("spatstat")
library("xtable")

# Data prep ----

# Training data, CV folds
traindata <- st_read("data/temp/temp_train.gpkg", quiet = T)
rCVfolds <- read_rds("results/case/temp_rfolds.rds")
rCVfolds <- ifelse(rCVfolds==0, 10, rCVfolds)
rCVfolds <- factor(rCVfolds, levels = 1:10)
sCVfolds <- read_rds("results/case/temp_sfolds.rds")$cluster
sCVfolds <- factor(sCVfolds, levels = 1:10)
spain <- read_sf("data/boundaries/spain.gpkg")

# tab results naive
tab_naive <- bind_rows(
  read_csv("results/case/temp_naive_baseline.csv") |> mutate(model = "Baseline"),
  read_csv("results/case/temp_naive_coord.csv") |> mutate(model = "Coordinates"),
  read_csv("results/case/temp_naive_EDF.csv") |> mutate(model = "EDF"),
  read_csv("results/case/temp_naive_RFsp.csv") |> mutate(model = "RFsp")) |>
  mutate(predictors = "naive")

# tab results complete
tab_complete <- bind_rows(
  read_csv("results/case/temp_complete_baseline.csv") |> mutate(model = "Baseline"),
  read_csv("results/case/temp_complete_coord.csv") |> mutate(model = "Coordinates"),
  read_csv("results/case/temp_complete_EDF.csv") |> mutate(model = "EDF"),
  read_csv("results/case/temp_complete_RFsp.csv") |> mutate(model = "RFsp")) |>
  mutate(predictors = "complete")

# Predictions
preds_naive <- c(rast("results/case/temp_naive_baseline.tif"),
                 rast("results/case/temp_naive_coord.tif"),
                 rast("results/case/temp_naive_EDF.tif"),
                 rast("results/case/temp_naive_RFsp.tif"))
names(preds_naive) <- c("Baseline","Coordinates", "EDF","RFsp")
preds_complete <- c(rast("results/case/temp_complete_baseline.tif"),
                 rast("results/case/temp_complete_coord.tif"),
                 rast("results/case/temp_complete_EDF.tif"),
                 rast("results/case/temp_complete_RFsp.tif"))
names(preds_complete) <- c("Baseline","Coordinates", "EDF","RFsp")


# Figure: points, folds, distances ----
p1 <- ggplot() +
  geom_sf(data = spain, alpha = 0) +
  geom_sf(data = traindata, aes(col = rCVfolds), size = 1) +
  theme_bw() + theme(legend.position = "none") +
  scale_color_brewer(palette = "Paired") +
  labs(col = "CV folds") +
  ggtitle("Random 10-fold CV station assignment")
p2 <- ggplot() +
  geom_sf(data = spain, alpha = 0) +
  geom_sf(data = traindata, aes(col = sCVfolds), size = 1) +
  theme_bw() + theme(legend.position = "none") +
  scale_color_brewer(palette = "Paired") +
  labs(col = "CV folds") +
  ggtitle("kNNDM 10-fold CV station assignment")
p3 <- plot(geodist(traindata, spain, cvfolds = rCVfolds), stat="ecdf") +
  ggtitle("Random 10-fold CV nearest neighbour distances") +
  theme(legend.title = element_text(size = 11),
        legend.text = element_text(size = 8))
p4 <- plot(geodist(traindata, spain, cvfolds = sCVfolds), stat="ecdf") +
  ggtitle("kNNDM 10-fold CV nearest neighbour distances") +
  theme(legend.title = element_text(size = 11),
        legend.text = element_text(size = 8))
pall <- plot_grid(p1, p2, p3, p4, ncol=2)
pall
# save_plot("figures/temp_CVmap.pdf", pall, base_height = 7.5, base_asp = 1.45)


# Figure: variograms ----

# Raw outcome
empvar1 <- variogram(as.formula("temp~1"), cutoff = 800000, traindata, width = 20000)
fitvar1 <- fit.variogram(empvar1, vgm(model="Sph", nugget = 0), fit.method = 1)

# Naive baseline
mod <- readRDS("results/case/temp_naive_baseline.rds")
traindata$res <- residuals(mod)
empvar2 <- variogram(as.formula("res~1"), cutoff = 800000, traindata, width = 20000)
fitvar2 <- fit.variogram(empvar2, vgm(model="Gau", nugget = 0), fit.method = 1)

# Naive coords
mod <- readRDS("results/case/temp_naive_coord.rds")
traindata$res <- residuals(mod)
empvar3 <- variogram(as.formula("res~1"), cutoff = 800000, traindata, width = 20000)

# Naive EDF
mod <- readRDS("results/case/temp_naive_EDF.rds")
traindata$res <- residuals(mod)
empvar4 <- variogram(as.formula("res~1"), cutoff = 800000, traindata, width = 20000)

# Naive RFsp
mod <- readRDS("results/case/temp_naive_RFsp.rds")
traindata$res <- residuals(mod)
empvar5 <- variogram(as.formula("res~1"), cutoff = 800000, traindata, width = 20000)

# Complete baseline
mod <- readRDS("results/case/temp_complete_baseline.rds")
traindata$res <- residuals(mod)
empvar6 <- variogram(as.formula("res~1"), cutoff = 800000, traindata, width = 20000)

# Complete coords
mod <- readRDS("results/case/temp_complete_coord.rds")
traindata$res <- residuals(mod)
empvar7 <- variogram(as.formula("res~1"), cutoff = 800000, traindata, width = 20000)

# Complete EDF
mod <- readRDS("results/case/temp_complete_EDF.rds")
traindata$res <- residuals(mod)
empvar8 <- variogram(as.formula("res~1"), cutoff = 800000, traindata, width = 20000)

# Complete RFsp
mod <- readRDS("results/case/temp_complete_RFsp.rds")
traindata$res <- residuals(mod)
empvar9 <- variogram(as.formula("res~1"), cutoff = 800000, traindata, width = 20000)

# # Figure
# pdf("figures/temp_variograms.pdf", width = 9, height = 9)
# grid.arrange(plot(empvar1, fitvar1, main = "Outcome"),
#              plot(empvar2, fitvar2, main = "Naive: Baseline"),
#              plot(empvar3, main = "Naive: Coordinates"),
#              plot(empvar4, main = "Naive: EDF"),
#              plot(empvar5, main = "Naive: RFsp"),
#              plot(empvar6, main = "Complete: Baseline"),
#              plot(empvar7, main = "Complete: Coordinates"),
#              plot(empvar8, main = "Complete: EDF"),
#              plot(empvar9, main = "Complete: RFsp"),
#              ncol = 3)
# dev.off()

# Figure: sample distribution ----
samples_ppp <- as.ppp(st_coordinates(traindata), W = as.owin(spain))
Gppp <- envelope(samples_ppp, Gest, fix.n = TRUE, global = TRUE)
Fppp <- envelope(samples_ppp, Fest, fix.n = TRUE, global = TRUE)
Kppp <- envelope(samples_ppp, Kest, fix.n = TRUE, global = TRUE)

# # Plot (in m)
# pdf("figures/temp_ppp.pdf", width = 9, height = 4)
# par(mfrow = c(1, 3))
# plot(Gppp, xlab = "r (m)", main = "A)")
# plot(Fppp, xlab = "r (m)", main = "B)")
# plot(Kppp, xlab = "r (m)", main = "C)")
# par(mfrow = c(1, 1))
# dev.off()

# Table: results ----
tabboth <- rbind(tab_naive, tab_complete) |>
  mutate(random_RMSE = paste0(round(random_RMSE_mean, 2),
                              " (", round(random_RMSE_sd, 2), ")"),
         random_Rsquared = paste0(round(random_Rsquared_mean, 2),
                              " (", round(random_Rsquared_sd, 2), ")")) |>
  dplyr::select(-random_RMSE_mean, -random_RMSE_sd,
                -random_Rsquared_mean, -random_Rsquared_sd) |>
  mutate(improxy = 100-impfeat,
         nonAOA = 100-AOA*100) |>
  select(predictors, model,
         random_RMSE, random_Rsquared,
         kNNDM_RMSE, kNNDM_Rsquared,
         nonAOA, improxy)
kable(tabboth, digits = 2)
print(xtable(tabboth), include.rownames = FALSE)

# Figure: Prediction maps ----
tempbreaks <- seq(6, 20, 2)
predmap1 <- tm_shape(preds_naive) +
  tm_raster(palette = hcl.colors(9, "Blue-Red 3"),
            style = "fixed", breaks = tempbreaks,
            title = "Temperature (ºC)", midpoint = NA,
            legend.show = FALSE) +
  tm_facets(nrow = 1) +
  tm_layout(main.title = "A) Naive model", main.title.size = 1.2)
predmap2 <- tm_shape(preds_complete) +
  tm_raster(palette = hcl.colors(9, "Blue-Red 3"),
            style = "fixed", breaks = tempbreaks,
            title = "Temperature (ºC)", midpoint = NA,
            legend.show = FALSE) +
  tm_facets(nrow = 1) +
  tm_layout(main.title = "B) Complete model", main.title.size = 1.2)
predlegend <- tm_shape(preds_complete[[1]]) +
  tm_raster(palette = hcl.colors(9, "Blue-Red 3"),
            style = "fixed", breaks = tempbreaks,
            midpoint = NA, title = "     Predicted air temperature (ºC)",
            legend.format = list(digits=0),
            legend.is.portrait = FALSE) +
  tm_layout(legend.position=c("center", "bottom"), legend.only = TRUE,
            legend.width = 1, legend.height = 0.1, legend.text.size = 5,
            legend.title.size = 14)
predmap_both <- tmap_arrange(predmap1, predmap2, predlegend,
                             nrow = 3, heights = c(0.425, 0.425, 0.15))
tmap_save(predmap_both, "figures/temp_predictions.pdf", width = 8, height = 5)
