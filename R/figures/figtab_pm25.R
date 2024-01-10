#------------------------------------------------------------------------------#
#                               Figures: PM2.5                                 #
#------------------------------------------------------------------------------#

library("tidyverse")
library("knitr")
library("CAST")
library("caret")
library("gridExtra")
library("gstat")
library("sf")
library("terra")
library("tmap")
library("spatstat")
library("xtable")

# Data prep ----

# Training data, CV folds
traindata <- st_read("data/AP/PM25_train.gpkg", quiet = T)
rCVfolds <- read_rds("results/case/pm25_rfolds.rds")
rCVfolds <- ifelse(rCVfolds==0, 10, rCVfolds)
rCVfolds <- factor(rCVfolds, levels = 1:10)
sCVfolds <- read_rds("results/case/pm25_sfolds.rds")$cluster
sCVfolds <- factor(sCVfolds, levels = 1:10)
spain <- st_read("data/boundaries/spain.gpkg", quiet = T)

# tab results naive
tab_naive <- bind_rows(
  read_csv("results/case/pm25_naive_baseline.csv") |> mutate(model = "Baseline"),
  read_csv("results/case/pm25_naive_coord.csv") |> mutate(model = "Coordinates"),
  read_csv("results/case/pm25_naive_EDF.csv") |> mutate(model = "EDF"),
  read_csv("results/case/pm25_naive_RFsp.csv") |> mutate(model = "RFsp")) |>
  mutate(predictors = "naive")

# tab results complete
tab_complete <- bind_rows(
  read_csv("results/case/pm25_complete_baseline.csv") |> mutate(model = "Baseline"),
  read_csv("results/case/pm25_complete_coord.csv") |> mutate(model = "Coordinates"),
  read_csv("results/case/pm25_complete_EDF.csv") |> mutate(model = "EDF"),
  read_csv("results/case/pm25_complete_RFsp.csv") |> mutate(model = "RFsp")) |>
  mutate(predictors = "complete")

# Predictions
preds_naive <- c(rast("results/case/pm25_naive_baseline.tif"),
                 rast("results/case/pm25_naive_coord.tif"),
                 rast("results/case/pm25_naive_EDF.tif"),
                 rast("results/case/pm25_naive_RFsp.tif"))
names(preds_naive) <- c("Baseline","Coordinates", "EDF","RFsp")
preds_complete <- c(rast("results/case/pm25_complete_baseline.tif"),
                    rast("results/case/pm25_complete_coord.tif"),
                    rast("results/case/pm25_complete_EDF.tif"),
                    rast("results/case/pm25_complete_RFsp.tif"))
names(preds_complete) <- c("Baseline","Coordinates", "EDF","RFsp")

# Figure: points, folds, distances ----
p1 <- ggplot() +
  geom_sf(data = spain, alpha = 0) +
  geom_sf(data = traindata, aes(col = rCVfolds), size = 1) +
  theme_bw() +
  scale_color_brewer(palette = "Paired") +
  labs(col = "CV folds") +
  ggtitle("Random 10-fold CV station assignment")
p2 <- ggplot() +
  geom_sf(data = spain, alpha = 0) +
  geom_sf(data = traindata, aes(col = sCVfolds), size = 1) +
  theme_bw() +
  scale_color_brewer(palette = "Paired") +
  labs(col = "CV folds") +
  ggtitle("kNNDM 10-fold CV station assignment")
p3 <- plot_geodist(traindata, spain, cvfolds = rCVfolds, showPlot = F, stat = "ecdf")$plot +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  ggtitle("Random 10-fold CV nearest neighbour distances")
p4 <- plot_geodist(traindata, spain, cvfolds = sCVfolds, showPlot = F, stat = "ecdf")$plot +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  ggtitle("kNNDM 10-fold CV nearest neighbour distances")
grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
# png("figures/pm25_CVmap.png", width = 3000, height = 2000, res = 300)
# grid.arrange(p1, p2, p3, p4, ncol = 2)
# dev.off()


# Figure: variograms ----

# Raw outcome
empvar1 <- variogram(as.formula("PM25~1"), cutoff = 700000, traindata, width = 20000)
fitvar1 <- fit.variogram(empvar1, vgm(model="Exp", nugget = 0), fit.method = 1)

# Naive baseline
mod <- readRDS("results/case/pm25_naive_baseline.rds")
traindata$res <- residuals(mod)
empvar2 <- variogram(as.formula("res~1"), cutoff = 700000, traindata, width = 20000)
fitvar2 <- fit.variogram(empvar2, vgm(model="Exp", nugget = 0), fit.method = 1)

# Naive coords
mod <- readRDS("results/case/pm25_naive_coord.rds")
traindata$res <- residuals(mod)
empvar3 <- variogram(as.formula("res~1"), cutoff = 700000, traindata, width = 20000)

# Naive EDF
mod <- readRDS("results/case/pm25_naive_EDF.rds")
traindata$res <- residuals(mod)
empvar4 <- variogram(as.formula("res~1"), cutoff = 700000, traindata, width = 20000)

# Naive RFsp
mod <- readRDS("results/case/pm25_naive_RFsp.rds")
traindata$res <- residuals(mod)
empvar5 <- variogram(as.formula("res~1"), cutoff = 700000, traindata, width = 20000)

# Complete baseline
mod <- readRDS("results/case/pm25_complete_baseline.rds")
traindata$res <- residuals(mod)
empvar6 <- variogram(as.formula("res~1"), cutoff = 700000, traindata, width = 20000)

# Complete coords
mod <- readRDS("results/case/pm25_complete_coord.rds")
traindata$res <- residuals(mod)
empvar7 <- variogram(as.formula("res~1"), cutoff = 700000, traindata, width = 20000)

# Complete EDF
mod <- readRDS("results/case/pm25_complete_EDF.rds")
traindata$res <- residuals(mod)
empvar8 <- variogram(as.formula("res~1"), cutoff = 700000, traindata, width = 20000)

# Complete RFsp
mod <- readRDS("results/case/pm25_complete_RFsp.rds")
traindata$res <- residuals(mod)
empvar9 <- variogram(as.formula("res~1"), cutoff = 700000, traindata, width = 20000)

# # Figure
# png("figures/pm25_variograms.png", width = 3500, height = 3600, res = 300)
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
# png("figures/pm25_ppp.png", width = 3200, height = 1200, res = 300)
# par(mfrow = c(1, 3))
# plot(Gppp, xlab = "r (m)", main = "A)")
# plot(Fppp, xlab = "r (m)", main = "B)")
# plot(Kppp, xlab = "r (m)", main = "C)")
# par(mfrow = c(1, 1))
# dev.off()

# Table: results ----
tabboth <- rbind(tab_naive, tab_complete) |>
  mutate(improxy = 100-impfeat,
         nonAOA = 100-AOA*100) |>
  select(predictors, model,
         random_RMSE, random_Rsquared,
         kNNDM_RMSE, kNNDM_Rsquared,
         nonAOA, improxy)
kable(tabboth, digits = 2)
print(xtable(tabboth), include.rownames = FALSE)


# Figure: Prediction maps ----
colbreaks <- seq(4,20,2)
predmap1 <- tm_shape(preds_naive) +
  tm_raster(palette = hcl.colors(10, "cividis"),
            style = "fixed", breaks = colbreaks,
            title = "PM\U02082.\U02085 (\u03bcg/m\U000b3)", midpoint = NA,
            legend.show = FALSE) +
  tm_facets(nrow = 1) +
  tm_layout(main.title = "A) Naive model", main.title.size = 1.2)
predmap2 <- tm_shape(preds_complete) +
  tm_raster(palette = hcl.colors(10, "cividis"),
            style = "fixed", breaks = colbreaks,
            title = "PM\U02082.\U02085 (\u03bcg/m\U000b3)", midpoint = NA,
            legend.show = FALSE) +
  tm_facets(nrow = 1) +
  tm_layout(main.title = "B) Complete model", main.title.size = 1.2)
predlegend <- tm_shape(preds_complete[[1]]) +
  tm_raster(palette = hcl.colors(10, "cividis"),
            style = "fixed", breaks = colbreaks,
            midpoint = NA, title = "               Predicted PM\U02082.\U02085 (\u03bcg/m\U000b3)",
            legend.format = list(digits=0),
            legend.is.portrait = FALSE) +
  tm_layout(legend.position=c("center", "bottom"), legend.only = TRUE,
            legend.width = 1, legend.height = 0.1, legend.text.size = 5,
            legend.title.size = 14)
predmap_both <- tmap_arrange(predmap1, predmap2, predlegend,
                             nrow = 3, heights = c(0.425, 0.425, 0.15))
# tmap_save(predmap_both, "figures/pm25_predictions.png", width = 8, height = 5)
