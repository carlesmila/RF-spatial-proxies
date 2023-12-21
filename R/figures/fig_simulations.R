#------------------------------------------------------------------------------#
#                        Figures: Simulation results                           #
#------------------------------------------------------------------------------#

library("cowplot")
library("tidyverse")

# Preparation ----
res10_partial <- read_csv("results/simulations/sim10_partial.csv")
res10_complete <- read_csv("results/simulations/sim10_complete.csv")
res40_partial <- read_csv("results/simulations/sim40_partial.csv")
res40_complete <- read_csv("results/simulations/sim40_complete.csv")
res <- bind_rows(res10_partial, res10_complete, res40_partial, res40_complete)
res$nsample <- NULL
names(res) <- gsub("coord", "Coordinates", names(res))
names(res) <- gsub("baseline", "Baseline", names(res))
res <- res |>
  mutate(dsample = case_when(
    dsample == "sregular" ~ "Regular",
    dsample == "random" ~ "Random",
    dsample == "wclust" ~ "Weakly\nclustered",
    dsample == "sclust" ~ "Strongly\nclustered"),
    simtype = case_when(
      simtype == "partial" ~ "Partial",
      simtype == "complete" ~ "Complete"))
res$range <- paste0("Range: ", res$range)
res$simtype <- paste0("Scenario: ", res$simtype)
res$dsample <- fct_inorder(res$dsample)
res$range <- fct_inorder(res$range)
res$simtype <- fct_inorder(res$simtype)


# True RMSE ----
accdata <- dplyr::select(res, range, dsample, simtype, contains("_surface_")) |>
  select(-RFGLS_surface_RMSE) |>
  pivot_longer(contains("_RMSE"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup()
p1 <- ggplot(accdata) +
  geom_boxplot(aes(y = RMSE, x = dsample, col = Model)) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(range ~ simtype) +
  xlab("Sampling pattern") + ylab("RMSE") +
  theme_bw() +
  theme(legend.position = "bottom")
# ggsave("figures/sims_trueRMSE.png", p1, width = 7, height = 5)


# Variable importance ----
varimpdata <- dplyr::select(res, range, dsample, simtype, contains("imp")) |>
  pivot_longer(-c(range, dsample, simtype), names_to = "metric", values_to = "featimp") |>
  group_by(metric) |>
  mutate(Model = strsplit(metric, "_")[[1]][1]) |>
  ungroup() |>
  mutate(featproxies = 100-featimp)
# with(varimpdata[varimpdata$dsample=="Random" & varimpdata$simtype == "Scenario: Partial" & varimpdata$range == "Range: 40",],
#      tapply(featproxies, Model, median))
# with(varimpdata[varimpdata$dsample=="Random" & varimpdata$simtype == "Scenario: Partial" & varimpdata$range == "Range: 40",],
#      tapply(featproxies, Model, IQR))
p2 <- ggplot(data = varimpdata) +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_hline(yintercept = 100, alpha = 0.2) +
  geom_boxplot(aes(x = dsample, y = featproxies, col = Model)) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(range ~ simtype) +
  ylim(0, 100) +
  ylab("% of total mean impurity decrease\nattributable to spatial proxies") + xlab("Sampling pattern") +
  theme_bw()+
  theme(legend.position = "bottom")
# ggsave("figures/sims_varimp.png", p2, width = 7, height = 5)


# Extrapolation ----
extradata <- dplyr::select(res, range, dsample, simtype, contains("AOA")) |>
   mutate(across(ends_with("AOA"), function(x) 100-x))
extradata <- extradata |>
  pivot_longer(-c(range, dsample, simtype), names_to = "stat", values_to = "value") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup()
with(extradata[extradata$Model=="EDF" & extradata$simtype == "Scenario: Partial" & extradata$range == "Range: 10",],
     tapply(value, dsample, median))
with(extradata[extradata$Model=="EDF" & extradata$simtype == "Scenario: Partial" & extradata$range == "Range: 10",],
     tapply(value, dsample, IQR))
p3 <- ggplot(extradata) +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_hline(yintercept = 100, alpha = 0.2) +
  geom_boxplot(aes(x = dsample, y = value, col = Model)) +
  facet_grid(range ~ simtype) +
  scale_color_brewer(palette = "Dark2") +
  ylab("% of the prediction area") + xlab("Metric") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("figures/sims_applicability.png", p3, width = 7, height = 5)


# Cross-validation ----
cvdata <- dplyr::select(res, range, dsample, simtype, contains("_RMSE")) |>
  select(-RFGLS_surface_RMSE)
cvdata <- cvdata |>
  pivot_longer(contains("_RMSE"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(model = strsplit(stat, "_")[[1]][1]) |>
  mutate(CV = strsplit(stat, "_")[[1]][2],
         CV = case_when(CV == "random" ~ "Random\n5-fold CV",
                        CV == "spatial" ~ "kNNDM\n5-fold CV",
                        CV == "surface" ~ "True\naccuracy")) |>
  ungroup() |>
  mutate(CV = fct_relevel(CV, c("True\naccuracy", "Random\n5-fold CV", "kNNDM\n5-fold CV")))
p4 <- cvdata[cvdata$simtype == "Scenario: Partial",] |>
  ggplot() +
  geom_boxplot(aes(y = RMSE, x = CV, colour = model)) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(range ~ dsample) +
  theme_bw() + ylab("RMSE") + xlab("Accuracy estimation method") +
  labs(col = "Model") +
  theme(legend.position = "bottom")
p5 <-  cvdata[cvdata$simtype == "Scenario: Complete",] |>
  ggplot() +
  geom_boxplot(aes(y = RMSE, x = CV, colour = model)) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(range ~ dsample) +
  theme_bw() + ylab("RMSE") + xlab("Accuracy estimation method") +
  labs(col = "Model") +
  theme(legend.position = "bottom")
ggsave("figures/sims_CV_partial.png", p4, width = 9, height = 5)
ggsave("figures/sims_CV_complete.png", p5, width = 9, height = 5)

# RF-GLS ----
accdata <- dplyr::select(res, range, dsample, simtype, contains("_surface_")) |>
  pivot_longer(contains("_RMSE"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1],
         Model = as.factor(ifelse(Model == "RFGLS", "RF-GLS", Model)),
         Model = fct_relevel(Model, "RF-GLS", after = 5)) |>
  ungroup()
accdata2 <- accdata[!grepl("RFGLS", accdata$stat),] |>
  group_by(range, dsample, simtype, Model) |>
  summarise(RMSE = median(RMSE)) |>
  group_by(range, dsample, simtype) |>
  filter(RMSE == min(RMSE)) |>
  ungroup()
accdata3 <- data.frame()
for(i in 1:nrow(accdata2)){
  accdata3 <- rbind(accdata3,
                    filter(accdata,
                           range == accdata2$range[i],
                           dsample == accdata2$dsample[i],
                           simtype == accdata2$simtype[i],
                           Model == accdata2$Model[i]))
}
accdata3$Model <- "Best-performing standard RF"
accdata3 <- rbind(accdata3, filter(accdata, Model == "RF-GLS"))

p6 <- ggplot(accdata3) +
  geom_boxplot(aes(y = RMSE, x = dsample, col = Model)) +
  scale_color_manual(values = c("grey60", "black")) +
  facet_grid(range ~ simtype) +
  labs(colour = "Model") +
  xlab("Sampling pattern") + ylab("RMSE") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("figures/sims_RF-GLS.png", p6, width = 7, height = 5)
