#------------------------------------------------------------------------------#
#                        Figures: Simulation results                           #
#------------------------------------------------------------------------------#

library("cowplot")
library("tidyverse")
library("RColorBrewer")

# Preparation ----
res10_proxyonly <- read_csv("results/simulations/sim10_proxyonly.csv")
res10_partial <- read_csv("results/simulations/sim10_partial.csv")
res10_complete <- read_csv("results/simulations/sim10_complete.csv")
res10_autocor <- read_csv("results/simulations/sim10_autocor.csv")
res40_proxyonly <- read_csv("results/simulations/sim40_proxyonly.csv")
res40_partial <- read_csv("results/simulations/sim40_partial.csv")
res40_complete <- read_csv("results/simulations/sim40_complete.csv")
res40_autocor <- read_csv("results/simulations/sim40_autocor.csv")
res <- bind_rows(res10_autocor, res10_partial, res10_complete, res10_proxyonly,
                 res40_autocor, res40_partial, res40_complete, res40_proxyonly)
res$nsample <- NULL
names(res) <- gsub("coord", "Coordinates", names(res))
names(res) <- gsub("baseline", "Baseline", names(res))
res <- res |>
  mutate(dsample = case_when(
    dsample == "sregular" ~ "Regular",
    dsample == "random" ~ "Random",
    dsample == "wclust" ~ "Weakly\nclustered",
    dsample == "sclust" ~ "Strongly\nclustered"),
    scenario = case_when(
      scenario == "partial" ~ "Scenario:\nMissing predictors",
      scenario == "complete" ~ "Scenario:\nComplete",
      scenario == "autocor" ~ "Scenario:\nAutocorrelated error",
      scenario == "proxyonly" ~ "Scenario:\nProxies only"
      ))
res$range <- paste0("Range: ", res$range)
res$dsample <- fct_inorder(res$dsample)
res$range <- fct_inorder(res$range)


# Main: True extrapolation ----
min_error_line <- rbind(
  data.frame(scenario = c("Scenario:\nComplete",
                          "Scenario:\nMissing predictors",
                          "Scenario:\nProxies only"), int=1),
  data.frame(scenario = c("Scenario:\nAutocorrelated error"), int=0))
figdata <- res |>
  dplyr::select(scenario, dsample, range, contains("surface_extra")) |>
  dplyr::select(!contains("RFGLS")) |>
  pivot_longer(-c("scenario", "dsample", "range"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup()
p <- ggplot() +
  geom_hline(data = min_error_line, aes(yintercept = int), lty=2, col="grey20") +
  geom_boxplot(data = figdata,
               aes(y=RMSE, x=dsample, col=Model), outlier.size=0.25, outlier.alpha=0.75) +
  facet_grid(range ~ scenario) +
  ylim(0, 5) +
  scale_color_brewer(palette = "Dark2") +
  xlab("Sampling pattern") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 8))
ggsave(p, filename = "figures/sims_TrueExtra.pdf", width = 9, height = 5)
rm("figdata", "p")

# Main: True interpolation ----
min_error_line <- rbind(
  data.frame(scenario = c("Scenario:\nComplete",
                          "Scenario:\nMissing predictors",
                          "Scenario:\nProxies only"), int=1),
  data.frame(scenario = c("Scenario:\nAutocorrelated error"), int=0))
figdata <- res |>
  dplyr::select(scenario, dsample, range, contains("surface_inter")) |>
  dplyr::select(!contains("RFGLS")) |>
  pivot_longer(-c("scenario", "dsample", "range"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup()
p <- ggplot() +
  geom_hline(data = min_error_line, aes(yintercept = int), lty=2, col="grey20") +
  geom_boxplot(data = figdata,
               aes(y=RMSE, x=dsample, col=Model), outlier.size=0.25, outlier.alpha=0.75) +
  facet_grid(range ~ scenario, scales = "free") +
  ylim(0,3) +
  scale_color_brewer(palette = "Dark2") +
  xlab("Sampling pattern") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 8))
ggsave(p, filename = "figures/sims_TrueInter.pdf", width = 9, height = 5)
rm("figdata", "p")

# Appendix: Extrapolation AOA ----
figdata <- res |>
  dplyr::select(scenario, dsample, range, contains("AOA_extra")) |>
  pivot_longer(-c("scenario", "dsample", "range"), names_to = "stat", values_to = "AOA") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup() |>
  mutate(extrap = 100 - AOA)
p <- ggplot(data = figdata) +
  geom_hline(aes(yintercept = 0), col = "grey50", alpha = 0.5) +
  geom_hline(aes(yintercept = 100), col = "grey50", alpha = 0.5) +
  geom_boxplot(aes(y=extrap, x=dsample, col=Model), outlier.size=0.25, outlier.alpha=0.75) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(range ~ scenario) +
  ylab("Feature extrapolation (%)") + xlab("Sampling pattern")  +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 8))
ggsave(p, filename = "figures/sims_AOAExtra.pdf", width = 9, height = 5)
rm("figdata", "p")

# Appendix: Interpolation AOA ----
figdata <- res |>
  dplyr::select(scenario, dsample, range, contains("AOA_inter")) |>
  pivot_longer(-c("scenario", "dsample", "range"), names_to = "stat", values_to = "AOA") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup() |>
  mutate(extrap = 100 - AOA)
p <- ggplot(data = figdata) +
  geom_hline(aes(yintercept = 0), col = "grey50", alpha = 0.5) +
  geom_hline(aes(yintercept = 100), col = "grey50", alpha = 0.5) +
  geom_boxplot(aes(y=extrap, x=dsample, col=Model), outlier.size=0.25, outlier.alpha=0.75) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(range ~ scenario) +
  ylab("Feature extrapolation (%)") + xlab("Sampling pattern")  +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 8))
ggsave(p, filename = "figures/sims_AOAInter.pdf", width = 9, height = 5)
rm("figdata", "p")


# Appendix: varimp ----
figdata <- res |>
  dplyr::select(scenario, dsample, range, contains("impfeat")) |>
  pivot_longer(-c("scenario", "dsample", "range"), names_to = "stat", values_to = "impfeat") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1],
         proxyfeat = 100-impfeat) |>
  ungroup()
p <- ggplot(data = figdata) +
  geom_hline(aes(yintercept = 0), col = "grey50", alpha = 0.5) +
  geom_hline(aes(yintercept = 100), col = "grey50", alpha = 0.5) +
  geom_boxplot(aes(y=proxyfeat, x=dsample, col=Model), outlier.size=0.5, outlier.alpha=0.75) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(aes(yintercept = 0), lty = 2) +
  facet_grid(range ~ scenario) +
  ylab("Proxy importance (%)") + xlab("Sampling pattern")  +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 8))
ggsave(p, filename = "figures/sims_varimp.pdf", width = 9, height = 5)
rm("figdata", "p")



# Main: validation extra autocor ----
figdata <- res |>
  dplyr::filter(scenario == "Scenario:\nAutocorrelated error") |>
  dplyr::mutate(Baseline_random_extra = Baseline_random_inter, # Random CV is the same for inter and extra
                Coordinates_random_extra = Coordinates_random_inter,
                EDF_random_extra = EDF_random_inter,
                RFsp_random_extra = RFsp_random_inter) |>
  dplyr::select(!contains("RFGLS")) |>
  dplyr::select(dsample, range, contains("_surface_extra"), contains("_test_extra"),
                contains("_random_extra"), contains("spatial_extra")) |>
  pivot_longer(-c("dsample", "range"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1],
         Evaluation = strsplit(stat, "_")[[1]][2],
         Area = strsplit(stat, "_")[[1]][3]) |>
  mutate(Evaluation = case_when(
    Evaluation == "random" ~ "Random\n5-fold CV",
    Evaluation == "spatial" ~ "kNNDM\n5-fold CV",
    Evaluation == "test" ~ "Prob. test\nsamples",
    Evaluation == "surface" ~ "True\naccuracy"),
    Evaluation = as.factor(Evaluation),
    Evaluation = fct_relevel(Evaluation, "Prob. test\nsamples"),
    Evaluation = fct_relevel(Evaluation, "True\naccuracy")) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup() |>
  mutate(dsample = gsub("\n", " ", dsample, fixed=TRUE),
         dsample = paste0("Sampling pattern:\n", dsample),
         dsample = fct_inorder(dsample))
p <- ggplot(data = figdata) +
  geom_boxplot(aes(y=RMSE, x=Evaluation, col=Model), outlier.size=0.5, outlier.alpha=0.75) +
  ylim(0.5, 5) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(range ~ dsample) +
  xlab("Evaluation method") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 7))
ggsave(p, filename = "figures/sims_valextra_autocor.pdf", width = 9, height = 5)
rm("figdata", "p")

# Appendix: validation extra complete ----
figdata <- res |>
  dplyr::filter(scenario == "Scenario:\nComplete") |>
  dplyr::mutate(Baseline_random_extra = Baseline_random_inter, # Random CV is the same for inter and extra
                Coordinates_random_extra = Coordinates_random_inter,
                EDF_random_extra = EDF_random_inter,
                RFsp_random_extra = RFsp_random_inter) |>
  dplyr::select(!contains("RFGLS")) |>
  dplyr::select(dsample, range, contains("_surface_extra"), contains("_test_extra"),
                contains("_random_extra"), contains("spatial_extra")) |>
  pivot_longer(-c("dsample", "range"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1],
         Evaluation = strsplit(stat, "_")[[1]][2],
         Area = strsplit(stat, "_")[[1]][3]) |>
  mutate(Evaluation = case_when(
    Evaluation == "random" ~ "Random\n5-fold CV",
    Evaluation == "spatial" ~ "kNNDM\n5-fold CV",
    Evaluation == "test" ~ "Prob. test\nsamples",
    Evaluation == "surface" ~ "True\naccuracy"),
    Evaluation = as.factor(Evaluation),
    Evaluation = fct_relevel(Evaluation, "Prob. test\nsamples"),
    Evaluation = fct_relevel(Evaluation, "True\naccuracy")) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup() |>
  mutate(dsample = gsub("\n", " ", dsample, fixed=TRUE),
         dsample = paste0("Sampling pattern:\n", dsample),
         dsample = fct_inorder(dsample))
p <- ggplot(data = figdata) +
  geom_boxplot(aes(y=RMSE, x=Evaluation, col=Model), outlier.size=0.5, outlier.alpha=0.75) +
  ylim(1, 4) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(range ~ dsample) +
  xlab("Evaluation method") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 7))
ggsave(p, filename = "figures/sims_valextra_complete.pdf", width = 9, height = 5)
rm("figdata", "p")

# Appendix: validation extra missing ----
figdata <- res |>
  dplyr::filter(scenario == "Scenario:\nMissing predictors") |>
  dplyr::mutate(Baseline_random_extra = Baseline_random_inter, # Random CV is the same for inter and extra
                Coordinates_random_extra = Coordinates_random_inter,
                EDF_random_extra = EDF_random_inter,
                RFsp_random_extra = RFsp_random_inter) |>
  dplyr::select(!contains("RFGLS")) |>
  dplyr::select(dsample, range, contains("_surface_extra"), contains("_test_extra"),
                contains("_random_extra"), contains("spatial_extra")) |>
  pivot_longer(-c("dsample", "range"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1],
         Evaluation = strsplit(stat, "_")[[1]][2],
         Area = strsplit(stat, "_")[[1]][3]) |>
  mutate(Evaluation = case_when(
    Evaluation == "random" ~ "Random\n5-fold CV",
    Evaluation == "spatial" ~ "kNNDM\n5-fold CV",
    Evaluation == "test" ~ "Prob. test\nsamples",
    Evaluation == "surface" ~ "True\naccuracy"),
    Evaluation = as.factor(Evaluation),
    Evaluation = fct_relevel(Evaluation, "Prob. test\nsamples"),
    Evaluation = fct_relevel(Evaluation, "True\naccuracy")) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup() |>
  mutate(dsample = gsub("\n", " ", dsample, fixed=TRUE),
         dsample = paste0("Sampling pattern:\n", dsample),
         dsample = fct_inorder(dsample))
p <- ggplot(data = figdata) +
  geom_boxplot(aes(y=RMSE, x=Evaluation, col=Model), outlier.size=0.5, outlier.alpha=0.75) +
  ylim(1, 5) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(range ~ dsample) +
  xlab("Evaluation method") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 7))
ggsave(p, filename = "figures/sims_valextra_missing.pdf", width = 9, height = 5)
rm("figdata", "p")

# Appendix: validation extra proxyonly ----
figdata <- res |>
  dplyr::filter(scenario == "Scenario:\nProxies only") |>
  dplyr::mutate(Baseline_random_extra = Baseline_random_inter, # Random CV is the same for inter and extra
                Coordinates_random_extra = Coordinates_random_inter,
                EDF_random_extra = EDF_random_inter,
                RFsp_random_extra = RFsp_random_inter) |>
  dplyr::select(!contains("RFGLS")) |>
  dplyr::select(dsample, range, contains("_surface_extra"), contains("_test_extra"),
                contains("_random_extra"), contains("spatial_extra")) |>
  pivot_longer(-c("dsample", "range"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1],
         Evaluation = strsplit(stat, "_")[[1]][2],
         Area = strsplit(stat, "_")[[1]][3]) |>
  mutate(Evaluation = case_when(
    Evaluation == "random" ~ "Random\n5-fold CV",
    Evaluation == "spatial" ~ "kNNDM\n5-fold CV",
    Evaluation == "test" ~ "Prob. test\nsamples",
    Evaluation == "surface" ~ "True\naccuracy"),
    Evaluation = as.factor(Evaluation),
    Evaluation = fct_relevel(Evaluation, "Prob. test\nsamples"),
    Evaluation = fct_relevel(Evaluation, "True\naccuracy")) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup() |>
  mutate(dsample = gsub("\n", " ", dsample, fixed=TRUE),
         dsample = paste0("Sampling pattern:\n", dsample),
         dsample = fct_inorder(dsample)) |>
  filter(Model != "Baseline")
cols <- c("Coordinates"="#D95F02",
          "EDF"="#7570B3",
          "RFsp"= "#E7298A")
p <- ggplot(data = figdata) +
  geom_boxplot(aes(y=RMSE, x=Evaluation, col=Model), outlier.size=0.5, outlier.alpha=0.75) +
  scale_color_manual(values = cols) +
  facet_grid(range ~ dsample) +
  ylim(1, 6) +
  xlab("Evaluation method") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 7))
ggsave(p, filename = "figures/sims_valextra_proxyonly.pdf", width = 9, height = 5)
rm("figdata", "p")



# Main: validation inter autocor ----
figdata <- res |>
  dplyr::filter(scenario == "Scenario:\nAutocorrelated error") |>
  dplyr::mutate(Baseline_random_extra = Baseline_random_inter, # Random CV is the same for inter and extra
                Coordinates_random_extra = Coordinates_random_inter,
                EDF_random_extra = EDF_random_inter,
                RFsp_random_extra = RFsp_random_inter) |>
  dplyr::select(!contains("RFGLS")) |>
  dplyr::select(dsample, range, contains("_surface_inter"), contains("_test_inter"),
                contains("_random_inter"), contains("spatial_inter")) |>
  pivot_longer(-c("dsample", "range"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1],
         Evaluation = strsplit(stat, "_")[[1]][2],
         Area = strsplit(stat, "_")[[1]][3]) |>
  mutate(Evaluation = case_when(
    Evaluation == "random" ~ "Random\n5-fold CV",
    Evaluation == "spatial" ~ "kNNDM\n5-fold CV",
    Evaluation == "test" ~ "Prob. test\nsamples",
    Evaluation == "surface" ~ "True\naccuracy"),
    Evaluation = as.factor(Evaluation),
    Evaluation = fct_relevel(Evaluation, "Prob. test\nsamples"),
    Evaluation = fct_relevel(Evaluation, "True\naccuracy")) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup() |>
  mutate(dsample = gsub("\n", " ", dsample, fixed=TRUE),
         dsample = paste0("Sampling pattern:\n", dsample),
         dsample = fct_inorder(dsample))
p <- ggplot(data = figdata) +
  geom_boxplot(aes(y=RMSE, x=Evaluation, col=Model), outlier.size=0.5, outlier.alpha=0.75) +
  ylim(0.5, 3.5) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(range ~ dsample) +
  xlab("Evaluation method") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 7))
ggsave(p, filename = "figures/sims_valinter_autocor.pdf", width = 9, height = 5)
rm("figdata", "p")

# Appendix: validation inter complete ----
figdata <- res |>
  dplyr::filter(scenario == "Scenario:\nComplete") |>
  dplyr::mutate(Baseline_random_extra = Baseline_random_inter, # Random CV is the same for inter and extra
                Coordinates_random_extra = Coordinates_random_inter,
                EDF_random_extra = EDF_random_inter,
                RFsp_random_extra = RFsp_random_inter) |>
  dplyr::select(!contains("RFGLS")) |>
  dplyr::select(dsample, range, contains("_surface_inter"), contains("_test_inter"),
                contains("_random_inter"), contains("spatial_inter")) |>
  pivot_longer(-c("dsample", "range"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1],
         Evaluation = strsplit(stat, "_")[[1]][2],
         Area = strsplit(stat, "_")[[1]][3]) |>
  mutate(Evaluation = case_when(
    Evaluation == "random" ~ "Random\n5-fold CV",
    Evaluation == "spatial" ~ "kNNDM\n5-fold CV",
    Evaluation == "test" ~ "Prob. test\nsamples",
    Evaluation == "surface" ~ "True\naccuracy"),
    Evaluation = as.factor(Evaluation),
    Evaluation = fct_relevel(Evaluation, "Prob. test\nsamples"),
    Evaluation = fct_relevel(Evaluation, "True\naccuracy")) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup() |>
  mutate(dsample = gsub("\n", " ", dsample, fixed=TRUE),
         dsample = paste0("Sampling pattern:\n", dsample),
         dsample = fct_inorder(dsample))
p <- ggplot(data = figdata) +
  geom_boxplot(aes(y=RMSE, x=Evaluation, col=Model), outlier.size=0.5, outlier.alpha=0.75) +
  ylim(1, 3) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(range ~ dsample) +
  xlab("Evaluation method") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 7))
ggsave(p, filename = "figures/sims_valinter_complete.pdf", width = 9, height = 5)
rm("figdata", "p")

# Appendix: validation inter missing ----
figdata <- res |>
  dplyr::filter(scenario == "Scenario:\nMissing predictors") |>
  dplyr::mutate(Baseline_random_extra = Baseline_random_inter, # Random CV is the same for inter and extra
                Coordinates_random_extra = Coordinates_random_inter,
                EDF_random_extra = EDF_random_inter,
                RFsp_random_extra = RFsp_random_inter) |>
  dplyr::select(!contains("RFGLS")) |>
  dplyr::select(dsample, range, contains("_surface_inter"), contains("_test_inter"),
                contains("_random_inter"), contains("spatial_inter")) |>
  pivot_longer(-c("dsample", "range"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1],
         Evaluation = strsplit(stat, "_")[[1]][2],
         Area = strsplit(stat, "_")[[1]][3]) |>
  mutate(Evaluation = case_when(
    Evaluation == "random" ~ "Random\n5-fold CV",
    Evaluation == "spatial" ~ "kNNDM\n5-fold CV",
    Evaluation == "test" ~ "Prob. test\nsamples",
    Evaluation == "surface" ~ "True\naccuracy"),
    Evaluation = as.factor(Evaluation),
    Evaluation = fct_relevel(Evaluation, "Prob. test\nsamples"),
    Evaluation = fct_relevel(Evaluation, "True\naccuracy")) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup() |>
  mutate(dsample = gsub("\n", " ", dsample, fixed=TRUE),
         dsample = paste0("Sampling pattern:\n", dsample),
         dsample = fct_inorder(dsample))
p <- ggplot(data = figdata) +
  geom_boxplot(aes(y=RMSE, x=Evaluation, col=Model), outlier.size=0.5, outlier.alpha=0.75) +
  ylim(1, 3.5) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(range ~ dsample) +
  xlab("Evaluation method") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 7))
ggsave(p, filename = "figures/sims_valinter_missing.pdf", width = 9, height = 5)
rm("figdata", "p")

# Appendix: validation inter proxyonly ----
figdata <- res |>
  dplyr::filter(scenario == "Scenario:\nProxies only") |>
  dplyr::mutate(Baseline_random_extra = Baseline_random_inter, # Random CV is the same for inter and extra
                Coordinates_random_extra = Coordinates_random_inter,
                EDF_random_extra = EDF_random_inter,
                RFsp_random_extra = RFsp_random_inter) |>
  dplyr::select(!contains("RFGLS")) |>
  dplyr::select(dsample, range, contains("_surface_inter"), contains("_test_inter"),
                contains("_random_inter"), contains("spatial_inter")) |>
  pivot_longer(-c("dsample", "range"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1],
         Evaluation = strsplit(stat, "_")[[1]][2],
         Area = strsplit(stat, "_")[[1]][3]) |>
  mutate(Evaluation = case_when(
    Evaluation == "random" ~ "Random\n5-fold CV",
    Evaluation == "spatial" ~ "kNNDM\n5-fold CV",
    Evaluation == "test" ~ "Prob. test\nsamples",
    Evaluation == "surface" ~ "True\naccuracy"),
    Evaluation = as.factor(Evaluation),
    Evaluation = fct_relevel(Evaluation, "Prob. test\nsamples"),
    Evaluation = fct_relevel(Evaluation, "True\naccuracy")) |>
  mutate(Model = strsplit(stat, "_")[[1]][1]) |>
  ungroup() |>
  mutate(dsample = gsub("\n", " ", dsample, fixed=TRUE),
         dsample = paste0("Sampling pattern:\n", dsample),
         dsample = fct_inorder(dsample)) |>
  filter(Model != "Baseline")
cols <- c("Coordinates"="#D95F02",
          "EDF"="#7570B3",
          "RFsp"= "#E7298A")
p <- ggplot(data = figdata) +
  geom_boxplot(aes(y=RMSE, x=Evaluation, col=Model), outlier.size=0.5, outlier.alpha=0.75) +
  scale_color_manual(values = cols) +
  facet_grid(range ~ dsample) +
  ylim(1, 4) +
  xlab("Evaluation method") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 7))
ggsave(p, filename = "figures/sims_valinter_proxyonly.pdf", width = 9, height = 5)
rm("figdata", "p")


# Main: RFGLS interpolation ----
figdata <- dplyr::select(res, range, dsample, scenario, contains("_surface_inter")) |>
  pivot_longer(contains("_surface_"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1],
         Model = as.factor(ifelse(Model == "RFGLS", "RF-GLS", Model)),
         Model = fct_relevel(Model, "RF-GLS", after = 5)) |>
  ungroup() |>
  filter(!(Model=="Baseline"&scenario=="Scenario:\nProxies only"))
figdata_others <- figdata[!grepl("RFGLS", figdata$stat),] |>
  group_by(range, dsample, scenario, Model) |>
  summarise(RMSE = median(RMSE)) |>
  group_by(range, dsample, scenario) |>
  filter(RMSE == min(RMSE)) |>
  ungroup()
figdata_RFGLS <- data.frame()
for(i in 1:nrow(figdata_others)){
  figdata_RFGLS <- rbind(figdata_RFGLS,
                         filter(figdata,
                                range == figdata_others$range[i],
                                dsample == figdata_others$dsample[i],
                                scenario == figdata_others$scenario[i],
                                Model == figdata_others$Model[i]))
}
figdata_RFGLS$Model <- "Best-performing standard RF"
figdata_RFGLS <- rbind(figdata_RFGLS, filter(figdata, Model == "RF-GLS"))

# Plot
min_error_line <- rbind(
  data.frame(scenario = c("Scenario:\nComplete",
                          "Scenario:\nMissing predictors",
                          "Scenario:\nProxies only"), int=1),
  data.frame(scenario = c("Scenario:\nAutocorrelated error"), int=0))
p <- ggplot(figdata_RFGLS) +
  geom_hline(data = min_error_line, aes(yintercept = int), lty=2, col="grey20") +
  geom_boxplot(aes(y = RMSE, x = dsample, col = Model), outlier.size=0.5, outlier.alpha=0.75) +
  scale_color_manual(values = c("grey60", "black")) +
  facet_grid(range ~ scenario) +
  labs(colour = "Model") +
  xlab("Sampling pattern") + ylab("RMSE") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 7))
ggsave(p, filename = "figures/sims_RFGLSinter.pdf", width = 9, height = 5)
rm("figdata", "p")

# Appendix: RFGLS extrapolation ----
figdata <- dplyr::select(res, range, dsample, scenario, contains("_surface_extra")) |>
  pivot_longer(contains("_surface_"), names_to = "stat", values_to = "RMSE") |>
  group_by(stat) |>
  mutate(Model = strsplit(stat, "_")[[1]][1],
         Model = as.factor(ifelse(Model == "RFGLS", "RF-GLS", Model)),
         Model = fct_relevel(Model, "RF-GLS", after = 5)) |>
  ungroup() |>
  filter(!(Model=="Baseline"&scenario=="Scenario:\nProxies only"))
figdata_others <- figdata[!grepl("RFGLS", figdata$stat),] |>
  group_by(range, dsample, scenario, Model) |>
  summarise(RMSE = median(RMSE)) |>
  group_by(range, dsample, scenario) |>
  filter(RMSE == min(RMSE)) |>
  ungroup()
figdata_RFGLS <- data.frame()
for(i in 1:nrow(figdata_others)){
  figdata_RFGLS <- rbind(figdata_RFGLS,
                         filter(figdata,
                                range == figdata_others$range[i],
                                dsample == figdata_others$dsample[i],
                                scenario == figdata_others$scenario[i],
                                Model == figdata_others$Model[i]))
}
figdata_RFGLS$Model <- "Best-performing standard RF"
figdata_RFGLS <- rbind(figdata_RFGLS, filter(figdata, Model == "RF-GLS"))

# Plot
min_error_line <- rbind(
  data.frame(scenario = c("Scenario:\nComplete",
                          "Scenario:\nMissing predictors",
                          "Scenario:\nProxies only"), int=1),
  data.frame(scenario = c("Scenario:\nAutocorrelated error"), int=0))
p <- ggplot(figdata_RFGLS) +
  geom_hline(data = min_error_line, aes(yintercept = int), lty=2, col="grey20") +
  geom_boxplot(aes(y = RMSE, x = dsample, col = Model), outlier.size=0.5, outlier.alpha=0.75) +
  scale_color_manual(values = c("grey60", "black")) +
  facet_grid(range ~ scenario) +
  ylim(0, 5) +
  labs(colour = "Model") +
  xlab("Sampling pattern") + ylab("RMSE") +
  theme_bw() +
  theme(legend.position = "bottom", axis.text = element_text(size = 7))
ggsave(p, filename = "figures/sims_RFGLSextra.pdf", width = 9, height = 5)
rm("figdata", "p")
