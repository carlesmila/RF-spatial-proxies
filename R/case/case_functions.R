#------------------------------------------------------------------------------#
#                      Function to fit case study models                       #
#------------------------------------------------------------------------------#

fitvalpred_rf <- function(covariates, proxies=NULL,
                          random_ctrl, spatial_ctrl,
                          traindf, rstack){

  # 1. AOA
  tune_grid <- round(seq(1, length(c(covariates, proxies)), length.out=10))
  tune_grid <- data.frame(mtry = tune_grid[!duplicated(tune_grid)],
                          splitrule="variance",
                          min.node.size=5)
  tune_ctrl <- trainControl(method="oob")
  tune_mod <- caret::train(traindf[c(covariates, proxies)], traindf[,"outcome"],
                           method="ranger", importance='impurity',
                           trControl=tune_ctrl, num.trees=300, tuneGrid=tune_grid)
  AOA <- suppressMessages(aoa(rstack, tune_mod))
  AOA <- as.numeric(global(AOA$AOA, na.rm=TRUE))
  names(AOA) <- "AOA"

  # 2. random CV
  random_grid <- data.frame(mtry = tune_mod$bestTune$mtry,
                            splitrule="variance",
                            min.node.size=5)
  random_mod <- train(traindf[c(covariates, proxies)], traindf[,"outcome"],
                      method="ranger", importance='impurity',
                      trControl=random_ctrl, num.trees = 300,  tuneGrid=random_grid)
  random_stats_mean <- sapply(random_mod$resample[c("RMSE", "Rsquared")], mean)
  names(random_stats_mean) <- paste0(names(random_stats_mean), "_mean")
  random_stats_sd <- sapply(random_mod$resample[c("RMSE", "Rsquared")], sd)
  names(random_stats_sd) <- paste0(names(random_stats_sd), "_sd")
  random_stats <- c(random_stats_mean, random_stats_sd)
  names(random_stats) <- paste0("random_", names(random_stats))

  # 3. kNNDM CV
  spatial_grid <- data.frame(mtry = tune_mod$bestTune$mtry,
                             splitrule="variance",
                             min.node.size=5)
  spatial_mod <- train(traindf[c(covariates, proxies)], traindf[,"outcome"],
                       method="ranger", importance='none',
                       trControl=spatial_ctrl, num.trees = 300, tuneGrid=spatial_grid)
  spatial_stats <- global_validation(spatial_mod)[c("RMSE", "Rsquared")]
  names(spatial_stats) <- paste0("kNNDM_", names(spatial_stats))

  # 4. Surface predictions
  preds <- predict(rstack, spatial_mod, na.rm=TRUE)

  # 5. Variable importance
  impfeat <- importance(random_mod$finalModel, type = 2)
  impfeat <- sum(impfeat[names(impfeat) %in% covariates])/sum(impfeat)*100
  names(impfeat) <- "impfeat"

  # Tidy and return results
  tabres <- as.data.frame(t(c(random_stats, spatial_stats, AOA, impfeat)))
  names(preds) <- c("prediction")
  list(tab = tabres, preds = preds, mod = tune_mod)
}
