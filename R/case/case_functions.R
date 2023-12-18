#------------------------------------------------------------------------------#
#                      Function to fit case study models                       #
#------------------------------------------------------------------------------#

# covariates <- basecovs
# proxies <- proxycovs
# random_ctrl <- rctrl
# spatial_ctrl <- sctrl
# traindf <- temp_df
# rstack <- tempstack[[c(basecovs, proxycovs)]]

fitvalpred_rf <- function(covariates, proxies=NULL,
                          random_ctrl, spatial_ctrl,
                          traindf, rstack){

  # 1. AOA
  tune_grid <- round(seq(2, length(c(covariates, proxies)), length.out=10))
  tune_grid <- data.frame(mtry = tune_grid[!duplicated(tune_grid)])
  tune_ctrl <- trainControl(method="oob")
  tune_mod <- caret::train(traindf[c(covariates, proxies)], traindf[,"outcome"],
                           method="rf", importance=TRUE,
                           trControl=tune_ctrl, ntree=300, tuneGrid=tune_grid)
  AOA <- suppressMessages(aoa(rstack, tune_mod))
  AOA <- as.numeric(global(AOA$AOA, na.rm=TRUE))
  names(AOA) <- "AOA"

  # 2. random CV
  random_grid <- data.frame(mtry = tune_mod$bestTune$mtry)
  random_mod <- train(traindf[c(covariates, proxies)], traindf[,"outcome"],
                      method="rf", importance=FALSE,
                      trControl=random_ctrl, ntree = 300,  tuneGrid=random_grid)
  random_stats <- global_validation(random_mod)[c("RMSE", "Rsquared")]
  names(random_stats) <- paste0("random_", names(random_stats))

  # 3. kNNDM CV
  spatial_grid <- data.frame(mtry = tune_mod$bestTune$mtry)
  spatial_mod <- train(traindf[c(covariates, proxies)], traindf[,"outcome"],
                       method="rf", importance=FALSE,
                       trControl=spatial_ctrl, ntree = 300, tuneGrid=spatial_grid)
  spatial_stats <- global_validation(spatial_mod)[c("RMSE", "Rsquared")]
  names(spatial_stats) <- paste0("kNNDM_", names(spatial_stats))

  # 4. Surface predictions
  preds <- predict(rstack, spatial_mod, na.rm=TRUE)

  # 5. Variable importance
  impfeat <- importance(random_mod$finalModel, type = 2)
  impfeat <- sum(impfeat[row.names(impfeat) %in% covariates, 1])/sum(impfeat[,1])*100
  names(impfeat) <- "impfeat"

  # Tidy and return results
  tabres <- as.data.frame(t(c(random_stats, spatial_stats, AOA, impfeat)))
  names(preds) <- c("prediction")
  list(tab = tabres, preds = preds, mod = tune_mod)
}
