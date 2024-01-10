fitval_rf <- function(covariates,
                      randomfolds, spatialfolds,
                      traindf, griddf){


  # 1. tune model + surface accuracy
  tune_grid <- round(seq(2, length(covariates), length.out=5))
  tune_ctrl <- trainControl(method="oob")
  tune_grid <- data.frame(mtry = tune_grid[!duplicated(tune_grid)])
  tune_mod <- caret::train(traindf[covariates], traindf[,"outcome"], method="rf", importance=TRUE,
                             trControl=tune_ctrl, ntree=100, tuneGrid=tune_grid)
  surface_RMSE <- postResample(pred = predict(tune_mod, newdata=griddf), obs = griddf$outcome)
  surface_RMSE <- surface_RMSE["RMSE"]

  # 2. random CV
  random_indx <- CreateSpacetimeFolds(data.frame(ID = randomfolds), spacevar = "ID", k = 5)
  random_grid <- data.frame(mtry = tune_mod$bestTune$mtry)
  random_ctrl <- trainControl(method="cv", index = random_indx$index, savePredictions='final')
  random_mod <- caret::train(traindf[covariates], traindf[,"outcome"], method="rf", importance=FALSE,
                             trControl=random_ctrl, ntree=100, tuneGrid=random_grid)
  random_RMSE <- global_validation(random_mod)["RMSE"]

  # 3. spatial CV
  spatial_indx <- CreateSpacetimeFolds(data.frame(ID = spatialfolds), spacevar = "ID", k = 5)
  spatial_grid <- data.frame(mtry = tune_mod$bestTune$mtry)
  spatial_ctrl <- trainControl(method="cv", index = spatial_indx$index, savePredictions='final')
  spatial_mod <- caret::train(traindf[covariates], traindf[,"outcome"], method="rf", importance=FALSE,
                              trControl=spatial_ctrl, ntree=100, tuneGrid=spatial_grid)
  spatial_RMSE <- global_validation(spatial_mod)["RMSE"]

  # 4. AOA
  AOA <- suppressMessages(aoa(griddf, tune_mod))
  AOA <- sum(AOA$AOA==1)/length(AOA$AOA)*100

  # 5. varImp
  if(all(covariates %in% paste0("cov", 1:6))){
    impfeat <- 100
  }else{
    impfeat <- importance(tune_mod$finalModel, type = 2)
    impfeat <- sum(impfeat[row.names(impfeat) %in% paste0("cov", 1:6), 1])/sum(impfeat[,1])*100
  }

  # Tidy and return results
  resfit <- data.frame(random_RMSE=random_RMSE,
                       spatial_RMSE=spatial_RMSE,
                       surface_RMSE=surface_RMSE,
                       AOA=AOA,
                       impfeat=impfeat,
                       mtry=tune_mod$bestTune$mtry)
  return(resfit)
}

sim_fields <- function(range, n_train=200,
                       sample_dist=c("sregular", "random", "wclust", "sclust"),
                       simtype){

  # Create an empty results object
  res <- data.frame()

  # Create grids (raster and point format) and sampling area
  rast_grid <- raster(ncols=100, nrows=100, xmn=0, xmx=100, ymn=0, ymx=100)
  point_grid <- st_as_sf(rasterToPoints(rast_grid, spatial = TRUE))
  study_area <- matrix(c(0,0,100,0,100,100,0,100,0,0), ncol=2, byrow=TRUE)
  study_area <- st_sfc(st_polygon(list(study_area)))

  # Simulate 6 covariates from a semivariogram and stack
  cov_mod <- vgm(model="Sph", psill=1, range=range, nugget=0)
  cov_mod <- gstat(formula=z~1, dummy=TRUE, beta=0, model=cov_mod, nmax=100)
  cov_points <- quiet(predict(cov_mod, point_grid, nsim=6))
  cov_stack <- rasterise_stack(cov_points, 1:6, paste0("cov", 1:6))
  rm("cov_mod")

  # Generate outcome
  out_rast <- cov_stack$cov1 + cov_stack$cov2*cov_stack$cov3 +
    cov_stack$cov4 + cov_stack$cov5*cov_stack$cov6

  # Prepare and add random noise
  rnoise <- raster(ncols=100, nrows=100, xmn=0, xmx=100, ymn=0, ymx=100)
  vals <- rnorm(100*100, sd=1)
  rnoise <- setValues(rnoise, vals)
  out_rast <- out_rast + rnoise
  names(out_rast)  <- "outcome"
  all_stack <- stack(out_rast, cov_stack)

  # Generate spatial proxies: coordinates
  proxies <- point_grid
  proxies$x <- sf::st_coordinates(point_grid)[,1]
  proxies$y <- sf::st_coordinates(point_grid)[,2]

  # Generate spatial proxies: EDF
  proxies$EDF1 <- st_distance(proxies, st_sfc(st_point(c(0,100))))
  proxies$EDF2 <- st_distance(proxies, st_sfc(st_point(c(100,100))))
  proxies$EDF3 <- st_distance(proxies, st_sfc(st_point(c(0,0))))
  proxies$EDF4 <- st_distance(proxies, st_sfc(st_point(c(100,0))))
  proxies$EDF5 <- st_distance(proxies, st_sfc(st_point(c(50,50))))
  proxies_stack <- rasterise_stack(proxies, 2:8, names(proxies)[2:8])
  all_stack <- stack(all_stack, proxies_stack)

  # Covariates to be included
  if(simtype == "complete"){
    covlist <- paste0("cov", 1:6)
  }else if(simtype == "partial"){
    covlist <- paste0("cov", 1:3)
  }

  # Go to next level - combinations of samcovlistple number and distribution
  train_grid <- expand.grid(n_train=n_train, sample_dist=sample_dist,
                            stringsAsFactors = FALSE)
  for(i in 1:nrow(train_grid)){

    # Fetch the indicators of sample number and distribution for the iteration
    dist_it <- train_grid$sample_dist[i]
    n_it <- train_grid$n_train[i]

    # Simulate sampling points according to parameters and constraints
    train_points <- sim1_samples(n_it, dist_it, study_area)

    # Generate spatial proxies: RFsp
    proxies <- point_grid
    proxies_dist <- st_distance(point_grid, train_points)
    proxies <- cbind(proxies_dist, proxies)
    proxies_stack <- rasterise_stack(proxies, which(sapply(proxies, is.numeric)),
                                     paste0("RFsp", 1:(ncol(proxies)-1)))
    all_stack2 <- stack(all_stack, proxies_stack)

    # Get training (delete duplicates) and surface data for modelling
    train_data <- as.data.frame(raster::extract(all_stack2, train_points))
    train_points <- train_points[!duplicated(train_data),]
    train_data <- as.data.frame(raster::extract(all_stack2, train_points))
    grid_data <- as.data.frame(raster::extract(all_stack2, point_grid))

    # Define folds - suppress CRS and randomCV messages in knndm
    rfolds <- sample(1:nrow(train_points)%%5, nrow(train_points))
    sfolds <- suppressMessages(knndm(train_points, ppoints = point_grid, k = 5, maxp = 0.5)$clusters)

    # Model 1: No spatial proxies
    mod_baseline <- fitval_rf(covlist,
                              rfolds, sfolds,
                              train_data, grid_data)
    names(mod_baseline) <- paste0("baseline_", names(mod_baseline))

    # Model 2: Coordinates
    coordlist <- c(covlist, "x", "y")
    mod_coord <- fitval_rf(coordlist,
                           rfolds, sfolds,
                           train_data, grid_data)
    names(mod_coord) <- paste0("coord_", names(mod_coord))

    # Model 3: EDF
    EDFlist <- c(covlist, "x", "y", paste0("EDF", as.character(1:5)))
    mod_EDF <- fitval_rf(EDFlist,
                         rfolds, sfolds,
                         train_data, grid_data)
    names(mod_EDF) <- paste0("EDF_", names(mod_EDF))

    # Model 4: RFsp
    RFsplist <- c(covlist, paste0("RFsp", as.character(1:nrow(train_data))))
    mod_RFsp <- fitval_rf(RFsplist,
                          rfolds, sfolds,
                          train_data, grid_data)
    names(mod_RFsp) <- paste0("RFsp_", names(mod_RFsp))

    # Model 5: RF-GLS
    suppressWarnings(
      mod_RFGLS <- RFGLS_estimate_spatial(coords = as.matrix(train_data[c("x","y")]),
                                          y = train_data$outcome,
                                          X = as.matrix(train_data[covlist]),
                                          param_estimate = TRUE,
                                          ntree = 100,
                                          mtry = mod_baseline$baseline_mtry)
    )
    suppressWarnings(
      pred_RFGLS <- RFGLS_predict_spatial(mod_RFGLS,
                                          as.matrix(grid_data[,c("x", "y")]),
                                          as.matrix(grid_data[,covlist]))$prediction
    )
    mod_RFGLS <- data.frame("RFGLS_surface_RMSE" = sqrt(mean((pred_RFGLS-grid_data$outcome)^2)))

    # Store results of the iteration
    res_it <- cbind(data.frame(range = range,
                               nsample = n_it,
                               dsample = dist_it,
                               simtype = simtype,
                               stringsAsFactors = FALSE),
                    mod_baseline, mod_coord, mod_EDF, mod_RFsp, mod_RFGLS)
    res <- bind_rows(res, res_it)

  } # End loop sampling

  row.names(res) <- NULL
  res
}
