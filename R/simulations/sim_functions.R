fitval_rf <- function(covariates,
                      randomfolds, interfolds, extrafolds,
                      traindf, testinterdf, testextradf,
                      gridinter, gridextra){


  # 1. tune model + surface accuracy + test data
  tune_grid <- round(seq(2, length(covariates), length.out=5))
  tune_ctrl <- trainControl(method="oob")
  tune_grid <- data.frame(mtry = tune_grid[!duplicated(tune_grid)],
                          splitrule="variance",
                          min.node.size=5)
  tune_mod <- caret::train(traindf[covariates], traindf[,"outcome"], method="ranger", importance='impurity',
                           trControl=tune_ctrl, num.trees=100, tuneGrid=tune_grid)
  surface_inter <- postResample(pred = predict(tune_mod, newdata=gridinter), obs = gridinter$outcome)
  surface_inter <- surface_inter["RMSE"]
  surface_extra <- postResample(pred = predict(tune_mod, newdata=gridextra), obs = gridextra$outcome)
  surface_extra <- surface_extra["RMSE"]
  test_inter <- postResample(pred = predict(tune_mod, newdata=testinterdf), obs = testinterdf$outcome)
  test_inter <- test_inter["RMSE"]
  test_extra <- postResample(pred = predict(tune_mod, newdata=testextradf), obs = testextradf$outcome)
  test_extra <- test_extra["RMSE"]

  # 2. random CV
  random_indx <- CreateSpacetimeFolds(data.frame(ID = randomfolds), spacevar = "ID", k = 5)
  random_grid <- data.frame(mtry = tune_mod$bestTune$mtry,
                            splitrule="variance",
                            min.node.size=5)
  random_ctrl <- trainControl(method="cv", index = random_indx$index, savePredictions='final')
  random_mod <- caret::train(traindf[covariates], traindf[,"outcome"], method="ranger", importance='none',
                             trControl=random_ctrl, num.trees=100, tuneGrid=random_grid)
  random_inter <- global_validation(random_mod)["RMSE"]

  # 3. kNNDM CV - interpolation
  spatial_indx1 <- CreateSpacetimeFolds(data.frame(ID = interfolds), spacevar = "ID", k = 5)
  spatial_grid1 <- data.frame(mtry = tune_mod$bestTune$mtry,
                              splitrule="variance",
                              min.node.size=5)
  spatial_ctrl1 <- trainControl(method="cv", index = spatial_indx1$index, savePredictions='final')
  spatial_mod1 <- caret::train(traindf[covariates], traindf[,"outcome"], method="ranger", importance='none',
                               trControl=spatial_ctrl1, num.trees=100, tuneGrid=spatial_grid1)
  spatial_inter <- global_validation(spatial_mod1)["RMSE"]

  # 4. kNNDM CV - extrapolation
  spatial_indx2 <- CreateSpacetimeFolds(data.frame(ID = extrafolds), spacevar = "ID", k = 5)
  spatial_grid2 <- data.frame(mtry = tune_mod$bestTune$mtry,
                              splitrule="variance",
                              min.node.size=5)
  spatial_ctrl2 <- trainControl(method="cv", index = spatial_indx2$index, savePredictions='final')
  spatial_mod2 <- caret::train(traindf[covariates], traindf[,"outcome"], method="ranger", importance='none',
                               trControl=spatial_ctrl2, num.trees=100, tuneGrid=spatial_grid2)
  spatial_extra <- global_validation(spatial_mod2)["RMSE"]

  # 5. AOA
  AOA_inter <- suppressMessages(aoa(gridinter, tune_mod, verbose = FALSE))
  AOA_inter <- sum(AOA_inter$AOA==1)/length(AOA_inter$AOA)*100
  AOA_extra <- suppressMessages(aoa(gridextra, tune_mod, verbose = FALSE))
  AOA_extra <- sum(AOA_extra$AOA==1)/length(AOA_extra$AOA)*100

  # 6. varImp
  if(all(covariates %in% paste0("cov", 1:6))){
    impfeat <- 100
  }else{
    impfeat <- importance(tune_mod$finalModel, type = 2)
    impfeat <- sum(impfeat[names(impfeat) %in% paste0("cov", 1:6)])/sum(impfeat)*100
  }

  # 7. Tidy and return results
  resfit <- data.frame(surface_inter=surface_inter,
                       test_inter=test_inter,
                       random_inter=random_inter,
                       spatial_inter=spatial_inter,
                       AOA_inter=AOA_inter,
                       surface_extra=surface_extra,
                       test_extra=test_extra,
                       spatial_extra=spatial_extra,
                       AOA_extra=AOA_extra,
                       impfeat=impfeat,
                       mtry=tune_mod$bestTune$mtry)
  return(resfit)
}




sim_fields <- function(seed,
                       range,
                       n_train=200,
                       sample_dist=c("sregular", "random", "wclust", "sclust"),
                       scenario){

  # Print progress
  print(seed)
  set.seed(seed)

  # Create an empty results object
  res <- data.frame()

  # Create grids (raster and point format) and sampling and extrapolation area
  rast_grid <- raster(ncols=300, nrows=100, xmn=0, xmx=300, ymn=0, ymx=100)
  point_grid <- st_as_sf(rasterToPoints(rast_grid, spatial = TRUE))
  sampling_area <- matrix(c(0,0,100,0,100,100,0,100,0,0), ncol=2, byrow=TRUE)
  sampling_area <- st_sfc(st_polygon(list(sampling_area)))
  extra_area <- matrix(c(200,0,300,0,300,100,200,100,200,0), ncol=2, byrow=TRUE)
  extra_area <- st_sfc(st_polygon(list(extra_area)))

  # Simulate 6 covariates from a semivariogram and stack
  cov_mod <- vgm(model="Sph", psill=1, range=range, nugget=0)
  cov_mod <- gstat(formula=z~1, dummy=TRUE, beta=0, model=cov_mod, nmax=100)
  cov_points <- quiet(predict(cov_mod, point_grid, nsim=6))
  cov_stack <- rasterise_stack(cov_points, 1:6, paste0("cov", 1:6))
  rm("cov_mod")

  # Generate outcome
  out_rast <- cov_stack$cov1 + cov_stack$cov2*cov_stack$cov3 +
    cov_stack$cov4 + cov_stack$cov5*cov_stack$cov6

  # If error is autocorrelated
  if(scenario == "autocor"){
    snoise_mod <- vgm(model="Sph", psill=1, range=25, nugget=0)
    snoise_mod <- gstat(formula=z~1, dummy=TRUE, beta=0, model=snoise_mod, nmax=100)
    snoise_points <- quiet(predict(snoise_mod, point_grid, nsim=1))
    snoise <- rasterise_stack(snoise_points, 1, "error")
    out_rast <- out_rast + snoise
    names(out_rast)  <- "outcome"
    all_stack <- stack(out_rast, cov_stack)
  }else{
    # otherwise add random noise
    rnoise <- raster(ncols=300, nrows=100, xmn=0, xmx=300, ymn=0, ymx=100)
    vals <- rnorm(300*100, sd=1)
    rnoise <- setValues(rnoise, vals)
    out_rast <- out_rast + rnoise
    names(out_rast)  <- "outcome"
    all_stack <- stack(out_rast, cov_stack)
  }

  # Generate spatial proxies: coordinates
  proxies <- point_grid
  proxies$x <- sf::st_coordinates(point_grid)[,1]
  proxies$y <- sf::st_coordinates(point_grid)[,2]

  # Generate spatial proxies: EDF
  proxies$EDF1 <- st_distance(proxies, st_sfc(st_point(c(0,100))))
  proxies$EDF2 <- st_distance(proxies, st_sfc(st_point(c(300,100))))
  proxies$EDF3 <- st_distance(proxies, st_sfc(st_point(c(0,0))))
  proxies$EDF4 <- st_distance(proxies, st_sfc(st_point(c(300,0))))
  proxies$EDF5 <- st_distance(proxies, st_sfc(st_point(c(150,50))))
  proxies_stack <- rasterise_stack(proxies, 2:8, names(proxies)[2:8])
  all_stack <- stack(all_stack, proxies_stack)

  # Covariates to be included
  if(scenario %in% c("complete", "autocor")){
    covlist <- paste0("cov", 1:6)
  }else if(scenario == "partial"){
    covlist <- paste0("cov", 1:3)
  }else if(scenario == "proxyonly"){
    covlist <- NULL
  }

  # Probability sampling for testing
  test_inter <- sim1_samples(100, "random", sampling_area)
  test_extra <- sim1_samples(100, "random", extra_area)

  # Go to next level - by sampling distribution
  for(dist_it in sample_dist){

    # Simulate sampling points according to parameters and constraints
    train_points <- sim1_samples(n_train, dist_it, sampling_area)

    # Generate spatial proxies: RFsp
    proxies <- point_grid
    proxies_dist <- st_distance(point_grid, train_points)
    proxies <- cbind(proxies_dist, proxies)
    proxies_stack <- rasterise_stack(proxies, which(sapply(proxies, is.numeric)),
                                     paste0("RFsp", 1:(ncol(proxies)-1)))
    all_stack2 <- stack(all_stack, proxies_stack)

    # Get training (delete duplicates), surface and test data
    train_data <- as.data.frame(raster::extract(all_stack2, train_points))
    train_points <- train_points[!duplicated(train_data),]
    train_data <- as.data.frame(raster::extract(all_stack2, train_points))
    gridinter_data <- as.data.frame(raster::extract(all_stack2, st_sf(geom=sampling_area)))
    gridextra_data <- as.data.frame(raster::extract(all_stack2, st_sf(geom=extra_area)))
    testinter_data <- as.data.frame(raster::extract(all_stack2, test_inter))
    testextra_data <- as.data.frame(raster::extract(all_stack2, test_extra))

    # Define folds - suppress CRS and randomCV messages in knndm (for sampling area)
    rfolds <- sample(1:nrow(train_points)%%5, nrow(train_points))
    sfolds_inter <- suppressMessages(suppressWarnings( # Random assignment + missing CRS messages
      knndm(train_points,
            predpoints = st_intersection(point_grid, sampling_area),
            k = 5, maxp = 0.5)$clusters
    ))
    sfolds_extra <- suppressMessages(suppressWarnings( # Random assignment + missing CRS messages
      knndm(train_points,
            predpoints = st_intersection(point_grid, extra_area),
            k = 5, maxp = 0.5)$clusters
    ))

    # Model 1: baseline
    if(scenario != "proxyonly"){ # only when there are covariates
      mod_baseline <- fitval_rf(covlist,
                                rfolds, sfolds_inter, sfolds_extra,
                                train_data, testinter_data, testextra_data,
                                gridinter_data, gridextra_data)
      names(mod_baseline) <- paste0("baseline_", names(mod_baseline))
      baseline_mtry <- mod_baseline$baseline_mtry

    }else{
      meanpred <- mean(train_data$outcome)
      mod_baseline <- data.frame(surface_inter=sqrt(mean((meanpred-gridinter_data$outcome)^2)),
                                 test_inter=NA,
                                 random_inter=NA,
                                 spatial_inter=NA,
                                 AOA_inter=NA,
                                 surface_extra=sqrt(mean((meanpred-gridextra_data$outcome)^2)),
                                 test_extra=NA,
                                 spatial_extra=NA,
                                 AOA_extra=NA,
                                 impfeat=NA,
                                 mtry=NA)
      names(mod_baseline) <- paste0("baseline_", names(mod_baseline))
      baseline_mtry <- 1
    }

    # Model 2: Coordinates
    coordlist <- c(covlist, "x", "y")
    mod_coord <- fitval_rf(coordlist,
                           rfolds, sfolds_inter, sfolds_extra,
                           train_data, testinter_data, testextra_data,
                           gridinter_data, gridextra_data)
    names(mod_coord) <- paste0("coord_", names(mod_coord))

    # Model 3: EDF
    EDFlist <- c(covlist, "x", "y", paste0("EDF", as.character(1:5)))
    mod_EDF <- fitval_rf(EDFlist,
                         rfolds, sfolds_inter, sfolds_extra,
                         train_data, testinter_data, testextra_data,
                         gridinter_data, gridextra_data)
    names(mod_EDF) <- paste0("EDF_", names(mod_EDF))

    # Model 4: RFsp
    RFsplist <- c(covlist, paste0("RFsp", as.character(1:nrow(train_data))))
    mod_RFsp <- fitval_rf(RFsplist,
                          rfolds, sfolds_inter, sfolds_extra,
                          train_data, testinter_data, testextra_data,
                          gridinter_data, gridextra_data)
    names(mod_RFsp) <- paste0("RFsp_", names(mod_RFsp))

    # Model 5: RF-GLS
    if(scenario=="proxyonly"){
      X_RFGLS <- matrix(rep(1, nrow(train_data)), ncol = 1)
      predinter_RFGLS <- matrix(rep(1, nrow(gridinter_data)), ncol = 1)
      predextra_RFGLS <- matrix(rep(1, nrow(gridextra_data)), ncol = 1)
    }else{
      X_RFGLS <- as.matrix(train_data[covlist])
      predinter_RFGLS <- as.matrix(gridinter_data[,covlist])
      predextra_RFGLS <- as.matrix(gridextra_data[,covlist])
    }

    suppressWarnings(
      mod_RFGLS <- RFGLS_estimate_spatial(coords = as.matrix(train_data[c("x","y")]),
                                          y = train_data$outcome,
                                          X = X_RFGLS,
                                          ntree = 100,
                                          mtry = baseline_mtry,
                                          param_estimate = TRUE)
    )
    suppressWarnings(
      RFGLS_surface_inter <- RFGLS_predict_spatial(mod_RFGLS,
                                                   as.matrix(gridinter_data[,c("x", "y")]),
                                                   predinter_RFGLS)$prediction
    )
    suppressWarnings(
      RFGLS_surface_extra <- RFGLS_predict_spatial(mod_RFGLS,
                                                   as.matrix(gridextra_data[,c("x", "y")]),
                                                   predextra_RFGLS)$prediction
    )
    mod_RFGLS <- data.frame(
      RFGLS_surface_inter = sqrt(mean((RFGLS_surface_inter-gridinter_data$outcome)^2)),
      RFGLS_surface_extra = sqrt(mean((RFGLS_surface_extra-gridextra_data$outcome)^2))
    )

    # Store results of the iteration
    res_it <- cbind(data.frame(range = range,
                               nsample = n_train,
                               dsample = dist_it,
                               scenario = scenario,
                               stringsAsFactors = FALSE),
                    mod_baseline, mod_coord, mod_EDF, mod_RFsp, mod_RFGLS)
    res <- bind_rows(res, res_it)

    # Cleaning
    rm("train_points", "train_data", "gridinter_data", "gridextra_data",
       "testinter_data", "testextra_data",
       "coordlist", "EDFlist", "RFsplist",
       "mod_baseline", "mod_coord", "mod_EDF", "mod_RFsp", "mod_RFGLS",
       "X_RFGLS", "predinter_RFGLS", "predextra_RFGLS",
       "res_it")

  } # End loop sampling

  row.names(res) <- NULL
  res
}
