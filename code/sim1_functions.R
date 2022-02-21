# Load packages
library("gstat")
library("caret")
library("NNDM")
library("sf")
library("raster")
library("dplyr")
library("tidyr")

# Source utils and
source("code/sim_utils.R")
options(dplyr.summarise.inform = FALSE)

#' Create an outcome according to van der Laan (2007) simulation framework
#' adapted to spatial prediction problems.
#' @details
#' 20 300x100 spherical stationary random fields with mean 0, sill 1,
#' nugget 0 and varying range are simulated and added according to van der Laan
#' (2007) equation. Random (var=1) and spatial (same parameters as covariates)
#' noise are added to the outcome.
#' @param cov_stack Raster stack of 20 elements names cov1, cov2...
#' @param snoise Raster with the spatial autocorrelation noise to be added.
van_der_laan <- function(cov_stack, snoise){

  # Continuous outcome using van der Laan's formula
  out <- cov_stack$cov1*cov_stack$cov2 + cov_stack$cov10^2 -
    cov_stack$cov3*cov_stack$cov17 - cov_stack$cov15*cov_stack$cov4 +
    cov_stack$cov9*cov_stack$cov5 + cov_stack$cov19 -
    cov_stack$cov20^2 + cov_stack$cov9 * cov_stack$cov8

  # Prepare random noise
  rnoise <- raster(ncols=300, nrows=100, xmn=0, xmx=300, ymn=0, ymx=100)
  vals <- rnorm(100*300, sd=1)
  rnoise <- setValues(rnoise, vals)

  # Add noise variables
  out <- out + rnoise + snoise
  names(out)  <- "outcome"

  return(out)
}


#' Sample simulation 1: random fields.
#' @details
#' Simulates a series of sampling points for simulation problem 1.
#' @param nsamples Integer. Number of samples to simulate.
#' @param dsamples Character. Spatial distribution of the samples. 5 are
#' possible: sregular, wregular, random, wclust, sclust.
#' @param sarea sf/sfc polygon where samples will be simulated.
sim1_samples <- function(nsamples, dsamples, sarea){

  if(dsamples=="sregular"){
    simpoints <- jitterreg_sample(sarea, nsamples, 2)
  }else if(dsamples=="wregular"){
    simpoints <- jitterreg_sample(sarea, nsamples, 5)
  }else if(dsamples=="random"){
    simpoints <- st_sample(sarea, nsamples)
  }else if(dsamples=="wclust"){
    simpoints <- clustered_sample(sarea, nsamples, 25, 5)
  }else if(dsamples=="sclust"){
    simpoints <- clustered_sample(sarea, nsamples, 10, 5)
  }

  simpoints <- st_sf(geometry=simpoints)
  simpoints
}


#' Fits a RF model and evaluates (random field simulation) it using several
#' methods.
#' @details
#' Fits a RF model and evaluates it using LOO CV, bLOO CV, NNDM LOO CV (interpolation
#' and extrapolation), and also returns the true interpolation and extrapolation
#' errors.
#' @param form String. Model formula.
#' @param range_indexTrain list. Indices for bLOO CV training data.
#' @param range_indexTest list. Indices for bLOO CV test data.
#' @param ndminter_indexTrain list. Indices for NNDM LOO CV (interpolation) training
#' data.
#' @param ndminter_indexTest list. Indices for NNDM LOO CV (interpolation) test data.
#' @param ndmextra_indexTrain list. Indices for NNDM LOO CV (extrapolation) training
#' data.
#' @param ndmextra_indexTest list. Indices for NNDM LOO CV (extrapolation) test data.
#' @param pgrid Data frame. Parameter grid of the model.
#' @param traindf Data frame. Training data to fit the model.
#' @param interdf Data frame. Interpolation surface data.
#' @param extradf Data frame. Extrapolation surface data.
fitval_rf <- function(form,
                      range_indexTrain, range_indexTest,
                      ndminter_indexTrain, ndminter_indexTest,
                      ndmextra_indexTrain, ndmextra_indexTest,
                      pgrid, traindf, interdf, extradf){

  # Validate with LOO and compute metrics
  loo_cntrl <- trainControl(method="LOOCV", savePredictions=TRUE,
                            allowParallel = TRUE)
  loo_mod <- train(form, data=traindf, method="rf",
                   trControl=loo_cntrl, tuneGrid=pgrid, ntree=100)
  loo_stats <- loo_mod$pred %>%
    summarise(RMSE = sqrt(mean((obs-pred)^2)),
              MAE = mean(abs(obs-pred)),
              R2 = cor(obs, pred)^2)
  names(loo_stats) <- paste0(names(loo_stats), "_LOO")

  # Compute CV statistics in interpolation surface
  interdf$preds <- predict(loo_mod, newdata=interdf)
  inter_stats <- interdf %>%
    summarise(RMSE = sqrt(mean((outcome-preds)^2)),
              MAE = mean(abs(outcome-preds)),
              R2 = cor(outcome, preds)^2)
  names(inter_stats) <- paste0(names(inter_stats), "_inter")

  # Compute CV statistics in extrapolation surface
  extradf$preds <- predict(loo_mod, newdata=extradf)
  extra_stats <- extradf %>%
    summarise(RMSE = sqrt(mean((outcome-preds)^2)),
              MAE = mean(abs(outcome-preds)),
              R2 = cor(outcome, preds)^2)
  names(extra_stats) <- paste0(names(extra_stats), "_extra")

  # Validate with bLOO and compute CV statistics
  bloo_cntrl <- trainControl(method="cv", index=range_indexTrain,
                             indexOut=range_indexTest,
                             savePredictions=TRUE, allowParallel = TRUE)
  bloo_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=bloo_cntrl, tuneGrid=pgrid, ntree=100))
  bloo_stats <- bloo_mod$pred %>%
    summarise(RMSE = sqrt(mean((obs-pred)^2)),
              MAE = mean(abs(obs-pred)),
              R2 = cor(obs, pred)^2)
  names(bloo_stats) <- paste0(names(bloo_stats), "_bLOO")

  # Validate with NDMinter and compute CV statistics
  ndminter_cntrl <- trainControl(method="cv",
                                 index=ndminter_indexTrain,
                                 indexOut=ndminter_indexTest,
                                 savePredictions=TRUE,
                                 allowParallel = TRUE)
  ndminter_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=ndminter_cntrl, tuneGrid=pgrid, ntree=100))
  ndminter_stats <- ndminter_mod$pred %>%
    summarise(RMSE = sqrt(mean((obs-pred)^2)),
              MAE = mean(abs(obs-pred)),
              R2 = cor(obs, pred)^2)
  names(ndminter_stats) <- paste0(names(ndminter_stats), "_NDMinter")

  # Validate with NDMextra and compute CV statistics
  ndmextra_cntrl <- trainControl(method="cv",
                                 index=ndmextra_indexTrain,
                                 indexOut=ndmextra_indexTest,
                                 savePredictions=TRUE,
                                 allowParallel = TRUE)
  ndmextra_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=ndmextra_cntrl, tuneGrid=pgrid, ntree=100))
  ndmextra_stats <- ndmextra_mod$pred %>%
    summarise(RMSE = sqrt(mean((obs-pred)^2)),
              MAE = mean(abs(obs-pred)),
              R2 = cor(obs, pred)^2)
  names(ndmextra_stats) <- paste0(names(ndmextra_stats), "_NDMextra")

  # Tidy and return results
  data.frame(loo_stats, bloo_stats, ndminter_stats, ndmextra_stats,
             inter_stats, extra_stats)
}


#' Simulation function 1: random fields.
#' @details
#' Simulates a landscape of 300x100 cells and sampling points in the training
#' area 100x100, fits a RF and evaluates it using LOO, bLOO, NNDM LOO CV; and surface
#' error (interpolation and extrapolation).
#' @param range Numeric. Autocorrelation range.
#' @param n_train Integer or vector of integers. Number of samples to simulate.
#' @param sample_dist String or vector of strings. Distribution of the sampling
#' points. 5 are possible: "sregular", "wregular", "random", "wclust","sclust".
sim_fields <- function(range,
                       n_train=c(100,200,300),
                       sample_dist=c("sregular", "wregular", "random",
                                     "wclust","sclust")){

  # Create an empty results object
  res<-data.frame()

  # Create grids (raster and point format) and sampling area
  rast_grid <- raster(ncols=300, nrows=100, xmn=0, xmx=300, ymn=0, ymx=100)
  point_grid <- st_as_sf(rasterToPoints(rast_grid, spatial = TRUE))
  inter_area <- matrix(c(0,0,100,0,100,100,0,100,0,0), ncol=2, byrow=TRUE)
  inter_area <- st_sfc(st_polygon(list(inter_area)))
  inter_grid <- point_grid[st_intersects(point_grid, inter_area, sparse=FALSE),]
  extra_area <- matrix(c(200,0,300,0,300,100,200,100,200,0), ncol=2, byrow=TRUE)
  extra_area <- st_sfc(st_polygon(list(extra_area)))
  extra_grid <- point_grid[st_intersects(point_grid, extra_area, sparse=FALSE),]

  # Simulate 20 covariates and a noise field from a semivariogram and stack
  cov_mod <- vgm(model="Sph", psill=1, range=range, nugget=0)
  cov_mod <- gstat(formula=z~1, dummy=TRUE, beta=0, model=cov_mod, nmax=100)
  cov_points <- quiet(predict(cov_mod, point_grid, nsim=21))
  cov_stack <- rasterise_stack(cov_points, 1:20, paste0("cov", 1:20))
  noise_stack <- rasterise_stack(cov_points, 21, "snoise")

  # Generate continuous outcome, add to stack
  out_rast <- van_der_laan(cov_stack, noise_stack)
  all_stack <- stack(out_rast, cov_stack)

  # Go to next level - combinations of sample number and distribution
  train_grid <- expand.grid(n_train=n_train, sample_dist=sample_dist,
                            stringsAsFactors = FALSE)
  for(i in 1:nrow(train_grid)){

    # Fetch the indicators of sample number and distribution for the iteration
    dist_it <- train_grid$sample_dist[i]
    n_it <- train_grid$n_train[i]

    # Simulate sampling points according to parameters and constraints
    train_points <- sim1_samples(n_it, dist_it, inter_area)

    # Define folds based on landscape autocorrelation range (bLOO_range)
    folds_bLOO <- bLOO(train_points, range, 0.5)

    # Define folds based on NNDM - interpolation
    folds_ndminter <- suppressMessages(
      nndm(train_points, inter_grid, range, 0.5))

    # Define folds based on NNDM - extrapolation
    folds_ndmextra <- suppressMessages(
      nndm(train_points, extra_grid, range, 0.5))

    # Get training and surface data for modelling
    train_data <- as.data.frame(raster::extract(all_stack, train_points))
    inter_data <- as.data.frame(raster::extract(all_stack, inter_grid))
    extra_data <- as.data.frame(raster::extract(all_stack, extra_grid))

    # Define formula and hyperparameter grid for model tuning
    form <- as.formula(paste0("outcome~", paste0("cov", 1:20, collapse="+")))
    pgrid <- data.frame(mtry=round(20/3, 0))

    # Model fitting and validation
    mod <- fitval_rf(form,
                     folds_bLOO$indx_train, folds_bLOO$indx_test,
                     folds_ndminter$indx_train, folds_ndminter$indx_test,
                     folds_ndmextra$indx_train, folds_ndmextra$indx_test,
                     pgrid, train_data, inter_data, extra_data)

    # Store results of the iteration
    res_it <- cbind(data.frame(range=range, nsample=n_it, dsample=dist_it,
                               stringsAsFactors = FALSE),
                    mod)
    res <- bind_rows(res, res_it)

  } # End loop sampling

  row.names(res) <- NULL
  res
}
