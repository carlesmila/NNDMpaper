# Libraries and utils ----
library("tidyverse")
library("raster")
library("NNDM")
library("sf")
library("caret")
library("gstat")
library("virtualspecies")
library("parallel")
library("doParallel")
library("pbapply")
source("code/sim_utils.R")

# No need for proj4 warnings
options("rgdal_show_exportToProj4_warnings"="none")


#' Sample simulation 2: virtual species.
#' @details
#' Simulates a series of sampling points for simulation problem 2.
#' @param nsamples Integer. Number of samples to simulate.
#' @param dsamples Character. Spatial distribution of the samples. 5 are
#' possible: sregular, wregular, random, wclust, sclust.
#' @param sarea sf/sfc polygon where samples will be simulated.
sim2_samples <- function(nsamples, dsamples, sarea){

  if(dsamples=="sregular"){
    simpoints <- jitterreg_sample(sarea, nsamples, 40000)
  }else if(dsamples=="wregular"){
    simpoints <- jitterreg_sample(sarea, nsamples, 80000)
  }else if(dsamples=="random"){
    simpoints <- st_sample(sarea, nsamples)
  }else if(dsamples=="wclust"){
    simpoints <- clustered_sample(sarea, nsamples, 25, 80000)
  }else if(dsamples=="sclust"){
    simpoints <- clustered_sample(sarea, nsamples, 10, 80000)
  }

  simpoints <- st_sf(geometry=simpoints)
  simpoints
}


#' Fits a RF model and evaluates it using several methods (species simulation).
#' @details
#' Fits a RF model and evaluates it using LOO CV, bLOO CV, NNDM LOO CV and true errors.
#' @param form String. Model formula.
#' @param rangeout_indexTrain list. Indices for bLOO CV (outcome range)
#' training data.
#' @param rangeout_indexTest list. Indices for bLOO CV (outcome range) test
#' data.
#' @param rangeres_indexTrain list. Indices for bLOO CV (residual range)
#' training data.
#' @param rangeres_indexTest list. Indices for bLOO CV (residual range) test
#' data.
#' @param ndmout_indexTrain list. Indices for NNDM LOO CV (outcome range) training
#' data.
#' @param ndmout_indexTest list. Indices for NNDM LOO CV (outcome range) test data.
#' @param ndmres_indexTrain list. Indices for NNDM LOO CV (residual range) training
#' data.
#' @param ndmres_indexTest list. Indices for NNDM LOO CV (residual range) test data.
#' @param pgrid Data frame. Parameter grid of the model.
#' @param traindf Data frame. Training data to fit the model.
#' @param surfdf Data frame. Surface data.
fitval_rf_species <- function(form,
                              rangeout_indexTrain, rangeout_indexTest,
                              rangeres_indexTrain, rangeres_indexTest,
                              ndmout_indexTrain, ndmout_indexTest,
                              ndmres_indexTrain, ndmres_indexTest,
                              pgrid, traindf, surfdf){

  # Validate with LOO and compute metrics
  loo_cntrl <- trainControl(method="LOOCV", savePredictions=TRUE)
  loo_mod <- train(form, data=traindf, method="rf",
                   trControl=loo_cntrl, tuneGrid=pgrid, ntree=100)
  loo_stats <- loo_mod$pred %>%
    summarise(RMSE = sqrt(mean((obs-pred)^2)),
              MAE = mean(abs(obs-pred)),
              R2 = cor(obs, pred)^2)
  names(loo_stats) <- paste0(names(loo_stats), "_LOO")

  # Compute CV statistics in surface
  surfdf$preds <- predict(loo_mod, newdata=surfdf)
  surf_stats <- surfdf %>%
    summarise(RMSE = sqrt(mean((outcome-preds)^2)),
              MAE = mean(abs(outcome-preds)),
              R2 = cor(outcome, preds)^2)
  names(surf_stats) <- paste0(names(surf_stats), "_surf")

  # Validate with bLOO and compute CV statistics (outcome range)
  bloo_out_cntrl <- trainControl(method="cv", index=rangeout_indexTrain,
                                 indexOut=rangeout_indexTest,
                                 savePredictions=TRUE)
  bloo_out_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=bloo_out_cntrl, tuneGrid=pgrid, ntree=100))
  bloo_out_stats <- bloo_out_mod$pred %>%
    summarise(RMSE = sqrt(mean((obs-pred)^2)),
              MAE = mean(abs(obs-pred)),
              R2 = cor(obs, pred)^2)
  names(bloo_out_stats) <- paste0(names(bloo_out_stats), "_bLOO_out")

  # Validate with bLOO and compute CV statistics (residual range)
  bloo_res_cntrl <- trainControl(method="cv", index=rangeres_indexTrain,
                                 indexOut=rangeres_indexTest,
                                 savePredictions=TRUE, allowParallel = TRUE)
  bloo_res_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=bloo_res_cntrl, tuneGrid=pgrid, ntree=100))
  bloo_res_stats <- bloo_res_mod$pred %>%
    summarise(RMSE = sqrt(mean((obs-pred)^2)),
              MAE = mean(abs(obs-pred)),
              R2 = cor(obs, pred)^2)
  names(bloo_res_stats) <- paste0(names(bloo_res_stats), "_bLOO_res")

  # Validate with NNDM LOO and compute CV statistics (outcome range)
  ndm_out_cntrl <- trainControl(method="cv",
                                index=ndmout_indexTrain,
                                indexOut=ndmout_indexTest,
                                savePredictions=TRUE)
  ndm_out_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=ndm_out_cntrl, tuneGrid=pgrid, ntree=100))
  ndm_out_stats <- ndm_out_mod$pred %>%
    summarise(RMSE = sqrt(mean((obs-pred)^2)),
              MAE = mean(abs(obs-pred)),
              R2 = cor(obs, pred)^2)
  names(ndm_out_stats) <- paste0(names(ndm_out_stats), "_NDM_out")

  # Validate with NNDM LOO and compute CV statistics (residual range)
  ndm_res_cntrl <- trainControl(method="cv",
                                index=ndmres_indexTrain,
                                indexOut=ndmres_indexTest,
                                savePredictions=TRUE,
                                allowParallel = TRUE)
  ndm_res_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=ndm_res_cntrl, tuneGrid=pgrid, ntree=100))
  ndm_res_stats <- ndm_res_mod$pred %>%
    summarise(RMSE = sqrt(mean((obs-pred)^2)),
              MAE = mean(abs(obs-pred)),
              R2 = cor(obs, pred)^2)
  names(ndm_res_stats) <- paste0(names(ndm_res_stats), "_NDM_res")

  # Tidy and return results
  data.frame(loo_stats, surf_stats,
             bloo_out_stats, bloo_res_stats,
             ndm_out_stats, ndm_res_stats)
}


#' Simulation function 2: virtual species.
#' @details
#' The function takes a virtual species simulated landscape for the Iberian
#' peninsula using bioclim data, simulates sampling points, computes outcome
#' and residual autocorrelation range, and fits a RF and evaluates it using LOO,
#' bLOO, and NNDM LOO CV; as well as the true error.
#' @param rgrid sf or sfc point object. Prediction grid.
#' @param rstack stack raster object. Contains the landscape data and the
#' simulated outcome.
#' @param sampling_area sf or sfc polygon object. Sampling area.
#' @param sample_dist String or vector or string. Distribution of the sampling
#' points. 5 are possible: "sregular", "wregular", "random", "wclust","sclust".
sim_species <- function(rgrid, rstack, sampling_area,
                        sample_dist=c("sregular", "wregular", "random",
                                      "wclust","sclust")){

  # Initiate results object and fixed information for all models
  res <- data.frame()
  grid_data <- as.data.frame(raster::extract(rstack, rgrid))
  form <- as.formula(paste0("outcome~", paste0("bio", 1:19, collapse="+")))
  pgrid <- data.frame(mtry=6)

  # Start sampling loop
  for(dist_it in sample_dist){

    # Simulate sampling points according to parameters and constraints
    train_points <- sim2_samples(100, dist_it, sampling_area)

    # Get training and surface data for modelling and validation
    train_data <- as.data.frame(raster::extract(rstack, train_points))

    #### Folds based on outcome range
    # Estimate outcome range
    train_points$outcome <- train_data$outcome
    empvar_out <- variogram(as.formula("outcome~1"), cutoff=600000, train_points)
    fitvar_out <- quiet(suppressWarnings( # Convergence warnings
      fit.variogram(empvar_out, vgm(model="Sph", nugget=0), fit.sills=c(F,T), fit.method = 1)))
    lrange_out <- fitvar_out$range[fitvar_out$model=="Sph"]
    # Define folds based on bLOO
    folds_bLOO_out <- bLOO(train_points, lrange_out, 0.5)
    # Define folds based on NNDM
    folds_ndm_out <- nndm(train_points, rgrid, lrange_out, 0.5)

    #### Folds based on residual range
    # Estimate residual range
    base_cntrl <- trainControl(method="none")
    base_mod <- train(form, data=train_data, method="rf",
                      trControl=base_cntrl, tuneGrid=pgrid, ntree=100)
    train_points$res <- train_data$outcome - predict(base_mod)
    empvar_res <- variogram(as.formula("res~1"), cutoff=600000, train_points)
    fitvar_res <- quiet(suppressWarnings(
      fit.variogram(empvar_res, vgm(model="Sph", nugget=0), fit.sills=c(F,T), fit.method = 1)))
    lrange_res <- fitvar_res$range[fitvar_res$model=="Sph"]
    # Define folds based on bLOO
    folds_bLOO_res <- bLOO(train_points, lrange_res, 0.5)
    # Define folds based on NNDM
    folds_ndm_res <- nndm(train_points, rgrid, lrange_res, 0.5)

    #### Model fitting and validation
    mod <- fitval_rf_species(form,
                             folds_bLOO_out$indx_train,
                             folds_bLOO_out$indx_test,
                             folds_bLOO_res$indx_train,
                             folds_bLOO_res$indx_test,
                             folds_ndm_out$indx_train,
                             folds_ndm_out$indx_test,
                             folds_ndm_res$indx_train,
                             folds_ndm_res$indx_test,
                             pgrid, train_data, grid_data)
    mod_all <- cbind(mod, data.frame(outrange=lrange_out, resrange=lrange_res))

    # Store results of the iteration
    res_it <- cbind(data.frame(dsample=dist_it, stringsAsFactors = FALSE),
                    mod_all)
    res <- bind_rows(res, res_it)
  }

  row.names(res) <- NULL
  res
}
