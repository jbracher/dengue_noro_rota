# This file serves to remove the sts objects from all noro and rota model fits stored for forecasting purposes.
# This is to avoid storing too large files in the github repository.
# Note that these sts objects are not used anywhere anyway as in all obtain_forecast files
# they are replaced by the respective sts objects from the hhh4addon package (this is done
# because the sts objects in the model fits contain NAs which served to steer the subset.)

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
source("auxiliary_functions.R")
source("basic_settings.R")

tps_dengue <- (19*52 - max_horizon):(23*52 - 1)

# dengue:
setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/dengue/model_fits")
for(lag_structure in c("ar1", "pois", "lin", "geom", "siraj", "unres")){
  folder_name_temp <- paste0("dengue_", lag_structure)
  file_names_temp <- list.files(folder_name_temp)
  for(f in file_names_temp){
    load(paste0(folder_name_temp, "/", f))
    fit_name_temp <- paste0("fit_dengue_", lag_structure, "_temp")
    fit_temp <- get(paste0(fit_name_temp))
    fit_temp$stsObj <- NULL
    assign(fit_name_temp, fit_temp)
    save(list = fit_name_temp, file = paste0(folder_name_temp, "/", f))
  }
}

library(hhh4addon)
data("noroBE")
a

prune_stsObj <- function(stsObj){
  stsObj@map <- SpatialPolygons(list())
  stsObj@state <- stsObj@alarm <- stsObj@upperbound <- matrix()
  stsObj
}

object_size(noroBE)
object_size(prune_stsObj(noroBE))

# noro and rota:
setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
for(disease in c("noro", "rota")){

  # # fits to full data sets:
  # folder_name_temp <- paste0(disease, "/model_fits/", disease, "_364")
  #
  # f <- paste0("/fits_", disease, "_364.rda")
  # load(paste0(folder_name_temp, f))
  # ls_name_temp <- paste0("fits_", disease)
  # ls_temp <- get(ls_name_temp)
  # for(i in seq_along(ls_temp)){
  #   for(j in seq_along(ls_temp[[i]])){
  #     ls_temp[[i]][[j]]$stsObj <- prune_stsObj(ls_temp[[i]][[j]]$stsObj)
  #   }
  # }
  # assign(paste0("fits_", disease), ls_temp)
  # save(list = ls_name_temp, file = paste0(folder_name_temp, "/", f))
  #
  #
  # for(model_version in c("full", "gravity")){
  #   f <- paste0("/fits_", disease, "_", model_version,
  #              "_vary_max_lag_364.rda")
  #   load(paste0(folder_name_temp, f))
  #   ls_name_temp <- paste0("fits_", disease, "_", model_version,
  #                          "_vary_max_lag")
  #   ls_temp <- get(ls_name_temp)
  #   for(i in seq_along(ls_temp)){
  #     for(j in seq_along(ls_temp[[i]])){
  #       ls_temp[[i]][[j]]$stsObj <- prune_stsObj(ls_temp[[i]][[j]]$stsObj)
  #     }
  #   }
  #   assign(paste0("fits_", disease), ls_name_temp)
  #   save(list = ls_name_temp, file = paste0(folder_name_temp, "/", f))
  # }

  # fits to subsets for forecasting:
  for(model_version in c("full", "gravity")){
    for(lag_structure in c("ar1", "pois", "geom", "lin", "unres")){
      folder_name_temp <- paste0(disease, "/model_fits/", disease,
                                 "_", model_version, "_", lag_structure)
      file_names_temp <- list.files(folder_name_temp)
      for(f in file_names_temp){
        load(paste0(folder_name_temp, "/", f))
        fit_name_temp <- paste0("fit_", disease, "_", model_version, "_", lag_structure, "_temp")
        fit_temp <- get(paste0(fit_name_temp))
        fit_temp$stsObj <- prune_stsObj(fit_temp$stsObj)
        assign(fit_name_temp, fit_temp)
        save(list = fit_name_temp, file = paste0(folder_name_temp, "/", f))
      }
    }
  }
}
