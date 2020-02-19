# Obtain forecasts and evaluate the multivariate logarithmic scores for norovirus forecasts
# at different horizons

# scp johannes@130.60.71.234:/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/noro/obtain_forecasts_noro.R obtain_forecasts_noro.R
# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
# setwd("noro")

library(surveillance)
library(hhh4addon)
source("../auxiliary_functions.R")

dir.create("logS")
dir.create("forecasts")

# get data:
data("noroBE")
# modify neighbourhood matrices to apply power law:
noroBE@neighbourhood <- noroBE@neighbourhood + 1

names_districts <- colnames(noroBE@observed)
n_units <- ncol(noroBE@observed)

# timepoints for which models were fitted (in fit_models_noro.R):
max_horizon <- 4 # evaluate up to horizon 4 weeks
tps <- (4*52 - max_horizon):(7*52 - 1)



# Evaluate for iterative forecasts (those are the ones used in the manuscript):
multiv_logS_noro <- results_detailed_noro <- list()
templ_multiv_logS <- matrix(ncol = max_horizon, nrow = length(tps),
                            dimnames = list(paste0("t=", tps),
                                            paste0("h", 1:max_horizon)))

templ_results_detailed <- data.frame(data_set = rep("noroBE", n_units*max_horizon*length(tps)),
                                     unit = NA,
                                     prediction_horizon = NA_integer_,
                                     lag_structure = "NA",
                                     prediction_time = NA_integer_,
                                     pred_mean = NA, pred_var = NA,
                                     lb50 = NA, ub50 = NA,
                                     lb95 = NA, ub95 = NA,
                                     obs = NA,
                                     unit_wise_log_score = NA,
                                     unit_wise_pit_l = NA,
                                     unit_wise_pit_u = NA,
                                     multiv_log_score = NA)
templ_results_detailed$lag_structure <- as.character(templ_results_detailed$lag_structure)


for(model_version in c("full", "gravity")){
  print(model_version)
  for(lag_structure in c("ar1", "geom", "pois", "lin", "unres")){
    print(lag_structure)

    # for the "end" model there is no distinction between "full" and "gravity":
    if(lag_structure == "end" & model_version == "gravity"){
        next()
    }

    # initialize matrices/fata.frames to store results:
    multiv_logS_noro[[model_version]][[lag_structure]] <- templ_multiv_logS
    results_detailed_noro[[model_version]][[lag_structure]] <- templ_results_detailed
    for(ind in tps){
      # load model:
      load(paste0("model_fits/noro_", model_version, "_", lag_structure,
                  "/fit_noro_", model_version, "_", lag_structure,  "_", ind, ".rda"))
      fit_noro_temp <- get(paste0("fit_noro_", model_version, "_", lag_structure, "_temp"))

      for(horizon in 1:max_horizon){
        if((horizon <= nrow(noroBE@observed) - ind) &
           (ind + horizon >= tps[1] + max_horizon)){
          capture.output(
            pred_noro_temp <- forecasting_nStepAhead(fit_noro_temp, stsObj = noroBE, tp = ind,
                                                horizon = horizon, n_sim = 1000)
          )
          # write out results:
          # find in which rows to store them:
          ind_rows0 <- min(which(is.na(
            results_detailed_noro[[model_version]][[lag_structure]]$prediction_horizon
            )))
          inds_rows <- seq(from = ind_rows0, length.out = n_units)

          results_detailed_noro[[model_version]][[lag_structure]][inds_rows, "unit"] <-
            1:n_units
          results_detailed_noro[[model_version]][[lag_structure]][inds_rows, "prediction_horizon"] <-
            horizon
          results_detailed_noro[[model_version]][[lag_structure]][inds_rows, "lag_structure"] <-
            lag_structure
          results_detailed_noro[[model_version]][[lag_structure]][inds_rows, "prediction_time"] <-
            ind
          results_detailed_noro[[model_version]][[lag_structure]][inds_rows, names(pred_noro_temp)] <-
            pred_noro_temp

          # store log scores separately:
          multiv_logS_noro[[model_version]][[lag_structure]][paste0("t=", ind), horizon] <-
            pred_noro_temp$multiv_log_score

          # print(horizon)
          # print(head(multiv_logS_noro[[model_version]][[lag_structure]]))
        }
      }
      print(ind)
      if(ind < 6) print(head(multiv_logS_noro[[model_version]][[lag_structure]]))
    }
    # write out results:
    write.csv(results_detailed_noro[[model_version]][[lag_structure]],
              file = paste0("forecasts/forecasts_noro_", model_version, "_",
                            lag_structure, ".csv"))
    write.csv(multiv_logS_noro[[model_version]][[lag_structure]],
              file = paste0("logS/multiv_logS_noro_", model_version, "_",
                            lag_structure, ".csv"))
  }
}
