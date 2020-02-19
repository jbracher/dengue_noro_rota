# Evaluate the logarithmic scores for dengue forecasts

# scp johannes@130.60.71.234:/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/dengue/obtain_forecasts_dengue.R obtain_forecasts_dengue.R
# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("dengue")

library(surveillance)
library(hhh4addon)
source("../auxiliary_functions.R")

# get data:
data("dengueSJ")

# timepoints for which models were fitted (in fit_models_dengue.R):
max_horizon <- 8
tps <- (19*52 - max_horizon):(23*52 - 1)
n_units <- ncol(dengueSJ@observed)
names_lag_structures <- c("ar1", "pois", "lin", "geom", "unres", "end", "siraj")

# create directory if necessary:
dir.create("logS")
dir.create("forecasts")

# Evaluate for iterative forecasts (those are the ones used in the manuscript):
templ_df <- data.frame(data_set = rep("dengue_sj", n_units*max_horizon*length(tps)),
                       unit = 1,
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
templ_df$lag_structure <- as.character(templ_df$lag_structure)

# lists to store results:
results_detailed_dengue <- logS_dengue <- list()

# run through different lag structures:
for(lag_structure in names_lag_structures){

  # matrices/data.frames to store results:
  results_detailed_dengue[[lag_structure]] <- templ_df
  logS_dengue[[lag_structure]] <-
    matrix(ncol = max_horizon, nrow = length(tps),
           dimnames = list(paste0("t_cond=", tps),
                           paste0("h", 1:max_horizon)))

  print(lag_structure)

  # run over prediction timepoints:
  for(ind in tps){
    # evaluation of log scores for order larger 1 involves simulation,
    # therefore set.seed
    set.seed(ind)

    # load model fit:
    file_name_fit <- paste0("model_fits/dengue_", lag_structure,
                            "/fit_dengue_", lag_structure, "_", ind, ".rda")
    load(file_name_fit)

    fit_dengue_temp <- get(paste0("fit_dengue_", lag_structure, "_temp"))
    fit_dengue_temp$stsObj <- dengueSJ # the stsObj of fit_dengue_temp contains NAs
    # which were introduced to steer to subset, replace by actual data:

    # evaluate for horizons 1 through 8:
    for(horizon in 1:max_horizon){
      if(horizon <= nrow(dengueSJ@observed) - ind){
        capture.output(
          pred_temp <- forecasting_nStepAhead(fit_dengue_temp, stsObj = dengueSJ, tp = ind,
                                              horizon = horizon, n_sim = 1000)
        )

        # write out results:
        ind_row <- min(which(is.na(results_detailed_dengue[[lag_structure]]$prediction_horizon)))
        results_detailed_dengue[[lag_structure]]$prediction_horizon[ind_row] <- horizon
        results_detailed_dengue[[lag_structure]]$lag_structure[ind_row] <- lag_structure
        results_detailed_dengue[[lag_structure]]$prediction_time[ind_row] <- ind
        results_detailed_dengue[[lag_structure]][ind_row, names(pred_temp)] <- pred_temp

        # also write out log score separately:
        logS_dengue[[lag_structure]][paste0("t_cond=", ind), horizon] <-
          pred_temp$multiv_log_score

      }
    }
    print(ind)
  }
  # write out results:
  write.csv(results_detailed_dengue[[lag_structure]], paste0("forecasts/forecasts_dengue_", lag_structure, ".csv"))
  # write out log scores:
  write.csv(logS_dengue[[lag_structure]], paste0("logS/logS_dengue_", lag_structure, ".csv"))
}
