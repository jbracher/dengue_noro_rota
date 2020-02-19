# Obtain forecasts and evaluate the multivariate logarithmic scores for norovirus forecasts
# at different horizons

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("noro")

library(surveillance)
library(hhh4addon)
source("../auxiliary_functions.R")

dir.create("dss")
dir.create("forecasts")

# get data:
data("noroBE")
# modify neighbourhood matrices to apply power law:
noroBE@neighbourhood <- noroBE@neighbourhood + 1

names_districts <- colnames(noroBE@observed)
n_units <- ncol(noroBE@observed)

# timepoints for which models were fitted (in fit_models_noro.R):
tps <- (4*52):(7*52 - 1)
max_horizon <- 4 # evaluate up to horizon 4 weeks



# Evaluate for iterative forecasts (those are the ones used in the manuscript):
multiv_dss_noro <- results_detailed_noro <- list()
templ_multiv_dss <- matrix(ncol = max_horizon, nrow = length(tps),
                           dimnames = list(paste0("t=", tps),
                                           paste0("h", 1:max_horizon)))

templ_results_moments <- data.frame(data_set = rep("noroBE", n_units*max_horizon*length(tps)),
                                    unit = NA,
                                    prediction_horizon = NA_integer_,
                                    lag_structure = "NA",
                                    prediction_time = NA_integer_,
                                    pred_mean = NA, pred_var = NA,
                                    cor1 = NA, cor2 = NA, cor3 = NA,
                                    cor4 = NA, cor5 = NA, cor6 = NA,
                                    cor7 = NA, cor8 = NA, cor9 = NA,
                                    cor10 = NA, cor11 = NA, cor12 = NA,
                                    obs = NA,
                                    unit_wise_dss = NA,
                                    multiv_dss = NA)
templ_results_moments$lag_structure <- as.character(templ_results_moments$lag_structure)


for(model_version in c("full", "gravity")){
  print(model_version)
  for(lag_structure in c("ar1", "geom", "pois", "lin", "unres")){
    print(lag_structure)

    # initialize matrices/fata.frames to store results:
    multiv_dss_noro[[model_version]][[lag_structure]] <- templ_multiv_dss
    results_detailed_noro[[model_version]][[lag_structure]] <- templ_results_moments

    for(ind in tps){
      # load model:
      load(paste0("model_fits/noro_", model_version, "_", lag_structure,
                  "/fit_noro_", model_version, "_", lag_structure,  "_", ind, ".rda"))
      fit_noro_temp <- get(paste0("fit_noro_", model_version, "_", lag_structure, "_temp"))
      fit_noro_temp$stsObj <- noroBE

      pred_noro_temp <-  predictive_moments(fit_noro_temp, t_condition = ind,
                                            lgt = min(max_horizon, nrow(noroBE@observed) - ind),
                                            return_cov_array = TRUE)

      for(horizon in 1:max_horizon){
        if(horizon <= nrow(noroBE@observed) - ind){
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
          results_detailed_noro[[model_version]][[lag_structure]][inds_rows, "pred_mean"] <-
            pred_noro_temp$mu_matrix[horizon, ]
          results_detailed_noro[[model_version]][[lag_structure]][inds_rows, "pred_var"] <-
            pred_noro_temp$var_matrix[horizon, ]
          results_detailed_noro[[model_version]][[lag_structure]][inds_rows, "obs"] <-
            pred_noro_temp$realizations_matrix[horizon, ]
          results_detailed_noro[[model_version]][[lag_structure]][inds_rows, paste0("cor", 1:12)] <-
            round(cov2cor(pred_noro_temp$cov_array[, , horizon]), 3)
          results_detailed_noro[[model_version]][[lag_structure]][inds_rows, "unit_wise_dss"] <-
            (pred_noro_temp$realizations_matrix[horizon, ] - pred_noro_temp$mu_matrix[horizon, ])/
            sqrt(pred_noro_temp$var_matrix[horizon, ])

          multiv_dss <- (t(pred_noro_temp$realizations_matrix[horizon, ] - pred_noro_temp$mu_matrix[horizon, ]) %*%
            solve(pred_noro_temp$cov_array[, , horizon]) %*%
            (pred_noro_temp$realizations_matrix[horizon, ] - pred_noro_temp$mu_matrix[horizon, ]) +
            determinant(pred_noro_temp$cov_array[, , horizon], logarithm = TRUE)$modulus[1])[1, 1]/(2*n_units)

          results_detailed_noro[[model_version]][[lag_structure]][inds_rows, "multiv_dss"] <-
            multiv_dss

          # store dss separately:
          multiv_dss_noro[[model_version]][[lag_structure]][paste0("t=", ind), horizon] <-
            multiv_dss

          print(horizon)
        }
      }
      print(ind)
    }

    # write out results:
    write.csv(results_detailed_noro[[model_version]][[lag_structure]],
              file = paste0("forecasts/moment_forecasts_noro_", model_version, "_",
                            lag_structure, ".csv"))
    write.csv(multiv_dss_noro[[model_version]][[lag_structure]],
              file = paste0("dss/multiv_dss_noro_", model_version, "_",
                            lag_structure, ".csv"))
  }
}

# summarize in data.frame:
summary_multiv_dss_noro <- data.frame(horizon = 1:4)

for(model_version in c("full", "gravity")){
  for(lag_structure in c("ar1", "geom", "pois", "lin", "unres")){
    summary_multiv_dss_noro[, paste0(model_version, "_", lag_structure)] <-
      colMeans(multiv_dss_noro[[model_version]][[lag_structure]], na.rm = TRUE)
  }
}

# add DSS for naive forecasts:
forecasts_naive_noro <- read.csv("forecasts/forecasts_noro_naive_glmnb.csv")

table(forecasts_naive_noro$prediction_time)
means_dss_naive <- means_dss_geom <- numeric(4)
dss_vals <- list()
for(hi in 1:4){
  subs_temp <- subset(forecasts_naive_noro, prediction_horizon == hi)
  dss_temp <- (log(subs_temp$pred_var) + (subs_temp$obs - subs_temp$pred_mean)^2/
                 (subs_temp$pred_var))/2
  dss_vals[[hi]] <- dss_temp
  means_dss_naive[hi] <- mean(dss_temp)
}

summary_multiv_dss_noro <- cbind(summary_multiv_dss_noro,
                                 naive_seasonal = means_dss_naive)


# store:
write.csv(summary_multiv_dss_noro, file = "dss/summary_multiv_dss_noro.csv", row.names = FALSE)
