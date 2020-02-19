# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("rota")

list.dirs()
dir.create("logS")
dir.create("forecasts")

library(surveillance)
library(hhh4addon)
library(MASS)
source("../auxiliary_functions_glmnb.R")


# get data:
data("rotaBE")

ts <- rotaBE@observed

# timepoints for which to fit models:
max_horizon <- 4
tps <- (4*52 - max_horizon):(7*52 - 1)
n_units <- ncol(rotaBE@observed)

multiv_logS_rota_naive <- matrix(ncol = max_horizon, nrow = length(tps),
                                 dimnames = list(paste0("t=", tps),
                                                 paste0("h", 1:max_horizon)))

# data.frame to store detailed results:
results_detailed_rota_naive <- data.frame(data_set = rep("rotaBE", n_units*max_horizon*length(tps)),
                                          unit = NA,
                                          prediction_horizon = NA_integer_,
                                          model = "naive",
                                          prediction_time = NA_integer_,
                                          pred_mean = NA, pred_var = NA,
                                          lb50 = NA, ub50 = NA,
                                          lb95 = NA, ub95 = NA,
                                          obs = NA,
                                          unit_wise_log_score = NA,
                                          multiv_log_score = NA)


for(ind in tps){
  max_horizon_temp <- min(max_horizon, max(tps) - ind + 1)
  # obtain 1 through 4-week-ahead forecasts:
  # undebug(naive_forecast_glmnb_multiv)
  pred_temp <- naive_forecast_glmnb_multiv(ts = ts,  t_cond = ind,
                                           max_horizon = max_horizon_temp)

  # determine in which rows to store
  inds_rows <- seq(from = min(which(is.na(results_detailed_rota_naive$prediction_time))),
                   length.out = n_units*max_horizon_temp)
  # insert results into df:
  results_detailed_rota_naive[inds_rows, names(pred_temp)] <- pred_temp
  # store log scores separately:
  multiv_logS_rota_naive[paste0("t=", ind), ] <-
    c(pred_temp$multiv_log_score[1:max_horizon_temp],
      rep(NA, max_horizon - max_horizon_temp))

  if(any(table(results_detailed_rota_naive$prediction_time) > n_units*max_horizon)) stop("Double entry in results_detailed")

  print(ind)
}

# add NAs at top for irrelevant forecasts:
for(i in 1:(max_horizon - 1)){
  multiv_logS_rota_naive[i, (1:(max_horizon - i))] <- NA
}

# write out results:
write.csv(multiv_logS_rota_naive, file = "logS/multiv_logS_rota_naive_glmnb.csv")
write.csv(results_detailed_rota_naive, file = "forecasts/forecasts_rota_naive_glmnb.csv")