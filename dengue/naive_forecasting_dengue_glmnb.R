library(surveillance)
library(hhh4addon)
library(MASS)

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("dengue")
source("../auxiliary_functions_glmnb.R")

data("dengueSJ")

tps <- (19*52 - 8):(23*52 - 1)

max_horizon <- 8
n_units <- ncol(dengueSJ@observed)



library(surveillance)
library(hhh4addon)
library(MASS)

# evaluation of log scores for order larger 1 involves simulation,
# therefore set.seed
seed <- 0

# get data:
data("dengueSJ")

ts <- dengueSJ@observed
unit <- 1
t_cond <- 988
max_horizon <- 8
freq <- 52

logS_dengue_naive <- matrix(ncol = max_horizon, nrow = length(tps),
                            dimnames = list(paste0("t=", tps),
                                            paste0("h", 1:max_horizon)))
results_detailed_dengue_naive <- data.frame(data_set = rep("dengueSJ", n_units*max_horizon*length(tps)),
                                            unit = 1,
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
  ind_rows0 <- min(which(is.na(results_detailed_dengue_naive$prediction_horizon)))
  inds_rows <- seq(from = ind_rows0, length.out = max_horizon)

  pred_naive_dengue_temp <- naive_forecast_glmnb(ts = ts, unit = 1, t_cond = ind,
                                               max_horizon = max_horizon, freq = 52)
  results_detailed_dengue_naive[inds_rows, names(pred_naive_dengue_temp)] <- pred_naive_dengue_temp
  # store log scores separately:
  logS_dengue_naive[paste0("t=", ind), ] <- pred_naive_dengue_temp$unit_wise_log_score

  print(ind)
}

# add NAs at top for irrelevant forecasts:
for(i in 1:(max_horizon - 1)){
  logS_dengue_naive[i, (1:(max_horizon - i))] <- NA
}

tail(logS_dengue_naive)

colMeans(logS_dengue_naive, na.rm = TRUE)

# write out results:
write.csv(results_detailed_dengue_naive,
          file = "forecasts/forecasts_dengue_naive_glmnb.csv")
write.csv(logS_dengue_naive,
          file = "logS/logS_dengue_naive_glmnb.csv")
