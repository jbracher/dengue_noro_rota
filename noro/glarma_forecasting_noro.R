# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("noro")

list.dirs()
dir.create("logS")
dir.create("forecasts")

library(surveillance)
library(hhh4addon)
library(glarma)
source("../auxiliary_functions_glarma.R")

# get data:
data("noroBE")

y <- noroBE@observed

christmas <- numeric(nrow(noroBE@observed))
christmas[(seq_along(christmas) %% 52) %in% c(0, 1)] <- 1

X <- t(sapply(2*pi*seq_along(y)/52,
              function (x) c(sin = sin(x), cos = cos(x))))
X <- cbind(intercept = 1, X, christmas = christmas)

# select orders:
glarma_fits_training_ar <-  glarma_fits_training_ma <- list()
for(max_lag in 1:12){
  glarma_fits_training_ar[[max_lag]] <- glarma_fits_training_ma[[max_lag]] <- list()
  for(unit in 1:12){
    glarma_fits_training_ar[[max_lag]][[unit]] <-
      glarma(y = y[1:208, unit], X = X[1:208, ], type = "NegBin", phiLags = 1:max_lag)
  }

  try({
    glarma_fits_training_ma[[max_lag]][[unit]] <-
      glarma(y = y[1:208, unit], X = X[1:208, ], type = "NegBin", thetaLags = 1:max_lag)
  })
  print(max_lag)
}

get_mean_aic <- function(li) mean(sapply(li, extractAIC))

plot(sapply(glarma_fits_training_ar, get_mean_aic))
get_mean_aic(glarma_fits_training_ma[[2]])


# choose order 3

# timepoints for which to fit models:
max_horizon <- 4
tps <- (4*52 - max_horizon):(7*52)
n_units <- ncol(noroBE@observed)

# data.frame to store detailed results:
results_detailed <- list()
results_detailed <- data.frame(data_set = rep("noroBE", n_units*max_horizon*length(tps)),
                               unit = NA,
                               prediction_horizon = NA_integer_,
                               model = "glarma",
                               prediction_time = NA_integer_,
                               pred_mean = NA, pred_var = NA,
                               lb50 = NA, ub50 = NA,
                               lb95 = NA, ub95 = NA,
                               obs = NA,
                               unit_wise_log_score = NA,
                               unit_wise_pit_l = NA,
                               unit_wise_pit_u = NA,
                               multiv_log_score = NA)

# matrix to store log scores:
templ_log_scores_n_step_ahead <- matrix(ncol = max_horizon, nrow = length(tps))
rownames(templ_log_scores_n_step_ahead) <- paste0("t_cond=", tps)
colnames(templ_log_scores_n_step_ahead) <- paste0("h", 1:max_horizon)

log_scores_n_step_ahead <- list()

for(unit in 1:n_units){
  print(paste("Starting unit", unit))
  log_scores_n_step_ahead[[unit]] <- templ_log_scores_n_step_ahead
  for(ind in tps){
    # obtain 1 through 8-week-ahead forecasts:
    pred_temp <- forecasting_nStepAhead_glarma(y = y[, unit], X = X, t_condition = ind,
                                               max_horizon = max_horizon, phiLags = 1:3,
                                               niter = 1000)

    # determine in which rows to store
    inds_rows <- seq(from = min(which(is.na(results_detailed$prediction_time))), length.out = max_horizon)
    # insert results into df:
    results_detailed$unit[inds_rows] <- unit
    results_detailed$prediction_horizon[inds_rows] <- 1:max_horizon
    results_detailed$prediction_time[inds_rows] <- ind
    results_detailed[inds_rows, names(pred_temp)] <- pred_temp

    # store log scores separately:
    log_scores_n_step_ahead[[unit]][paste0("t_cond=", ind), ] <- pred_temp$unit_wise_log_score

    if(any(table(results_detailed$prediction_time) > n_units*max_horizon)) stop("Double entry inresults_detailed")

    print(ind)
  }
}


# aggregate:
multiv_log_scores_n_step_ahead <- Reduce("+", log_scores_n_step_ahead)/n_units

# add NAs at top for irrelevant forecasts:
for(i in 1:(max_horizon - 1)){
  multiv_log_scores_n_step_ahead[i, (1:(max_horizon - i))] <- NA
}

# write out results:
write.csv(multiv_log_scores_n_step_ahead, file = "logS/multiv_logS_noro_glarma.csv")
write.csv(results_detailed, file = "forecasts/forecasts_noro_glarma.csv")

colMeans(multiv_log_scores_n_step_ahead, na.rm = TRUE)
