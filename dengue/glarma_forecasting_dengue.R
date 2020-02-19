# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("dengue")

# load packages:
library(surveillance)
library(hhh4addon)
library(glarma)
source("../auxiliary_functions_glarma.R")

# create folder structure if necessary:
list.dirs()
dir.create("model_fits")
dir.create(paste0("model_fits/dengue_glarma"))

# get and plot data:
data(dengueSJ)
plot(dengueSJ)

# bring into glarma format:
y <- as.vector(dengueSJ@observed)
X <- t(sapply(2*pi*seq_along(y)/52,
              function (x) c(sin = sin(x), cos = cos(x))))
X <- cbind(intercept = 1, X)

# select orders:
glarma_fits_training_ar <-  glarma_fits_training_ma <- list()
for(max_lag in 1:12){
  glarma_fits_training_ar[[max_lag]] <-
    glarma(y = y[1:980], X = X[1:980, ], type = "NegBin", phiLags = 1:max_lag)
  try({
    glarma_fits_training_ma[[max_lag]] <-
      glarma(y = y[1:980], X = X[1:980, ], type = "NegBin", thetaLags = 1:max_lag)
  })
  print(max_lag)
}
# best fit is achieved with order 10

plot(sapply(glarma_fits_training_ar, extractAIC))
sapply(glarma_fits_training_ma, function(X ) if(!is.null(X)){extractAIC(X)} else NA)

# one-step-ahead forecasts (can be obtained without simulation):
max_horizon <- 8
tps <- (19*52 - max_horizon):(23*52 - 1)

preds_glarma <- matrix(ncol = 3, nrow = length(tps))
colnames(preds_glarma) <- c("mu", "size", "obs")
rownames(preds_glarma) <- paste0("t", tps)

for(ind in tps){
  # fit model
  glarmafit_temp <- glarma(y = y[1:ind], X = X[1:ind, ], type = "NegBin", phiLags = 1:10)
  # obtain predictive mean and size
  mu_temp <- forecast(glarmafit_temp, n.ahead = 1, newdata = X[ind + 1, , drop = FALSE])$mu
  size_temp <- coef(glarmafit_temp, type = "NB")
  obs_temp <- y[ind + 1]
  # store:
  preds_glarma[paste0("t", ind), ] <- c(mu_temp, size_temp, obs_temp)
  print(ind)
}

# obtain log scores:
osa_logS <- dnbinom(x = preds_glarma[, "obs"],
                    mu = preds_glarma[, "mu"],
                    size = preds_glarma[, "size"], log = TRUE)
mean(osa_logS[-(1:7)]) # have to remove first observations as these are not yet in the validation period.


# data.frame to store detailed results:
results_detailed <- data.frame(data_set = rep("dengue_sj", max_horizon*length(tps)),
                                     unit = 1,
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
log_scores_n_step_ahead <- matrix(ncol = 8, nrow = length(tps))
rownames(log_scores_n_step_ahead) <- paste0("t_cond=", tps)
colnames(log_scores_n_step_ahead) <- paste0("h", 1:max_horizon)

for(ind in tps){
  # obtain 1 through 8-week-ahead forecasts:
  pred_temp <- forecasting_nStepAhead_glarma(y = y, X = X, t_condition = ind,
                                                 max_horizon = max_horizon, phiLags = 1:10,
                                                 niter = 1000)

  # determine in which rows to store
  inds_rows <- seq(from = min(which(is.na(results_detailed$prediction_time))), length.out = max_horizon)
  # insert results into df:
  results_detailed$prediction_horizon[inds_rows] <- 1:max_horizon
  results_detailed$prediction_time[inds_rows] <- ind
  results_detailed[inds_rows, names(pred_temp)] <- pred_temp

  # store log scores separately:
  log_scores_n_step_ahead[paste0("t_cond=", ind), ] <- pred_temp$unit_wise_log_score

  if(any(table(results_detailed$prediction_time) > max_horizon)) stop("Double entry inresults_detailed")
  print(ind)
}


colMeans(log_scores_n_step_ahead, na.rm = TRUE)

# write out:
write.csv(log_scores_n_step_ahead, file = "logS/logS_dengue_glarma.csv")
write.csv(results_detailed, file = "forecasts/forecasts_dengue_glarma.csv")
