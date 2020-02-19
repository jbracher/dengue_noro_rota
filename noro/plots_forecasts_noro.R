# create plots to describe norovirus forecasts:

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("noro")

library(surveillance)
library(hhh4addon)
data("noroBE")

# get helper functions:
source("../auxiliary_functions.R")
source("../basic_settings.R")

# get results for geometric lags:
forecasts_noro_full_geom <- read.csv("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/noro/forecasts/forecasts_noro_full_geom.csv")
forecasts_noro_full_ar1 <- read.csv("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/noro/forecasts/forecasts_noro_full_ar1.csv")

# make subsets according to horizons (subs_noro_full_geom) and combination
# of horizon and unit (subs_noro_full_geom_u)
subs_noro_full_geom <- subs_noro_full_geom_u <-
  subs_noro_full_ar1 <- subs_noro_full_ar1_u <- list()
for(horizon in 1:4){
  subs_noro_full_geom[[horizon]] <- subset(forecasts_noro_full_geom,
                                           forecasts_noro_full_geom$prediction_horizon == horizon)
  subs_noro_full_ar1[[horizon]] <- subset(forecasts_noro_full_ar1,
                                           forecasts_noro_full_ar1$prediction_horizon == horizon)

  subs_noro_full_geom_u[[horizon]] <- list()
  for(un in 1:12){
    subs_noro_full_geom_u[[horizon]][[un]] <- subset(forecasts_noro_full_geom,
                                                     forecasts_noro_full_geom$prediction_horizon == horizon &
                                                       forecasts_noro_full_geom$unit == un)
  }
}

# timepoints in calendar time (for plotting)
timepoints <- seq(from = 2015, by = 1/52, length.out = 3*52 + 1)


par(mfrow = c(4, 2), mar = c(3.2, 4.5, 3, 1), font.main = 1, family = "serif", las = 1)

cex.obs <- 0.6
h_left <- 1 # choose horizon shown in left column of plot
h_right <- 4 # choose horizon shown in right column of plot
fact_ylim <- 1.7 # factor to steer ylim

par(mfrow = c(4, 2), mar = c(3.2, 4.5, 3, 1), font.main = 1, family = "serif", las = 1)

cex.obs <- 0.6
h_left <- 1 # choose horizon shown in left column of plot
h_right <- 4 # choose horizon shown in right column of plot
fact_ylim <- 1.7 # factor to steer ylim

# plot horizons 1 and 4 in two columns:
for(unit in 1:12){
  # left column: horizon h_left weeks
  # start with empty plot
  plot(NULL, xlim = range(timepoints),
       ylim = c(0, fact_ylim*max(subs_noro_full_geom_u[[h_left]][[unit]]$obs)),
       xlab = "", ylab = "",
       axes = FALSE, main = paste0(full_names_districts[unit], ", h = 1"))
  axis(1, at = 2013:2018); axis(2); box()
  if(unit %% 4 == 1){
    legend("topleft",
           legend = c("95% PI", "50% PI", "pred.mean"),
           col = c(cols_models_nr_transp["geom"],
                   cols_models_nr["geom"],
                   cols_models_nr["geom"]),
           pch = c(15, 15, 0), bty = "n")
  }

  # add 95% PI
  shade_forecast(timepoints = tail(timepoints, length(subs_noro_full_geom_u[[h_left]][[unit]]$lb95)),
                 lb = subs_noro_full_geom_u[[h_left]][[unit]]$lb95,
                 ub = subs_noro_full_geom_u[[h_left]][[unit]]$ub95,
                 col = cols_models_nr_transp["geom"])
  # overplot with 50% PI
  shade_forecast(timepoints = tail(timepoints, nrow(subs_noro_full_geom_u[[h_left]][[unit]])),
                 lb = subs_noro_full_geom_u[[h_left]][[unit]]$lb50,
                 ub = subs_noro_full_geom_u[[h_left]][[unit]]$ub50,
                 col = cols_models_nr["geom"])
  lines(tail(timepoints, nrow(subs_noro_full_geom_u[[h_left]][[unit]])),
        subs_noro_full_geom_u[[h_left]][[unit]]$pred_mean,
        col = "white")

  # add observations
  points(tail(timepoints, length(subs_noro_full_geom_u[[h_left]][[unit]]$obs)),
         subs_noro_full_geom_u[[h_left]][[unit]]$obs, pch = 16, cex = cex.obs)

  # right column: horizon h_right weeks
  plot(NULL, xlim = range(timepoints),
       ylim = c(0, fact_ylim*max(subs_noro_full_geom_u[[h_right]][[unit]]$obs)),
       xlab = "", ylab = "",
       axes = FALSE, main = paste0(full_names_districts[unit], ",  h = 3"))
  axis(1, at = 2013:2018); axis(2); box()
  shade_forecast(timepoints = tail(timepoints, length(subs_noro_full_geom_u[[h_right]][[unit]]$lb95)),
                 lb = subs_noro_full_geom_u[[h_right]][[unit]]$lb95,
                 ub = subs_noro_full_geom_u[[h_right]][[unit]]$ub95,
                 col = cols_models_nr_transp["geom"])
  shade_forecast(timepoints = tail(timepoints, nrow(subs_noro_full_geom_u[[h_right]][[unit]])),
                 lb = subs_noro_full_geom_u[[h_right]][[unit]]$lb50,
                 ub = subs_noro_full_geom_u[[h_right]][[unit]]$ub50,
                 col = cols_models_nr["geom"])

  lines(tail(timepoints, nrow(subs_noro_full_geom_u[[h_right]][[unit]])),
        subs_noro_full_geom_u[[h_right]][[unit]]$pred_mean,
        col = "white")

  points(tail(timepoints, length(subs_noro_full_geom_u[[h_right]][[unit]]$obs)),
         subs_noro_full_geom_u[[h_right]][[unit]]$obs, pch = 16, cex = cex.obs)
}



###
# fans starting at a given time point to illustrate how forecasts get wider:

# timepoints at which ans are to start:
t_cond_plot <- 300 + 0:6*7
# unit for which to show plot:
unit_plot <- 2

par(mfrow = c(1, 1))
# start with empty plot:
plot(NULL, xlim = range(timepoints[t_cond_plot - 4*52]), ylim = c(0, 30),
     xlab = "", ylab = "")

# add the different fanplots:
for(i in seq_along(t_cond_plot)){
  # get subset
  subs_temp <- subset(forecasts_noro_full_geom,
                      unit == unit_plot & prediction_time == t_cond_plot[i])
  # add observation at time of forecast
  obs_i <- noroBE@observed[t_cond_plot[i], unit_plot]; names(obs_i) <- NULL
  subs_temp <- rbind(NA, subs_temp)
  subs_temp[1, c("pred_mean", "lb50", "ub50", "lb95", "ub95", "obs")] <- obs_i

  # add to plot
  t_plot_temp <- t_cond_plot[i] + 0:4
  shade_forecast(timepoints[t_plot_temp - 4*52], lb = subs_temp$lb95, ub = subs_temp$ub95,
                 col = cols_models_nr["geom"])
  shade_forecast(timepoints[t_plot_temp - 4*52], lb = subs_temp$lb50, ub = subs_temp$ub50,
                 col = cols_models_nr["geom"])
  points(timepoints[t_plot_temp - 4*52], subs_temp$obs, pch = 16, cex = cex.obs)
}

# add observations
points(timepoints,
       tail(noroBE@observed[, unit_plot], 3*52),
       pch = 16, cex = cex.obs)

# dip corresponds to Christmas break


# plausibility check against predictive_moments:
t_plausibility <- 310:313 # choose a timepoint to perform plausibility check:

# load corresponding models
load(paste0("model_fits/noro_full_geom/fit_noro_full_geom_", t_plausibility[1], ".rda"))
fit1 <- fit_noro_full_geom_temp
load(paste0("model_fits/noro_full_geom/fit_noro_full_geom_", t_plausibility[2], ".rda"))
fit2 <- fit_noro_full_geom_temp
load(paste0("model_fits/noro_full_geom/fit_noro_full_geom_", t_plausibility[3], ".rda"))
fit3 <- fit_noro_full_geom_temp
load(paste0("model_fits/noro_full_geom/fit_noro_full_geom_", t_plausibility[4], ".rda"))
fit4 <- fit_noro_full_geom_temp

library(hhh4addon)
# obtain moment forecasts
p1 <- predictive_moments(fit1, t_condition = t_plausibility[1], lgt = 4)
p2 <- predictive_moments(fit2, t_condition = t_plausibility[2], lgt = 4)
p3 <- predictive_moments(fit3, t_condition = t_plausibility[3], lgt = 4)
p4 <- predictive_moments(fit4, t_condition = t_plausibility[4], lgt = 4)

# compare to moments in forecast files
# means
subset(forecasts_noro_full_geom, prediction_horizon == 3 &
         unit == 1 & prediction_time %in% t_plausibility)$pred_mean
c(p1$mu_matrix[3, 1],
  p2$mu_matrix[3, 1],
  p3$mu_matrix[3, 1],
  p4$mu_matrix[3, 1])

# variances
subset(forecasts_noro_full_geom, prediction_horizon == 3 &
         unit == 1 & prediction_time %in% t_plausibility)$pred_var
c(p1$var_matrix[3, 1],
  p2$var_matrix[3, 1],
  p3$var_matrix[3, 1],
  p4$var_matrix[3, 1])
# works

# PIT plots:
par(mfrow = c(2, 2), las = 1)

# horizon 1 wk
plot_pit_hist(pit_l = subs_noro_full_geom[[1]]$unit_wise_pit_l,
              pit_u = subs_noro_full_geom[[1]]$unit_wise_pit_u,
              main = "Horizon 1 week")

# horizon 2 wk
plot_pit_hist(pit_l = subs_noro_full_geom[[2]]$unit_wise_pit_l,
              pit_u = subs_noro_full_geom[[2]]$unit_wise_pit_u,
              main = "Horizon 2 weeks")

# horizon 3 wk
plot_pit_hist(pit_l = subs_noro_full_geom[[3]]$unit_wise_pit_l,
              pit_u = subs_noro_full_geom[[3]]$unit_wise_pit_u,
              main = "Horizon 3 weeks")

# horizon 4 wk
plot_pit_hist(pit_l = subs_noro_full_geom[[4]]$unit_wise_pit_l,
              pit_u = subs_noro_full_geom[[4]]$unit_wise_pit_u,
              main = "Horizon 4 weeks")

# coverage of 50% and 95% prediction intervals:
coverage_noro_full_geom <- sapply(subs_noro_full_geom, eval_coverage)
coverage_noro_full_ar1 <- sapply(subs_noro_full_ar1, eval_coverage)

# plot:
plot(NULL, xlim = c(1, 4), ylim = c(0, 1),
     xlab = "horizon", ylab = "coverage")

abline(h = 0.5, lty = "dashed")
lines(coverage_noro_full_geom["coverage50", ], cex = 0.5, type = "b",
      col = cols_models_nr["geom"], pch = 0, lty = "dashed")
lines(coverage_noro_full_ar1["coverage50", ], cex = 0.5, type = "b",
      col = cols_models_nr["ar1"], pch = 0, lty = "dashed")

abline(h = 0.95, lty = "solid")
lines(coverage_noro_full_geom["coverage95", ], cex = 0.5, type = "b",
      col = cols_models_nr["geom"], pch = 15, lty = "solid")
lines(coverage_noro_full_ar1["coverage95", ], cex = 0.5, type = "b",
      col = cols_models_nr["ar1"], pch = 15, lty = "solid")

legend("bottom", legend = c("fixed", "geometric", "Poisson", "nominal level", "50% PI", "95% PI"),
       ncol = 3,
       col = c(cols_models_nr[c("ar1", "geom", "pois")], "black", "black", "black"),
       lty = c(NA, NA, NA, NA, "dashed", "solid"),
       pch = c(15, 15, 15, 15, 0, 15), bty = "n")

