# create plots to describe dengue forecasts:

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("dengue")

source("../auxiliary_functions.R")
source("../basic_settings.R")

max_horizon <- 8
tps <- (19*52 - max_horizon):(23*52 - 1)


# get KCDE results:
# get results from KCDE methods:
# (requires downloading repository from https://github.com/reichlab/article-disease-pred-with-kcde):
forecasts_dengue_kcde <- readRDS("forecasts/kcde-predictions.rds")

# get our results for ar1 models:
forecasts_dengue_ar1 <- read.csv("forecasts/forecasts_dengue_ar1.csv")

# get our results for higher-order lags:
forecasts_dengue_geom <- read.csv("forecasts/forecasts_dengue_geom.csv")
forecasts_dengue_pois <- read.csv("forecasts/forecasts_dengue_pois.csv")
forecasts_dengue_lin <- read.csv("forecasts/forecasts_dengue_lin.csv")
forecasts_dengue_unres <- read.csv("forecasts/forecasts_dengue_unres.csv")

max_horizon <- 8
tps <- (19*52 - max_horizon):(23*52 - 1)


# make subsets by horizon:
subs_kcde <- subs_dengue_geom <- subs_dengue_ar1 <-
  subs_dengue_pois <- subs_denge_unres <- subs_dengue_lin <- list()


for(horizon in 1:8){
  subs_kcde[[horizon]] <- subset(forecasts_dengue_kcde, prediction_horizon == horizon &
                                   seasonality == TRUE &
                                   bw_parameterization == "full"
  )

  subs_dengue_ar1[[horizon]] <- subset(forecasts_dengue_ar1, prediction_horizon == horizon &
                                         prediction_time + horizon >= tps[1] + max_horizon + 1 &
                                         prediction_time + horizon <= tail(tps, 1) + 1)

  subs_dengue_geom[[horizon]] <- subset(forecasts_dengue_geom, prediction_horizon == horizon &
                                          prediction_time + horizon >= tps[1] + max_horizon + 1 &
                                          prediction_time + horizon <= tail(tps, 1) + 1)

  subs_dengue_pois[[horizon]] <- subset(forecasts_dengue_pois, prediction_horizon == horizon &
                                          prediction_time + horizon >= tps[1] + max_horizon + 1 &
                                          prediction_time + horizon <= tail(tps, 1) + 1)

  subs_dengue_lin[[horizon]] <- subset(forecasts_dengue_lin, prediction_horizon == horizon &
                                         prediction_time + horizon >= tps[1] + max_horizon + 1 &
                                         prediction_time + horizon <= tail(tps, 1) + 1)

  subs_denge_unres[[horizon]] <- subset(forecasts_dengue_unres, prediction_horizon == horizon &
                                          prediction_time + horizon >= tps[1] + max_horizon + 1 &
                                          prediction_time + horizon <= tail(tps, 1) + 1)
}

library(hhh4addon)
data("dengueSJ")
# timepoints in calendar time:
timepoints <- seq(from = dengueSJ@start[1] +
                    (dengueSJ@start[2] + tps[1] + max_horizon)/dengueSJ@freq,
                  by = 1/dengueSJ@freq, length.out = nrow(subs_dengue_ar1[[1]]))



par(las = 1, mfrow = c(4, 1))

# plot 50% PI:

for(horizon in c(1, 2, 4, 8)){
  plot(NULL, xlim = range(timepoints), ylim = c(0, 500), xlab = "", ylab = "",
       main = paste0("Horizon ", horizon, "wk"))

  shade_forecast(timepoints,
                 subs_dengue_geom[[horizon]]$lb95,
                 subs_dengue_geom[[horizon]]$ub95,
                 col = cols_models_dengue_transp["geom"])

  shade_forecast(timepoints,
                 subs_dengue_geom[[horizon]]$lb50,
                 subs_dengue_geom[[horizon]]$ub50,
                 col = cols_models_dengue["geom"])

  lines(timepoints,
        subs_dengue_geom[[horizon]]$pred_mean,
        col = "white"
  )

  points(timepoints, subs_dengue_geom[[1]]$obs, pch = 16, cex = 0.5)

  if(horizon == 1){
    legend("topleft",
           legend = c("95% PI", "50% PI", "pred.mean"),
           col = c(cols_models_dengue_transp["geom"],
                   cols_models_dengue["geom"],
                   cols_models_dengue["geom"]),
           pch = c(15, 15, 0), bty = "n")
  }

}



par(las = 1, mfrow = c(4, 1))

# plot 50% PI:

rgb_mix_geom_ar1 <- (col2rgb(cols_models_dengue["ar1"]) +
                       col2rgb(cols_models_dengue["geom"]))/2/260
col_mix_geom_ar1 <- rgb(rgb_mix_geom_ar1[1], rgb_mix_geom_ar1[2], rgb_mix_geom_ar1[3])


for(horizon in c(1, 2, 4, 8)){
  plot(NULL, xlim = range(timepoints), ylim = c(0, 500), xlab = "", ylab = "",
       main = paste0("Horizon ", horizon, "wk"))

  shade_forecast(timepoints,
                 subs_dengue_geom[[horizon]]$lb50,
                 subs_dengue_geom[[horizon]]$ub50,
                 col = cols_models_dengue["geom"])

  shade_forecast(timepoints,
                 subs_dengue_ar1[[horizon]]$lb50,
                 subs_dengue_ar1[[horizon]]$ub50,
                 col = cols_models_dengue["ar1"])

  shade_forecast(timepoints,
                 pmax(subs_dengue_ar1[[horizon]]$lb50, subs_dengue_geom[[horizon]]$lb50),
                 pmin(subs_dengue_ar1[[horizon]]$ub50, subs_dengue_geom[[horizon]]$ub50),
                 col = col_mix_geom_ar1)

  points(timepoints, subs_dengue_geom[[1]]$obs, pch = 16, cex = 0.5)

  if(horizon == 1){
    legend("topleft",
           legend = c("geometric", "fixed"),
           col = c(cols_models_dengue["geom"],
                   cols_models_dengue["ar1"]),
           pch = c(15, 15, 0), bty = "n")
  }

}


par(las = 1, mfrow = c(4, 1))

rgb_mix_geom_kcde <- (col2rgb(cols_models_dengue["kcde"]) +
                        col2rgb(cols_models_dengue["geom"]))/2/260
col_mix_geom_kcde <- rgb(rgb_mix_geom_kcde[1], rgb_mix_geom_kcde[2], rgb_mix_geom_kcde[3])

# plot 50% PI:

for(horizon in c(1, 2, 4, 8)){
  plot(NULL, xlim = range(timepoints), ylim = c(0, 500), xlab = "", ylab = "",
       main = paste0("Horizon ", horizon, "wk"))

  shade_forecast(timepoints,
                 subs_dengue_geom[[horizon]]$lb50,
                 subs_dengue_geom[[horizon]]$ub50,
                 col = cols_models_dengue["geom"])

  shade_forecast(timepoints,
                 subs_kcde[[horizon]]$interval_pred_lb_50,
                 subs_kcde[[horizon]]$interval_pred_ub_50,
                 col = cols_models_dengue["kcde"])

  shade_forecast(timepoints,
                 pmax(subs_dengue_geom[[horizon]]$lb50,
                      subs_kcde[[horizon]]$interval_pred_lb_50),
                 pmin(subs_dengue_geom[[horizon]]$ub50,
                      subs_kcde[[horizon]]$interval_pred_ub_50),
                 col = col_mix_geom_kcde)

  points(timepoints, subs_dengue_geom[[1]]$obs, pch = 16, cex = 0.5)

  if(horizon == 1){
    legend("topleft",
           legend = c("geometric", "KCDE"),
           col = c(cols_models_dengue["geom"],
                   cols_models_dengue["kcde"]),
           pch = c(15, 15, 0), bty = "n")
  }

}


# PIT plots:
par(mfrow = c(2, 2), las = 1, font.main = 1, family = "serif", las = 1,
    mar = c(3, 3, 1, 1))

# horizon 1 wk
plot_pit_hist(pit_l = subs_dengue_geom[[1]]$unit_wise_pit_l,
              pit_u = subs_dengue_geom[[1]]$unit_wise_pit_u,
              main = "Horizon 1 week")

# horizon 2
plot_pit_hist(pit_l = subs_dengue_geom[[2]]$unit_wise_pit_l,
              pit_u = subs_dengue_geom[[2]]$unit_wise_pit_u,
              main = "Horizon 2 weeks")

# horizon 4 wk
plot_pit_hist(pit_l = subs_dengue_geom[[4]]$unit_wise_pit_l,
              pit_u = subs_dengue_geom[[4]]$unit_wise_pit_u,
              main = "Horizon 4 weeks")

# horizon 8 wk
plot_pit_hist(pit_l = subs_dengue_geom[[8]]$unit_wise_pit_l,
              pit_u = subs_dengue_geom[[8]]$unit_wise_pit_u,
              main = "Horizon 8 weeks")

