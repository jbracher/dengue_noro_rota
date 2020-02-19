# Fitting and storing dengue models for all timepoints at which forecasts are generated

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("dengue")
# scp johannes@130.60.71.234:/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/dengue/fit_models_dengue.R fit_models_dengue.R

# load packages:
library(surveillance)
library(hhh4addon)

# create folder structure if necessary:
list.dirs()
dir.create("model_fits")
for(lag_structure in c("ar1", "geom", "pois", "lin", "unres", "end", "siraj")){
  dir.create(paste0("model_fits/dengue_", lag_structure))
}

# get some additional functions:
# source("../auxiliary_functions.R")

# get and plot data:
data(dengueSJ)
plot(dengueSJ)

names_lag_structures <- c("end", "ar1", "pois", "lin", "geom", "unres", "siraj")


# define controls for different lag weighting schemes:
max_lag_param <- 5
max_lag_unres <- 4
ctrl_dengue <- list(
  ar = list(f = addSeason2formula(f = ~ 1, S = 2), lag = 1),
  end = list(f = addSeason2formula(f = ~ 1, S = 1, period = 52)),
  family = "NegBin1",
  max_lag  = max_lag_param
)

ctrl_dengue_end <- list(
  ar = list(f = ~-1),
  end = list(f = addSeason2formula(f = ~ 1, S = 1, period = 52)),
  family = "NegBin1"
)

# timepoints for which to fit models:
tps <- (19*52 - 8):(23*52 - 1)

# for the linear lag we adopt a grid-based optimization
vals_kappa_lin <- 1:99/100
vals_par_lag_lin <- log(vals_kappa_lin/(1 - vals_kappa_lin))

# weighting scheme based on serial intervals taken from Siraj:
wgts_siraj <- c(0.001, 0.999*c(0.2, 0.425, 0.25, 0.125)) # cannot asign probability 0 to lag 1.
par_lag_siraj <- log(wgts_siraj[-1]/(1 - sum(wgts_siraj[-1])))
# unrestricted_lag(par_lag_siraj, 1, 5) # works

# fit models for all time points during the evaluation period
# we originally also experimented with models which have a min_lag larger than 1,
# i.e. force the first couple of weights to 0, but this did not yield improvements.
# therefore only run min_lag = 1
for(min_lag in 1:1){

  # adapt control settings to min_lag:
  ctrl_dengue_temp <- ctrl_dengue
  ctrl_dengue_temp$subset <- (min_lag + max_lag_param):nrow(dengueSJ)
  ctrl_dengue_temp$ar$lag <- ctrl_dengue_temp$ne$lag <- min_lag
  ctrl_dengue_temp$min_lag <- min_lag
  ctrl_dengue_temp$max_lag <- min_lag + max_lag_param - 1

  # define controls with different weighting schemes:
  ctrl_dengue_pois_temp <- ctrl_dengue_temp; ctrl_dengue_pois_temp$funct_lag <- poisson_lag
  ctrl_dengue_lin_temp <- ctrl_dengue_temp; ctrl_dengue_lin_temp$funct_lag <- linear_lag

  ctrl_dengue_unres_temp <- ctrl_dengue_temp; ctrl_dengue_unres_temp$funct_lag <- unrestricted_lag
  ctrl_dengue_unres_temp$max_lag <- min_lag + max_lag_unres - 1

  ctrl_dengue_siraj_temp <- ctrl_dengue_temp; ctrl_dengue_siraj_temp$funct_lag <- unrestricted_lag
  ctrl_dengue_siraj_temp$par_lag <- par_lag_siraj

  # run over timepoints in validation period:
  for(ind in tps){

    # steer subset via NAs:
    dengueSJ_temp <- dengueSJ
    if(ind < nrow(dengueSJ)){
      dengueSJ_temp@observed[(ind + 1):nrow(dengueSJ_temp), ] <- NA
    }

    # fit endemic-only model
    fit_dengue_end_temp <- hhh4(dengueSJ_temp, ctrl_dengue_end)
    save(fit_dengue_end_temp, file = paste0("model_fits/dengue_end",
                                            "/fit_dengue_end_", ind, ".rda"))

    # fit ar1 model:
    fit_dengue_ar1_temp <- hhh4(dengueSJ_temp, ctrl_dengue_temp)
    save(fit_dengue_ar1_temp, file = paste0("model_fits/dengue_ar1",
                                          "/fit_dengue_ar1_", ind, ".rda"))

    # fit model with geometric lags:
    start_par_lag <- ifelse(ind == tps[1], 0.5, fit_dengue_geom_temp$par_lag)
    fit_dengue_geom_temp <- profile_par_lag(dengueSJ_temp, ctrl_dengue_temp,
                                          start_par_lag = start_par_lag)
    save(fit_dengue_geom_temp, file = paste0("model_fits/dengue_geom",
                                           "/fit_dengue_geom_", ind, ".rda"))

    # fit model with Poisson lags:
    start_par_lag <- ifelse(ind == tps[1], 0.5, fit_dengue_pois_temp$par_lag)
    fit_dengue_pois_temp <- profile_par_lag(dengueSJ_temp, ctrl_dengue_pois_temp,
                                          start_par_lag = start_par_lag)
    save(fit_dengue_pois_temp, file = paste0("model_fits/dengue_pois",
                                           "/fit_dengue_pois_", ind, ".rda"))

    # fit model with linear lags using grid-based optimization:
    fit_dengue_lin_temp <- fit_par_lag(dengueSJ_temp, ctrl_dengue_lin_temp,
                                       range_par = vals_par_lag_lin, use_update = FALSE)$best_mo
    start_par_lag <- ifelse(ind == tps[1], 0.5, fit_dengue_lin_temp$par_lag)
    fit_dengue_lin_temp <- profile_par_lag(dengueSJ_temp, ctrl_dengue_lin_temp,
                                           start_par_lag = start_par_lag)
    save(fit_dengue_lin_temp, file = paste0("model_fits/dengue_lin",
                                            "/fit_dengue_lin_", ind, ".rda"))

    # fit model with unconstrained weights:
    start_par_lag <- if(ind == tps[1]){
      rep(0.5, max_lag_unres - 1)
    }else{
      fit_dengue_unres_temp$par_lag
    }
    fit_dengue_unres_temp <- profile_par_lag(dengueSJ_temp, ctrl_dengue_unres_temp,
                                           start_par_lag = start_par_lag)
    save(fit_dengue_unres_temp, file = paste0("model_fits/dengue_unres",
                                            "/fit_dengue_unres_", ind, ".rda"))

    # fit model with weights taken from Siraj 2017:
    fit_dengue_siraj_temp <- hhh4_lag(dengueSJ_temp, ctrl_dengue_siraj_temp)
    save(fit_dengue_siraj_temp, file = paste0("model_fits/dengue_siraj",
                                              "/fit_dengue_siraj_", ind, ".rda"))

    print(ind)
  }
}
