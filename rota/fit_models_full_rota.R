# Fitting and storing rotavirus models for all timepoints at which forecasts are generated
# Full models (gravity models done in separate file)

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("rota")
# scp johannes@130.60.71.234://home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/rota/fit_models_full_rota.R fit_models_full_rota.R

# create folder structure if necessary:
list.dirs()
dir.create("model_fits")
for(lag_structure in c("ar1", "geom", "pois", "lin", "unres")){
  dir.create(paste0("model_fits/rota_full_", lag_structure))
}

library(surveillance)
library(hhh4addon)

# get data:
data("rotaBE")
# modify neighbourhood matrices to apply power law:
rotaBE@neighbourhood <- rotaBE@neighbourhood + 1

# define variable for Christmas breaks:
christmas <- numeric(nrow(rotaBE@observed))
christmas[(seq_along(christmas) %% 52) %in% c(0, 1)] <- 1

# number of sine/cosine waves to include in the endemic and epidemic components:
S_end <- 1
S_epid <- 1
# min_lag and max_lag:
min_lag <- 1
max_lag <- 5
max_horizon <- 4

# subset:
subs <- (max_lag + 1):nrow(rotaBE@observed)

# controls for full model
ctrl_full <- list(end = list(f = addSeason2formula(~0 + fe(1, unitSpecific = TRUE) + christmas, S = S_end)),
                  ar = list(f = ~-1),
                  ne = list(f = addSeason2formula(~0 + fe(1, unitSpecific = TRUE), S = S_epid),
                            weights = W_powerlaw(maxlag = 5, normalize = TRUE, log = TRUE)),
                  lag = min_lag, # for AR(1) version
                  min_lag = min_lag,
                  max_lag = max_lag,
                  family = "NegBin1",
                  subset = subs,
                  data = list(christmas = christmas))
ctrl_full_pois <- ctrl_full; ctrl_full_pois$funct_lag <- poisson_lag
ctrl_full_lin <- ctrl_full; ctrl_full_lin$funct_lag <- linear_lag
ctrl_full_unres <- ctrl_full; ctrl_full_unres$funct_lag <- unrestricted_lag

# timepoints for which to fit models:
tps <- (4*52 - max_horizon):(7*52)

# for the linear lag we adopt a grid-based optimization
vals_kappa_lin <- 1:99/100
vals_par_lag_lin <- log(vals_kappa_lin/(1 - vals_kappa_lin))

# fit models for all time points during the evaluation period
for(ind in tps){

  # steer subset via NAs:
  rotaBE_temp <- rotaBE
  if(ind < nrow(rotaBE)){
    rotaBE_temp@observed[(ind + 1):nrow(rotaBE_temp), ] <- NA
  }

  # fit ar1 model:
  fit_rota_full_ar1_temp <- hhh4(rotaBE_temp, ctrl_full)
  save(fit_rota_full_ar1_temp, file = paste0("model_fits/rota_full_ar1",
                                             "/fit_rota_full_ar1_", ind, ".rda"))

  # fit model with geometric lags:
  start_par_lag <- ifelse(ind == tps[1], 0.5, fit_rota_full_geom_temp$par_lag)
  fit_rota_full_geom_temp <- profile_par_lag(rotaBE_temp, ctrl_full,
                                             start_par_lag = start_par_lag)
  save(fit_rota_full_geom_temp, file = paste0("model_fits/rota_full_geom",
                                              "/fit_rota_full_geom_", ind, ".rda"))

  # fit model with Poisson lags:
  start_par_lag <- ifelse(ind == tps[1], 0.5, fit_rota_full_pois_temp$par_lag)
  fit_rota_full_pois_temp <- profile_par_lag(rotaBE_temp, ctrl_full_pois,
                                             start_par_lag = start_par_lag)
  save(fit_rota_full_pois_temp, file = paste0("model_fits/rota_full_pois",
                                              "/fit_rota_full_pois_", ind, ".rda"))

  # fit model with linear lags using grid-based optimization:
  # the linear lag version is fitted with optimization over a grid to avoid getting stuck in local optima:
  # start_par_lag <- ifelse(ind == tps[1], 0.5, fit_rota_full_lin_temp$par_lag)
  # fit_rota_full_lin_temp <- profile_par_lag(rotaBE_temp, ctrl_full_lin,
  #                                           start_par_lag = start_par_lag)
  fit_rota_full_lin_temp <- fit_par_lag(rotaBE_temp, ctrl_full_lin,
                                        range_par = vals_par_lag_lin, use_update = FALSE)$best_mod
  save(fit_rota_full_lin_temp, file = paste0("model_fits/rota_full_lin",
                                             "/fit_rota_full_lin_", ind, ".rda"))

  # fit model with unconstrained weights:
  start_par_lag <- if(ind == tps[1]){
    rep(0.5, max_lag - 1)
  }else{
    fit_rota_full_unres_temp$par_lag
  }
  fit_rota_full_unres_temp <- profile_par_lag(rotaBE_temp, ctrl_full_unres,
                                              start_par_lag = start_par_lag)
  save(fit_rota_full_unres_temp, file = paste0("model_fits/rota_full_unres",
                                               "/fit_rota_full_unres_", ind, ".rda"))

  print(ind)
}
