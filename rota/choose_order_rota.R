# compare fits for different values of the order p:

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting//")
setwd("rota")

library(surveillance)
library(hhh4addon)

data("rotaBE")
# modify neighbourhood matrices to apply power law:
rotaBE@neighbourhood <- rotaBE@neighbourhood + 1

# define variable for Christmas breaks:
christmas <- numeric(nrow(rotaBE@observed))
christmas[(seq_along(christmas) %% 52) %in% c(0, 1)] <- 1

# number of sine/cosine waves to include in the endemic and epidemic components:
S_end <- 1
S_epid <- 1
# max_lag:
max_lag <- 5
# subset:
n_seas <- 7
subs <- 8:(n_seas*52)

# for the linear lag we adopt a grid-based optimization
vals_kappa_lin <- 1:99/100
vals_par_lag_lin <- log(vals_kappa_lin/(1 - vals_kappa_lin))

#### Full model:

# define controls:
ctrl_full <- list(end = list(f = addSeason2formula(~0 + fe(1, unitSpecific = TRUE) + christmas, S = S_end)),
                  ar = list(f = ~-1),
                  ne = list(f = addSeason2formula(~0 + fe(1, unitSpecific = TRUE), S = S_epid),
                            weights = W_powerlaw(maxlag=5, normalize=TRUE, log=TRUE)),
                  max_lag = max_lag,
                  family = "NegBin1",
                  subset = subs,
                  data = list(christmas = christmas))
ctrl_full_pois <- ctrl_full; ctrl_full_pois$funct_lag <- poisson_lag
ctrl_full_lin <- ctrl_full; ctrl_full_lin$funct_lag <- linear_lag
ctrl_full_unres <- ctrl_full; ctrl_full_unres$funct_lag <- unrestricted_lag

# initialize lists to store models:
fits_rota_full_vary_max_lag <- list()

fits_rota_full_vary_max_lag$geom[[1]] <-
  fits_rota_full_vary_max_lag$pois[[1]] <-
  fits_rota_full_vary_max_lag$lin[[1]] <-
  fits_rota_full_vary_max_lag$unres[[1]] <- hhh4(rotaBE, ctrl_full)

# fit models with different orders:
for(max_lag in 2:7){
  # adapt max_lag in controls:
  ctrl_full$max_lag <- ctrl_full_pois$max_lag <-
    ctrl_full_lin$max_lag <- ctrl_full_unres$max_lag <- max_lag

  # fit models:
  fits_rota_full_vary_max_lag$geom[[max_lag]] <- profile_par_lag(rotaBE, ctrl_full)
  fits_rota_full_vary_max_lag$pois[[max_lag]] <- profile_par_lag(rotaBE, ctrl_full_pois)

  # use grid-based optimization for linear lag:
  fits_rota_full_vary_max_lag$lin[[max_lag]] <- fit_par_lag(rotaBE, ctrl_full_lin,
                                                            range_par = vals_par_lag_lin,
                                                            use_update = FALSE)$best_mod

  if(max_lag > 2){
    fits_rota_full_vary_max_lag$unres[[max_lag]] <-
      profile_par_lag(rotaBE, ctrl_full_unres,
                      start_par_lag = c(fits_rota_full_vary_max_lag$unres[[max_lag - 1]]$par_lag, -2))
  }else{
    fits_rota_full_vary_max_lag$unres[[max_lag]] <-
      profile_par_lag(rotaBE, ctrl_full_unres)
  }

  print(max_lag)
}

# evaluate AICs:
AICs_vary_max_lag <- matrix(NA, ncol = 4, nrow = 7,
                            dimnames = list(NULL, c("pois", "lin", "geom", "unres")))

for(lag_structure in c("pois", "lin", "geom", "unres")){
  for(max_lag in 1:7){
    AICs_vary_max_lag[max_lag, lag_structure] <-
      AIC(fits_rota_full_vary_max_lag[[lag_structure]][[max_lag]])
  }
}

# store:
save(fits_rota_full_vary_max_lag, file = paste0("model_fits/rota_", n_seas*52 , "/fits_rota_full_vary_max_lag_", n_seas*52, ".rda"))
write.csv(AICs_vary_max_lag, file = paste0("AIC/AIC_rota_full_vary_max_lag_", n_seas*52, ".csv"))



### gravity model

# define controls:
ctrl_gravity <- list(end = list(f = addSeason2formula(~0 + fe(1, unitSpecific = TRUE) + christmas, S = S_end)),
                     ar = list(f = ~-1),
                     ne = list(f = addSeason2formula(~1 + log(pop), S = S_epid),
                               weights = W_powerlaw(maxlag=5, normalize=TRUE, log=TRUE)),
                     max_lag = max_lag,
                     family = "NegBin1",
                     subset = subs,
                     data = list(christmas = christmas, pop = rotaBE@populationFrac))
ctrl_gravity_pois <- ctrl_gravity; ctrl_gravity_pois$funct_lag <- poisson_lag
ctrl_gravity_lin <- ctrl_gravity; ctrl_gravity_lin$funct_lag <- linear_lag
ctrl_gravity_unres <- ctrl_gravity; ctrl_gravity_unres$funct_lag <- unrestricted_lag


# initialize lists to store models:
fits_rota_gravity_vary_max_lag <- list()
fits_rota_gravity_vary_max_lag$geom[[1]] <-
  fits_rota_gravity_vary_max_lag$pois[[1]] <-
  fits_rota_gravity_vary_max_lag$lin[[1]] <-
  fits_rota_gravity_vary_max_lag$unres[[1]] <- hhh4(rotaBE, ctrl_gravity)

# fit models for different orders:
for(max_lag in 2:7){
  # adapt max_lag in controls:
  ctrl_gravity$max_lag <- ctrl_gravity_pois$max_lag <-
    ctrl_gravity_lin$max_lag <- ctrl_gravity_unres$max_lag <- max_lag

  # fit models:
  fits_rota_gravity_vary_max_lag$geom[[max_lag]] <- profile_par_lag(rotaBE, ctrl_gravity)
  fits_rota_gravity_vary_max_lag$pois[[max_lag]] <- profile_par_lag(rotaBE, ctrl_gravity_pois)

  # use grid-based optimization for linear lag:
  fits_rota_gravity_vary_max_lag$lin[[max_lag]] <- fit_par_lag(rotaBE, ctrl_gravity_lin,
                                                               range_par = vals_par_lag_lin,
                                                               use_update = FALSE)$best_mod

  if(max_lag > 2){
    fits_rota_gravity_vary_max_lag$unres[[max_lag]] <-
      profile_par_lag(rotaBE, ctrl_gravity_unres,
                      start_par_lag = c(fits_rota_gravity_vary_max_lag$unres[[max_lag - 1]]$par_lag, -2))
  }else{
    fits_rota_gravity_vary_max_lag$unres[[max_lag]] <-
      profile_par_lag(rotaBE, ctrl_gravity_unres)
  }

  print(max_lag)
}

# evaluate AICs:
AICs_gravity_vary_max_lag <- matrix(NA, ncol = 4, nrow = 7,
                                    dimnames = list(NULL, c("pois", "lin", "geom", "unres")))

for(lag_structure in c("pois", "lin", "geom", "unres")){
  for(max_lag in 1:7){
    AICs_gravity_vary_max_lag[max_lag, lag_structure] <-
      AIC(fits_rota_gravity_vary_max_lag[[lag_structure]][[max_lag]])
  }
}

# store:
save(fits_rota_gravity_vary_max_lag, file = paste0("model_fits/rota_", n_seas*52 ,"/fits_rota_gravity_vary_max_lag_", n_seas*52, ".rda"))
write.csv(AICs_gravity_vary_max_lag, file = paste0("AIC/AIC_rota_gravity_vary_max_lag_", n_seas*52, ".csv"))
