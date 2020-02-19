# compare the AICs of the different fitted models with the respective chosen order
# and make some descriptive plots

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("rota")

# load libraries, set colours and get some auxiliary functions:
source("../auxiliary_functions.R")
source("../basic_settings.R")


# get full models from fits to choose order:
load("model_fits/rota_364/fits_rota_full_vary_max_lag_364.rda")
load("model_fits/rota_364/fits_rota_gravity_vary_max_lag_364.rda")

# choose models with respective orders:
fits_rota <- list(
  full = list(ar1 = fits_rota_full_vary_max_lag$geom[[1]],
              geom = fits_rota_full_vary_max_lag$geom[[5]],
              pois = fits_rota_full_vary_max_lag$pois[[5]],
              lin = fits_rota_full_vary_max_lag$lin[[5]],
              unres = fits_rota_full_vary_max_lag$unres[[5]]),
  gravity = list(ar1 = fits_rota_gravity_vary_max_lag$geom[[1]],
                 geom = fits_rota_gravity_vary_max_lag$geom[[5]],
                 pois = fits_rota_gravity_vary_max_lag$pois[[5]],
                 lin = fits_rota_gravity_vary_max_lag$lin[[5]],
                 unres = fits_rota_gravity_vary_max_lag$unres[[5]])
)
fits <- list(rota = fits_rota) # code taken from Rnw where rota is stored in the same list

# save(fits_rota, file = "model_fits/rota_364/fits_rota_364.rda")

### compare AICs:
AICs <- matrix(NA, nrow = 5, ncol = 2,
               dimnames = list(c("ar1", "pois", "lin", "geom", "unres"),
                               c("gravity", "full")))

for(model_version in c("gravity", "full")){
  for(lag_structure in rownames(AICs)){
    AICs[lag_structure, model_version] <- AIC(fits_rota[[model_version]][[lag_structure]])
  }
}

# store:
# write.csv(AICs, file = "AIC/AICs_rota_364.csv")

### Descriptive plot as in manuscript:

nu_seas <- overdisp <- power_param <- list()

disease <- "rota"
power_param[[disease]] <- overdisp[[disease]] <-
  matrix(nrow = length(names_hhh4_lag_structures_nr), ncol = length(names_hhh4_model_versions_nr),
         dimnames = list(names_hhh4_lag_structures_nr, names_hhh4_model_versions_nr))
for(model_version in names_hhh4_model_versions_nr){
  for(lag_structure in names_hhh4_lag_structures_nr){
    nu_seas[[disease]][[model_version]][[lag_structure]] <-
      get_nu_seas(fits[[disease]][[model_version]][[lag_structure]], regions = names_districts,
                  christmas = c(1, rep(0, 50), 1))
    overdisp[[disease]][lag_structure, model_version] <-
      exp(-fits[[disease]][[model_version]][[lag_structure]]$coefficients["-log(overdisp)"])
    power_param[[disease]][lag_structure, model_version] <-
      fits[[disease]][[model_version]][[lag_structure]]$coefficients["neweights.d"]
  }
}

# Plot:

par(font.main = 1, family = "sans", mar = c(4, 5, 3, 1), las = 1)
layout_matr <- matrix(c(1, 1, 2, 2, 3, 4), byrow = TRUE, ncol = 6)
layout(layout_matr)

### rota

# Lag weights
lwd_weights <- 2
plot(1 - 0.2, 1, type = "h", ylim = 0:1, xlim = c(0.9, 5.1),
     # ylab = expression(lag~weight~group(lfloor, hat(u)[d], rfloor)),
     ylab = "weight",
     xlab = "lag", main = "", col = cols_models_nr["ar1"], lwd = lwd_weights)
points(seq_along(fits$rota$full$pois$distr_lag) - 0.1,
       fits$rota$full$pois$distr_lag, type = "h", col = cols_models_nr["pois"], lwd = lwd_weights)
points(1:10 + 0.0, fits$rota$full$lin$distr_lag[1:10], type = "h", col = cols_models_nr["lin"], lwd = lwd_weights)
points(seq_along(fits$rota$full$geom$distr_lag) + 0.1,
       fits$rota$full$geom$distr_lag, type = "h", col = cols_models_nr["geom"], lwd = lwd_weights)
points(seq_along(fits$rota$full$unres$distr_lag) + 0.2,
       fits$rota$full$unres$distr_lag, type = "h", col = cols_models_nr["unres"], lwd = lwd_weights)
legend("topright", col = cols_models_nr, lty = 1,
       legend = c("first-order", "Poisson", "triangular", "geometric", "unrestricted"),
       bty = "n", lwd = 2)

# total endemic component
order_seas_rota <- c(27:52, 1:26)

plot(colSums(nu_seas$rota$full$ar1)[order_seas_rota], type = "l", ylim = c(0, 47), col = cols_models_nr["ar1"],
     xlab = "calendar week", ylab = "total endemic component", main = "(a) rotavirus",
     axes = FALSE)
axis(1, at = c(1, 13, 26, 39, 52), labels = order_seas_rota[c(1, 13, 26, 39, 52)])
axis(2); box()
# text(-18, 25, expression(sum(nu[it], g, 1)), xpd=TRUE, las = 2, srt = 90)
lines(colSums(nu_seas$rota$full$lin)[order_seas_rota], col = cols_models_nr["lin"])
lines(colSums(nu_seas$rota$full$geom)[order_seas_rota], col = cols_models_nr["geom"])
lines(colSums(nu_seas$rota$full$pois)[order_seas_rota], col = cols_models_nr["pois"], lty = "dashed")
lines(colSums(nu_seas$rota$full$unres)[order_seas_rota], col = cols_models_nr["unres"], lty = "dashed")
# or add weights? for 1--4

# power parameter
plot(c(1, 0.9, 1.1, 0.9, 1.1), power_param$rota[, "full"], axes = FALSE,
     xlim = c(0.7, 1.3), ylim = c(0.9, 1.25),
     col = cols_models_nr[rownames(power_param$rota)], pch = 15,
     # ylab = expression(power~law~decay~hat(rho)),
     ylab = "power law decay",
     xlab = "")
box()
axis(2)

# overdispersion
plot(c(1, 0.9, 1.1, 0.9, 1.1), overdisp$rota[, "full"], axes = FALSE,
     xlim = c(0.7, 1.3), ylim = c(0.170, 0.183),
     col = cols_models_nr[rownames(overdisp$rota)], pch = 15,
     # ylab = expression(overdisp.~parameter~hat(psi)),
     ylab = "overdispersion",
     xlab = "")
box()
axis(2)

