# Generate descriptive graphics for dengue analysis:

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("dengue")

# get some helper functions:
source("../auxiliary_functions.R")
source("../basic_settings.R")

library(surveillance)
library(hhh4addon)

# get and plot data:
data("dengueSJ")
plot(dengueSJ)

# fit models:

# define controls
ctrls_dengue <- list()
ctrls_dengue$ar1 <- list(
  ar = list(f = addSeason2formula(f = ~ 1, S = 2), lag = 1),
  end = list(f = addSeason2formula(f = ~ 1, S = 1, period = 52)),
  subset = 11:(19*52),
  family = "NegBin1"
)
# for different weighting schemes:
# the max_lag values come from choose_order_dengue.R
ctrls_dengue$geom <- ctrls_dengue$ar1; ctrls_dengue$pois$funct_lag <- geometric_lag; ctrls_dengue$geom$max_lag <- 5
ctrls_dengue$pois <- ctrls_dengue$ar1; ctrls_dengue$pois$funct_lag <- poisson_lag; ctrls_dengue$pois$max_lag <- 5
ctrls_dengue$ar2 <- ctrls_dengue$ar1; ctrls_dengue$ar2$funct_lag <- ar2_lag; ctrls_dengue$pois$max_lag <- 5
ctrls_dengue$lin <- ctrls_dengue$ar1; ctrls_dengue$lin$funct_lag <- linear_lag; ctrls_dengue$lin$max_lag <- 5
ctrls_dengue$unres <- ctrls_dengue$ar1; ctrls_dengue$unres$funct_lag <- unrestricted_lag; ctrls_dengue$unres$max_lag <- 4

# based on serial interval:
wgts_siraj <- c(0.001, 0.999*c(0.2, 0.425, 0.25, 0.125)) # cannot asign probability 0 to lag 1
par_lag_siraj <- log(wgts_siraj[-1]/(1 - sum(wgts_siraj[-1])))
ctrls_dengue$siraj <- ctrls_dengue$geom; ctrls_dengue$siraj$funct_lag <- unrestricted_lag
ctrls_dengue$siraj$par_lag <- par_lag_siraj

# fit models:
fits_dengue <- list()

fits_dengue$ar1 <- hhh4(dengueSJ,
                        control = ctrls_dengue$ar1)

fits_dengue$ar2 <- profile_par_lag(dengueSJ,
                                   control = ctrls_dengue$ar2)

fits_dengue$geom <- profile_par_lag(dengueSJ,
                                    control = ctrls_dengue$geom)

fits_dengue$pois <- profile_par_lag(dengueSJ,
                                    control = ctrls_dengue$pois)

fits_dengue$lin <- profile_par_lag(dengueSJ,
                                   control = ctrls_dengue$lin)

fits_dengue$unres <- profile_par_lag(dengueSJ,
                                     control = ctrls_dengue$unres)

fits_dengue$siraj <- hhh4_lag(dengueSJ, control = ctrls_dengue$siraj)

# compute AICs:
AICs_dengue <- lapply(fits_dengue, AIC)

# read in the AICs with varying orders as computed in choose_order_dengue:
AICs_vary_max_lag_dengue <- read.csv("AIC/AICs_dengue_vary_max_lag.csv")
ref <- AICs_vary_max_lag_dengue[1, "geom"]
AICs_vary_max_lag_dengue <- AICs_vary_max_lag_dengue - ref

# create figure:

par(mfrow = c(1, 3), las = 1, mar = c(4, 4, 0.5, 1))

# plot AICs for different values of p:
plot(2:7, AICs_vary_max_lag_dengue[2:7, "geom"], type = "b", col = cols_models_dengue["geom"],
     xlab = "p", ylab  = "improvement in AIC", pch = 15, cex = 0.9, ylim = c(6555, 6605)-ref-3)
lines(2:7, AICs_vary_max_lag_dengue[2:7, "pois"], col = cols_models_dengue["pois"],
      type = "b", pch = 15, cex = 0.9)
lines(2:7, AICs_vary_max_lag_dengue[2:7, "lin"], col = cols_models_dengue["lin"],
      type = "b", pch = 15, cex = 0.9)
lines(2:7, AICs_vary_max_lag_dengue[2:7, "unres"], col = cols_models_dengue["unres"],
      type = "b", pch = 15, cex = 0.9)


# plot lag weights:
lwd_weights <- 2
plot(1 - 0.2, 1, type = "h", ylim = 0:1, xlim = c(0.9, 5.1),
     # ylab = expression(lag~weight~group(lfloor, hat(u)[d], rfloor)),
     ylab = "weight",
     xlab = "lag / serial interval", main = "", col = cols_models_dengue["ar1"], lwd = lwd_weights)
points(1:6 - 0.1, fits_dengue$pois$distr_lag[1:6], type = "h", col = cols_models_dengue["pois"], lwd = lwd_weights)
points(1:6 + 0, fits_dengue$lin$distr_lag[1:6], type = "h", col = cols_models_dengue["lin"], lwd = lwd_weights)
points(1:6 + 0.1, fits_dengue$geom$distr_lag[1:6], type = "h", col = cols_models_dengue["geom"], lwd = lwd_weights)
points(1:6 + 0.2, fits_dengue$unres$distr_lag[1:6], type = "h", col = cols_models_dengue["unres"], lwd = lwd_weights)
points(1:6 + 0.2, fits_dengue$siraj$distr_lag[1:6], type = "h", col = cols_models_dengue["siraj"], lwd = lwd_weights)


legend("topright", col = cols_models_dengue[c("fixed", "pois", "lin", "geom", "unres", "siraj")],
       lty = 1, legend = c("fixed", "Poisson",
                           "triangular", "geometric",
                           "unrestricted", "siraj"),
       bty = "n", lwd = 2)

# plot acf of Pearson residuals:
cond_pearson_resids_dengue <- acfs_cond_pr_dengue <-
  means_cond_pr_dengue <- vars_cond_pr_dengue <- list()
for(lag_structure in names(fits_dengue)){
  size_temp <- exp(-fits_dengue[[lag_structure]]$coefficients["-log(overdisp)"])
  fitted_values_temp <- fits_dengue[[lag_structure]]$fitted.values
  cond_sds_temp <- sqrt(fitted_values_temp + size_temp*fitted_values_temp^2)
  obs_values_temp <- fits_dengue[[lag_structure]]$stsObj@observed[ctrls_dengue$ar1$subset, ]
  pearson_resids_temp <- (fitted_values_temp - obs_values_temp)/cond_sds_temp
  cond_pearson_resids_dengue[[lag_structure]] <- pearson_resids_temp
  acfs_cond_pr_dengue[[lag_structure]] <- my_acf(pearson_resids_temp)
  means_cond_pr_dengue[[lag_structure]] <- mean(pearson_resids_temp)
  vars_cond_pr_dengue[[lag_structure]] <- var(pearson_resids_temp)
}

myplot_acf(acfs_cond_pr_dengue$ar1, ylim = c(-0.3, 0.5), nobs = length(ctrls_dengue$ar1$subset),
           shift_x = -0.2, col = cols_models_dengue[1])
myplot_acf(acfs_cond_pr_dengue$pois, shift_x = -0.1, col = cols_models_dengue["pois"], add = TRUE)
myplot_acf(acfs_cond_pr_dengue$lin, shift_x = 0, col = cols_models_dengue["lin"], add = TRUE)
myplot_acf(acfs_cond_pr_dengue$geom, shift_x = 0.1, col = cols_models_dengue["geom"], add = TRUE)
myplot_acf(acfs_cond_pr_dengue$unres, shift_x = 0.2, col = cols_models_dengue["unres"], add = TRUE)
myplot_acf(acfs_cond_pr_dengue$siraj, shift_x = 0.2, col = cols_models_dengue["siraj"], add = TRUE)

