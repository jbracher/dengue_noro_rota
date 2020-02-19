# compare different orders p for dengue models

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("dengue")

# load packages:
library(RColorBrewer)
library(surveillance)
library(hhh4addon)

# get some additional functions:
# source("../auxiliary_functions.R")

# get and plot data:
data(dengueSJ)

names_lag_structures <- c("ar1", "pois", "lin", "geom", "unres")
cols_lag_structures <- brewer.pal(8, "Dark2")[c(8, 4:7)]
names(cols_lag_structures) <- names_lag_structures

plot(dengueSJ)
abline(v = 19*52, col = "red", lty = 2)

# define controls for different lag weighting schemes:
ctrls_dengue <- list()
ctrls_dengue$ar1 <- list(
  ar = list(f = addSeason2formula(f = ~ 1, S = 2), lag = 1),
  end = list(f = addSeason2formula(f = ~ 1, S = 1, period = 52)),
  subset = 11:(19*52),
  family = "NegBin1"
)
ctrls_dengue$geom <- ctrls_dengue$ar1; ctrls_dengue$pois$funct_lag = geometric_lag; ctrls_dengue$geom$max_lag <- 10
ctrls_dengue$pois <- ctrls_dengue$ar1; ctrls_dengue$pois$funct_lag = poisson_lag; ctrls_dengue$pois$max_lag <- 10
ctrls_dengue$ar2 <- ctrls_dengue$ar1; ctrls_dengue$ar2$funct_lag = ar2_lag; ctrls_dengue$ar2$max_lag <- 2
ctrls_dengue$lin <- ctrls_dengue$ar1; ctrls_dengue$lin$funct_lag = linear_lag; ctrls_dengue$lin$max_lag <- 10
ctrls_dengue$unres <- ctrls_dengue$ar1; ctrls_dengue$unres$funct_lag = unrestricted_lag; ctrls_dengue$unres$max_lag <- 10


# fit models varying order p:

# p = 1:
fit_dengue_ar1 <- hhh4(dengueSJ, ctrls_dengue$ar1)

fits_dengue_vary_max_lag <- list()

fits_dengue_vary_max_lag$geom[[1]] <-
  fits_dengue_vary_max_lag$pois[[1]] <-
  fits_dengue_vary_max_lag$lin[[1]] <-
  fits_dengue_vary_max_lag$unres[[1]] <- fit_dengue_ar1

for(max_lag in 2:10){
  ctrls_dengue$geom$max_lag <- max_lag
  fits_dengue_vary_max_lag$geom[[max_lag]] <- profile_par_lag(dengueSJ, ctrls_dengue$geom)

  ctrls_dengue$pois$max_lag <- max_lag
  fits_dengue_vary_max_lag$pois[[max_lag]] <- profile_par_lag(dengueSJ, ctrls_dengue$pois)

  ctrls_dengue$lin$max_lag <- max_lag
  fits_dengue_vary_max_lag$lin[[max_lag]] <- profile_par_lag(dengueSJ, ctrls_dengue$lin)

  ctrls_dengue$unres$max_lag <- max_lag
  start_par_lag <- if(max_lag == 2){
    NULL
  }else{
    c(fits_dengue_vary_max_lag$unres[[max_lag - 1]]$par_lag, -3)
  }
  fits_dengue_vary_max_lag$unres[[max_lag]] <- profile_par_lag(dengueSJ, ctrls_dengue$unres, start_par_lag = start_par_lag)

  print(max_lag)
}

lapply(fits_dengue_vary_max_lag$unres, function(x) x$convergence_profile)
# warnings stem from unrestricted with high lags.

save(fits_dengue_vary_max_lag, file = "model_fits/fits_denge_vary_max_lag.rda")
load("model_fits/fits_denge_vary_max_lag.rda")

AICs_vary_max_lag_dengue <- matrix(ncol = 5, nrow = 10,
                                   dimnames = list(NULL, c("max_lag", "geom", "pois", "lin", "unres")))
AICs_vary_max_lag_dengue[, "max_lag"] <- 1:10
AICs_vary_max_lag_dengue[, "geom"] <- unlist(lapply(fits_dengue_vary_max_lag$geom, AIC))
AICs_vary_max_lag_dengue[, "pois"] <- unlist(lapply(fits_dengue_vary_max_lag$pois, AIC))
AICs_vary_max_lag_dengue[, "lin"] <- unlist(lapply(fits_dengue_vary_max_lag$lin, AIC))
AICs_vary_max_lag_dengue[, "unres"] <- unlist(lapply(fits_dengue_vary_max_lag$unres, AIC))



write.csv(AICs_vary_max_lag_dengue, file = "AIC/AICs_dengue_vary_max_lag.csv", row.names = FALSE)
