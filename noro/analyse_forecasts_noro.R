# Analyse the log scores for the norovirus forecasts

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("noro")

# load libraries, set colours and get some auxiliary functions:
source("../basic_settings.R")
source("../auxiliary_functions.R")

max_horizon <- 4

# read in log scores:
multiv_logS_full_ar1 <- read.csv("logS/multiv_logS_noro_full_ar1.csv")
multiv_logS_full_geom <- read.csv("logS/multiv_logS_noro_full_geom.csv")

multiv_logS_gravity_ar1 <- read.csv("logS/multiv_logS_noro_gravity_ar1.csv")
multiv_logS_gravity_geom <- read.csv("logS/multiv_logS_noro_gravity_geom.csv")

# multiv_logS_end <- read.csv("logS/multiv_logS_noro_full_end.csv")

multiv_logS_naive <- read.csv("logS/multiv_logS_noro_naive_glmnb.csv")

multiv_logS_glarma <- read.csv("logS/multiv_logS_noro_glarma.csv")

# put together in data.frame:
logS_horizons_noro <- data.frame(horizon = 1:4,
                   ar1 = colMeans(multiv_logS_full_ar1[, -1], na.rm = TRUE),
                   geom = colMeans(multiv_logS_full_geom[, -1], na.rm = TRUE),
                   gravity_ar1 = colMeans(multiv_logS_gravity_ar1[, -1], na.rm = TRUE),
                   gravity_geom = colMeans(multiv_logS_gravity_geom[, -1], na.rm = TRUE),
                   glarma = colMeans(multiv_logS_glarma[, -1], na.rm = TRUE),
                   naive = colMeans(multiv_logS_naive[, -1], na.rm = TRUE)
                   )

# store:
# write.csv(logS_horizons_noro, file = "logS/summary_multiv_logS_noro.csv", row.names = FALSE)

# apply permutation tests ar1 vs geom:
library(surveillance)

pvals_ar1_vs_geom <- matrix(ncol = 2, nrow = 4,
                            dimnames = list(c("h1", "h2", "h3", "h4"),
                                            c("full", "gravity")))

set.seed(124)

for(h in 1:4){
  inds_temp <- which(!is.na(multiv_logS_full_ar1[, paste0("h", h)]))
  pvals_ar1_vs_geom[h, "full"] <-
    permutationTest(multiv_logS_full_ar1[inds_temp, paste0("h", h)],
                    multiv_logS_full_geom[inds_temp, paste0("h", h)])$pVal.permut

  pvals_ar1_vs_geom[h, "gravity"] <-
    permutationTest(multiv_logS_gravity_ar1[inds_temp, paste0("h", h)],
                    multiv_logS_gravity_geom[inds_temp, paste0("h", h)])$pVal.permut

}

write.csv(pvals_ar1_vs_geom,
         file = "logS/pvals_logS_ar1_vs_geom_noro_horizons.csv", row.names = TRUE)



# apply permutation tests full vs gravity:
library(surveillance)

pvals_gravity_vs_full <- matrix(ncol = 2, nrow = 4,
                            dimnames = list(c("h1", "h2", "h3", "h4"),
                                            c("ar1", "geom")))

set.seed(123)

for(h in 1:4){
  inds_temp <- which(!is.na(multiv_logS_full_ar1[, paste0("h", h)]))

  pvals_gravity_vs_full[h, "ar1"] <-
    permutationTest(multiv_logS_full_ar1[inds_temp, paste0("h", h)],
                    multiv_logS_gravity_ar1[inds_temp, paste0("h", h)])$pVal.permut

  pvals_gravity_vs_full[h, "geom"] <-
    permutationTest(multiv_logS_full_geom[inds_temp, paste0("h", h)],
                    multiv_logS_gravity_geom[inds_temp, paste0("h", h)])$pVal.permut

}

# store:
write.csv(pvals_gravity_vs_full,
          file = "logS/pvals_logS_gravity_vs_full_noro_horizons.csv", row.names = TRUE)


# Plot as in manuscript:
plot(1:4, logS_horizons_noro$ar1, type = "b", col = cols_models_nr["ar1"], pch = 16, cex = 0.9, lwd=1.5,
     ylim = c(2.33, 2.5), xlab = "forecast horizon", ylab = "mean multivariate logS", axes = FALSE,
     main = "(a) Norovirus")
axis(1, at = 1:4); axis(2); box()
lines(1:4, logS_horizons_noro$geom, type = "b", col = cols_models_nr["geom"], pch = 16, cex = 0.9, lwd=1.5)

lines(1:4, logS_horizons_noro$gravity_ar1, type = "b", col = cols_models_nr["ar1"],
      pch = 17, cex = 0.9, lty = 2, lwd=1.5)
lines(1:4, logS_horizons_noro$gravity_geom, type = "b", col = cols_models_nr["geom"],
      pch = 17, cex = 0.9, lty = 2, lwd=1.5)

lines(1:4, logS_horizons_noro$naive, type = "b", col = cols_models_nr["naive"], cex = 0.7, pch = 15)

lines(1:4, logS_horizons_noro$glarma, type = "b", col = cols_models_nr["glarma"], cex = 0.7, pch = 15)


legend("bottomright", legend = c("first-order, gravity model",
                                 "first-order, full model",
                                 "geometric weights, gravity model",
                                 "geometric weights, full model",
                                 "GLARMA",
                                 "naive"),
       col = c(cols_models_nr[c("ar1", "ar1", "geom", "geom", "glarma", "naive")]),
       lty = c(2, 1, 2, 1, 1, 1), pch = c(17, 16, 17, 16, 15, 15), cex = 0.8, bty = "n")