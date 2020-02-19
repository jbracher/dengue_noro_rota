# Analyse the log scores for the dengue forecasts

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("dengue")

library(surveillance)
names_lag_structures <- c("ar1", "pois", "lin", "geom", "unres", "siraj")

source("../basic_settings.R")


# compute the mean log scores for our forecasts:
logS_dengue <- list()
mean_logS_dengue <- matrix(NA, ncol = 6, nrow = 8,
                           dimnames = list(paste0("h", 1:8),
                                           c("ar1", "pois", "lin", "geom",
                                             "unres", "siraj")))

# run through different lag structures:
for(lag_structure in colnames(mean_logS_dengue)){
  # read in results (generated in evaluate_logS_dengue.R)
  logS_dengue_temp <- read.csv(paste0("logS/logS_dengue_", lag_structure, ".csv"))

  logS_dengue[[lag_structure]] <- matrix(nrow = 4*52, ncol = 8,
                                         dimnames = list(paste("t=", (19*52 + 1):(23*52)),
                                                         paste0("h", 1:8)))
  # extract scores and for different horizons and store in list logS_dengue
  for(horizon in 1:8){
    logS_dengue[[lag_structure]][, paste0("h", horizon)] <-
      logS_dengue_temp[seq(from = 10 - horizon, length.out = 4*52), paste0("h", horizon)]
  }

  # evaluate and store mean logS
  mean_logS_dengue[, lag_structure] <- colMeans(logS_dengue[[lag_structure]])
}

# meaningful column names:
colnames(mean_logS_dengue) <- paste0("hhh4_", colnames(mean_logS_dengue))

# get Ray's results for hhh4 for comparison
# (requires downloading repository from https://github.com/reichlab/article-disease-pred-with-kcde):
results_ray_hhh4 <- readRDS("/home/johannes/Documents/Ray_article/article-disease-pred-with-kcde-master/inst/results/dengue_sj/prediction-results/surveillance-predictions.rds")
# extract and store in list:
scores_ray_hhh4 <- NA*logS_dengue$ar1
for(i in 1:8){
  scores_ray_hhh4[, i] <- subset(results_ray_hhh4, prediction_horizon == i)$log_score
}

# plot to see agreement:
plot(logS_dengue$ar1[, 1], type = "l")
lines(-scores_ray_hhh4[, 1], col = "red")

plot(logS_dengue$ar1[, 2], type = "l")
lines(-scores_ray_hhh4[, 2], col = "red")

plot(logS_dengue$ar1[, 3], type = "l")
lines(-scores_ray_hhh4[, 3], col = "red")

plot(colMeans(logS_dengue$ar1), type = "l")
lines(-colMeans(scores_ray_hhh4), col = "red")
lines(colMeans(logS_dengue$geom))


# get results from KCDE methods:
# (requires downloading repository from https://github.com/reichlab/article-disease-pred-with-kcde):
results_kcde <- readRDS("/home/johannes/Documents/Ray_article/article-disease-pred-with-kcde-master/inst/results/dengue_sj/prediction-results/kcde-predictions.rds")

# get mean scores:
mean_logS_kcde.0 <- aggregate(results_kcde$log_score,
                              by = list(bw_parameterization = results_kcde$bw_parameterization,
                                        seasonality = results_kcde$seasonality,
                                        prediction_horizon = results_kcde$prediction_horizon),
                              FUN = mean)
mean_logS_kcde.0$specification <- paste0(mean_logS_kcde.0$bw_parameterization, ".",
                                         ifelse(mean_logS_kcde.0$seasonality, "seas", "non_seas"))
mean_logS_kcde <- reshape(mean_logS_kcde.0[, c("prediction_horizon", "specification", "x")], timevar = "specification",
                          idvar = "prediction_horizon", direction = "wide")
colnames(mean_logS_kcde) <- gsub(colnames(mean_logS_kcde), pattern = "x.", replacement = "kcde.")
mean_logS_kcde <- subset(mean_logS_kcde, prediction_horizon <= 8)

# get scores of naive forecasting method and GLARMA:
results_naive <- read.csv("logS/logS_dengue_naive_glmnb.csv")
row.names(results_naive) <- results_naive[, 1]; results_naive <- results_naive[, -1]

results_glarma <- read.csv("logS/logS_dengue_glarma.csv")
row.names(results_glarma) <- results_glarma[, 1]; results_glarma <- results_glarma[, -1]

# select the same weeks as in hhh4 and KCDE:
logS_dengue_glarma <- logS_dengue_naive <-
  matrix(nrow = 4*52, ncol = 8,
         dimnames = list(paste("t=", (19*52 + 1):(23*52)),
                         paste0("h", 1:8)))
for(horizon in 1:8){
  logS_dengue_naive[, paste0("h", horizon)] <-
    results_naive[seq(from = 10 - horizon, length.out = 4*52), paste0("h", horizon)]
  logS_dengue_glarma[, paste0("h", horizon)] <-
    results_glarma[seq(from = 10 - horizon, length.out = 4*52), paste0("h", horizon)]
}

# check agreement:
plot(logS_dengue_glarma[, 1], type = "l")
lines(logS_dengue$ar1[, 1], col = "red")
# work together

mean_logS_naive <- colMeans(logS_dengue_naive, na.rm = TRUE)
mean_logS_glarma <- colMeans(logS_dengue_glarma, na.rm = TRUE)


# put everything together:
summary_logS_dengue <- cbind(prediction_horizon = 1:8,
                      -1*mean_logS_kcde[, -1],
                      mean_logS_dengue,
                      naive = mean_logS_naive,
                      glarma = mean_logS_glarma)
# and write out:
write.csv(summary_logS_dengue, file = "logS/summary_logS_dengue.csv", row.names = FALSE)


# compute p-values vs KCDE:

set.seed(333) # set.seed as this is a random procedure

pvals_vs_kcde_dengue <- matrix(nrow = 8, ncol = 8,
                        dimnames = list(paste0("h.", 1:8),
                                        c("prediction_horizon", "hhh4_ar1", "hhh4_geom",
                                          "hhh4_pois", "hhh4_lin", "hhh4_unres",
                                          "hhh4_siraj", "glarma")))

pvals_vs_kcde_dengue[, "prediction_horizon"] <- 1:8

# run over hhh4 methods:
for (lag_structure in c("ar1", "geom", "pois", "lin", "unres", "siraj")) {
  for(h in 1:8){
    pvals_vs_kcde_dengue[h, paste0("hhh4_", lag_structure)] <-
      permutationTest(logS_dengue[[lag_structure]][, h],
                      -results_kcde$log_score[results_kcde$seasonality == TRUE &
                                                results_kcde$bw_parameterization == "full" &
                                                results_kcde$prediction_horizon == h])$pVal.permut
  }
}

# glarma:
for(h in 1:8){
  pvals_vs_kcde_dengue[h, "glarma"] <-
    permutationTest(logS_dengue_glarma[, h],
                    -results_kcde$log_score[results_kcde$seasonality == TRUE &
                                              results_kcde$bw_parameterization == "full" &
                                              results_kcde$prediction_horizon == h])$pVal.permut
}


# write out results
write.csv(pvals_vs_kcde_dengue, file = "logS/pvals_vs_kcde_dengue.csv", row.names = FALSE)



# plot as in manuscript:

# ned to adapt naming a bit:
summary_logS_dengue <- summary_logS_dengue[, c("kcde.full.seas",
                                               "hhh4_ar1", "hhh4_pois", "hhh4_lin",
                                               "hhh4_geom", "hhh4_unres", "hhh4_siraj",
                                               "naive", "glarma")]
# adapt column names:
colnames(summary_logS_dengue)[1] <- "kcde"
colnames(summary_logS_dengue) <- gsub("hhh4_", "", colnames(summary_logS_dengue))

delta_logS_dengue <- summary_logS_dengue - summary_logS_dengue[, "kcde"]

# adapt column names of p-values:
colnames(pvals_vs_kcde_dengue) <- gsub("hhh4_", "", colnames(pvals_vs_kcde_dengue))

par(mfrow = c(1, 2), font.main = 1, family = "sans", mar = c(3.8, 5, 2.5, 1), las=1)
matplot(summary_logS_dengue,
        type = "b",
        col = cols_models_dengue[colnames(summary_logS_dengue)], lty = 1,
        pch = 16, cex = 0.5, ylim = c(3.8, 6.4), xlab = "forecast horizon",
        ylab = "mean logS")

legend("topleft", legend = c("first-order", "Poisson", "triangular",
                             "geometric", "unrestricted", "literature",
                             "KCDE", "GLARMA", "naive"),
       col = cols_models_dengue,
       lty = 1, bty = "n", ncol = 2, cex = 0.75, lwd = 2)


delta_logS_dengue <- summary_logS_dengue - summary_logS_dengue[, "kcde"]

matplot(delta_logS_dengue,
        type = "b", col = cols_models_dengue[colnames(delta_logS_dengue)], lty = 1,
        pch = c(16, rep(NA, 8)), cex = 0.5, ylim = c(-0.1, 0.3), xlab = "forecast horizon",
        ylab = "difference in mean logS", main = " ")

# add points with variable sizes (not fo naive and kcde):
for(lag_structure in names_models_dengue[-which(names_models_dengue %in% c("kcde", "naive"))]){
  point_sizes_temp <- rep(0.5, 8)
  point_sizes_temp[pvals_vs_kcde_dengue[, lag_structure] < 0.1] <- 0.85
  point_sizes_temp[pvals_vs_kcde_dengue[, lag_structure] < 0.01] <- 1

  pch_temp <- rep(16, 8)
  pch_temp[pvals_vs_kcde_dengue[, lag_structure] < 0.1] <- 1
  pch_temp[pvals_vs_kcde_dengue[, lag_structure] < 0.01] <- 16

  points(1:8,
         delta_logS_dengue[, paste0(lag_structure)],
         col = cols_models_dengue[lag_structure],
         cex = point_sizes_temp,
         pch = pch_temp)
}

legend("bottom", legend = c("p > 0.1", "p < 0.1", "p < 0.01"), pch = c(16, 1, 16),
       pt.cex = c(0.5, 0.85, 1), ncol = 3, bty = "n", cex = 0.8)
text("p-values against KCDE:", x = 3.9, y = -0.07, cex = 0.9)

