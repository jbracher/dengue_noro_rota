# evaluate the log scores of one-step-ahead forecasts of rotavirus
# this can be done without simulation

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("rota")

library(surveillance)
library(hhh4addon)

# get data:
data("rotaBE")
# modify neighbourhood matrices to apply power law:
rotaBE@neighbourhood <- rotaBE@neighbourhood + 1

names_model_versions <- c("full", "gravity")
names_lag_structures <- c("ar1", "pois", "lin", "geom", "unres")
names_districts <- colnames(rotaBE@observed)

# timepoints for which models were fitted (in fit_models_dengue.R):
tps <- (4*52 - 1):(7*52 - 1)

# initialize list for storing:
osa_rota <- list()
template_pred <- matrix(NA, ncol = ncol(rotaBE@observed), nrow = length(tps),
                        dimnames = list(tps, colnames(rotaBE@observed)))
template_psi <- matrix(NA, ncol = 1, nrow = length(tps),
                       dimnames = list(tps, NULL))


# evaluate log scores:
for(model_version in names_model_versions){
  for(lag_structure in names_lag_structures){

    if(lag_structure == "end" & model_version != "full") next() # end model only available for "full"

    osa_rota[[model_version]][[lag_structure]] <-
      list(pred = template_pred,
           observed = rotaBE@observed[tps + 1, ],
           psi = template_psi)
    for(ind in tps){
      # get fit:
      file_name_fit <- paste0("model_fits/rota_", model_version, "_", lag_structure,
                              "/fit_rota_", model_version, "_", lag_structure,"_", ind, ".rda")
      load(file_name_fit)
      fit_temp <- get(paste0("fit_rota_", model_version, "_", lag_structure, "_temp"))

      # one-step-ahead forecast:
      if(lag_structure == "ar1"){
        osa_temp <- oneStepAhead(fit_temp, tp = rep(ind, 2), type = "final")
      }else{
        suppressMessages(
          osa_temp <- oneStepAhead_hhh4lag(fit_temp, tp = rep(ind, 2), type = "final")
        )
      }

      # store:
      osa_rota[[model_version]][[lag_structure]]$pred[as.character(ind), ] <-
        osa_temp$pred
      osa_rota[[model_version]][[lag_structure]]$psi[as.character(ind), ] <-
        osa_temp$psi

    }

    # write out scores:
    # write.csv(scores(osa_rota[[model_version]][[lag_structure]]),
    #           file = paste0("logS/logS_rota_", model_version, "_", lag_structure, "_h1.csv"))
  }
}



# create summary:
logS_rota <- matrix(NA, ncol = length(names_model_versions), nrow = length(names_lag_structures),
                    dimnames = list(names_lag_structures, names_model_versions))
for(model_version in names_model_versions){
  for(lag_structure in names_lag_structures){
    logS_rota[lag_structure, model_version] <-
      mean(scores(osa_rota[[model_version]][[lag_structure]])[, "logs"])
  }
}

# add results for GLARMA models:
logS_glarma_rota <- read.csv("logS/multiv_logS_rota_glarma.csv")
logS_rota <- rbind(logS_rota, glarma = NA)
logS_rota["glarma", "full"] <- mean(logS_glarma_rota[, "h1"], na.rm = TRUE)


# add result for naive seasonal approach:
logS_naive_rota <- read.csv("logS/multiv_logS_rota_naive_glmnb.csv")
logS_rota <- rbind(logS_rota, naive = NA)
logS_rota["naive", "full"] <- mean(logS_naive_rota[, "h1"], na.rm = TRUE)

# write.csv(logS_rota, file = "logS/summary_osa_logS_rota.csv")


# perform permutation tests gravity vs. full model:
pvals_permut_rota_gravity_vs_full <- matrix(NA, nrow = length(names_lag_structures), ncol = 1,
                                            dimnames = list(names_lag_structures,
                                                            "h1"))

set.seed(123)

for(lag_structure in rownames(pvals_permut_rota_gravity_vs_full)){
  pvals_permut_rota_gravity_vs_full[lag_structure, ] <-
    permutationTest(scores(osa_rota$gravity[[lag_structure]])[, "logs"],
                    scores(osa_rota$full[[lag_structure]])[, "logs"])$pVal.permut
}

# write.csv(pvals_permut_rota_gravity_vs_full, file = "logS/pvals_logS_gravity_vs_full_rota.csv")



# and perform permutation tests vs. the geometric lag model:
pvals_permut_rota_vs_geom <- matrix(NA, nrow = 4, ncol = 2,
                                    dimnames = list(c("ar1", "pois", "lin", "unres"),
                                                    names_model_versions))

set.seed(123)

for(lag_structure in rownames(pvals_permut_rota_vs_geom)){
  for (model_version in colnames(pvals_permut_rota_vs_geom)){
    pvals_permut_rota_vs_geom[lag_structure, model_version] <-
      permutationTest(scores(osa_rota[[model_version]][[lag_structure]])[, "logs"],
                      scores(osa_rota[[model_version]]$geom)[, "logs"])$pVal.permut
  }
}

# add p-value for glarma:
# get glarma results:
mutliv_logS_glarma_rota <- read.csv("logS/multiv_logS_rota_glarma.csv")
multiv_logS_glarma_h1_rota <- mutliv_logS_glarma_rota$h1[!is.na(mutliv_logS_glarma_rota$h1)]

# check that aligned:
plot(multiv_logS_glarma_h1_rota, type = "l")
lines(scores(osa_rota[[model_version]][[lag_structure]])[, "logs"], col = "red")
# works

# compute and add p-values:
pvals_permut_rota_vs_geom <- rbind(pvals_permut_rota_vs_geom, glarma = NA)
pvals_permut_rota_vs_geom["glarma", "full"] <-
  permutationTest(multiv_logS_glarma_h1_rota,
                  scores(osa_rota[["full"]]$geom)[, "logs"])$pVal.permut

pvals_permut_rota_vs_geom["glarma", "gravity"] <-
  permutationTest(multiv_logS_glarma_h1_rota,
                  scores(osa_rota[["gravity"]]$geom)[, "logs"])$pVal.permut

# write.csv(pvals_permut_rota_vs_geom, file = "logS/pvals_logS_vs_geom_rota.csv")
