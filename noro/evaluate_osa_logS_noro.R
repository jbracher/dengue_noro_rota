# evaluate the log scores of one-step-ahead forecasts of norovirus
# this can be done without simulation

# setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/")
setwd("noro")

library(surveillance)
library(hhh4addon)

# get data:
data("noroBE")
# modify neighbourhood matrices to apply power law:
noroBE@neighbourhood <- noroBE@neighbourhood + 1

names_model_versions <- c("full", "gravity")
names_lag_structures <- c("ar1", "pois", "lin", "geom", "unres")
names_districts <- colnames(noroBE@observed)

# timepoints for which models were fitted (in fit_models_dengue.R):
tps <- (4*52 - 1):(7*52 - 1)

# initialize list for storing:
osa_noro <- list()
template_pred <- matrix(NA, ncol = ncol(noroBE@observed), nrow = length(tps),
                        dimnames = list(tps, colnames(noroBE@observed)))
template_psi <- matrix(NA, ncol = 1, nrow = length(tps),
                       dimnames = list(tps, NULL))


# evaluate log scores:
for(model_version in names_model_versions){
  for(lag_structure in names_lag_structures){

    if(lag_structure == "end" & model_version != "full") next() # end model only available for "full"

    osa_noro[[model_version]][[lag_structure]] <-
      list(pred = template_pred,
           observed = noroBE@observed[tps + 1, ],
           psi = template_psi)
    for(ind in tps){
      # get fit:
      file_name_fit <- paste0("model_fits/noro_", model_version, "_", lag_structure,
                              "/fit_noro_", model_version, "_", lag_structure,"_", ind, ".rda")
      load(file_name_fit)
      fit_temp <- get(paste0("fit_noro_", model_version, "_", lag_structure, "_temp"))

      # one-step-ahead forecast:
      if(lag_structure == "ar1"){
        osa_temp <- oneStepAhead(fit_temp, tp = rep(ind, 2), type = "final")
      }else{
        suppressMessages(
          osa_temp <- oneStepAhead_hhh4lag(fit_temp, tp = rep(ind, 2), type = "final")
        )
      }

      # store:
      osa_noro[[model_version]][[lag_structure]]$pred[as.character(ind), ] <-
        osa_temp$pred
      osa_noro[[model_version]][[lag_structure]]$psi[as.character(ind), ] <-
        osa_temp$psi

    }

    # write out scores:
    # write.csv(scores(osa_noro[[model_version]][[lag_structure]]),
    #           file = paste0("logS/logS_noro_", model_version, "_", lag_structure, "_h1.csv"))
  }
}



# create summary:
logS_noro <- matrix(NA, ncol = length(names_model_versions), nrow = length(names_lag_structures),
                    dimnames = list(names_lag_structures, names_model_versions))
for(model_version in names_model_versions){
  for(lag_structure in names_lag_structures){
    logS_noro[lag_structure, model_version] <-
      mean(scores(osa_noro[[model_version]][[lag_structure]])[, "logs"])
  }
}

# add results for GLARMA models:
logS_glarma_noro <- read.csv("logS/multiv_logS_noro_glarma.csv")
logS_noro <- rbind(logS_noro, glarma = NA)
logS_noro["glarma", "full"] <- mean(logS_glarma_noro[, "h1"], na.rm = TRUE)


# add result for naive seasonal approach:
logS_naive_noro <- read.csv("logS/multiv_logS_noro_naive_glmnb.csv")
logS_noro <- rbind(logS_noro, naive = NA)
logS_noro["naive", "full"] <- mean(logS_naive_noro[, "h1"], na.rm = TRUE)

# write.csv(logS_noro, file = "logS/summary_osa_logS_noro.csv")


# perform permutation tests gravity vs. full model:
pvals_permut_noro_gravity_vs_full <- matrix(NA, nrow = length(names_lag_structures), ncol = 1,
                                            dimnames = list(names_lag_structures,
                                                            "h1"))

set.seed(123)

for(lag_structure in rownames(pvals_permut_noro_gravity_vs_full)){
  pvals_permut_noro_gravity_vs_full[lag_structure, ] <-
    permutationTest(scores(osa_noro$gravity[[lag_structure]])[, "logs"],
                    scores(osa_noro$full[[lag_structure]])[, "logs"])$pVal.permut
}

# write.csv(pvals_permut_noro_gravity_vs_full, file = "logS/pvals_logS_gravity_vs_full_noro.csv")



# and perform permutation tests vs. the geometric lag model:
pvals_permut_noro_vs_geom <- matrix(NA, nrow = 4, ncol = 2,
                                    dimnames = list(c("ar1", "pois", "lin", "unres"),
                                                    names_model_versions))

set.seed(123)

for(lag_structure in rownames(pvals_permut_noro_vs_geom)){
  for (model_version in colnames(pvals_permut_noro_vs_geom)){
    pvals_permut_noro_vs_geom[lag_structure, model_version] <-
      permutationTest(scores(osa_noro[[model_version]][[lag_structure]])[, "logs"],
                      scores(osa_noro[[model_version]]$geom)[, "logs"])$pVal.permut
  }
}

# add p-value for glarma:
# get glarma results:
mutliv_logS_glarma_noro <- read.csv("logS/multiv_logS_noro_glarma.csv")
multiv_logS_glarma_h1_noro <- mutliv_logS_glarma_noro$h1[!is.na(mutliv_logS_glarma_noro$h1)]

# check that aligned:
plot(multiv_logS_glarma_h1_noro, type = "l")
lines(scores(osa_noro[[model_version]][[lag_structure]])[, "logs"], col = "red")
# works

# compute and add p-values:
pvals_permut_noro_vs_geom <- rbind(pvals_permut_noro_vs_geom, glarma = NA)
pvals_permut_noro_vs_geom["glarma", "full"] <-
  permutationTest(multiv_logS_glarma_h1_noro,
                  scores(osa_noro[["full"]]$geom)[, "logs"])$pVal.permut

pvals_permut_noro_vs_geom["glarma", "gravity"] <-
  permutationTest(multiv_logS_glarma_h1_noro,
                  scores(osa_noro[["gravity"]]$geom)[, "logs"])$pVal.permut

# write.csv(pvals_permut_noro_vs_geom, file = "logS/pvals_logS_vs_geom_noro.csv")
