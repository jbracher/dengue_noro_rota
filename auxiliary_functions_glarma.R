# function to obtain forecasts at higher horizons using the glarma pakage:
forecasting_nStepAhead_glarma <- function(y, X, t_condition, max_horizon,
                                          phiLags, type = "NegBin", niter = 100){
  # fit glarma model:
  glarmafit <- glarma(y = y[1:t_condition], X = X[1:t_condition, ], type = type, phiLags = phiLags)
  fitted_size <- coef(glarmafit, type = "NB")
  # glarma forecasting function generates sample paths.
  # generate many of these and retain the sampled conditional expectations mu:
  sampled_mu <- sampled_Y <- matrix(nrow = niter, ncol = max_horizon)
  for(i in 1:niter){
    forecast_temp <- forecast(glarmafit, n.ahead = max_horizon,
                              newdata = X[t_condition + (1:max_horizon), , drop = FALSE],
                              newoffset = matrix(0, nrow = max_horizon))
    sampled_mu[i, ] <- forecast_temp$mu
    sampled_Y[i, ] <- forecast_temp$Y
  }
  # then perform Rao-Blackwellization, i.e. average over predictive probability masses
  # associated with the different sample paths:
  pred_mean <- pred_var <-
    lb50 <- ub50 <- lb95 <- ub95 <-
    pit_l <- pit_u <-
    log_scores <- unit_wise_log_score <- obs <- numeric(max_horizon)

  for(h in 1:max_horizon){ # run over different horizons
    # important: average probability masses *before* log transformation:
    obs[h] <- y[t_condition + h]
    log_scores[h] <- -log(mean(dnbinom(obs[h],
                                       mu = sampled_mu[, h],
                                       size = fitted_size)))


    # choose appropriate support:
    support_temp <- 0:qnbinom(p = 0.995, mu = max(c(obs[h], sampled_mu[, h]), na.rm = TRUE),
                              size = fitted_size)
    pred_densities_temp <- numeric(length(support_temp))

    for(i in 1:niter){
      pred_densities_temp <- pred_densities_temp +
        1/niter*dnbinom(support_temp, mu = sampled_mu[i, h],
                        size = fitted_size)
    }


    pred_mean[h] <- sum(support_temp*pred_densities_temp)
    pred_var[h] <- sum(support_temp^2*pred_densities_temp) - pred_mean[h]^2
    unit_wise_log_score[h] <- -log(pred_densities_temp[
      obs[h] + 1
      ])

    pred_cumul_distr_temp <- cumsum(pred_densities_temp)
    lb50[h] <- min(which(pred_cumul_distr_temp >= 0.25)) - 1 # -1 bc support includes 0
    ub50[h] <- max(which(pred_cumul_distr_temp <= 0.75)) # no -1  as we want to be slightly conservative
    lb95[h] <- min(which(pred_cumul_distr_temp >= 0.025)) - 1
    ub95[h] <- max(which(pred_cumul_distr_temp <= 0.975))

    if(!is.na(obs[h])){
      pit_l[h] <- if(obs[h] == 0){
        0
      }else{
        pred_cumul_distr_temp[obs[h]] # + 1 bc support includes 0
      }
      pit_u[h] <- pred_cumul_distr_temp[obs[h] + 1]
    }else{
      pit_l[h] <- pit_u[h] <- NA
    }
  }

  # return object:
  return(list(pred_mean = pred_mean, pred_var = pred_var,
              lb50 = lb50, ub50 = ub50,
              lb95 = lb95, ub95 = ub95,
              obs = obs,
              unit_wise_log_score = log_scores,
              unit_wise_pit_l = pit_l,
              unit_wise_pit_u = pit_u,
              multiv_log_score = NA))
}
