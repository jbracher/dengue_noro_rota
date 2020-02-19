# function to obtain naive seasonal forecasts as suggesed by Leo:
naive_forecast_glmnb <- function(ts, t_cond, max_horizon = 1, freq = 52){

  dummies_season <- as.factor(rep(1:freq, length.out = length(ts)))

  # get training data:
  subset_training <- 1:t_cond
  dat_training <- ts[subset_training]
  dummies_season_training <- dummies_season[subset_training]

  # fit model:
  fit_nb <- glm.nb(dat_training ~ dummies_season_training)

  # subset test:
  subset_test <- t_cond + (1:max_horizon)
  dummies_season_test <- dummies_season[subset_test]

  # predict:
  mu <- predict.glm(fit_nb,
                    newdata = data.frame(dummies_season_training = dummies_season_test),
                    type = "response")

  size <- fit_nb$theta

  pred_lb50 <- qnbinom(0.25, mu = mu, size = size)
  pred_ub50 <- qnbinom(0.75, mu = mu, size = size)
  pred_lb95 <- qnbinom(0.025, mu = mu, size = size)
  pred_ub95 <- qnbinom(0.975, mu = mu, size = size)

  obs <- ts[t_cond + (1:max_horizon)]

  unit_wise_log_score <- -dnbinom(x = obs, mu = mu, size = size, log = TRUE)
  multiv_log_score <- NA

  return(list(prediction_time = t_cond,
              prediction_horizon = 1:max_horizon,
              pred_mean = mu, pred_var = mu + mu^2/size,
              lb50 = pred_lb50, ub50 = pred_ub50,
              lb95 = pred_lb95, ub95 = pred_ub95,
              obs = obs,
              unit_wise_log_score = unit_wise_log_score,
              multiv_log_score = multiv_log_score))
}

# function to obtain naive seasonal forecasts as suggesed by Leo:
naive_forecast_glmnb_multiv <- function(ts, t_cond, max_horizon = 1, freq = 52){

  n_units <- ncol(ts)
  n_timepoints <- nrow(ts)

  dummies_season <- matrix(as.factor(rep(1:freq, length.out = length(ts))),
                           ncol = n_units, nrow = n_timepoints)
  dummies_regions <- matrix(as.factor(1:n_units),
                            ncol = n_units, nrow = n_timepoints,
                            byrow = TRUE)

  # get training data:
  subset_training <- 1:t_cond

  dat_training <- as.vector(ts[subset_training, ])
  dummies_season_training <- as.vector(dummies_season[subset_training, ])
  dummies_regions_training <- as.vector(dummies_regions[subset_training, ])

  # fit model:
  fit_nb <- glm.nb(dat_training ~ dummies_season_training + dummies_regions_training)

  # subset test:
  subset_test <- t_cond + (1:max_horizon)
  dummies_season_test <- as.vector(dummies_season[subset_test, ])
  dummies_regions_test <- as.vector(dummies_regions[subset_test, ])

  # predict:
  mu <- predict.glm(fit_nb,
                    newdata = data.frame(dummies_season_training = dummies_season_test,
                                         dummies_regions_training = dummies_regions_test),
                    type = "response")

  size <- fit_nb$theta

  pred_lb50 <- qnbinom(0.25, mu = mu, size = size)
  pred_ub50 <- qnbinom(0.75, mu = mu, size = size)
  pred_lb95 <- qnbinom(0.025, mu = mu, size = size)
  pred_ub95 <- qnbinom(0.975, mu = mu, size = size)

  obs <- as.vector(ts[subset_test, ])

  unit_wise_log_score <- -dnbinom(x = obs, mu = mu, size = size, log = TRUE)
  matr_unit_wise_log_score <- matrix(unit_wise_log_score, ncol = n_units)
  multiv_log_score <- rep(rowMeans(matr_unit_wise_log_score), n_units)

  return(list(prediction_time = t_cond,
              unit = rep(1:n_units, each = max_horizon),
              prediction_horizon = rep(1:max_horizon, n_units),
              pred_mean = mu, pred_var = mu + mu^2/size,
              lb50 = pred_lb50, ub50 = pred_ub50,
              lb95 = pred_lb95, ub95 = pred_ub95,
              obs = obs,
              unit_wise_log_score = unit_wise_log_score,
              multiv_log_score = multiv_log_score))
}
