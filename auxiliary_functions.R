# scp johannes@130.60.71.234:/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis_forecasting/auxiliary_functions.R auxiliary_functions.R

csv_to_sts <- function(file, names, start, end, ...){
  # read data:
  dat <- read.csv2(file,
                   sep = ",", header = TRUE, stringsAsFactors = FALSE)
  dat$year <- dat$week <- NA
  # handle week 53:
  is_week53 <- which(grepl("w53", dat$time))
  dat <- dat[-is_week53, ]
  # handle time variable:
  for(i in 1:nrow(dat)){
    temp <- as.numeric(strsplit(dat$time[i], "w", fixed = TRUE)[[1]])
    dat$year[i] <- temp[1]
    dat$week[i] <- temp[2]
  }
  dat$time <- NULL
  colnames(dat)[1:length(names)] <- names
  # restrict to selected range
  if( (tail(dat$year, 1) < end[1]) ||
      (tail(dat$year, 1) == end[1]) & tail(dat$week, 1) < end[2]){
    stop("Either start or end is outside of the range of the provided data.")
  }

  dat <- subset(dat, year >= start[1])
  dat <- subset(dat, !(year == start[1] & week < start[2]))

  dat <- subset(dat, year <= end[1])
  dat <- subset(dat, !(year == end[1] & week > end[2]))

  dat$week <- NULL
  dat$year <- NULL

  # to sts object:
  stsObj <- new("sts", observed = dat, start = start,...)
  return(stsObj)
}

# obtain some quantities from the fits: endemic components:
get_nu_seas <- function(hhh4Obj, regions, christmas){
  pars <- hhh4Obj$coefficients
  nu <- matrix(ncol = 52, nrow = length(regions), dimnames = list(regions))
  for(region in regions){
    nu[region, ] <- exp(pars[paste0("end.1.", region)] +
                          pars[paste0("end.sin(2 * pi * t/52)")]*sin(2*pi/52*(1:52)) +
                          pars[paste0("end.cos(2 * pi * t/52)")]*cos(2*pi/52*(1:52)) +
                          pars[paste0("end.christmas")] * christmas)
  }

  return(nu)
}

# for comparison: empirical moments:
get_emp_moments <- function(sts, start = 1){
  obs <- sts@observed
  start_week <- sts@start[2]
  freq <- sts@freq
  end_week <- (start_week + nrow(obs) - 1) %% freq

  obs <- rbind(matrix(nrow = start_week - 1, ncol = ncol(obs)),
               obs,
               matrix(nrow = freq - end_week, ncol = ncol(obs)))
  weekwise_means <- weekwise_vars <- matrix(ncol = ncol(obs), nrow = freq)
  for(i in 1:ncol(sts)){
    matr_temp <- matrix(obs[, i], ncol = freq, byrow = TRUE)
    weekwise_means[, i] <- colMeans(matr_temp, na.rm = TRUE)
    weekwise_vars[, i] <- apply(matr_temp, 2, var, na.rm = TRUE)
  }
  colnames(weekwise_means) <- colnames(weekwise_vars) <- colnames(obs)
  reorder <- c(seq(from = start, to = freq),
               seq(to = start - 1, length.out = start - 1))
  weekwise_means <- weekwise_means[reorder, , drop = FALSE]
  weekwise_vars <- weekwise_vars[reorder, , drop = FALSE]
  return(list(mean = weekwise_means, var = weekwise_vars))
}

### THIS VERSION DOES NOT HANDLE THE END OF THE YEAR CORRECTLY!!!
# # compute the stationary moments of calendar week-wise averages
# compute_sm_of_means <- function(fit, n_seasons = 7, start = 1, return_Sigma = FALSE){
#   # obtain dimensions:
#   n_units <- ncol(fit$stsObj@observed)
#   freq <- fit$stsObj@freq
#   # get stationary moments for n_seasons seasons:
#   sm <- stationary_moments(fit, start = start, n_seasons = 3, return_Sigma = TRUE)
#   inds_first <- seq(from = 1, length.out = n_units*freq)
#   inds_last <- seq(to = 3*n_units*freq, length.out = n_units*freq)
#   inds_middle <- seq(to = 2*n_units*freq, length.out = n_units*freq)
#   # we avoid spanning up the matrix for the full n_seasons seasons by just computing
#   # everything for 3 seasons and re-using the middle part of the matrix n_seasons - 2 times.
#
#   new_mu_vector <- sm$mu_vector[inds_first]
#   new_Sigma <- (sm$Sigma[inds_first, inds_first] +
#                   (n_seasons - 2)*sm$Sigma[inds_middle, inds_middle] +
#                   sm$Sigma[inds_last, inds_last])/n_seasons^2
#
#   new_sm <- list()
#
#   if(return_Sigma){
#     new_sm$mu_vector <- new_mu_vector
#     new_sm$Sigma <- new_Sigma
#   }
#
#   new_sm$mu_matrix <- matrix(new_mu_vector, ncol = n_units, byrow = TRUE)
#   new_sm$var_matrix <- matrix(diag(new_Sigma), ncol = n_units, byrow = TRUE)
#
#   rownames(new_sm$mu_matrix) <- rownames(new_sm$var_matrix) <- 1:freq
#   colnames(new_sm$mu_matrix) <- colnames(new_sm$var_matrix) <- colnames(sm$mu_matrix)
#
#   return(new_sm)
# }

#########################################
# compute de-correlated Pearson residuals:
compute_decorr_pearson_residuals <- function(stat_mom_of_means, emp_mom){
  diffs <- (as.vector(t(emp_mom$mean)) - stat_mom_of_means$mu_vector)
  resid_uncorr <- t(chol(solve(stat_mom_of_means$Sigma)))%*%t(t(diffs))
  return(resid_uncorr)
}

# more computing intensive version, only kept for comparison:
compute_sm_of_means_detailed <- function(fit, n_seasons, start = 1, return_Sigma = FALSE){
  # obtain dimensions:
  n_units <- ncol(fit$stsObj@observed)
  freq <- fit$stsObj@freq
  # get stationary moments for n_seasons seasons:
  sm <- stationary_moments(fit, start = start, n_seasons = n_seasons, return_Sigma = TRUE)
  # transformation matrix to move to calendar week-wise means:
  trafo_matr <- matrix(0, nrow = n_units*freq, ncol = n_units*freq*n_seasons)
  for(i in 1:n_seasons){
    diag(trafo_matr[, ((i - 1)*(n_units*freq) + 1:(n_units*freq))]) <- 1/n_seasons
  }
  # aggregate moments:
  new_sm <- list()

  new_mu_vector <- trafo_matr%*%sm$mu_vector
  new_Sigma <- trafo_matr%*%sm$Sigma%*%t(trafo_matr)

  if(return_Sigma){
    new_sm$mu_vector <- new_mu_vector
    new_sm$Sigma <- new_Sigma
  }

  new_sm$mu_matrix <- matrix(new_mu_vector, ncol = n_units, byrow = TRUE)
  new_sm$var_matrix <- matrix(diag(new_Sigma), ncol = n_units, byrow = TRUE)

  rownames(new_sm$mu_matrix) <- rownames(new_sm$var_matrix) <- 1:freq
  colnames(new_sm$mu_matrix) <- colnames(new_sm$var_matrix) <- colnames(sm$mu_matrix)

  return(new_sm)
}

compute_sm_of_means <- function(fit, n_seasons, start = 1, return_Sigma = FALSE){
  sm2 <- stationary_moments(fit, return_Sigma = TRUE, n_seasons = 2, start = start)
  wdt_block <- ncol(sm2$Sigma)/2

  block1 <- sm2$Sigma[1:wdt_block, 1:wdt_block]
  block2 <- (n_seasons - 1)/n_seasons*sm2$Sigma[1:wdt_block, wdt_block + 1:wdt_block]
  block3 <- (n_seasons - 1)/n_seasons*sm2$Sigma[wdt_block + 1:wdt_block, 1:wdt_block]

  sm_of_means <- list()
  sm_of_means$mu_matrix <- sm2$mu_matrix[1:(nrow(sm2$mu_matrix)/2), , drop = FALSE]
  sm_of_means$mu_vector <- sm2$mu_vector[1:wdt_block]
  sm_of_means$var_matrix <- sm2$var_matrix[1:(nrow(sm2$mu_matrix)/2), , drop = FALSE]/n_seasons
  sm_of_means$Sigma <- (block1 + block2 + block3)/n_seasons

  return(sm_of_means)
}

# idea: fill only block diagonal and
compute_sm_of_means2 <- function(fit, n_seasons, start = 1, return_Sigma = FALSE){
  # obtain dimensions:
  n_units <- ncol(fit$stsObj@observed)
  freq <- fit$stsObj@freq
  # get stationary moments for n_seasons seasons:
  sm2 <- stationary_moments(fit, start = start, n_seasons = 2, return_Sigma = TRUE)

  lgt_n_seas <- n_units*freq*n_seasons

  long_mu_vector <- rep(sm2$mu_vector, length.out = lgt_n_seas)
  large_Sigma <- matrix(0, nrow = lgt_n_seas, ncol = lgt_n_seas)
  for(i in 1:(n_seasons - 1)){
    large_Sigma[(i - 1)*n_units*freq + (1:(2*n_units*freq)),
          (i - 1)*n_units*freq + (1:(2*n_units*freq))] <-
      sm2$Sigma
  }

  # transformation matrix to move to calendar week-wise means:
  trafo_matr <- matrix(0, nrow = n_units*freq, ncol = n_units*freq*n_seasons)
  for(i in 1:n_seasons){
    diag(trafo_matr[, ((i - 1)*(n_units*freq) + 1:(n_units*freq))]) <- 1/n_seasons
  }
  # aggregate moments:
  new_sm <- list()

  new_mu_vector <- trafo_matr%*%long_mu_vector
  new_Sigma <- trafo_matr%*%large_Sigma%*%t(trafo_matr)

  if(return_Sigma){
    new_sm$mu_vector <- new_mu_vector
    new_sm$Sigma <- new_Sigma
  }

  new_sm$mu_matrix <- matrix(new_mu_vector, ncol = n_units, byrow = TRUE)
  new_sm$var_matrix <- matrix(diag(new_Sigma), ncol = n_units, byrow = TRUE)

  rownames(new_sm$mu_matrix) <- rownames(new_sm$var_matrix) <- 1:freq
  colnames(new_sm$mu_matrix) <- colnames(new_sm$var_matrix) <- colnames(sm$mu_matrix)

  return(new_sm)
}

##################
# p-values for region-wise versions:
get_pval <- function(stat_mom_of_means, emp_mom, correction_df = 0){
  diffs <- as.vector(t(emp_mom$mean)) - stat_mom_of_means$mu_vector
  test_stat <- (t(diffs)%*%solve(stat_mom_of_means$Sigma)%*%diffs)[1, 1]
  df <- length(diffs) - correction_df
  pval <- pchisq(test_stat, df, lower.tail = FALSE)
  return(list(test_stat = test_stat, df = df, pval = pval))
}

# p-values for region-wise versions with 2/3 correction:
get_pval2 <- function(stat_mom_of_means, emp_mom, correction_df = 0){
  h <- function(x) x^(2/3)
  diffs <- h(as.vector(t(emp_mom$mean))) - h(stat_mom_of_means$mu_vector)
  cov_diffs <- diag(2/3*as.vector(stat_mom_of_means$mu_vector)^(-1/3)) %*% stat_mom_of_means$Sigma %*% diag(2/3*as.vector(stat_mom_of_means$mu_vector)^(-1/3))
  test_stat <- t(diffs)%*%solve(cov_diffs)%*%diffs
  df <- length(diffs) - correction_df
  pval <- pchisq(test_stat, df, lower.tail = FALSE)
  return(list(test_stat = test_stat, df = df, pval = pval))
}

### one-step-ahead forecasts also updating the lag decay parameter
# arguments: fit: an hhh4_lag object from which to generate one-step-ahead-forecasts
#             tp: range of timepoints for which to generate forecasts, i.e. tp[1] + 1, ..., tp[2] + 1.
osa_refitting_par_lag <- function(fit, tp){
  # extract info from arguments:
  control <- fit$control
  sts <- fit$stsObj
  start_subset <- control$subset[1]
  tp_vect <- tp[1]:tp[2]

  # initialize osa list to store results:
  templ <- matrix(ncol = ncol(sts@observed), nrow = length(tp_vect),
                  dimnames = list(tp_vect + 1, colnames(sts@observed)))
  osa <- list(pred = templ, observed = templ, psi = templ,
              par_lag = matrix(NA, nrow = length(tp_vect), ncol = length(fit$par_lag)))
  if(control$family == "Poisson") osa$psi <- NULL

  mod_temp <- fit

  # run through time points for which to obtain one-step-ahead-predictions
  for(i in seq_along(tp_vect)){
    # adapt subset:
    control$subset <- start_subset:(tp_vect[i])
    # adopt grid for lag decay parameter
    # new_grid <- get_grid_around_previous_par_lag(mod_temp)
    # update model:
    mod_temp <- profile_par_lag(sts, control = control, start_par_lag = mod_temp$par_lag)
    # get one-step-ahead prediction
    osa_temp <- suppressMessages(oneStepAhead_hhh4lag(mod_temp, tp = c(tp_vect[i], tp_vect[i]), type = "final"))
    # store result:
    osa$pred[i, ] <- osa_temp$pred
    osa$observed[i, ] <- osa_temp$observed
    if(control$family != "Poisson"){
      osa$psi[i, ] <- osa_temp$psi
    }
    osa$par_lag[i, ] <- mod_temp$par_lag
    print(i)
  }

  return(osa)
}

# Functions for decomposing the fitted values (used in plotting)
decompose_epidemic_component <- function(fit){
  # extract info:
  sts <- fit$stsObj
  max_lag <- if(class(fit)[1] == "hhh4lag") fit$max_lag else 1
  subset <- fit$control$subset
  n_units <- ncol(sts@observed)
  param <- hhh4addon:::lambda_tilde_complex_neighbourhood(fit, periodic = FALSE,
                                                          subset = 1:max(subset))

  # initialize:
  contributions <- array(dim = c(max(subset),
                                 n_units,
                                 n_units,
                                 max_lag),
                         dimnames = list(1:max(subset),
                                         paste0("from.", colnames(sts@observed)),
                                         paste0("to.", colnames(sts@observed)),
                                         paste0("lag", 1:max_lag)))
  # fill:
  for(t in subset){
    phi_temp <- param$lambda[,,t]
    obs_preceding <- sts@observed[t - max_lag:1, , drop = FALSE]
    for(lag in 1:max_lag){
      inds <- seq(to = n_units*(max_lag - lag + 1), length.out = n_units)
      phi_this_lag <- phi_temp[, inds]
      contributions[t, , , lag] <- t(phi_this_lag)*matrix(obs_preceding[max_lag - lag + 1, ], ncol = n_units, nrow = n_units)
    }
  }
  return(contributions)
}

decompose_coarse <- function(fit){
  sts <- fit$stsObj
  max_lag <- if(class(fit)[1] == "hhh4lag") fit$max_lag else 1
  subset <- fit$control$subset
  n_units <- ncol(sts@observed)
  param <- hhh4addon:::lambda_tilde_complex_neighbourhood(fit, periodic = FALSE,
                                                          subset = 1:max(subset))
  decomposition_epidemic <- decompose_epidemic_component(fit)
  contributions_coarse <- array(dim = c(max(subset),
                                        5,
                                        n_units),
                                dimnames = list(1:max(subset),
                                                c("endemic",
                                                  "epidemic.self.lag1", "epidemic.self.higher_lags",
                                                  "epidemic.other.lag1", "epidemic.other.higher_lags"),
                                                paste0("to.", colnames(sts@observed))))
  for(t in subset){
    contributions_coarse[t, "endemic", ] <- param$nu[t, ]
    contributions_coarse[t, "epidemic.self.lag1", ] <- diag(decomposition_epidemic[t, , , 1])
    contributions_coarse[t, "epidemic.other.lag1", ] <-
      colSums(decomposition_epidemic[t, , , 1]) - contributions_coarse[t, "epidemic.self.lag1", ]
    if(max_lag > 1){
      contributions_higher_lags_temp <- apply(decomposition_epidemic[t, , , 2:max_lag], 1:2, sum)
      contributions_coarse[t, "epidemic.self.higher_lags", ] <- diag(contributions_higher_lags_temp)
      contributions_coarse[t, "epidemic.other.higher_lags", ] <-
        colSums(contributions_higher_lags_temp) - contributions_coarse[t, "epidemic.self.higher_lags", ]
    }else{
      contributions_coarse[t, "epidemic.self.higher_lags", ] <- contributions_coarse[t, "epidemic.other.higher_lags", ] <- 0
    }
  }
  return(contributions_coarse)
}

# fancy plot of fitted values:
plot_fit_strat <- function(fit, unit = 1, col = c("lightgrey", brewer.pal(n = 4, name = 'RdBu')[c(2, 1, 3, 4)]),
                           ylim = NULL, cex.points = 0.7, draw_yax = TRUE,...){
  polygon_from_curve <- function(tp, coords, ...){
    tp <- tp[!is.na(coords)]; tp <- c(tp[1], tp, tail(tp, 1))
    coords <- coords[!is.na(coords)]; coords <- c(0, coords, 0)
    polygon(tp, coords, ...)
  }
  obs <- fit$stsObj@observed[, unit]
  timepoints_calendar <- seq(from = fit$stsObj@start[1] + fit$stsObj@start[2]/fit$stsObj@freq,
                             length.out = nrow(fit$stsObj@observed), by  =1/fit$stsObj@freq)
  fitted_vals <- fitted.values(fit)[, unit]
  decomp <- decompose_coarse(fit)
  cumsums <- decomp
  if(is.null(ylim)) ylim <- c(0, max(c(cumsums[, , unit], obs), na.rm = TRUE))
  for(i in 2:5) cumsums[, i, ] <- cumsums[, i, ] + cumsums[, i - 1, ]
  plot(NULL, xlim = range(timepoints_calendar), ylim = ylim,
       xlab = NA, ylab = "no. of reported cases", axes = FALSE, ...)
  axis(1, labels = NA)
  axis(1, lty = 0, at = 2011:2017 + 0.5, labels = 2011:2017)
  if(draw_yax) axis(2)
  box()
  polygon_from_curve(timepoints_calendar, cumsums[, "epidemic.other.higher_lags", unit], col = col[5], border = col[5])
  polygon_from_curve(timepoints_calendar, cumsums[, "epidemic.other.lag1", unit], col = col[4], border = col[4])
  polygon_from_curve(timepoints_calendar, cumsums[, "epidemic.self.higher_lags", unit], col = col[3], border = col[3])
  polygon_from_curve(timepoints_calendar, cumsums[, "epidemic.self.lag1", unit], col = col[2], border = col[2])
  polygon_from_curve(timepoints_calendar, cumsums[, "endemic", unit], col = col[1], border = col[1])
  lines(timepoints_calendar[fit$control$subset], fitted_vals)
  points(timepoints_calendar, obs, pch = 16, cex = cex.points)
}

# adding the (within-sample) predictive interval:
add_pi <- function(fit, unit, upper, lower,...){
  timepoints_calendar <- seq(from = fit$stsObj@start[1] + fit$stsObj@start[2]/fit$stsObj@freq,
                             length.out = nrow(fit$stsObj@observed), by  =1/fit$stsObj@freq)
  subset <- fit$control$subset
  fitted_vals <- fitted.values(fit)[, unit]
  size <- exp(fit$coefficients["-log(overdisp)"])
  upper <- qnbinom(upper, mu = fitted_vals, size = size)
  lower <- qnbinom(lower, mu = fitted_vals, size = size)

  lines(timepoints_calendar[subset], upper, lty = "dotted")
  lines(timepoints_calendar[subset], lower, lty = "dotted")
}

# Functions for residual analyses:
# function to compute within-regions ACFs:
my_acf <- function(pearson_resids){
  acf_matr <- matrix(ncol = ncol(pearson_resids), nrow = 52,
                     dimnames = list(NULL, colnames(pearson_resids)))
  for(i in 1:ncol(pearson_resids)){
    acf_matr[, i] <- acf(pearson_resids[, i], lag.max = 52, plot = FALSE)$acf[-1,,1]
  }
  return(acf_matr)
}

# a customized plotting function:
myplot_acf <- function(macf, unit = 1, nobs = 358, shift_x = 0, lwd = 2, ylim = c(-0.2, 0.2), add = FALSE, ...){
  if(!add){
    plot(c(1:5, 7) + shift_x,
         macf[c(1:5, 52), unit],
         xlim = c(0.5, 7.5), ylim = ylim, type = "h",
         axes = FALSE, ylab = "residual ACF",
         lwd = lwd, xlab = "lag", ...)
    axis(1, at = c(1:5, 7), labels = c(1:5, 52))
    axis(2)
    box()
    abline(h = 0)
    x <- 6; d <- 0.03; y_ax <- 1.1*min(ylim)
    rect(x - d, y_ax -d, x + d, y_ax + d, col = "white", border = NA, xpd = TRUE)
    lines(c(x - 2*d, x), y_ax + c(d/2, -d/2), xpd = TRUE)
    lines(c(x, x + 2*d), y_ax + c(d/2, -d/2), xpd = TRUE)
    abline(h = c(-1.96, 1.96)/sqrt(nobs), lty = "dotted")
  }else{
    lines(c(1:5, 7) + shift_x,
          macf[c(1:5, 52), unit], type = "h", lwd = lwd, ...)
  }
}

# Plotting functions for periodically stationary moments:
# Means:
plot_stat_means <- function(stat_mom, emp_mom, envelopes = NULL, disease, unit = 1,...){
  plot(stat_mom[[disease]]$full$geom$mu_matrix[, unit], type = "l",
       xlab = "calendar week", ylab = "mean", col = cols_model_versions[1],
       axes = FALSE,...)
  axis(1, at = c(1, 13, 26, 39, 52)); axis(2); box()
  lines(stat_mom[[disease]]$simple_seas$geom$mu_matrix[, unit], col = cols_model_versions[2])
  # lines(stat_mom[[disease]]$no_cross$geom$mu_matrix[, unit], col = cols_model_versions[3])
  # lines(stat_mom[[disease]]$neither$geom$mu_matrix[, unit], col = cols_model_versions[4])
  # lines(stat_mom[[disease]]$full$ar1$mu_matrix[, unit], col = cols_model_versions[5])
  # lines(stat_mom_unstrat[[disease]]$geom$mu_matrix, col = cols_model_versions[6], lty = "dashed")

  if(!is.null(envelopes)){
    lines(envelopes[[disease]]$full$geom$quantiles_means$`p=0.97`[, unit],
          col = cols_model_versions[1], lty = 2)
    lines(envelopes[[disease]]$full$geom$quantiles_means$`p=0.03`[, unit],
          col = cols_model_versions[1], lty = 2)

    lines(envelopes[[disease]]$simple_seas$geom$quantiles_means$`p=0.97`[, unit],
          col = cols_model_versions[2], lty = 2)
    lines(envelopes[[disease]]$simple_seas$geom$quantiles_means$`p=0.03`[, unit],
          col = cols_model_versions[2], lty = 2)
  }

  points(emp_mom[[disease]]$mean[, unit], pch = 1, cex = 0.6)
}

# Standard deviations:
plot_stat_sds <- function(stat_mom, emp_mom, envelopes = NULL, disease, unit = 1,...){
  plot(sqrt(stat_mom[[disease]]$full$geom$var_matrix[, unit]), type = "l",
       xlab = "calendar week", ylab = "standard deviation", col = cols_model_versions[1],
       axes = FALSE, ...)
  axis(1, at = c(1, 13, 26, 39, 52)); axis(2); box()
  lines(sqrt(stat_mom[[disease]]$simple_seas$geom$var_matrix[, unit]), col = cols_model_versions[2])
  # lines(sqrt(stat_mom[[disease]]$no_cross$geom$var_matrix[, unit]), col = cols_model_versions[3])
  # lines(sqrt(stat_mom[[disease]]$neither$geom$var_matrix[, unit]), col = cols_model_versions[4])
  # lines(sqrt(stat_mom[[disease]]$full$ar1$var_matrix[, unit]), col = cols_model_versions[5])
  # lines(sqrt(stat_mom_unstrat[[disease]]$geom$var_matrix), col = cols_model_versions[6], lty = "dashed")

  if(!is.null(envelopes)){
    lines(sqrt(envelopes[[disease]]$full$geom$quantiles_vars$`p=0.97`[, unit]),
          col = cols_model_versions[1], lty = 2)
    lines(sqrt(envelopes[[disease]]$full$geom$quantiles_vars$`p=0.03`[, unit]),
          col = cols_model_versions[1], lty = 2)

    lines(sqrt(envelopes[[disease]]$simple_seas$geom$quantiles_vars$`p=0.97`[, unit]),
          col = cols_model_versions[2], lty = 2)
    lines(sqrt(envelopes[[disease]]$simple_seas$geom$quantiles_vars$`p=0.03`[, unit]),
          col = cols_model_versions[2], lty = 2)
  }

  points(sqrt(emp_mom[[disease]]$var[, unit]), pch = 1, cex = 0.6)
}

# means, with "confidence" bands:
plot_sm_bands <- function(sm_of_means, emp_mom,
                          disease, model_version, lag_structure,
                          unit, add = FALSE, col = "black", ylim = c(0, 20), ...){
  means <- sm_of_means[[disease]][[model_version]][[lag_structure]]$mu_matrix[, unit]
  sds <- sqrt(sm_of_means[[disease]][[model_version]][[lag_structure]]$var_matrix[, unit])
  upper <- means + 1.96*sds
  lower <- means - 1.96*sds

  if(!add){
    par(mar = c(1, 4.2, 4, 1))
    plot(NULL, ylim = ylim, xlim = c(1, 52),
         cex = 0.7, xlab = "", ylab = expression(hat(mu)[it]),
         axes = FALSE, ...)
    axis(2)
    box()
  }
  lines(sm_of_means[[disease]][[model_version]][[lag_structure]]$mu_matrix[, unit],
        col = col)
  lines(upper, lty = "dotted", col = col)
  lines(lower, lty = "dotted", col = col)
  points(emp_mom[[disease]]$mean[, unit], cex = 0.6)
}

# residuals relative to periodically stationary moments:
plot_stat_resids <- function(sm_of_means, emp_mom,
                             disease,
                             model_version1, lag_structure1,
                             model_version2, lag_structure2,
                             unit, add = FALSE, col = "black", ylim = c(-5, 5), ...){

  my_lines <- function(v1, v2, col1 = cols_model_versions[1], col2 = cols_model_versions[2]){
    inds1larger <- which(abs(v1) > abs(v2))
    inds2larger <- which(abs(v2) > abs(v1))
    lines(inds1larger, v1[inds1larger], type = "h", col = col1)
    lines(inds2larger, v2[inds2larger], type = "h", col = col2, lwd = 1)
    lines(inds2larger, v1[inds2larger], type = "h", col = col1)
    lines(inds1larger, v2[inds1larger], type = "h", col = col2, lwd = 1)
  }

  emp_means <- emp_mom[[disease]]$mean[, unit]

  means1 <- sm_of_means[[disease]][[model_version1]][[lag_structure1]]$mu_matrix[, unit]
  sds1 <- sqrt(sm_of_means[[disease]][[model_version1]][[lag_structure1]]$var_matrix[, unit])
  pr1 <- (emp_means - means1)/sds1

  means2 <- sm_of_means[[disease]][[model_version2]][[lag_structure2]]$mu_matrix[, unit]
  sds2 <- sqrt(sm_of_means[[disease]][[model_version2]][[lag_structure2]]$var_matrix[, unit])
  pr2 <- (emp_means - means2)/sds2

  plot(NULL, xlim = c(1, 52), axes = FALSE,
       type = "h", ylim = ylim, col = cols_model_versions[2],
       xlab = "calendar week", ylab = "Pearson residual",...)
  my_lines(pr1, pr2)
  axis(1, at = c(1, 13, 26, 39, 52)); axis(2); box()
  abline(h = 0, col = "black")
  abline(h = c(-1.96, 1.96), col = "black", lty = "dotted")
}

plot_stat_resids2 <- function(stat_resids, disease,
                              model_version1, lag_structure1,
                              model_version2, lag_structure2,
                              unit = 1, ylim = c(-5, 5), ...){

  my_lines <- function(v1, v2, col1 = cols_model_versions[1], col2 = cols_model_versions[2]){
    inds1larger <- which(abs(v1) > abs(v2))
    inds2larger <- which(abs(v2) > abs(v1))
    lines(inds1larger, v1[inds1larger], type = "h", col = col1)
    lines(inds2larger, v2[inds2larger], type = "h", col = col2, lwd = 1)
    lines(inds2larger, v1[inds2larger], type = "h", col = col1)
    lines(inds1larger, v2[inds1larger], type = "h", col = col2, lwd = 1)
  }

  pr1 <- stat_resids[[disease]][[model_version1]][[lag_structure1]][, unit]
  pr2 <- stat_resids[[disease]][[model_version2]][[lag_structure2]][, unit]

  plot(NULL, xlim = c(1, 52), axes = FALSE,
       type = "h", ylim = ylim, col = cols_model_versions[2],
       xlab = "calendar week", ylab = "Pearson residual",...)
  my_lines(pr1, pr2)
  axis(1, at = c(1, 13, 26, 39, 52)); axis(2); box()
  abline(h = 0, col = "black")
  abline(h = c(-1.96, 1.96), col = "black", lty = "dotted")
}


# function to simulate p-values
sim_teststats <- function(fit, n_sim, n_seasons, seed = 123){
  set.seed(seed)
  # sm_of_means_sim <- compute_sm_of_means(fit, n_seasons = n_seasons, return_Sigma = TRUE) # only for testing
  teststats_sim <- teststats_trafo_sim <- numeric(n_sim)
  for(i in 1:n_sim){
    # generate data and fit model to simulated data::
    if(is.null(fit$control$funct_lag)){
      max_lag <- 1
      dat_sim0 <- simulate(fit, y.start = fit$stsObj@observed[1, ]) # run once to go through "burn-in" period (independence of starting values)
      dat_sim <- simulate(fit, y.start = dat_sim0@observed[52 - 1, ]) # feed values from the first run into a second as starting values
      fit_sim <- hhh4(dat_sim, fit$control)
    }else{
      max_lag <- fit$control$max_lag
      dat_sim0 <- hhh4addon:::simulate.hhh4lag(fit, y.start = fit$stsObj@observed[1:max_lag, , drop = FALSE])
      dat_sim <- hhh4addon:::simulate.hhh4lag(fit, y.start = dat_sim0@observed[n_seasons*52 - (max_lag - 1):0, , drop = FALSE])
      fit_sim <- hhh4addon::profile_par_lag(dat_sim, fit$control)
    }

    # compute stationary moments of week-wise averages:
    sm_of_means_sim <- compute_sm_of_means(fit_sim, n_seasons = n_seasons, return_Sigma = TRUE)
    # compute empirical moments:
    subset_emp_mom <- seq(from = min(fit$control$subset) - max_lag, to = max(fit$control$subset))
    emp_mom_sim <- get_emp_moments(dat_sim[subset_emp_mom, ])
    # get test statistic:
    teststats_sim[i] <- get_pval(stat_mom_of_means = sm_of_means_sim, emp_mom = emp_mom_sim)$test_stat
    teststats_trafo_sim[i] <- get_pval2(stat_mom_of_means = sm_of_means_sim, emp_mom = emp_mom_sim)$test_stat
    print(i)
  }
  return(list(teststats_sim = teststats_sim, teststats_trafo_sim = teststats_trafo_sim))
}


# same simulation procedure, but to retain envelopes:
sim_envelopes_sm_of_means <- function(fit, n_seasons, n_sim = 100,
                                      probs = c(0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99),
                                      start = 1, return_samples = FALSE, seed = 123){
  n_units <- ncol(fit$stsObj@observed)
  freq <- fit$stsObj@freq
  set.seed(seed)

  samples_weekwise_means <- samples_weekwise_vars <-
      array(dim = c(freq, n_units, n_sim))

  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  for(i in 1:n_sim){
    # generate data and fit model to simulated data::
    if(is.null(fit$control$funct_lag)){
      max_lag <- 1
      dat_sim0 <- simulate(fit, y.start = fit$stsObj@observed[1, ]) # run once to go through "burn-in" period (independence of starting values)
      dat_sim <- simulate(fit, y.start = dat_sim0@observed[n_seasons*freq - 1, ]) # feed values from the first run into a second as starting values
    }else{
      max_lag <- fit$control$max_lag
      dat_sim0 <- hhh4addon:::simulate.hhh4lag(fit, y.start = fit$stsObj@observed[1:max_lag, , drop = FALSE])
      dat_sim <- hhh4addon:::simulate.hhh4lag(fit, y.start = dat_sim0@observed[n_seasons*freq - (max_lag - 1):0, , drop = FALSE])
    }

    # compute empirical moments:
    subset_emp_mom <- seq(from = min(fit$control$subset) - max_lag, to = max(fit$control$subset))
    emp_mom_sim <- get_emp_moments(dat_sim[subset_emp_mom, ], start = start)

    samples_weekwise_means[,,i] <- emp_mom_sim$mean
    samples_weekwise_vars[,,i] <- emp_mom_sim$var

    setTxtProgressBar(pb, i)
  }

  quantiles_weekwise_means <- list()
  quantiles_weekwise_vars <- list()

  for(j in seq_along(probs)){
    quantiles_weekwise_means[[paste0("p=", probs[j])]] <-
      apply(samples_weekwise_means, 1:2, quantile, probs = probs[j])
    quantiles_weekwise_vars[[paste0("p=", probs[j])]] <-
      apply(samples_weekwise_vars, 1:2, quantile, probs = probs[j])
  }

  return(list(samples_mean = if(return_samples) samples_weekwise_means else NULL,
              samples_var = if(return_samples) samples_weekwise_vars else NULL,
              quantiles_means = quantiles_weekwise_means,
              quantiles_vars = quantiles_weekwise_vars))
}


# function to evaluate logS via Rao-Blackwellization
# and return district-wise summaries of predictive distributions
forecasting_nStepAhead <- function(fit, stsObj, tp_cond, horizon, n_sim){

  is_hhh4_lag <- (class(fit)[1] == "hhh4lag") # check if hhh4 or hhh4lag object
  n_units <- ncol(stsObj@observed)

  fit$stsObj <- stsObj # replace stsObj in fit by original stsObj (as subsets steered via NA)
  max_lag <- ifelse(is_hhh4_lag, fit$max_lag, 1)
  # wrapper function to handle apply commands
  dnb <- function(mu, size, x){
    dnbinom(x = x, mu = mu, size = size)
  }

  # for horizon h = 1 no simulation is necessary:
  if(horizon <= 1){
    forecast <-
      if(is_hhh4_lag){
        suppressMessages(
          oneStepAhead_hhh4lag(fit, tp = rep(tp_cond + horizon - 1, 2), type = "final")
        )
      }else{
        oneStepAhead(fit, tp = rep(tp_cond + horizon - 1, 2), type = "final")
      }
    log_score <- scores(forecast, which = "logs")
  }else{ # for horizons >= 2 we need to simulate:
    tp_to_simulate <- tp_cond + 1:(horizon)
    # need to do one more time point due to bug in surveillance

    # generate sample paths:
    sims <- simulate(fit, subset = tp_to_simulate,
                     y.start = fit$stsObj@observed[tp_cond - (max_lag:1) + 1, , drop = FALSE],
                     simplify = TRUE, nsim = n_sim)
    sims[length(tp_to_simulate),,] <- stsObj@observed[tp_cond + horizon, ] # put true observation back

    # obtain log scores for each sample path:
    sim_log_scores <- sim_mu <- matrix(NA, ncol = n_units, nrow = n_sim)
    for(i in 1:n_sim){
      fit_temp <- fit
      fit_temp$stsObj@observed[tp_to_simulate, ] <- sims[,,i] # plug simulated path into fit object

      # do one-step ahead forecast given the simulated path:
      forecast_temp <-
        if(is_hhh4_lag){
          suppressMessages(
            oneStepAhead_hhh4lag(fit_temp, tp = rep(tp_cond + horizon - 1, 2), type = "final")
          )
        }else{
          oneStepAhead(fit_temp, tp = rep(tp_cond + horizon - 1, 2), type = "final")
        }

      # store the predictive means and size (in the NB distribution) under the respective
      # simulated path:
      sim_mu[i, ] <- forecast_temp$pred
      if(i == 1) size <- rep(exp(forecast_temp$psi), length.out = n_units) # stays the same in all iterations

      # store the log score obtained under the respective simulated path:
      scores_temp <- scores(forecast_temp, individual = TRUE)
      if(ncol(stsObj@observed) == 1){
        sim_log_scores[i, ] <- scores_temp["logs"]
      }else{
        sim_log_scores[i, ] <- scores_temp[, "logs"]
      }
    }
    # average over log scores obtained with different sample paths:
    log_score <- -log(mean(exp(-rowSums(sim_log_scores))))/ncol(sim_log_scores)
  }

  # extract characteristics of predictive distributions and the scores:
  templ_vect <- numeric(n_units); names(templ_vect) <- colnames(stsObj@observed)
  unit_wise_log_score <- pred_mean <- pred_var <- pred_lb50 <- pred_ub50 <-
    pred_lb95 <- pred_ub95 <- unit_wise_pit_l <- unit_wise_pit_u <- templ_vect

  for(unit in 1:n_units){
    # get observed value:
    obs_unit_temp <- stsObj@observed[tp_cond + horizon, unit]

    if(horizon == 1){ # for horizon 1 can extract directly from return of oneStepAhead
      mu_unit <- forecast$pred[unit]
      size_unit <- rep(exp(forecast$psi), length.out = n_units)[unit] # catch case of NegBin1
      support_temp <- 0:max(qnbinom(p = 0.99,
                                    mu = mu_unit,
                                    size = size_unit),
                            obs_unit_temp + 2, na.rm = TRUE)
      pred_densities_temp <- dnbinom(support_temp, size = size_unit, mu = mu_unit)

      pred_mean[unit] <- mu_unit
      pred_var[unit] <- mu_unit + 1/size_unit*mu_unit^2
    }else{ # otherwise need to average over samples:
      support_temp <- 0:max(qnbinom(p = 0.99,
                                    mu = max(sim_mu[, unit]),
                                    size = size[unit]),
                            obs_unit_temp + 2, na.rm = TRUE)

      pred_densities_temp <- rowMeans(sapply(sim_mu[, unit], dnb,
                                             x = support_temp,
                                             size = size[unit]))
      pred_mean[unit] <- sum(support_temp*pred_densities_temp)
      pred_var[unit] <- sum(support_temp^2*pred_densities_temp) - pred_mean[unit]^2
    }

    # compute limits of prediction intervals:
    pred_cumul_distr_temp <- cumsum(pred_densities_temp)
    pred_lb50[unit] <- min(which(pred_cumul_distr_temp >= 0.25)) - 1 # -1 bc support includes 0
    pred_ub50[unit] <- max(which(pred_cumul_distr_temp <= 0.75)) # no -1  as we want to be slightly conservative
    pred_lb95[unit] <- min(which(pred_cumul_distr_temp >= 0.025)) - 1
    pred_ub95[unit] <- max(which(pred_cumul_distr_temp <= 0.975))

    # compute unit-wise log scores:
    unit_wise_log_score[unit] <- -log(pred_densities_temp[
      obs_unit_temp + 1
      ])

    # and unit-wise PIT value:
    if(!is.na(obs_unit_temp)){
      unit_wise_pit_l[unit] <- if(obs_unit_temp == 0){
        0
      }else{
        pred_cumul_distr_temp[obs_unit_temp] # + 1 bc support includes 0
      }
      unit_wise_pit_u[unit] <- pred_cumul_distr_temp[obs_unit_temp + 1]
    }else{
      unit_wise_pit_l[unit] <- unit_wise_pit_u[unit] <- NA
    }
  }

  # return object:
  return(list(pred_mean = pred_mean, pred_var = pred_var,
              lb50 = pred_lb50, ub50 = pred_ub50,
              lb95 = pred_lb95, ub95 = pred_ub95,
              obs = stsObj@observed[tp_cond + horizon, ],
              unit_wise_log_score = unit_wise_log_score,
              unit_wise_pit_l = unit_wise_pit_l,
              unit_wise_pit_u = unit_wise_pit_u,
              multiv_log_score = log_score))
}

# function to obtain naive seasonal forecasts, similar to Farrington/Noufaily
naive_forecast <- function(ts, t, horizon = 1, half_window_width = 0, freq = 52){
  n_units <- ncol(ts)
  seas <- floor(t/freq) # number of past seasons available
  # determine indices of past observations within
  # a certain window around the calendar week to be forecasted
  subset_training0 <- (t + horizon) - (1:seas)*freq
  subset_training <- subset_training0
  if(half_window_width > 0){
    for(i in 1:half_window_width){
      subset_training <- c(subset_training0 - i,
                           subset_training,
                           subset_training0 + i
      )
    }
  }
  subset_training <- subset_training[subset_training > 0]
  dat_training <- ts[subset_training, , drop = FALSE]

  mu <- size <- numeric(n_units)

  for(unit in 1:n_units){
    dat_unit <- dat_training[, unit]
    # fit Poisson model if approximately equidispersed, NB otherwise
    if(var(dat_unit) > 1.1*mean(dat_unit)){
      mod <- glm.nb(dat_unit ~ 1)
      mu[unit] <- exp(mod$coefficients)
      size[unit] <- mod$theta
    }else{
      mod <- glm(dat_unit ~ 1, family = poisson)
      mu[unit] <- exp(mod$coefficients)
      size[unit] <- 100
    }
  }

  pred_lb50 <- qnbinom(0.25, mu = mu, size = size)
  pred_ub50 <- qnbinom(0.75, mu = mu, size = size)
  pred_lb95 <- qnbinom(0.025, mu = mu, size = size)
  pred_ub95 <- qnbinom(0.975, mu = mu, size = size)

  obs <- ts[t + horizon, ]

  unit_wise_log_score <- -dnbinom(x = obs, mu = mu, size = size, log = TRUE)
  multiv_log_score <- rep(mean(unit_wise_log_score), n_units)

  return(list(pred_mean = mu, pred_var = mu + mu^2/size,
              lb50 = pred_lb50, ub50 = pred_ub50,
              lb95 = pred_lb95, ub95 = pred_ub95,
              obs = obs,
              unit_wise_log_score = unit_wise_log_score,
              multiv_log_score = multiv_log_score,
              half_window_width = half_window_width))
}

# Helper function to modify alpha value of a given colour
modify_alpha <- function(col, alpha) {
  col_rgb <- col2rgb(col, alpha = TRUE)
  col_rgb[4,] <- alpha*255
  col_rgb <- col_rgb/255.0
  return(rgb(col_rgb[1,], col_rgb[2,], col_rgb[3,], col_rgb[4,]))
}

# visualize forecasts:
# plotting function:
shade_forecast <- function(timepoints, lb, ub, col = rgb(0, 0, 0, alpha = 0.5), ...){
  polygon(x = c(timepoints, rev(timepoints)),
          y = c(lb, rev(ub)),
          col = col, border = col, ...)
}

# compute coverage of 50% and 95% prediction intervals:
eval_coverage <- function(subs){
  c(
    coverage50 = mean(subs$obs >= subs$lb50 &
                        subs$obs <= subs$ub50),
    coverage95 = mean(subs$obs >= subs$lb95 &
                        subs$obs <= subs$ub95)
  )
}

# two functions for plotting randomized PIT histograms
plot_pit_hist <- function(pit_l, pit_u, ...){
  pit_vals <- pit_l + runif(20*length(pit_l), 0, 1)*(pit_u - pit_l)
  hist(pit_vals, freq = FALSE,  xlab = "PIT", breaks = 0:10/10, ...)
}

add_pit_hist <- function(pit_l, pit_u, shift = 0.05, ...){
  pit_vals <- pit_l + runif(length(pit_l), 0, 1)*(pit_u - pit_l)
  hi <- hist(pit_vals, plot = FALSE)
  points(hi$breaks[-length(hi$breaks)] + shift, hi$density, type = "h",
         lwd = 4,...)
}