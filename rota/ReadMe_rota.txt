rotavirus analysis for the Bracher/Held: Endemic-epidemic models with discrete-time serial interval
distributions for infectious disease prediction

Role of the different files nad order in which they should be run:

- choose_order_rota.R: fits models to the whole data set with different values for the order p.
- analyse_fitted_models_rota.R: compare AICs and create descriptive plots for different model variants.
- fit_models_full_rota.R, fit_models_gravity_rota.R: fits models for each subset required for "rolling" forecasts, stores them.
- evaluate_osa_logS_rota.R: evaluate log scores for one-step-ahead-forecasts from the models fitted in fitmodels_rota.R
- obain_forecasts_rota: obtain forecast distributions and evaluate log scores for higher horizons
- obain_forecasts_moments_rota: analytically obtain forecast moments for higher horizons
- naive_forecasting_rota.R: apply naive forecasting method (using glm.nb) to rotavirus data
glarma_forecasting_rota: apply GLARMA models to rotavirus data and obtain forecasts
- analyse_forecasts_rota.R: comparison of mean logS obtained with different models
- plots_forecasts_rota.R: generate some plots visualizing the rotavirus forecasts (as in the Supplementary Material)
