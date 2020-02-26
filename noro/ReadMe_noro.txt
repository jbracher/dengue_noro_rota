Norovirus analysis for the Bracher/Held: Endemic-epidemic models with discrete-time serial interval
distributions for infectious disease prediction

Role of the different files nad order in which they should be run:

- choose_order_noro.R: fits models to the whole data set with different values for the order p.
- analyse_fitted_models_noro.R: compare AICs and create descriptive plots for different model variants.
- fit_models_full_noro.R, fit_models_gravity_noro.R: fits models for each subset required for "rolling" forecasts, stores them.
- evaluate_osa_logS_noro.R: evaluate log scores for one-step-ahead-forecasts from the models fitted in fitmodels_noro.R
- obain_forecasts_noro: obtain forecast distributions and evaluate log scores for higher horizons
- obain_forecasts_moments_noro: analytically obtain forecast moments for higher horizons
- naive_forecasting_noro.R: apply naive forecasting method (using glm.nb) to norovirus data
glarma_forecasting_noro: apply GLARMA models to norovirus data and obtain forecasts
- analyse_forecasts_noro.R: comparison of mean logS obtained with different models
- plots_forecasts_noro.R: generate some plots visualizing the norovirus forecasts (as in the Supplementary Material)

- AIC/ contains summaries of the AICs obtained by the different models
- forecasts/ contains summaries of the forecasts from the different models (predictive moments and intervals, PIT values etc)
- logS/ contains summaries of the obtained logarithmic scores
- model_fits/ contains all fitted models from which forecasts were obtained (as rda files)
