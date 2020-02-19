Dengue analysis for the Bracher/Held: Endemic-epidemic models with discrete-time serial interval
distributions for infectious disease prediction

Role of the different files and order in which they depend on each other:

- choose_order_dengue.R: fits models with different values for the order p, compares their AIC.
- describe_fits_dengue.R: fits models (to training data) with values for p chosen in choose_order_dengue.R, creates figure describing model fit.
- fit_models_dengue.R: fits models for each subset required for "rolling" forecasts, stores them
- obtain_forecasts_dengue.R: obtain forecast distributions from models fitted in fit_models_dengue and evaluate log scores for different horizons
- naive_forecasting_dengue.R: apply naive forecasting method (using glm.nb) to dengue counts
- glarma_forecasts_dengue.R: fit GLARMA
- analyse_logS_dengue.R: comparison of mean logS obtained with different methods

