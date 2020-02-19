# Handling some libraries and names/colour schemes needed in almost all analyses

library(surveillance)
library(hhh4addon)
library(RColorBrewer)

# For dengue:
data("dengueSJ")
names_models_dengue <- c("ar1", "pois", "lin", "geom",
                         "unres", "siraj", "kcde", "glarma", "naive")
cols_models_dengue <- c(brewer.pal(8, "Dark2")[c(4, 5, 7, 6, 3:1, 8)], "black")
names(cols_models_dengue) <- names_models_dengue

cols_models_dengue_transp <- modify_alpha(cols_models_dengue, 0.6)
names(cols_models_dengue_transp) <- names(cols_models_dengue)


# For noro/rota:
data("noroBE")
data("rotaBE")

# modify neighbourhood matrices to apply power law:
noroBE@neighbourhood <- noroBE@neighbourhood + 1
rotaBE@neighbourhood <- rotaBE@neighbourhood + 1

names_districts <- c("chwi", "frkr", "lich", "mahe", "mitt", "neuk",
                     "pank", "rein", "span", "zehl", "scho", "trko")
full_names_districts <- c("Charlottenburg-Wilmersdorf", "Friedrichshain-Kreuzberg", "Lichtenberg",
                          "Marzahn-Hellersdorf", "Mitte", "Neukölln", "Pankow", "Reinickendorf",
                          "Spandau", "Steglitz-Zehlendorf", "Tempelhof-Schöneberg", "Treptow-Köpenick")
names(full_names_districts) <-  names_districts
names_diseases <- c("noro", "rota")

names_hhh4_model_versions_nr <- c("full", "gravity")

names_hhh4_lag_structures_nr <- c("ar1", "pois", "lin", "geom", "unres")

names_cols_nr <- c("ar1", "pois", "lin", "geom",
                   "unres", "glarma", "naive")
cols_models_nr <- cols_models_dengue[names_cols_nr]

cols_models_nr_transp <- modify_alpha(cols_models_dengue, 0.6)
names(cols_models_nr_transp) <- names(cols_models_nr)