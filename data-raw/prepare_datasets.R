# data-raw/prepare_datasets.R
# Création des jeux de données intégrés au package mmrClustVar

library(usethis)

### ---------------------------
###  Jeux NUMÉRIQUES
### ---------------------------

iris_num <- iris[, 1:4]
mtcars_num <- mtcars
airquality_num <- airquality

### ---------------------------
###  Jeux CATÉGORIELS
### ---------------------------

# Arthritis (package vcd)
if (!requireNamespace("vcd", quietly = TRUE)) {
    install.packages("vcd")
}
library(vcd)
data("Arthritis")
arthritis_cat <- Arthritis

# Titanic (package titanic)
if (!requireNamespace("titanic", quietly = TRUE)) {
    install.packages("titanic")
}
library(titanic)
titanic_cat <- titanic::titanic_train

# HouseVotes84 (package mlbench)
if (!requireNamespace("mlbench", quietly = TRUE)) {
    install.packages("mlbench")
}
library(mlbench)
data("HouseVotes84", package = "mlbench")
housevotes_cat <- HouseVotes84

### ---------------------------
###  Jeux MIXTES
### ---------------------------

# iris mixte = Species en caractère
iris_mixed <- iris
iris_mixed$Species <- as.character(iris$Species)

# Credit (package ISLR) - jeu mixte
if (!requireNamespace("ISLR", quietly = TRUE)) {
    install.packages("ISLR")
}
library(ISLR)
credit_mix <- ISLR::Credit

# AdultUCI (package arules)
if (!requireNamespace("arules", quietly = TRUE)) {
    install.packages("arules")
}
library(arules)
data("AdultUCI")
adult <- AdultUCI

set.seed(1)
adult_small <- adult[sample(nrow(adult), 200), ]

### ---------------------------
###  Sauvegarde dans /data
### ---------------------------

usethis::use_data(
    iris_num, mtcars_num, airquality_num,
    arthritis_cat, titanic_cat, housevotes_cat,
    iris_mixed, credit_mix, adult_small,
    overwrite = TRUE
)