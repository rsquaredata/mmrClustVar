# data-raw/prepare_datasets.R
# Creation of the datasets bundled with the mmrClustVar package

library(usethis)

### ---------------------------
###  NUMERIC DATASETS
### ---------------------------

iris_num     <- iris[, 1:4]
mtcars_num   <- mtcars
airquality_num <- airquality

### ---------------------------
###  CATEGORICAL DATASETS
### ---------------------------

# Arthritis (from package vcd)
if (!requireNamespace("vcd", quietly = TRUE)) {
    install.packages("vcd")
}
library(vcd)
data("Arthritis")
arthritis_cat <- Arthritis

# Titanic (from package titanic)
if (!requireNamespace("titanic", quietly = TRUE)) {
    install.packages("titanic")
}
library(titanic)
titanic_cat <- titanic::titanic_train

# HouseVotes84 (from package mlbench)
if (!requireNamespace("mlbench", quietly = TRUE)) {
    install.packages("mlbench")
}
library(mlbench)
data("HouseVotes84", package = "mlbench")
housevotes_cat <- HouseVotes84

### ---------------------------
###  MIXED DATASETS
### ---------------------------

# iris_mixed = Species as character
iris_mixed <- iris
iris_mixed$Species <- as.character(iris$Species)

# Credit (from package ISLR) – mixed dataset
if (!requireNamespace("ISLR", quietly = TRUE)) {
    install.packages("ISLR")
}
library(ISLR)
credit_mix <- ISLR::Credit

# AdultUCI (from package arules)
if (!requireNamespace("arules", quietly = TRUE)) {
    install.packages("arules")
}
library(arules)
data("AdultUCI")
adult <- AdultUCI

# Sample of 200 rows from AdultUCI
set.seed(1)
adult_small <- adult[sample(nrow(adult), 200), ]

### ---------------------------
###  Save datasets into /data
### ---------------------------

usethis::use_data(
    iris_num, mtcars_num, airquality_num,
    arthritis_cat, titanic_cat, housevotes_cat,
    iris_mixed, credit_mix, adult_small,
    overwrite = TRUE
)