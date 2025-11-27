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
###  METAL UNIVERSE DATASET (mixed)
### ---------------------------

# Load the raw CSV file
metal_universe <- read.csv(
    "data-raw/metal_universe_raw.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# Fix potential trailing spaces in column names
names(metal_universe) <- sub("cat_production_preference\\s*$",
                             "cat_production_preference",
                             names(metal_universe))
names(metal_universe) <- sub("cat_mosh_preference\\s*$",
                             "cat_mosh_preference",
                             names(metal_universe))

# Convert core metadata
metal_universe$group_name       <- as.character(metal_universe$group_name)
metal_universe$subgenre         <- as.factor(metal_universe$subgenre)
metal_universe$country          <- as.factor(metal_universe$country)
metal_universe$cat_front_gender <- as.factor(metal_universe$cat_front_gender)

# Convert numeric variables
num_cols <- grep("^num_", names(metal_universe), value = TRUE)
for (col in num_cols) {
    metal_universe[[col]] <- as.numeric(metal_universe[[col]])
}

# Convert categorical variables
cat_cols <- grep("^cat_", names(metal_universe), value = TRUE)
for (col in cat_cols) {
    metal_universe[[col]] <- as.factor(metal_universe[[col]])
}


### ---------------------------
###  Save datasets into /data
### ---------------------------

usethis::use_data(
    iris_num, mtcars_num, airquality_num,
    arthritis_cat, titanic_cat, housevotes_cat,
    iris_mixed, credit_mix, adult_small,
    metal_universe,
    overwrite = TRUE
)