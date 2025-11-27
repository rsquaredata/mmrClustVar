# mmrClustVar

R package implementing multiple variable clustering algorithms (numeric, categorical, mixed), with an R6-based architecture, interpretability tools, and an interactive Shiny app.

------------------------------------------------------------------------

## Overview

**mmrClustVar** is an academic R package developed as part of the Master 2 SISE (Statistics & Computer Science for Data Science) at Université Lyon 2.

It provides:

- k-means (numeric)
- k-modes (categorical)
- k-prototypes (mixed)
- k-medoids (general)
- Automatic/assisted selection of K
- Interpretability tools (inertia, adhesion, profiles)
- Shiny application for exploration
- Full export functionalities

------------------------------------------------------------------------

## Repository Structure

```         
mmrClustVar/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── mmrClustVar.R
│   ├── mmrClustVarBase.R
│   ├── mmrClustVarKMeans.R
│   ├── mmrClustVarKModes.R
│   ├── mmrClustVarKPrototypes.R
│   ├── mmrClustVarKMedoids.R
│   ├── run_app.R
│   ├── utils_inertia.R
│   └── datasets_documentation.R
├── inst/shiny/mmrClustVar_app/
├── data/
├── data-raw/prepare_datasets.R
└── README.md
```

------------------------------------------------------------------------

## Installation

```r
install.packages("remotes")   # if needed
remotes::install_github("rsquaredata/mmrClustVar")

library(mmrClustVar)
```

------------------------------------------------------------------------

## Quick Start — Using the R6 Class

### 1. Fit

```r
df <- iris[,1:4]

obj <- mmrClustVar$new(
  method="kmeans",
  K=3,
  scale=TRUE
)

obj$fit(df)
```

### 2. Summary

```r
obj$summary()
```

### 3. Clusters & Centers

```r
obj$get_clusters()
obj$get_centers()
obj$get_inertia()
```

### 4. Predict New Variables

```r
obj$predict(df[,1,drop=FALSE])
```

### 5. Plots

```r
obj$plot("clusters")
obj$plot("inertia")
obj$plot("membership")
obj$plot("profiles")
```

### 6. Inertia Path

```r
obj$compute_inertia_path(K_seq=2:6, X=df)
obj$plot("inertia")
```

------------------------------------------------------------------------

# Using the Shiny App

```r
run_mmrClustVar_app()
```

Features:

- Import datasets or upload CSV/XLSX  
- Select active/supplementary variables  
- All algorithms available  
- Diagnostics + plots  
- Export clusters, summary, full ZIP bundle  

------------------------------------------------------------------------

# References

- Chavent et al. (2012), *ClustOfVar: An R Package for the Clustering of Variables*.
- Husson, Josse & Pagès (2017), *Exploratory Multivariate Analysis by Example using R*.
- MacQueen (1967), *Some methods for classification and analysis of multivariate observations*.
- Huang (1998), *Extensions to k-means for categorical values*.
- Kaufman & Rousseeuw (2005), *Finding Groups in Data*.
- Rakotomalala (2025), *R programming & classification lectures*.
- R Core Team (2025), *R Language*.
- Chang (2025), *R6: Encapsulated object-oriented programming for R*.

------------------------------------------------------------------------

# Authors

- Marin Nagy
- Mazilda Zehraoui
- Rina Razafimahefa

------------------------------------------------------------------------

# License

MIT License
