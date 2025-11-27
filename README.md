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

``` r
devtools::install_github("username/mmrClustVar")
library(mmrClustVar)
```

------------------------------------------------------------------------

## Quick Start — Using the R6 Class

### 1. Fit

``` r
df <- iris[,1:4]

obj <- mmrClustVar$new(
  method="kmeans",
  K=3,
  scale=TRUE
)

obj$fit(df)
```

### 2. Summary

``` r
obj$summary()
```

### 3. Clusters & Centers

``` r
obj$get_clusters()
obj$get_centers()
obj$get_inertia()
```

### 4. Predict New Variables

``` r
obj$predict(df[,1,drop=FALSE])
```

### 5. Plots

``` r
obj$plot("clusters")
obj$plot("inertia")
obj$plot("membership")
obj$plot("profiles")
```

### 6. Inertia Path

``` r
obj$compute_inertia_path(K_seq=2:6, X=df)
obj$plot("inertia")
```

------------------------------------------------------------------------

# Using the Shiny App

``` r
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

-   Chavent et al. (2012), [*ClustOfVar: An R Package for the Clustering of Variables*](https://arxiv.org/pdf/1112.0295), JSS.
-   Husson, Josse & Pagès (2017), [*Exploratory Multivariate Analysis by Example using R*](http://staff.ustc.edu.cn/~ynyang/vector/books/Husson-Le-Pages.pdf).
-   MacQueen (1967), [*Some methods for classification and analysis of multivariate observations*](https://www.cs.cmu.edu/~bhiksha/courses/mlsp.fall2010/class14/macqueen.pdf).
-   Huang (1998), [*Extensions to k-means for categorical values*](https://cse.hkust.edu.hk/~qyang/Teaching/537/Papers/huang98extensions.pdf).
-   Kaufman & Rousseeuw (2005), [*Finding Groups in Data*](https://www.researchgate.net/profile/Peter-Rousseeuw/publication/220695963_Finding_Groups_in_Data_An_Introduction_To_Cluster_Analysis/links/60fbbe85169a1a0103b20e91/Finding-Groups-in-Data-An-Introduction-To-Cluster-Analysis.pdf).
-   Rakotomalala (2025), *R programming & classification lectures*.
-   R Core Team (2025), [*R Language*](https://www.r-project.org/).
-   Chang (2025), [*R6: Encapsulated object-oriented programming for R*](https://r6.r-lib.org/)

------------------------------------------------------------------------

# Authors

-   Marin Nagy
-   Mazilda Zehraoui
-   Rina Razafimahefa

------------------------------------------------------------------------

# License

MIT License
