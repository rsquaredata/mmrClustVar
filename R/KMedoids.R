#' K-Medoids Engine for Variable Clustering (Cramér's V)
#'
#' Variable clustering by k-medoids on a dissimilarity matrix
#' based on Cramér's V between variables (columns).
#'
#' Dissimilarity between two variables j and l:
#'   d(j, l) = 1 - V(X[, j], X[, l])
#'
#' This class is intended to be used internally by mmrClustVar.
#'
mmrClustVarKMedoids <- R6::R6Class(
  "KMedoids",
  
  inherit = mmrClustVarBase,
  
  public = list(
    
    #' @description
    #' Constructor for the K-medoids engine.
    #' @param K Number of clusters (>= 2)
    #' @param scale Ignored here (Cramér's V works on categorical vars),
    #'   kept only for compatibility with the facade.
    #' @param lambda Unused here, kept for compatibility.
    initialize = function(K,
                          scale  = FALSE,
                          lambda = 1,
                          ...) {
      
      # For Cramér's V-based k-medoids, we do NOT scale
      super$initialize(
        K          = K,
        scale      = FALSE,
        lambda     = lambda,
        method_name = "kmedoids"
      )
    }
    
    # fit(), predict(), print(), summary(), plot() are inherited
    # from mmrClustVarBase and use our private methods below.
  ),
  
  private = list(
    
    # Indices (or names) of the medoid variables
    FMedoids = NULL,
    
    # ------------------------------------------------------------
    # Cramér's V between two variables (vectors) x and y
    # Both are converted to factors (with discretization for numeric).
    # ------------------------------------------------------------
    cramer_v = function(x, y, n_bins = 5L) {
      
      # Helper to convert numeric -> factor via discretization
      to_factor <- function(z) {
        if (is.numeric(z)) {
          # Discretization into n_bins categories
          # include.lowest = TRUE to keep all obs
          z_cut <- try(
            cut(z, breaks = n_bins, include.lowest = TRUE),
            silent = TRUE
          )
          if (inherits(z_cut, "try-error")) {
            # Fallback: all in one category if cut fails
            factor(rep("bin1", length(z)))
          } else {
            z_cut
          }
        } else if (!is.factor(z)) {
          factor(z)
        } else {
          z
        }
      }
      
      x_f <- to_factor(x)
      y_f <- to_factor(y)
      
      tbl <- table(x_f, y_f)
      n   <- sum(tbl)
      
      # If degenerate (1 modalité d’un côté), V = 0
      if (n == 0L || nrow(tbl) < 2L || ncol(tbl) < 2L) {
        return(0)
      }
      
      # Chi-square test without continuity correction
      suppressWarnings({
        chi2 <- stats::chisq.test(tbl, correct = FALSE)
      })
      
      chi2_stat <- as.numeric(chi2$statistic)
      if (!is.finite(chi2_stat) || chi2_stat < 0) {
        return(0)
      }
      
      min_dim <- min(nrow(tbl), ncol(tbl))
      if (min_dim <= 1L) {
        return(0)
      }
      
      v <- sqrt(chi2_stat / (n * (min_dim - 1)))
      
      if (!is.finite(v) || is.na(v)) {
        v <- 0
      }
      
      v
    },
    
    # ------------------------------------------------------------
    # Core K-medoids algorithm for variables
    # X: data.frame with columns = variables
    # Returns list(clusters, centers, inertia, converged)
    # ------------------------------------------------------------
    run_clustering = function(X) {
      
      if (!requireNamespace("cluster", quietly = TRUE)) {
        stop("[mmrClustVarKMedoids] Package 'cluster' is required (for pam).")
      }
      
      p <- ncol(X)
      var_names <- colnames(X)
      
      if (p < 2L) {
        stop("[mmrClustVarKMedoids] Need at least 2 variables (columns).")
      }
      
      # --- 1) Compute dissimilarity matrix D (p x p) using 1 - Cramér's V ----
      D <- matrix(0, nrow = p, ncol = p)
      colnames(D) <- rownames(D) <- var_names
      
      for (j in seq_len(p - 1L)) {
        for (l in seq.int(j + 1L, p)) {
          v_jl <- private$cramer_v(X[[j]], X[[l]])
          d_jl <- 1 - v_jl      # dissimilarité = 1 - V
          D[j, l] <- d_jl
          D[l, j] <- d_jl
        }
      }
      
      diag(D) <- 0
      
      # --- 2) Apply K-medoids (PAM) on the dissimilarity matrix --------------
      K <- private$FNbGroupes
      
      pam_res <- cluster::pam(
        x    = D,
        k    = K,
        diss = TRUE
      )
      
      # pam_res$clustering: vector of length p, cluster of each variable
      clusters_vec <- pam_res$clustering
      names(clusters_vec) <- var_names
      
      # pam_res$id.med: indices of medoids in 1..p
      medoid_indices <- pam_res$id.med
      medoid_names   <- var_names[medoid_indices]
      
      private$FMedoids <- medoid_names
      
      # --- 3) Define "centers" as a small data.frame of medoid variables -----
      centers <- data.frame(
        cluster = seq_along(medoid_names),
        medoid  = medoid_names,
        stringsAsFactors = FALSE
      )
      
      # --- 4) Inertia: sum of dissimilarities to each cluster's medoid -------
      inertia <- 0
      for (k in seq_along(medoid_indices)) {
        m_idx    <- medoid_indices[k]
        members  <- which(clusters_vec == k)
        if (length(members) > 0L) {
          inertia <- inertia + sum(D[members, m_idx])
        }
      }
      
      list(
        clusters  = clusters_vec,
        centers   = centers,
        inertia   = inertia,
        converged = TRUE
      )
    },
    
    # ------------------------------------------------------------
    # Assign a new variable to the closest medoid
    # x_new: vector (new variable)
    # var_name: name of the variable
    # ------------------------------------------------------------
    predict_one_variable = function(x_new, var_name) {
      
      if (is.null(private$FX_active) || is.null(private$FMedoids)) {
        stop("[mmrClustVarKMedoids] Model must be fitted before predict().")
      }
      
      medoids <- private$FMedoids
      X_ref   <- private$FX_active   # data.frame of active vars (columns)
      
      # Compute dissimilarity to each medoid
      dists <- sapply(medoids, function(m_name) {
        x_ref <- X_ref[[m_name]]
        v     <- private$cramer_v(x_new, x_ref)
        1 - v
      })
      
      best_k <- which.min(dists)
      
      out <- data.frame(
        variable            = var_name,
        cluster             = best_k,
        distance_to_medoid  = dists[best_k],
        stringsAsFactors    = FALSE
      )
      
      out
    }
  )
)
