
#' @export

Kprototypes <- R6::R6Class(
  "K-prototypes",
  
  inherit = .Cluster,
  
  private = list(
    
    # --- Private attributs ---
    
    FScale       = NULL,  # TRUE / FALSE (standardisation des variables quantitatives actives)
    FNumCols     = NULL,  # indices des variables quantitatives dans FX_active
    FCatCols     = NULL,  # indices des variables qualitatives dans FX_active
    FLambda      = NULL,  # pondération partie catégorielle (k-prototypes)
    
    # --- Helper methods ---
    
    prepare_X = function(X, update_structure = TRUE) {
      # Check X matrix
      private$check_X(X, update_structure)
      # Identification des colonnes numériques / catégorielles
      num_idx <- which(vapply(X, is.numeric, logical(1L)))
      cat_idx <- which(vapply(X, function(col) is.factor(col) || is.character(col), logical(1L)))
      
      # Conversion des qualitatives en facteur
      if (length(cat_idx) > 0L) {
        for (j in cat_idx) {
          if (!is.factor(X[[j]])) {
            X[[j]] <- factor(X[[j]])
          }
        }
      }
      
      if (update_structure) {
        private$FNumCols <- num_idx
        private$FCatCols <- cat_idx
      }
      
      return(X)
    },
    
    compute_num_profile = function(X, vars_k, num_idx) {
      # vars_k : indices de colonnes dans le cluster k
      num_in_cluster <- intersect(vars_k, num_idx)
      if (length(num_in_cluster) == 0L) return(NULL)
      
      Xk_num <- as.matrix(X[, num_in_cluster, drop = FALSE])
      rowMeans(Xk_num, na.rm = TRUE)  # profil numérique : moyenne par individu
    },
    
    compute_cat_profile = function(X, vars_k, cat_idx) {
      # vars_k : indices de colonnes dans le cluster k
      cat_in_cluster <- intersect(vars_k, cat_idx)
      if (length(cat_in_cluster) == 0L) return(NULL)
      
      Xk_cat <- X[, cat_in_cluster, drop = FALSE]
      # profil catégoriel : mode par individu
      apply(Xk_cat, 1L, function(row_i) {
        row_i <- as.character(row_i)
        tab <- table(row_i, useNA = "no")
        if (length(tab) == 0L) NA_character_ else names(tab)[which.max(tab)]
      })
    },
    
    # Rattachement d'une seule variable descriptive
    predict_one_variable = function(x_new, var_name) {
      # x_new : vecteur (numérique ou factor/character)
      # var_name : nom de la variable descriptive
      
      method <- private$FMethod
      n <- nrow(private$FX_active)
      K <- private$FNbGroupes
      
      if (length(x_new) != n) {
        stop("predict_one_variable() : variable '", var_name,
             "' doesn't have the same number of observations as the active data.")
      }
      
      centers    <- private$FCenters   # liste de K prototypes mixtes (num + cat)
      X_active   <- private$FX_active
      clusters   <- private$FClusters
      num_idx    <- private$FNumCols
      cat_idx    <- private$FCatCols
      lambda     <- private$FLambda
      
      if (is.null(centers) || is.null(clusters)) {
        stop("predict_one_variable() : prototypes are not available. Call fit() first.")
      }
      
      # --- Cas numérique : même logique que k-means, avec prototypes$num ---
      if (is.numeric(x_new)) {
        distances <- numeric(K)
        adhesions <- numeric(K)
        
        for (k in seq_len(K)) {
          zk_num <- centers[[k]]$num
          if (is.null(zk_num)) {
            distances[k] <- NA_real_
            adhesions[k] <- NA_real_
          } else {
            r <- suppressWarnings(stats::cor(x_new, zk_num, use = "pairwise.complete.obs"))
            if (is.na(r)) r <- 0
            d2 <- 1 - r^2
            distances[k] <- d2
            adhesions[k] <- r^2
          }
        }
        
        best_k  <- which.min(distances)
        best_d  <- distances[best_k]
        best_ad <- adhesions[best_k]
        
        res <- data.frame(
          variable      = var_name,
          cluster       = best_k,
          adhesion      = best_ad,
          distance_best = best_d,
          stringsAsFactors = FALSE
        )
        for (k in seq_len(K)) {
          col_name <- paste0("d", k)
          res[[col_name]] <- distances[k]
        }
        return(res)
      }
      
      # --- Cas qualitatif : simple matching sur profil catégoriel par cluster ---
      if (!(is.factor(x_new) || is.character(x_new))) {
        stop("predict_one_variable() : for 'kprototypes', variable '",
             var_name, "' must be either numerical or qualitative (factor/character).")
      }
      x_char <- as.character(x_new)
      
      distances <- numeric(K)
      adhesions <- numeric(K)
      
      for (k in seq_len(K)) {
        vars_k <- which(clusters == k)               # variables actives dans le cluster k
        cat_in_cluster <- intersect(vars_k, cat_idx) # indices globaux des quali dans ce cluster
        
        if (length(cat_in_cluster) == 0L) {
          # Aucun support qualitatif dans ce cluster → dissimilarité max
          distances[k] <- lambda * 1
          adhesions[k] <- 0
        } else {
          Xk_cat <- X_active[, cat_in_cluster, drop = FALSE]
          # Profil catégoriel par individu : modalité la plus fréquente dans le cluster
          zk_cat <- apply(Xk_cat, 1L, function(row_i) {
            row_i <- as.character(row_i)
            tab <- table(row_i, useNA = "no")
            if (length(tab) == 0L) NA_character_ else names(tab)[which.max(tab)]
          })
          
          mismatch <- x_char != zk_cat
          d_raw <- mean(mismatch, na.rm = TRUE)  # simple matching dans [0,1]
          if (is.na(d_raw)) d_raw <- 1
          distances[k] <- lambda * d_raw         # même échelle que dans l'algo
          adhesions[k] <- 1 - d_raw              # proportion de matches, non pondérée
        }
      }
      
      best_k  <- which.min(distances)
      best_d  <- distances[best_k]
      best_ad <- adhesions[best_k]
      
      res <- data.frame(
        variable      = var_name,
        cluster       = best_k,
        adhesion      = best_ad,
        distance_best = best_d,
        stringsAsFactors = FALSE
      )
      for (k in seq_len(K)) {
        col_name <- paste0("d", k)
        res[[col_name]] <- distances[k]
      }
      
      return(res)
    },
    
    compute_membership = function(X, clusters, centers) {
      p <- length(clusters)
      membership <- rep(NA_real_, p)
      
      # r^2 pour les variables numériques,
      # 1 - dissimilarité simple matching pour les variables qualitatives
      num_idx <- private$FNumCols
      cat_idx <- private$FCatCols
      X_used  <- X
      
      for (j in seq_len(p)) {
        kj      <- clusters[j]
        proto_k <- centers[[kj]]
        
        if (j %in% num_idx) {
          # Partie numérique : r^2(X_j, P_k^num)
          xj     <- X_used[[j]]
          zk_num <- proto_k$num
          if (is.null(zk_num)) {
            membership[j] <- NA_real_
          } else {
            r <- suppressWarnings(
              stats::cor(xj, zk_num, use = "pairwise.complete.obs")
            )
            if (is.na(r)) r <- 0
            membership[j] <- r^2
          }
          
        } else if (j %in% cat_idx) {
          # Partie catégorielle : 1 - d_raw, avec d_raw dans [0,1]
          xj_char <- as.character(X_used[[j]])
          zk_cat  <- proto_k$cat   # vecteur de longueur n
          
          if (is.null(zk_cat) || length(zk_cat) != length(xj_char)) {
            membership[j] <- 0
          } else {
            mismatch <- xj_char != zk_cat
            d_raw <- mean(mismatch, na.rm = TRUE)
            if (is.na(d_raw)) d_raw <- 1
            membership[j] <- 1 - d_raw
          }
          
        } else {
          membership[j] <- NA_real_
        }
      }
      
      membership_label <- "Adherence score (r^2 for numbers, 1 - dissimilarity for categories)"
      
      return(list(
        content = membership,
        label = membership_label
      ))
    },
    
    # --- Main algorithm method ---
    
    run_kprototypes = function(X) {
      # X : data.frame des variables actives
      n <- nrow(X)
      p <- ncol(X)
      
      num_idx <- private$FNumCols   # indices des variables quantitatives
      cat_idx <- private$FCatCols   # indices des variables qualitatives
      
      K       <- private$FNbGroupes
      lambda  <- private$FLambda
      
      max_iter <- 100L
      tol      <- 1e-6
      
      if (K < 1L || K > p) {
        stop("K must be between 1 and the number of active variables")
      }
      
      # --- helpers distances ---
      
      dist_num <- function(j, proto_num) {
        # j : indice global de la variable
        if (is.null(proto_num) || !(j %in% num_idx)) return(Inf)
        xj <- as.numeric(X[[j]])
        r  <- suppressWarnings(cor(xj, proto_num,
                                   use = "pairwise.complete.obs"))
        if (is.na(r)) 1 else 1 - r^2
      }
      
      dist_cat <- function(j, proto_cat) {
        # j : indice global de la variable
        if (is.null(proto_cat) || !(j %in% cat_idx)) return(Inf)
        xj <- as.character(X[[j]])
        mismatch <- xj != proto_cat
        d <- mean(mismatch, na.rm = TRUE)
        if (is.na(d)) d <- 1
        lambda * d
      }
      
      # --- initialisation des clusters ---
      
      # on choisit K variables différentes comme noyaux initiaux
      init_idx <- sample.int(p, size = K, replace = FALSE)
      clusters <- rep(NA_integer_, p)
      clusters[init_idx] <- seq_len(K)
      
      # prototypes initiaux (à partir des noyaux uniquement)
      prototypes <- vector("list", K)
      for (k in seq_len(K)) {
        vars_k <- which(clusters == k)
        prototypes[[k]] <- list(
          num = private$compute_num_profile(X, vars_k, num_idx),
          cat = private$compute_cat_profile(X, vars_k, cat_idx)
        )
      }
      
      # affectation initiale des variables non noyaux
      for (j in which(is.na(clusters))) {
        d_k <- numeric(K)
        if (j %in% num_idx) {
          for (k in seq_len(K)) d_k[k] <- dist_num(j, prototypes[[k]]$num)
        } else {
          for (k in seq_len(K)) d_k[k] <- dist_cat(j, prototypes[[k]]$cat)
        }
        clusters[j] <- which.min(d_k)
      }
      
      # --- boucle de réallocation ---
      
      inertia <- Inf
      converged <- FALSE
      iter <- 0L
      
      while (iter < max_iter && !converged) {
        iter <- iter + 1L
        
        # 1) recalcul des prototypes pour chaque cluster
        prototypes <- vector("list", K)
        for (k in seq_len(K)) {
          vars_k <- which(clusters == k)
          prototypes[[k]] <- list(
            num = private$compute_num_profile(X, vars_k, num_idx),
            cat = private$compute_cat_profile(X, vars_k, cat_idx)
          )
        }
        
        # 2) réaffectation de toutes les variables
        new_clusters <- clusters
        for (j in seq_len(p)) {
          d_k <- numeric(K)
          if (j %in% num_idx) {
            for (k in seq_len(K)) d_k[k] <- dist_num(j, prototypes[[k]]$num)
          } else {
            for (k in seq_len(K)) d_k[k] <- dist_cat(j, prototypes[[k]]$cat)
          }
          new_clusters[j] <- which.min(d_k)
        }
        
        # 3) calcul de l'inertie intra-cluster
        new_inertia <- 0
        for (j in seq_len(p)) {
          k <- new_clusters[j]
          if (j %in% num_idx) {
            new_inertia <- new_inertia + dist_num(j, prototypes[[k]]$num)
          } else {
            new_inertia <- new_inertia + dist_cat(j, prototypes[[k]]$cat)
          }
        }
        
        # 4) test de convergence
        if (all(new_clusters == clusters) ||
            abs(inertia - new_inertia) < tol) {
          converged <- TRUE
          clusters  <- new_clusters
          inertia   <- new_inertia
          break
        }
        
        clusters <- new_clusters
        inertia  <- new_inertia
      }
      
      res <- list(
        clusters   = clusters,
        centers = prototypes,
        inertia    = inertia,
        converged  = converged
      )
      
      return(res)
    }
  ),
  
  public = list(
    
    initialize = function(K = 2L, lambda = 1, scale = TRUE, random_state = NULL) {
      if (K < 2L) stop("K must be greater than or equal to 2")
      if (lambda <= 0) stop("lambda must be strictly positive.")
      
      private$FMethod <- "kprototypes"
      private$FNbGroupes <- as.integer(K)
      private$FScale     <- isTRUE(scale)
      private$FLambda    <- lambda
      
      set.seed(random_state)
    },
    
    fit = function(X) {
      X <- private$prepare_X(X, update_structure = TRUE)
      if (isTRUE(private$FScale)) {
        X <- private$scale_active_variables(X, private$FNumCols)
      }
      
      private$FX_active <- X
      
      res <- private$run_kprototypes(X)
      
      private$FClusters    <- res$clusters
      private$FCenters     <- res$centers
      private$FInertia     <- res$inertia
      private$FConvergence <- isTRUE(res$converged)
      
      invisible(self)
    },
    
    predict = function(X_new) {
      if (is.null(private$FX_active)) {
        stop("The model has not yet been trained. Call fit() first.")
      }
      private$check_X_new(X_new)
      
      res_list <- vector("list", length = ncol(X_new))
      for (j in seq_len(ncol(X_new))) {
        res_list[[j]] <- private$predict_one_variable(
          x_new    = X_new[[j]],
          var_name = colnames(X_new)[j]
        )
      }
      
      private$FX_descr <- X_new # mémorisation des variables descriptives
      
      res <- do.call(rbind, res_list)
      rownames(res) <- NULL
      return(res)
    }
  )
)