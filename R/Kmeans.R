


#' @export

Kmeans <- R6::R6Class(
  "K-means",
  
  inherit = .Cluster,
  
  private = list(
    
    # --- Private attributs ---
    
    FScale       = NULL,  # TRUE / FALSE (standardisation des variables quantitatives actives)
    FNumCols     = NULL,  # indices des variables quantitatives dans FX_active
    
    # --- Helper methods ---
    
    prepare_X = function(X, update_structure = TRUE) {
      # Check X matrix
      private$check_X(X, update_structure)
      
      if (update_structure) {
        # Identification des colonnes numériques / catégorielles
        num_idx <- which(vapply(X, is.numeric, logical(1L)))
        private$FNumCols <- num_idx
        
        if (length(num_idx) == 0L) {
          stop("No quantitative variables found.")
        }
      }
      
      return(X)
    },
    
    # Rattachement d'une seule variable descriptive
    predict_one_variable = function(x_new, var_name) {
      # x_new : vecteur (numérique ou factor/character)
      # var_name : nom de la variable descriptive
      
      method <- private$FMethod
      n <- nrow(private$FX_active)
      K <- private$FNbGroupes
      
      if (length(x_new) != n) {
        stop("predict_one_variable(): variable '", var_name,
             "' doesn't have the same number of observations as the active data.")
      }
      
      if (!is.numeric(x_new)) {
        stop("predict_one_variable(): for “kmeans”, variable '",
             var_name, "'must be numeric.")
      }
      
      centers <- private$FCenters  # liste de K composantes latentes Z_k
      if (is.null(centers)) {
        stop("predict_one_variable(): prototypes are not available. Call fit() first.")
      }
      
      distances <- numeric(K)
      adhesions <- numeric(K)
      
      for (k in seq_len(K)) {
        zk <- centers[[k]]
        if (is.null(zk)) {
          distances[k] <- NA_real_
          adhesions[k] <- NA_real_
        } else {
          r <- suppressWarnings(stats::cor(x_new, zk, use = "pairwise.complete.obs"))
          if (is.na(r)) r <- 0
          d2 <- 1 - r^2
          distances[k] <- d2
          adhesions[k] <- r^2
        }
      }
      
      best_k   <- which.min(distances)
      best_d2  <- distances[best_k]
      best_r2  <- adhesions[best_k]
      
      res <- data.frame(
        variable      = var_name,
        cluster       = best_k,
        adhesion      = best_r2,
        distance_best = best_d2,
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
      
      # r^2(X_j, Z_k(j))
      X_mat <- as.matrix(X)
      for (j in seq_len(p)) {
        kj <- clusters[j]
        zk <- centers[[kj]]
        xj <- X_mat[, j]
        r  <- suppressWarnings(stats::cor(xj, zk, use = "pairwise.complete.obs"))
        if (is.na(r)) r <- 0
        membership[j] <- r^2
      }
      membership_label <- "r^2 (correlation with the latent comp of the cluster)"
      
      return(list(
        content = membership,
        label = membership_label
      ))
    },
    
    # --- Main algorithm method ---
    
    run_kmeans = function(X) {
      # X : data.frame ou matrice n x p, uniquement des variables quantitatives
      K <- private$FNbGroupes
      n <- nrow(X)
      p <- ncol(X)
      
      # Sécurité : vérifier que tout est numérique
      is_num <- vapply(X, is.numeric, logical(1L))
      if (!all(is_num)) {
        stop("run_kmeans() : X must contain only quantitative variables.")
      }
      
      X_mat <- as.matrix(X)
      
      # Initialisation des clusters (partition des variables)
      # On essaie de répartir à peu près équitablement
      clusters <- rep(seq_len(K), length.out = p)
      clusters <- sample(clusters)  # randomisation
      
      max_iter    <- 100L
      tol         <- 1e-6
      inertia_old <- Inf
      converged   <- FALSE
      Z_list      <- vector("list", K)  # composantes latentes
      
      for (iter in seq_len(max_iter)) {
        # 1) Calcul des composantes latentes Z_k pour chaque cluster
        for (k in seq_len(K)) {
          vars_k <- which(clusters == k)
          
          if (length(vars_k) == 0L) {
            # Cluster vide : on choisit une variable au hasard et on la force dans ce cluster
            cand <- sample(seq_len(p), 1L)
            clusters[cand] <- k
            vars_k <- cand
          }
          
          Xk <- X_mat[, vars_k, drop = FALSE]
          
          if (ncol(Xk) == 1L) {
            # Une seule variable : la composante latente est le profil de cette variable
            Z_list[[k]] <- as.numeric(Xk[, 1L])
          } else {
            # Plusieurs variables : on prend la 1ère composante principale
            # (scores sur les individus)
            pc <- stats::prcomp(Xk, center = FALSE, scale. = FALSE)
            Z_list[[k]] <- as.numeric(pc$x[, 1L])
          }
        }
        
        # 2) Calcul des distances d^2(j,k) = 1 - r^2(X_j, Z_k)
        dist_mat <- matrix(NA_real_, nrow = p, ncol = K)
        for (j in seq_len(p)) {
          xj <- X_mat[, j]
          for (k in seq_len(K)) {
            zk <- Z_list[[k]]
            r  <- suppressWarnings(stats::cor(xj, zk, use = "pairwise.complete.obs"))
            if (is.na(r)) r <- 0
            d2 <- 1 - r^2
            dist_mat[j, k] <- d2
          }
        }
        
        # 3) Nouvelle affectation : pour chaque variable, le cluster de distance minimale
        new_clusters <- apply(dist_mat, 1L, which.min)
        
        # Inertie intra-cluster = somme des distances au meilleur cluster
        min_d2   <- dist_mat[cbind(seq_len(p), new_clusters)]
        inertia  <- sum(min_d2, na.rm = TRUE)
        
        # 4) Critère de convergence
        if (all(new_clusters == clusters) || abs(inertia_old - inertia) < tol) {
          converged   <- TRUE
          clusters    <- new_clusters
          inertia_old <- inertia
          break
        }
        
        clusters    <- new_clusters
        inertia_old <- inertia
      }
      
      # Résultat : on retourne la partition et les composantes latentes
      res <- list(
        clusters = clusters,       # vecteur de longueur p
        centers  = Z_list,         # liste de K vecteurs (longueur n)
        inertia  = inertia_old,    # inertie intra totale
        converged = converged
      )
      
      return(res)
    }
  ),
  
  public = list(
    
    initialize = function(K = 2L, scale = TRUE, random_state = NULL) {
      if (K < 2L) stop("K must be greater than or equal to 2.")
      
      private$FMethod <- "kmeans"
      private$FNbGroupes <- as.integer(K)
      private$FScale     <- isTRUE(scale)
      set.seed(random_state)
    },
    
    fit = function(X) {
      X <- private$prepare_X(X, update_structure = TRUE)
      if (isTRUE(private$FScale)) {
        X <- private$scale_active_variables(X, private$FNumCols)
      }
      private$FX_active <- X
      
      res <- private$run_kmeans(X)
      
      private$FClusters    <- res$clusters
      private$FCenters     <- res$centers
      private$FInertia     <- res$inertia
      private$FConvergence <- isTRUE(res$converged)
      
      invisible(res)
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