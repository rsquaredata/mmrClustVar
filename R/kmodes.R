library(R6)
# source("C:/Users/moulo/Desktop/SISE/Programmation R/git/clustering.R")
#install.packages("FactoMineR")
#install.packages("plotly")
library(FactoMineR)  # à mettre dans Depends/Imports du package, pas forcément ici
source("C:/Users/moulo/Desktop/SISE/Programmation R/git/cluster.R")
KModesVariables <- R6Class("KModesVariables",
                           
                           inherit = Clustering,
                           
                           private = list(
                             
                             .X                  = NULL,   # données complètes (toutes variables)
                             .X_active           = NULL,   # sous-ensemble des variables actives (pour le fit)
                             .X_descriptive      = NULL,   # sous-ensemble des variables descriptives (optionnel)
                             
                             .modes              = NULL,   # centres (variables fictives, K x n)
                             .clusters           = NULL,   # assignation des variables actives (longueur = nb var actives)
                             .max_iter           = 50,
                             
                             .cost               = NA_real_,        # coût final total (dissimilarité intra-cluster)
                             .cost_history       = NULL,           # trajectoire des coûts par itération
                             .cluster_costs      = NULL,           # coût par cluster
                             .n_iter             = NULL,           # nb d’itérations réellement effectuées
                             
                             .active_vars        = NULL,           # noms des variables actives
                             .descriptive_vars   = NULL,           # noms des variables descriptives
                             .K_auto             = FALSE,          # TRUE si K a été choisi automatiquement
                             
                             .mca_result         = NULL,           # résultat de l’ACM (FactoMineR::MCA)
                             
                             # ---------------------------
                             # Distance Hamming entre 2 variables (vecteur)
                             # (on la garde pour compatibilité, mais on vectorise ailleurs)
                             # ---------------------------
                             .col_dissimilarity = function(x, y) {
                               sum(x != y)
                             },
                             
                             # ---------------------------
                             # Initialisation CAO (corrigée) sur variables actives
                             # ---------------------------
                             .cao_init = function(X, K) {
                               
                               p <- ncol(X)   # nb variables
                               n <- nrow(X)
                               
                               modes <- matrix(NA, nrow = K, ncol = n)
                               colnames(modes) <- paste0("ind_", 1:n)
                               
                               # 1) Mode initial = variable avec la modalité majoritaire la plus forte
                               freq_score <- sapply(1:p, function(j) {
                                 max(prop.table(table(X[[j]])))
                               })
                               
                               first <- which.max(freq_score)
                               modes[1, ] <- X[[first]]
                               
                               # 2) Modes suivants : variable max distante des précédentes
                               for (k in 2:K) {
                                 
                                 dist_to_existing <- sapply(1:p, function(j) {
                                   sum(sapply(1:(k-1), function(m) {
                                     private$.col_dissimilarity(X[[j]], modes[m, ])
                                   }))
                                 })
                                 
                                 next_idx <- which.max(dist_to_existing)
                                 modes[k, ] <- X[[next_idx]]
                               }
                               
                               as.data.frame(modes, stringsAsFactors = FALSE)
                             },
                             
                             # ---------------------------
                             # Mise à jour des modes
                             # (mode catégoriel pour chaque individu)
                             # ---------------------------
                             .update_modes = function(X, clusters, K) {
                               
                               n <- nrow(X)
                               p <- ncol(X)
                               
                               new_modes <- matrix(NA, nrow = K, ncol = n)
                               
                               for (k in 1:K) {
                                 
                                 vars_in_cluster <- which(clusters == k)
                                 
                                 if (length(vars_in_cluster) == 0) {
                                   # cluster vide : on garde l'ancien mode
                                   new_modes[k, ] <- private$.modes[k, ]
                                   next
                                 }
                                 
                                 # Mode pour chaque individu i (ligne)
                                 for (i in 1:n) {
                                   vals <- X[i, vars_in_cluster]
                                   new_modes[k, i] <- names(which.max(table(vals)))
                                 }
                               }
                               
                               as.data.frame(new_modes, stringsAsFactors = FALSE)
                             },
                             
                             # ---------------------------
                             # Sélection automatique de K
                             # (simple mais correcte)
                             # ---------------------------
                             .auto_select_K = function(X) {
                               
                               K_values <- 2:6
                               costs <- numeric(length(K_values))
                               p <- ncol(X)
                               n <- nrow(X)
                               
                               X_mat <- as.matrix(X)
                               
                               for (i in seq_along(K_values)) {
                                 
                                 K <- K_values[i]
                                 
                                 modes_K <- private$.cao_init(X, K)
                                 modes_mat <- as.matrix(modes_K)
                                 
                                 # Assignation vectorisée
                                 dist_mat <- matrix(0, nrow = K, ncol = p)
                                 for (m in 1:K) {
                                   mode_vec <- as.vector(modes_mat[m, ])
                                   dist_mat[m, ] <- colSums(X_mat != mode_vec)
                                 }
                                 clusters_K <- apply(dist_mat, 2, which.min)
                                 
                                 # Coût total
                                 costs[i] <- sum(dist_mat[cbind(clusters_K, 1:p)])
                               }
                               
                               K_values[which.min(costs)]
                             },
                             
                             # ---------------------------
                             # Sélection automatique des variables actives
                             # (basée sur l’entropie des modalités)
                             # ---------------------------
                             .auto_select_active_vars = function(X) {
                               p <- ncol(X)
                               cols <- colnames(X)
                               
                               entropies <- sapply(X, function(col) {
                                 tab <- prop.table(table(col))
                                 -sum(tab * log(tab + 1e-12))
                               })
                               
                               thr <- median(entropies)
                               active <- cols[entropies >= thr]
                               if (length(active) == 0L || length(active) == p) {
                                 # Sécurité : au pire, on prend la moitié des variables
                                 active <- cols[order(entropies, decreasing = TRUE)][1:ceiling(p / 2)]
                               }
                               descriptive <- setdiff(cols, active)
                               
                               list(active = active, descriptive = descriptive)
                             },
                             
                             # ---------------------------
                             # Calcul complet des distances modes <-> variables (vectorisé)
                             # X_mat : n x p (caractères)
                             # modes : K x n (data.frame)
                             # retour : matrice K x p
                             # ---------------------------
                             .compute_mode_var_distances = function(X_mat, modes) {
                               K <- nrow(modes)
                               p <- ncol(X_mat)
                               modes_mat <- as.matrix(modes)
                               
                               dist_mat <- matrix(0, nrow = K, ncol = p)
                               for (m in 1:K) {
                                 mode_vec <- as.vector(modes_mat[m, ])
                                 # X_mat != mode_vec -> matrice n x p, colSums -> vecteur p
                                 dist_mat[m, ] <- colSums(X_mat != mode_vec)
                               }
                               dist_mat
                             }
                           ),
                           
                           public = list(
                             
                             # ---------------------------
                             # FIT 
                             # ---------------------------
                             fit = function(X,
                                            active_vars      = NULL,
                                            descriptive_vars = NULL,
                                            auto_select_actives = TRUE,
                                            compute_mca      = TRUE,
                                            mca_ncp          = 5,
                                            verbose          = TRUE) {
                               
                               X <- as.data.frame(X)
                               
                               # Tout en caractères
                               X[] <- lapply(X, as.character)
                               
                               if (!all(sapply(X, is.character))) {
                                 stop("Toutes les variables doivent être qualitatives (ou convertibles en caractères).")
                               }
                               
                               private$.X <- X
                               all_vars <- colnames(X)
                               p_total  <- ncol(X)
                               n        <- nrow(X)
                               
                               # ---------------------------
                               # Gestion variables actives / descriptives
                               # ---------------------------
                               if (!is.null(active_vars)) {
                                 # L'utilisateur a explicitement choisi les actives
                                 active_vars <- intersect(active_vars, all_vars)
                                 if (length(active_vars) == 0L) {
                                   stop("Aucune variable active valide parmi les noms fournis.")
                                 }
                                 
                                 if (!is.null(descriptive_vars)) {
                                   descriptive_vars <- intersect(descriptive_vars, all_vars)
                                 } else {
                                   descriptive_vars <- setdiff(all_vars, active_vars)
                                 }
                                 
                               } else if (auto_select_actives) {
                                 # Sélection automatique basée sur entropie
                                 sels <- private$.auto_select_active_vars(X)
                                 active_vars      <- sels$active
                                 descriptive_vars <- sels$descriptive
                                 
                                 if (verbose) {
                                   message("Variables actives sélectionnées automatiquement :")
                                   message(paste(active_vars, collapse = ", "))
                                 }
                                 
                               } else {
                                 # Si l'utilisateur ne veut pas de sélection auto et ne fournit pas d'actives,
                                 # on considère que toutes les variables sont actives.
                                 active_vars      <- all_vars
                                 descriptive_vars <- character(0)
                               }
                               
                               private$.active_vars      <- active_vars
                               private$.descriptive_vars <- descriptive_vars
                               
                               X_active <- X[active_vars]
                               private$.X_active <- X_active
                               private$.X_descriptive <- if (length(descriptive_vars) > 0) X[descriptive_vars] else NULL
                               
                               p <- ncol(X_active)
                               if (p < 2L) {
                                 stop("Il faut au moins 2 variables actives pour faire un K-Modes de variables.")
                               }
                               
                               # ---------------------------
                               # Déterminer K
                               # ---------------------------
                               K <- self$get_n_cluster()
                               if (is.null(K)) {
                                 K <- private$.auto_select_K(X_active)
                                 private$.K_auto <- TRUE
                                 if (verbose) message("K sélectionné automatiquement : ", K)
                                 self$set_n_cluster(K)
                               } else {
                                 private$.K_auto <- FALSE
                               }
                               
                               # ---------------------------
                               # Initialisation
                               # ---------------------------
                               private$.modes <- private$.cao_init(X_active, K)
                               
                               # Pour accélérer les calculs : matrice n x p
                               X_mat <- as.matrix(X_active)
                               
                               cost_history <- numeric(0)
                               clusters <- NULL
                               
                               # ---------------------------
                               # Boucle K-Modes (variables actives)
                               # ---------------------------
                               for (iter in 1:private$.max_iter) {
                                 
                                 # Distances modes <-> variables (K x p)
                                 dist_mat <- private$.compute_mode_var_distances(X_mat, private$.modes)
                                 
                                 # Assignation : pour chaque variable (colonne), on choisit le mode le plus proche
                                 clusters <- apply(dist_mat, 2, which.min)
                                 names(clusters) <- colnames(X_active)
                                 
                                 # Coût total et coûts par cluster
                                 total_cost <- sum(dist_mat[cbind(clusters, 1:p)])
                                 cost_history <- c(cost_history, total_cost)
                                 
                                 cluster_costs <- tapply(dist_mat[cbind(clusters, 1:p)],
                                                         INDEX = clusters,
                                                         FUN = sum)
                                 
                                 # Mise à jour des modes
                                 new_modes <- private$.update_modes(X_active, clusters, K)
                                 
                                 if (verbose) {
                                   message(sprintf("Itération %d - Coût total = %d", iter, total_cost))
                                 }
                                 
                                 # Convergence ?
                                 if (all(new_modes == private$.modes)) {
                                   if (verbose) message("Convergence atteinte en ", iter, " itérations.")
                                   
                                   private$.clusters      <- clusters
                                   private$.cost          <- total_cost
                                   private$.cost_history  <- cost_history
                                   private$.cluster_costs <- cluster_costs
                                   private$.n_iter        <- iter
                                   
                                   break
                                 }
                                 
                                 private$.modes <- new_modes
                                 private$.clusters <- clusters
                                 private$.cost          <- total_cost
                                 private$.cost_history  <- cost_history
                                 private$.cluster_costs <- cluster_costs
                                 private$.n_iter        <- iter
                                 
                                 if (iter == private$.max_iter && verbose) {
                                   warning("Max_iter atteint sans convergence.")
                                 }
                               }
                               
                               # ---------------------------
                               # ACM (MCA) pour représentation graphique
                               # ---------------------------
                               if (compute_mca) {
                                 if (requireNamespace("FactoMineR", quietly = TRUE)) {
                                   quali.sup.idx <- if (length(descriptive_vars) > 0) {
                                     match(descriptive_vars, all_vars)
                                   } else {
                                     NULL
                                   }
                                   
                                   mca_res <- FactoMineR::MCA(
                                     X,
                                     quali.sup = quali.sup.idx,
                                     ncp       = mca_ncp,
                                     graph     = FALSE
                                   )
                                   private$.mca_result <- mca_res
                                 } else if (verbose) {
                                   warning("FactoMineR n'est pas disponible. L'ACM n'a pas été calculée.")
                                 }
                               }
                               
                               invisible(self)
                             },
                             
                             # ---------------------------
                             # PREDICT: assigner nouvelles variables (ex : variables descriptives)
                             # newvars : data.frame (colonnes = variables à clusteriser)
                             # ---------------------------
                             predict = function(newvars) {
                               
                               if (is.null(private$.modes)) {
                                 stop("Le modèle doit être ajusté avec fit() avant d'utiliser predict().")
                               }
                               
                               newvars <- as.data.frame(newvars)
                               newvars[] <- lapply(newvars, as.character)
                               
                               X_ref <- private$.X_active
                               if (is.null(X_ref)) {
                                 stop("Les données d'origine ne sont plus disponibles pour vérifier la compatibilité.")
                               }
                               
                               # Vérif : même nombre d'individus
                               if (nrow(newvars) != nrow(X_ref)) {
                                 stop("newvars doit avoir le même nombre d'individus (lignes) que les données utilisées dans fit().")
                               }
                               
                               X_new <- as.matrix(newvars)
                               K <- self$get_n_cluster()
                               
                               # Distances modes <-> nouvelles variables (K x p_new)
                               dist_mat <- private$.compute_mode_var_distances(X_new, private$.modes)
                               clusters_new <- apply(dist_mat, 2, which.min)
                               names(clusters_new) <- colnames(newvars)
                               
                               return(clusters_new)
                             },
                             
                             # ---------------------------
                             # PRINT
                             # ---------------------------
                             print = function(...) {
                               cat("KModesVariables Object\n")
                               cat("----------------------\n")
                               cat("Nombre de clusters (K) :", self$get_n_cluster(), "\n")
                               cat("Variables totales     :", if (!is.null(private$.X)) ncol(private$.X) else NA, "\n")
                               cat("Variables actives     :", length(private$.active_vars), "\n")
                               cat("Variables descriptives:", length(private$.descriptive_vars), "\n")
                               if (!is.null(private$.clusters)) {
                                 cat("Variables actives clusterisées :", length(private$.clusters), "\n")
                               } else {
                                 cat("Modèle non encore ajusté.\n")
                               }
                               invisible(self)
                             },
                             
                             # ---------------------------
                             # SUMMARY
                             # ---------------------------
                             summary = function(...) {
                               
                               cat("=== Résumé KModesVariables ===\n\n")
                               
                               K <- self$get_n_cluster()
                               cat("Nombre de clusters (K)        :", K, "\n")
                               cat("K sélectionné automatiquement :", if (isTRUE(private$.K_auto)) "Oui" else "Non", "\n\n")
                               
                               if (is.null(private$.X)) {
                                 cat("Pas de données associées.\n")
                                 return(invisible(self))
                               }
                               
                               n <- nrow(private$.X)
                               p_tot <- ncol(private$.X)
                               
                               cat("Nombre d'individus (n)        :", n, "\n")
                               cat("Nombre de variables total (p) :", p_tot, "\n")
                               cat("Variables actives             :", paste(private$.active_vars, collapse = ", "), "\n")
                               if (length(private$.descriptive_vars) > 0) {
                                 cat("Variables descriptives        :", paste(private$.descriptive_vars, collapse = ", "), "\n")
                               } else {
                                 cat("Variables descriptives        : aucune\n")
                               }
                               cat("\n")
                               
                               if (is.null(private$.clusters)) {
                                 cat("Modèle non encore ajusté.\n")
                                 return(invisible(self))
                               }
                               
                               cat("Nombre d'itérations effectuées:", private$.n_iter, "\n")
                               cat("Coût total final              :", private$.cost, "\n")
                               if (!is.null(private$.cost_history)) {
                                 cat("Coût (1ère itération)         :", private$.cost_history[1], "\n")
                                 cat("Coût (dernière itération)     :", tail(private$.cost_history, 1), "\n")
                               }
                               cat("\n")
                               
                               # Résumé par cluster
                               cl <- private$.clusters
                               tab_cl <- table(cl)
                               cat("=== Résumé des clusters (variables actives) ===\n")
                               res_cluster <- data.frame(
                                 Cluster         = as.integer(names(tab_cl)),
                                 Taille          = as.integer(tab_cl),
                                 Proportion      = as.numeric(tab_cl) / length(cl),
                                 Cout_Cluster    = if (!is.null(private$.cluster_costs)) {
                                   private$.cluster_costs[as.character(names(tab_cl))]
                                 } else NA_real_
                               )
                               print(res_cluster, row.names = FALSE)
                               cat("\n")
                               
                               # Liste des variables par cluster
                               cat("=== Variables par cluster ===\n")
                               split_vars <- split(names(cl), cl)
                               for (k in sort(unique(cl))) {
                                 cat("Cluster", k, ":\n")
                                 cat("  ", paste(split_vars[[as.character(k)]], collapse = ", "), "\n")
                               }
                               cat("\n")
                               
                               # Résumé rapide de l'ACM si disponible
                               if (!is.null(private$.mca_result)) {
                                 eig <- private$.mca_result$eig
                                 cat("=== ACM (MCA) ===\n")
                                 cat("Inertie expliquée (axes 1 & 2) :\n")
                                 cat(sprintf(" Axe 1 : %.2f %%\n", eig[1, 2]))
                                 cat(sprintf(" Axe 2 : %.2f %%\n", eig[2, 2]))
                                 cat("\n")
                                 cat("Les coordonnées des variables et des clusters\n")
                                 cat("peuvent être utilisées pour des graphiques dans l'application.\n")
                               } else {
                                 cat("ACM non calculée (ou FactoMineR non disponible).\n")
                               }
                               
                               invisible(self)
                             },
                             
                             # ---------------------------
                             # Getters
                             # ---------------------------
                             get_modes            = function() private$.modes,
                             get_clusters         = function() private$.clusters,              # sur variables actives
                             get_cost             = function() private$.cost,
                             get_cost_history     = function() private$.cost_history,
                             get_cluster_costs    = function() private$.cluster_costs,
                             get_active_vars      = function() private$.active_vars,
                             get_descriptive_vars = function() private$.descriptive_vars,
                             get_mca              = function() private$.mca_result
                           )
)

