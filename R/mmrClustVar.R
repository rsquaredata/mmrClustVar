#' Clustering de variables (façade R6)
#'
#' Classe R6 "façade" qui encapsule les différentes implémentations
#' de clustering de variables (k-means, k-modes, k-prototypes, k-medoids)
#' et fournit une interface unique pour l'utilisateur.
#'
#' L'objet se crée avec \code{mmrClustVar$new(...)} puis expose
#' les méthodes principales \code{$fit()}, \code{$predict()},
#' \code{$print()}, \code{$summary()} et \code{$plot()}.
#'
#' @docType class
#' @name mmrClustVar
#' @export
#'
#' @section Méthodes principales :
#' \describe{
#'   \item{\code{$new(method = c("kmeans", "kmodes", "kprototypes",
#'                              "kmedoids", "auto"),
#'                  K,
#'                  scale = TRUE,
#'                  lambda = 1, ...)}}{
#'     Constructeur. Initialise la méthode, le nombre de clusters,
#'     l'option de standardisation et le paramètre \code{lambda}.
#'   }
#'
#'   \item{\code{$fit(X)}}{
#'     Apprend un modèle de clustering de variables à partir des
#'     variables actives fournies dans \code{X}.
#'   }
#'
#'   \item{\code{$predict(X_new)}}{
#'     Rattache des variables supplémentaires aux clusters existants.
#'   }
#'
#'   \item{\code{$print()}}{
#'     Affiche un résumé succinct du modèle (méthode, K, nb de variables,
#'     inertie, convergence).
#'   }
#'
#'   \item{\code{$summary()}}{
#'     Affiche un résumé détaillé (tailles de clusters, inerties,
#'     indicateurs d'adhésion, etc.).
#'   }
#'
#'   \item{\code{$plot(type = c("inertia", "clusters", "membership",
#'                              "profiles"))}}{
#'     Génère des graphiques pour l'analyse du modèle.
#'   }
#' }
NULL

mmrClustVar <- R6::R6Class(
    "mmrClustVar",
    
    public = list(
        
        initialize = function(
        method = c("kmeans", "kmodes", "kprototypes", "kmedoids", "auto"),
        K,
        scale  = TRUE,
        lambda = 1,
        ...
        ) {
            method <- match.arg(method)
            private$method  <- method
            private$K       <- K
            private$scale   <- scale
            private$lambda  <- lambda
            private$engine  <- NULL  # sera créé à la volée dans fit()
        },
        
        fit = function(X) {
            if (!is.data.frame(X)) {
                stop("[mmrClustVar] X doit être un data.frame.")
            }
            
            # si l'engine n'existe pas encore, ou method == "auto",
            # on choisit la classe fille ici
            if (is.null(private$engine) || identical(private$method, "auto")) {
                private$engine <- private$build_engine(X)
            }
            
            private$engine$fit(X)
            invisible(self)
        },
        
        predict = function(X_new) {
            if (is.null(private$engine)) {
                stop("[mmrClustVar] predict() appelé avant fit().")
            }
            private$engine$predict(X_new)
        },
        
        print = function(...) {
            cat("Classe 'mmrClustVar'\n")
            cat("  Méthode demandée :", private$method, "\n")
            cat("  K                :", private$K, "\n")
            cat("  scale            :", private$scale, "\n")
            cat("  lambda           :", private$lambda, "\n\n")
            
            if (!is.null(private$engine)) {
                cat("Résumé du moteur interne :\n")
                private$engine$print()
            } else {
                cat("(Aucun modèle appris pour l'instant : appelez $fit(X)).\n")
            }
            
            invisible(self)
        },
        
        summary = function(...) {
            if (is.null(private$engine)) {
                cat("[mmrClustVar] summary() : aucun modèle (fit() non appelé).\n")
                return(invisible(NULL))
            }
            private$engine$summary(...)
        },
        
        plot = function(type = c("inertia", "clusters", "membership", "profiles"), ...) {
            if (is.null(private$engine)) {
                stop("[mmrClustVar] plot() : aucun modèle (fit() non appelé).")
            }
            private$engine$plot(type = type, ...)
        },
        
        get_clusters = function() {
            if (is.null(private$engine)) return(NULL)
            private$engine$get_clusters()
        },
        
        get_centers = function() {
            if (is.null(private$engine)) return(NULL)
            private$engine$get_centers()
        },
        
        get_inertia = function() {
            if (is.null(private$engine)) return(NA_real_)
            private$engine$get_inertia()
        },
        
        get_convergence = function() {
            if (is.null(private$engine)) return(NA)
            private$engine$get_convergence()
        },
        
        get_method = function() {
            if (!is.null(private$engine)) {
                return(private$engine$get_method())
            }
            private$method
        },
        
        get_K = function() {
            private$K
        },
        
        get_lambda = function() {
            private$lambda
        }
    ),
    
    private = list(
        method = NULL,  # "kmeans", "kmodes", "kprototypes", "kmedoids" ou "auto"
        K      = NULL,
        scale  = NULL,
        lambda = NULL,
        
        engine = NULL,  # instance d'une classe fille mmrClustVar*
        
        build_engine = function(X) {
            method <- private$method
            K      <- private$K
            scale  <- private$scale
            lambda <- private$lambda
            
            # --- Cas "auto" : choix en fonction des types de colonnes ---
            if (identical(method, "auto")) {
                is_num <- vapply(X, is.numeric, logical(1L))
                is_cat <- (!is_num) & vapply(
                    X,
                    function(col) is.factor(col) || is.character(col),
                    logical(1L)
                )
                
                if (all(is_num)) {
                    method_effective <- "kmeans"
                } else if (all(is_cat)) {
                    method_effective <- "kmodes"
                } else {
                    method_effective <- "kprototypes"
                }
                
                cat("[mmrClustVar] method = 'auto' → méthode choisie :", method_effective, "\n")
            } else {
                method_effective <- method
            }
            
            # --- Instanciation de la classe fille correspondante ---
            engine <- switch(
                method_effective,
                "kmeans" = mmrClustVarKMeans$new(
                    K      = K,
                    scale  = scale,
                    lambda = lambda
                ),
                "kmodes" = mmrClustVarKModes$new(
                    K      = K,
                    scale  = scale,
                    lambda = lambda
                ),
                "kprototypes" = mmrClustVarKPrototypes$new(
                    K      = K,
                    scale  = scale,
                    lambda = lambda
                ),
                "kmedoids" = mmrClustVarKMedoids$new(
                    K      = K,
                    scale  = scale,
                    lambda = lambda
                ),
                stop("[mmrClustVar] Méthode non supportée : ", method_effective)
            )
            
            engine
        }
    )
)
