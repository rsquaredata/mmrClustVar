# mmrClustVar — FAQ

Inventaire de **questions susceptibles d'être posées** concernant :

- Théorie du clustering de variables
- Implémentation dans mmrClustVar
- Architecture R6 + moteurs spécialisés
- Shiny app
- Choix de design, arbitrages, limites
- Étude de cas *metal_universe*
- Perspectives d'amélioration

Chaque question contient :
- **Réponse synthétique**
- **Réponse longue** (en spoiler)

---

# 1. Questions générales sur le projet

---

## Pourquoi faire du clustering de *variables* ?
- Pour identifier des groupes de variables redondantes ou fortement liées.
- Pour réduire la dimension, simplifier l'interprétation des donnés, structurer le dataset en blocs thématiques.
- Équivalent d'une analyse de variables, complémentaire du clustering d'individus.
<details>
  <summary></summary>
Le clustering de variables permet de travailler sur la structure interne du dataset, non pas au niveau des individus mais au niveau des *dimensions*.

L'idée est de détecter des **blocs de variables cohérents**, càd des variables qui :
- sont fortement **corrélées** entre elles (bloc numérique),
- ou qui représentent des **profils de modalités** (co-occurrences) comparables (bloc qualitatif),
- ou un mélange des deux pour (bloc mixte).

Cela répond à plusieurs besoins méthodologiques :
- **réduction de dimension interprétable** : chaque cluster peut être résumé par un prototype (ACP locale, mode, médoïde).
- **détection de redondance** : les variables quasi dupliquées ou dérivées se retrouvent ensemble.
- **amélioration des modèles downstream** : éviter les multicolinéarité, structurer l'information, créer des méta-variables.
- **structuration thématique** : on découvre des "blocs logiques" dans les données, ce que l'ACP globale ne montre pas toujours.

C'est donc une approche complémentaire à l'ACP et au clustering d'individus pour qui permet une analyse exploratoire avancée.

Là où le clustering d'individus répond à : "quels individus se ressemblent ?", le clustering de varaibles répond à : "quelles dimensions racontent la même histoire ?".

Cette question est essentielle en data science, car bien souvent, **les modèles ne souffrent pas d'un manque de données... mais d'un excès de variables redondantes".
</details>

---

# 2. Questions sur les algorithmes implémentés

---

## Pourquoi avoir réimplémenté k-means, k-modes, k-prototypes et k-medoids ?

Parce qu'ils couvrent les trois types de données : numériques, qualitatives, mixtes.  
Et ils respectent le cahier des charges : techniques de réallocation + traitement catégoriel.

<details>
  <summary></summary>
Le cahier des charges imposait :
- **au moins 3 méthodes**,
- dont **au moins une à réallocation**,
- et **au moins une pour les variables qualitatives**.

Nous avons donc implémenté :

| Algorithme       | Type de données | Technique | Raisons du choix |
|------------------|------------------|-----------|------------------|
| k-means          | numérique        | réallocation | standard du clustering num |
| k-modes          | qualitatif       | réallocation | dédié aux profils modaux |
| k-prototypes     | mixte            | réallocation | méthode hybride num+cat |
| k-medoids        | tout type        | dissimilarités | robuste, interprétable |

Cela va **au-delà du sujet**, qui n'en exigeait que 3.
</details>
---

## Pourquoi utiliser l'ACP locale dans k-means (principe de ClustOfVar) ?

Parce qu'il faut un prototype qui résume un **cluster de variables**, et la 1e composante principale locale maximise la variance expliquée.

<details>
  <summary></summary>
Dans ClustOfVar et dans la littérature, un cluster de variables numériques se résume par :
- la **1e composante principale locale**,
- calculée uniquement sur les variables appartenant au cluster,
- elle maximise simultanément :
  - la variance capturée dans ce cluster,
  - la corrélation moyenne avec les variables qu'elle représente.

Donc, utiliser l'ACP locale comme prototype est :
- mathématiquement justifié,
- cohérent avec ClustOfVar,
- interprétable,
- efficace en pratique.
</details>

**variance** : mesure à quel point les valeurs d’une variable sont dispersées autour de leur moyenne → $Var(X) = \frac{1}{n} \sum_{i=1}^n (x_i - \bar{x}^2)$ où $x_i$ = valeurs individuelles, $\bar{x}$ = moyenne, $(x_i - \bar{x}^2$ = écart à la moyenne.

---

## Pourquoi ne pas utiliser directement une ACP globale ?

Une ACP globale ne capture pas des structures locales.

<details> 
<summary></summary>
Une ACP globale :
- impose des directions uniques valables pour tout le jeu de données,
- ne permet pas d'identifier des blocs distincts,
- mélange les thématiques.

Au contraire :
- une ACP locale capture une structure interne propre à chaque cluster,
- et définit un prototype spécifique → beaucoup plus précis.
- </details>

---

## Pourquoi avoir choisi la dissimilarité *simple matching* pour k-modes ?

Parce que c'est la dissimilarité standard pour ce type d'algorithme : elle compte le nombre de désaccords.

<details>
  <summary></summary>
Le simple matching mesure :
$$
d(x_j, x_{\ell} = \text{nombre de positions où les modalités diffèrent}.
$$

Ses avantages :
- cohérent avec le mode comme prototype,
- extrêmement simple,
- interprétable,
- rapide à calculer,
- largement utilisé dans la littérature (Huang, 1997).
</details>

---

## Pourquoi la distance mixte de k-prototypes est-elle une somme pondérée ?

Il faut pouvoir combiner deux espaces de nature différente : numérique et catégoriel.

<details>
  <summary></summary>
La formule :
$$
d_{tot} = d_{num} + \lambda d_{cat}
$$

permet :
- d’équilibrer les contributions,
- d’éviter que l’un des deux blocs domine,
- de régler la sensibilité de la méthode suivant les données.

$\lambda$ est un vrai hyperparamètre → nous avons laissé l’utilisateur en contrôle, ce qui était cohérent avec le sujet.
</details>

---

## Pourquoi avoir ajouté k-medoids alors qu’il n’était pas demandé ?

Parce qu’il est complémentaire, robuste, et très pédagogique pour illustrer l’approche dissimilarité.

<details>
  <summary></summary>
k-medoids :
- fonctionne sur *toutes* dissimilarités,
- produit des prototypes interprétables (médoïdes),
- plus robuste que k-means (moins sensibles aux outliers),
- permet une comparaison instructive avec  k-means et k-prototypes.

C'est une valeur ajoutée du projet.
</details>

---

# 3. Questions sur la stratégie de choix de K

---

## Pourquoi avoir choisi la courbe d’inertie (méthode du coude) ?

C’est l’outil le plus simple. Mais la méthode du coude est un point de départ, nous avons implémenté d'autres outils : adhésion moyenne, part d'inertie expliquée, indicateurs cohésion/séparation, et des visualisations d'interprétation.

Le choix de K ne repose jamais sur un seul critère, mais sur une convergence de signaux.

<details>
  <summary></summary>
La courbe d’inertie intra-cluster :
- est facile à lire,
- fonctionne pour tous les algorithmes à réallocation,
- est cohérente avec ClustOfVar,
- fait partie du programme du cours (Ricco en parle).

Elle est complétée par :
- l’adhésion moyenne,
- les inerties par cluster,
- les plots d’interprétation → le choix ne repose pas sur un seul indicateur.
</details>

## Quelles sont les autres stratégies possibles ?

Les stratégies pertinentes pour le choix du nombre de clusters K en clustering de variables sont :

- **adhésion moyenne** (degree of memebership) : plus l'adhésion moyenne est éléve, plus la partition moyenne est cohérente. Une chute buutale indique un mauvais choix de K.
- **part d'inertie expliquée** : on maximise ce critère pour trouver le K optimal.
- **stabilité des partitions** : si de pertites perturbations donnent des partitions radicalement différentes, K n'est pas bon.
- **comparaison de l'inertie inter vs intra cluster**.
- **analyse des prototypes (ACP locale/profil modal)** : quand les prototypes deviennent difficiles à interpréter ou trop proches, c'est que K est trop grand.
- **mesure de séparation géométrique (silhouette adaptée)** (méthode avancée) : on peut utiliser une silhouette modifiée pour les variables.
- **critères basés sur l'information** : on peut adapter AIC, BIC ou ICL, mais ça demande un modèle probabiliste.
- **visualisation des redondances et blocs** : les heatmaps de corrélation ou de dissimilarité peuvent suggérer un nombre naturel de blocs.
- **approche subjective : lisibilité + utilité** : le meilleur K est aussi celui qui a du sens pour l'analyse.

<sub>**Modèle probabiliste** : dans un clustering de variables, c'est un modèle qui dit que chaque variable est générée par un *mélange de distributions* correspondant aux clusters. Exemples : mélange de gaussiennes, mixture multinomiale, modèle latent pour catégories.  
`mmrClustVar` n'utilise pas de modèle probabiliste mais des méthodes algorithmiques (k-means-mile).

**Modèle latent** : modèle où les clusters ne sont pas directement observés, ils sont des *variables cachées* qui génèrent les données.  
Exemples : LCA (Latent Class Analysis), mixtures en EM, modèles factoriels avec classes.<su</sub>

<details>
  <summary></summary>

- l'**adhésion moyenne** mesure comment une variable est représentée par le prototype du cluster. Selon la méthode :
  - numérique → $r^2$
  - catégoriel → taux d'accord modal
  - mixte → somme pondérée
  - k-medoids : prototype = médoïde → $adhesion(x_j) = 1 - \frac{D(x_j, medoid_k)}{\max_{\ell}D(x_j, medoid_{\ell})}$
Graphiquemebt, on trace l'adhésion moyenne en fonction de K → on cherche un plateau.

- **part d'inertie expliquée** :
  - $\text{part_exp}(K) = 1 - \frac{\mathbf{I}_{intra}(K)}{\mathbf{I}_{tot}}$.
  - On choisit un K qui 1. explique suffisamment de variance, 2. sans explosion du nombre de K.

- **stabilité des partitions** :
  - on peut :
    - bootstrapper les variables (plusieurs runs aléatoires avec remise),
    - refaire le clustering plusieurs fois,
    - mesurer la similarité entre partitions (ARI, Rand Index, Jaccard):
      - Rand Index : proportion de paires de variables classées de la même façon dans les 2 partitions ($\in \[0,1\]$) ; $ARI = \frac{RI - \mathbb{E}(RI)}{1 - \mathbb{E}(RI)}$
      - ARI (adjusted rand Index) : rand index corrigé qui tient compte de l'aléatoire ($\in \[-\infty, 1\]$) : 1 = clusters identiques ; 0 = ce qu'on attend par hasard ; negatif = pure que hasard.
      - Jaccard : mesure de similarité plus stricte qui ne prend en compte que les paires réellement dans un même cluster.
  - un bon K = partition stable.

- **comparaison inertie inter vs intra** : l'objectif est de maximiser $\frac{inter}{intra}$

- **analyse des prototypes** :
  - l'ACP locale n'explique pas grand chose → cluster mal défini
  - modes trop proches d'un cluster à l'autre
  - prototype redondant → K trop grand

- **silhouette adaptée** :
  - silhouette classique = distance entre individus.
  - pour les **variables**, on remplace :
    - distance variable-cluster → $1-r^2$ ou simple matching.
    - distance variable-autre cluster → même mesure
  - calcul d'une silhouette moyenne pour chaque K.

- **critères basés sur l'information** : rarement utilisé pour le clusterinf de variables **non probabilistes** mais possible vi a
  - modèles de mélange sur les scores de prototypes
  - modèles latents
  
<su>bScore des prototype : à quel point une variable ressemble à son prototype, par exemple dans `mmrClustVar`
- k-means : prototype = composante principale **locale** ; score = $r^2$
- k-modes : prototype = profil modal ; score = taux d'accord
- k-prototypes : score = somme pondérée $r^2 + \lambda \times \text{simple matching}$
- k-medoids : prototype = variable réelle ; score = dissimilarité $\mathcal{D}$</sub>

- **visualisation des redondances et blocs** : on regarde
  - les "blocs clairs" dans la matrice $r^2$
  - les zones homogènes dans la matrice simple matching
  - les dendrogrammes en aval (non hiérarchiques mais indicatifs)

- **approche subjective** :
  - on équilibre bonne inertie / bonne adhésion / bonne interprétabilité / cohérence avec le domaine
  - approche terrain des statisticiens

</details>

---

# 4. Questions sur l'architecture R6

---

## Pourquoi avoir utilisé une architecture R6 plutôt qu’un ensemble de fonctions ?

Parce qu’un package de clustering moderne doit être orienté objet : un modèle = un objet avec des méthodes.

<details>
  <summary></summary>
R6 permet :
- encapsulation des données (X actif, clusters, prototypes),
- méthodes propres : `fit`, `predict`, `summary`, `plot`,
- moteurs spécialisés interchangeables,
- façade `Interface` uniforme → simplicité pour l’utilisateur,
- extensibilité → ajouter un algorithme ne casse rien.

C'est l’architecture utilisée dans **tidymodels**, dans **mlr3** et dans **scikit-learn** (philosophie identique).
</details>

---

## Pourquoi avoir fait une classe "façade" + engines spécialisés ?

Pour séparer la logique utilisateur (interface) et la logique algorithmique (moteurs).

<details>
  <summary></summary>
La façade :
- uniformise l’interface,
- gère le mode "auto",
- centralise `summary`, `plot`, `interpretation`.

Les moteurs :
- contiennent la logique algorithmique pure,
- définissent leurs propres prototypes,
- implémentent les hooks (`run_clustering`, `predict_one_variable`, etc.).

→ Cela rend le package modulaire, et extensible.
</details>

---

# 5. Questions sur Shiny

---

## Quelles sont les fonctionnalités principales de votre application Shiny ?

Importer un dataset, choisir les variables, lancer le clustering, afficher les résultats, les diagnostics, et exporter les clusters.

<details>
  <summary></summary>>
- Import CSV/XLSX  
- Choix variables actives / illustratives  
- Choix de la méthode (ou mode auto)  
- Lancement du clustering (via Interface)  
- Résultats :
  - clusters
  - prototypes
  - inerties
  - adhésion
- Graphiques :
  - inertie (K-path)
  - heatmap
  - barplots d’adhésion
  - dendrogramme
- Export ZIP :
  - clusters CSV
  - résumé
  - diagnostics

C’est **plus complet que l’exigence minimale du sujet**.
</details>

---

# 6. Étude de cas *metal_universe*

---

## Pourquoi avoir créé un dataset thématique fictif ?

Pour avoir un terrain de jeu contrôlé où la structure de redondance est connue à l’avance.

<details>
  <summary></summary>
*metal_universe* permet :
- d’avoir des corrélations maîtrisées,
- d’avoir des sous-genres cohérents,
- de tester la sensibilité aux variables de bruit,
- de comparer les algorithmes sur un cas réaliste,
- d’avoir un dataset "signature" qui va au-delà du projet standard.
</details>

---

## Comment avez-vous généré les variables numériques fictives ?

À partir de plages plausibles définies par sous-genre + ajout de bruit + création de variables dérivées.

<details>
  <summary></summary>
Étapes :
1. sous-genres crédibles  
2. plages numériques plausibles par style  
3. tirages dans ces plages (uniforme / normale tronquée)  
4. ajout de bruit pur  
5. création de variables dérivées linéaires  
6. recodage catégoriel secondaire

C’est exactement ce qu’on attend d'un dataset "didactique".
</details>

---

## Est-ce que les algorithmes récupèrent vraiment les blocs simulés ?

Oui : les blocs rythmiques, violents et socio-géographiques ressortent clairement.

<details>
  <summary></summary>
k-prototypes et k-medoids identifient :
- un bloc rythmique (BPM, speed, tempo variability)
- un bloc "violence du son"
- un bloc géographique / sociologique
- un cluster mélangeant modérément les deux
- un cluster pour le bruit explicite

→ donc la structure conçue est bien retrouvée.
</details>

---

# 7. Questions sur les limites et perspectives

---

## Quelles sont les limites principales du package ?

ACP locale, tuning de $\lambda$ manuel, coût de k-medoids, pas encore de méthode hiérarchique.

<details>
  <summary></summary>
- ACP locale = modèle linéaire  
- $\lambda$ fixe = tuning manuel  
- dissimilarité complète en O(p²) pour k-medoids  
- pas de vrai clustering hiérarchique de variables  
- initialisation perfectible
- pas de sélection automatique de $\lambda$

<sub>Concernant l'initialisation :
- actuellement, l'initialisation est heuristique
```r
init <- sample(1:K, p, replace = TRUE)
```
- meilleures heuristiques :
  - **k-means++ adapté aux variables** :
    - idée : choisir des prototypes initialement très différents
    - premier prototype = variable au hasard → prototypes suivants = variables éloignées des prototypes déjà choisis
    - améliore la stabilité des clusters
  - **initialisation "PCA seeds"** : on fait une ACP sur toutes les variables, on utilise les coordonnées comme cecteurs, et on regroupe en $K$ avec k-means, juste pour obtenur des seeds.
  - **initialisation par clustering hiérarchique** : hclust sur les corrélations entre variables, puis couper l'arbre en *k* → initialisation structurée
  - **modes fréquents pour k-modes** : prendre les variables les plus centrales comme prototypes initiaux

`hclust()` : fonction de R qui fait du **clustering hiérarchique**.  
L  fonction construit une **arbre hiérarchique** (dendrogramme) en fusionnant progressivement les objets les plus proches selon une méthode de liaison (ward, complete, average, single, etc.).  
Usage typique :
```r
d <- dist(X)
hc <- hclust(d, method = "ward.D2")
plot(hc)
```

**variable centrale** : variable qui a la modalité la plus fréquente pour la majorité de ses modalités, càd qui a **faible entropie / forte fréquence modale**.

Pour l'initalisation automatique de $\lambda$ :
- pas de réponse officielle mais plusieurs stratégies existent
- **équilibrage des variances** : $\lambda = \frac{\text{variance moyenne du bloc numérique}}{\text{dissimilarité moyenne du bloc catégoriel}}$
- **cross validation non supervisée** : on teste plusieurs $\lambda = 0.1, 0.5, 1, 2, 5, \ldots$, on clacule l'inertie intra, et on choisit le $\lambda$ qui minimise $\mathbf{I}(K, \lambda)$
- **stabilité des partitions** : $\lambda$ optimal = celui qui donne la partition la plus stable (ARI maximal en bootstrap)
- **équilibre adhésion numérique / adhésion catégorielle** : $\lambda$ tel que $adhesion_{num} = adhesion_{cat}$</sub>

</details>

---

## Quelles améliorations seraient prioritaires ?

Automatisation de $\lambda$, nouvelles similarités, méthodes hiérarchiques.

<details>
  <summary></summary>
1. calibration automatique de λ  
2. k-means++ ou spectral init  
3. corrélations non linéaires (distance de distance corrélationnelle)  
4. clustering hiérarchique de variables entièrement implémenté  
5. version étendue du dataset metal avec web scraping  
</details>

---

# 8. Questions "Pourquoi pas… ?"

---

## Pourquoi ne pas avoir implémenté une méthode hiérarchique ?

Parce que ce n’était pas obligatoire et aurait doublé la charge algorithmique.

<details>
  <summary></summary>
Un vrai clustering hiérarchique de variables nécessite de :
- définir une dissimilarité entre clusters,
- recalculer un prototype,
- garantir une cohérence avec les données mixtes.

La littérature sur les méthodes hiérarchiques mixtes est plus complexe → hors périmètre du projet, mais prévu dans les perspectives.
</details>

---

## Pourquoi ne pas avoir utilisé directement le package *ClustOfVar* ?

Pour réimplémenter nos propres algorithmes.

<details>
  <summary></summary>
ClustOfVar est une référence, mais :
- le sujet impose d’arbitrer ce qu’on réimplémente vs ce qu’on délègue,
- notre but était pédagogique,
- nos moteurs sont modulaires,
- l’interface R6 est plus moderne que les fonctions de ClustOfVar.
</details>

---

# 9. Questions de validation du cahier des charges

---

## Comment prouvez-vous que votre package respecte le cahier des charges ?

Nous dépassons tous les critères : méthodes, R6, Shiny, auto-K, documentation.

<details>
  <summary></summary>
- 4 algos (3 demandés)  
- plusieurs techniques de réallocation (1 demandé)
- 2 algos gèrent qualitatif (1 demandé)
- stratégie de choix du K (plusieurs proposées)  
- classe R6 complète  
- application Shiny complète  
- documentation anglaise 
- datasets internes dont un dataset thématique original
</details>

---

# 10. Questions de type « faille logique »

---

## Que se passe-t-il si une variable est constante ?

Elle est soit rejetée, soit son adhésion est 0 — c’est géré dans la base class.

<details>
  <summary></summary>
Dans `ClusterrBase` :
- détection des variances nulles,
- exclusion automatique ou avertissement,
- inertie/adhésion fixées de manière cohérente.
</details>

---

## Que se passe-t-il si toutes les variables sont catégorielles ou toutes numériques ?

Le moteur approprié est sélectionné automatiquement.

<details>
  <summary></summary>
En mode auto : 
- tout numérique → k-means  
- tout qualitatif → k-modes  
- mixte → k-prototypes
En mode manuel :
- message d'erreur si l'utilisateur choisit un moteur inadapté.

k-medoids reste accessible manuellement.
</details>
---

# 11. Questions sur les résultats

---

## Comment interpréter l’adhésion ?

C’est un score qui mesure à quel point une variable est proche de son cluster.

<details>
  <sumamry></summary>
- numérique : basé sur $r^2$ avec le prototype  
- catégoriel : basé sur le taux d'accord modal  
- mixte : combinaison pondérée  

Plus l’adhésion est haute, plus la variable "appartient vraiment" à son groupe.
</details>

---

## Que représente la part d’inertie expliquée ?

La qualité globale de la partition.

<details>
  <summary></summary>
Part d'inertie explique : $part_exp = 1 - \frac{\mathbf{I}_{intra}}{\mathbf{I}_{totale}}$

- plus elle est élevée, meilleure est la partition  
- permet de comparer K  
- score présent dans summary et plots  
</details>

---

# 12. Questions sur le mode auto

---

## Comment fonctionne le mode "auto” ?

Il détecte le type de variables et choisit l’algorithme adapté.

<details>
  <summary></summary>
Pipeline :
1. analyse du type de colonnes (num, cat, mixte)  
2. sélection de :
   - k-means pour numérique  
   - k-modes pour qualitatif  
   - k-prototypes pour mixte  
3. paramètres passés à l’engine  
4. appel du clustering  

C’est transparent pour l’utilisateur.
</details>

---

# 13. Bonus (questions difficiles)

---

## Comment comparer vos résultats à ceux d’autres packages ?

Comparer :
- inerties
- stabilité (bootstrap)
- silhouettes adaptées aux variables
- visualisations

---

## Complexité temporelle des algorithmes ?

- k-means / k-modes / k-prototypes : O(iter × K × p)  
- k-medoids : O(p²) + réallocations

<details>
  <summary></summary>
$O(\cdot)$ : complexité algorithmique, càd vitess à laquelle le temps de calcul augment quand
- $p$ = nombre de variables
- $K$ = nombre de clusters
- $iter$ = nombre d'itération

$O(iter \times K \times p)$ : pour chaque variable, à chaque itération, on compare aux $K$ prototypes

$O(p^2)$ : on calcule une **matrice complète de dissimilarité** entre **toutes les $p$ variables**. Plus $p$ est grand, plus la complexité explose.

</details>  

---

## Et si l’utilisateur met λ = 0 ?

Le bloc catégoriel est ignoré → on retombe sur un k-means local sur la partie numérique.

---