# mmrClustVar — FAQ

Inventory of **potential questions** regarding:

- Variable clustering therory
- Implementation in mmrClustVar
- R6 Architecture + specialized engines
- Shiny app
- Design choices, trade-offs, limitations
- *metal_universe* case study
- Future improvements

Eache question contains:
- **Short answer**
- **Long answer** (in spoiler)

---

# 1. General questions about the project

---

## Why perform variable clustering?
- To identify groups of redundant or strongly related variables.
- To reduce dimensionality and simplify interpretation.
- To structure the data into meaningful thematic blocks.
- It is the variable-oriented counterpart of individual clustering

<details><summary></summary>

Variable clustering focuses on the internal structure of the dataset at the *variable* level rather than on the *individual* level.

The goal is to detect **coherent blocks of varaibles*, that is, variables that:
- are strongly **correlated** with each other (numeric block),
- represent similar **modal profiles** (categorical blocks),
- or mix both types (mixed blocks).

This meets several methodological needs:
- **interpretable dimensionality reduction** via block prototypes (local PCA, modes, medoids);
- **redundance detection**: neear-duplicate or dervided variables naturally cluster together,
- **improved downstream modelling**: discovering multicollinearity, structuring information, creating meta-variables.
- **thematic structuring**:  discovering logical subspaces that global PCA cannot isolate.

While individual clustering answers *"whuch individuals look alike?"*, variable clustering answers *"which dimensions tell the same story?"*

This is essential becauese most models suffer not from *too few** data points, but from **to many redundant variables**.
</details>

---

# 2. Questions about the implemented algorithms

---

## Why reimplément k-means, k-modes, k-prototypes and k-medoids?

Because together they cover numerci, catargorical and mixed data, and meet the project requirements; relocation algorithms + support for categorical variaables.

<details><summary></summary>

The assignment required:
- **3 methods**,
- including **one relocation method**,
- and **one categorical method**.

We implemented:

| Algorithm       | Data type | Technique  | Rationale                   |
|--- -------------|-----------|------------|-----------------------------|
| k-means      | numeric     | relocation | standard numeric clustering |
| k-modes      | categorical | relocation | dedicated to modal profiles |
| k-prototypes | mixed       | relocation |  hybrid num+cat approach    |
| k-medoids    | any type    | dissimilarities | robust, interpretable  |

This goes **beyond** the requirement, which only asked for three.
</details>
---

## Why use local PCA in k-means (ClustOfVar principle)?

Because a cluster of variables needs a prototype summarising them, and the first local principal component maximises explained variance.

<details><summary></summary>
In ClustOfVar and in the literature, numeric variable clusters are summarised by:
- the **the first principal component**,
- computed only on variables belonging to the cluster,
- maximising both variance capture and average correlation with member variables.

Thus, using local PCA is:
- mathematically justified,
- consistent with ClustOfVar,
- interpretable,
- effective in practive.

**Variance** : dispersion around values around the means → $Var(X) = \frac{1}{n} \sum_{i=1}^n (x_i - \bar{x}^2)$ où $x_i$ = valeurs individuelles, $\bar{x}$ = moyenne, $(x_i - \bar{x}^2$ = écart à la moyenne.
</details>

---

## Wy not use a global PCA?

Because a global PCA cannot reveal local block structure.

<details><summary></summary>

A global PCA:
- imposes directions valid for the entire dataset,
- cannot isolate distinct blocks,
- mixes thematic structures.

Local PCA:
- captures cluster-specific internal structure,
- provides more precise prototypes.
- </details>

---

## Why use simple matching dissimilarity in k-modes?

Because it is the standard for this algorithm: it counts disagreements.

<details><summary></summary>

Simple matching measures:
$$
d(x_j, x_{\ell} = \text{number of positions where categories differ}.
$$

Advantages:
- consistent witht the modal prototype,
- simple and interpretable,
- fast to compute,
- widely used (Huang, 1997).
</details>

---

## Why is the mixed distance in k-prototypes a weigthed sum?

Because numeric and categorical spaces must be combined in a controlled way.

<details><summary></summary>
La formule :
$$
d_{tot} = d_{num} + \lambda d_{cat}
$$

This:
- balances contributions,
- prevents one block from dominating,
- lets users tune sensitivity.

$\lambda$ is a real hyperparameter → making it user-controllable is consistent with the topic.
</details>

---

## Wgy include k-medoids when it was not required?

Because it is complementary, robust, and pedagogically useful to illustrate dissimilarity-based clustering.

<details><summary></summary>
k-medoids:
- works with *any* dissimilarity,
- yields interpretable prototypes (actual variables),
- more robust than k-means,
- enables informative comparison with k-means et k-prototypes.

It is an added value.
</details>

---

# 3. Choosing the number of clusters K

---

## Why use the inertia curve?

Because it is simple and standard. But our implementation provides many additional tools: average adhesion, explained inertia, cohesion/separation indicators, interpretation plots.

The choice of K should never rely on a single criterion.

<details><summary></summary>

The intra-cluster inertia curve:
- is easy to read,
- works for relocation, algorithms,
- is consistent with ClustOfVar.

We complement it with:
- mean adhesion,
- per-cluster inertias,
- interpretation plots

→ K is chosen through convergence of multiple signals.
</details>

## What other strategies are relevant?

The most useful strrategies for variables clustering include:

- **meand adhesion** (degree of memebership)
- **explained inertia**
- **partitions stability**
- **inter vs intra inertia ratio**.
- **prototypes interpretation quality (ACP locale/profil modal)**
- **adapted silhouette scores**
- **information criteria (model based)**
- **visual inspection (correlation / dissimilarity heatmaps)**
- **subjective interpretability**

<sub>**Modèle probabiliste** : dans un clustering de variables, c'est un modèle qui dit que chaque variable est générée par un *mélange de distributions* correspondant aux clusters. Exemples : mélange de gaussiennes, mixture multinomiale, modèle latent pour catégories.  
`mmrClustVar` n'utilise pas de modèle probabiliste mais des méthodes algorithmiques (k-means-mile).  
**Modèle latent** : modèle où les clusters ne sont pas directement observés, ils sont des *variables cachées* qui génèrent les données.  
Exemples : LCA (Latent Class Analysis), mixtures en EM, modèles factoriels avec classes.</sub>

<details>
  <summary></summary>

- l'**mean adhesion** measyres how a variable is represented by the cluster prototype. Depending on the method:
  - numeric → $r^2$
  - categorical → proportion of mode similarity
  - mixed → weighted sum
  - k-medoids : prototype = medoid → $adhesion(x_j) = 1 - \frac{D(x_j, medoid_k)}{\max_{\ell}D(x_j, medoid_{\ell})}$
Visually, we plot the mean adhesion per K and look for a plateau.

- **explained inertia**:
  - $\text{exp_part}(K) = 1 - \frac{\mathbf{I}_{intra}(K)}{\mathbf{I}_{tot}}$.
  - we choose a value for K that 1. explains enough variance, 2. without an exploding number of K. 

- **clustering stability**:
  - available alternatives:
    - bootstraping the variables (several runs with random sampling andreplacement),
    - running several clusterings,
    - measuring the similarity between culsterings (ARI, Rand Index, Jaccard):
      - Rand Index: proportion of variable pairs with the same ordering in each clustering($\in \[0,1\]$) ; $ARI = \frac{RI - \mathbb{E}(RI)}{1 - \mathbb{E}(RI)}$
      - ARI (adjusted rand Index): corrected erand index corrigé which takes randomness into account ($\in \[-\infty, 1\]$) : 1 = identival clusters ; 0 = expected random result ; negatif = worse than randomness).
      - Jaccard: stricter similarity measure stricte that only takes into account pairs without are actually in the same cluster.
  - good K = stable clustering.

- **inter vs intra comparison**: the goal is to maximize $\frac{inter}{intra}$

- **prototype analysis**:
  - the local PCA provides little explanation → ill-defined cluster
  - modes that are too close from one one cluster to another
  - redundant prototype → K is too big

- **silhouette adaptation**:
  - usual silhouette for individuals = distance between individuals.
  - for **variables**, we replace:
    - variable-cluster distance → $1-r^2$ or simple matching.
    - variable-other cluster dustance → same measure.
  - mean silhouette computer for each K.

- **information-based criteria**: seldom used for **non probabilistic** variable clustering but can be employed with
  - models mixing prototype scores
  - latent models
  
<su>Prototype score: how much a variable resembles its prototype, for example in `mmrClustVar`
- k-means: prototype = local PC1 ; score = $r^2$
- k-modes: prototype = modal profil ; score = matching rate
- k-prototypes: score = weighted sum $r^2 + \lambda \times \text{simple matching}$
- k-medoids : prototype = real variable ; score = dissimilarity $\mathcal{D}$</sub>

- **visualisation of redundancy and blocks** : we look for
  - "clear blocs" in the $r^2$ matrix
  - homogenous areas in the simple matching matrix
  - dendrograms available in the app (non hierarchical but still informative)

- **subjective criteria** :
  - balancing good inertia / good adhesion / good interpretability / consistence with the domain
  - statistician approach

</details>

---

# 4. R6 architecture

---

## Wgy use R6 classes instead of a set of functions?

Because moderne clustering packages are object-oriented: a model = an object with methods.

<details><summary></summary>

R6 provides:
- encapsulation (active X, clusters, prototypes),
- dedicated methods: `fit`, `predict`, `summary`, `plot`,
- interchangeable engines,
- a uniform user-facing interface,
- extensibility.

This mirrors **tidymodels**, **mlr3** and the philosophy of **scikit-learn**.
</details>

---

## Why implement a façade class + specialized engines?

To separate user logic (interface) from algorithmic logic (engines).

<details><summary></summary>

The façade :
- unifies the API,
- manages "auto" mode,
- centralises `summary`, `plot`, `interpretation`.

The engines:
- contain the pure algorithmic logic,
- implement their own prototypes,
- exposes hooks (`run_clustering`, `predict_one_variable`, etc.).

→ modular and extensible design.
</details>

---

# 5. Shiny Application

---

## What are the main features of the Shiny app?

Impor dataset, choose variables, run clustering, display results and diagnostics, export clusters.

<details><summary></summary>>
- CSV/XLSX import
- Active/illustrative variable selection
- Method selection (or auto mode)  
- Clustering via the interface
- Outputs:
  - clusters
  - prototypes
  - inerties
  - adhesion
- Plots:
  - inertia (K-path)
  - heatmap
  - adhesion barplots
  - dendrogram
- ZIP export:
  - clusters CSV
  - summary
  - diagnostics

This exceeds the minimal project requirements.
</details>

---

# 6. *metal_universe* study case

---

## Why create a fictional thematice dataset?

To work on a controlled environment where redundancy structure is know.

<details><summary></summary>

*metal_universe* provided :
- dcontrolled correlations,
- coherents subgenres,
- noise variables,
- algorithm comparison on a realistic case,
- a "signature" dataset beyond a standard project.
</details>

---
How where the numeric variables generated?

Based on plausible ranges per subgenre + noise + derived variables.

<details><summary></summary>

Steps:
1. credible subgenres
2. plausible numeric ranges by style
3. sampling (uniform / truncated normal)
4. added noise
5. linear derived variables  
6. secondary categorical recoding

Exactly what is expected from a didactic dataset.
</details>

---

## Do the algorithms recover the simulated blocks?

Yes: rhythmic, aggressive, anf socio-geographical block clearly emerge.

<details>
  <summary></summary>

for k = 5, k-prototypes and k-medoids identify:
- a rhythmic block (BPM, speed, tempo variability)
- a "ound aggression" block
- a socio-geographical block
- a moderately mixed cluster
- one cluster for explicit noise

→ the intended structured is recovered.
</details>

---

# 7. Limitations and future directions

---

## What are the main limitations?

Local PCA, manual $\lambda$ tunin, k-medoids cost, no hierarchical method yet.

<details> <summary></summary>
- local PCA = linear modle  
- fixed $\lambda$ = manual tuning
- full dissimilary matrix in O(p²) for k-medoids  
- no true hierarchical variable clustering
- improvable initialization
- no automatic $\lambda$ selection

<sub>Regarding initialization:
- current initialization is random
```r
init <- sample(1:K, p, replace = TRUE)
```
- better heuristic :
  - **k-means++ adapted to the variables**:
    - idea: choose initial prototypes that differ greatly
    - first prototype = random variable → next prototypes = most distant variables
    - provides better stability to clusters
  - **"PCA seeds" initialization**: PCA on all variables, coordinates used as vectors on fait une ACP sur toutes les variables, grouping in $K$ with k-means to get seeds.
  - **initialization via hierarchical clustering** : hclust on variable correlations entre variables, then cutting the tree ein *k* → structured initialization
  - **frequent modes for k-modes** : choose most central variables as prototypes

`hclust()` : R function de R that implements **hierarchical clustering**.  
The function builds a **hierarchical tree** (dendrogram) by progressively merging the closest objects following a liaison method (ward, complete, average, single, etc.).  
Typical usage :
```r
d <- dist(X)
hc <- hclust(d, method = "ward.D2")
plot(hc)
```

**central variable** : variable with the most frequent category for the majority of its categories, i.e. with a **weak entropy / strong modal frequency**.

Automatic tuning of $\lambda$:
- no official answer but several strategies exist
- **balancing on variances**: $\lambda = \frac{\text{mean variance of the num block}}{\text{mean dissimilarity of the cat block}}$
- **unsupervised cross validation**: testing several $\lambda = 0.1, 0.5, 1, 2, 5, \ldots$, computing intra inertia, picking the $\lambda$ which minimizes $\mathbf{I}(K, \lambda)$
- **cluster stability** : optimal $\lambda$ = the one with provides the most stable clustering (maximal ARI maximal with bootstrapping)
- **balancing numeric adhesion / categorical adhesion** : choosing $\lambda$ so that que $adhesion_{num} = adhesion_{cat}$</sub>

</details>

---

## What imporvements have highest priority?

Automatic $\lambda$, new similarities, hierarchical methods.

<details><summary></summary>

1. automatic λ calibration
2. imprived initialization (k-means++, spectral init)
3. non-linear correlations
4. full hierarchical variable clustering
5. extended *metam_universe* dataset with web scraping  
</details>

---

# 8. "Why not…?" questions

---

## Why not implement a hierarchical method?

Because it was not required and would significantly increase algorithmic scope.

<details><summary></summary>

True hierarchical varaible clustering required:
- defining cluster-to-cluster dissimilarities,
- recomputing prototypes,
- ensuring compatiblity with mixed data.

The literature is complex → planned for future work, out of scope here.
</details>

---

## PWhy not use the ClustOfVar package?

Because the goal was to reimplement our own algorithms.

<details><summary></summary>

ClustOfVar is a reference, but:
- the project requires implementation decisions,
- pedagogical purpose,
- modular engines,
- modern R6 API vs legacy functions.
</details>

---

# 9. Compliance with assignment

---

## How do you prove that your package meets the requirements?

We exceed all criteria: methods, R6, Shiny, auto-K, English documentation.

<details><summary></summary>

- 4 algorithms (3 required)  
- multiple relecations methods (1 required)
- 2 categorical-capable algorithms (1 required)
- multiple K-selection strategies
- full R6 class architecture
- complete Shiny app
- English documentation
- internal datasets including a custom thematic one
</details>

---

# 10. Logical edge cases

---

## What happens if a variable is constant?

It is rejected or assigned zero adhesion — handled in the base class.

<details>
  <summary></summary>
In `ClusterBase` :
- detection of zero variance,
- automatic exclusion or warning,
- coherent handling of inertias/adhesion.
</details>

---

## What if all variables are categorical or all numeric?

The appropriate engine is selected automatically.

<details><summary></summary>

Auto mode: 
- all numeric → k-means  
- all categorical → k-modes  
- mixed → k-prototypes

Manual mode:
- error message if incompatible engine is chosen.
- k-medoids always available manually.
</details>
---

# 11. Interpreting results

---

## How to interpret adhesion?

It measures how well a variable fits its cluster.

<details><sumamry></summary>

- numeric: based on $r^2$
- categorical: modal agreement
- mixed: weighted combination

Higher adhesion: stronger membership.
</details>

---

## What does explained inertia represent?

Global quality of the partition.

<details><summary></summary>

Explained inertia: $\text{exp_part} = 1 - \frac{\mathbf{I}_{intra}}{\mathbf{I}_{total}}$

- higher means better partitioning
- helps compare K  
- reported summary and plots  
</details>

---

# 12. Auto mode

---

## How does the "auto" mode work?

it detects variable types and selects the appropriate algorithm.

<details><summary></summary>

Pipeline:
1. analyse column types (num, cat, mixed)  
2. choose:
   - k-means for numeric
   - k-modes for categorical
   - k-prototypes for mixed  
3. pass parameters to engine  
4. run clustering  

Fully transparent for the user.
</details>

---

# 13. More questions

---

## How to compare results with other packages?

Compare:
- inertias
- stability (bootstrap)
- adapted silhouettes
- visualisations

---

## Time complexity of the algorithms?

- k-means / k-modes / k-prototypes : O(iter × K × p)  
- k-medoids : O(p²) + reallocations

<details><summary></summary>

$O(\cdot)$ : algorithmic complexity, i.e. how fast computing time increases
- $p$ = number of variables
- $K$ = number of clusters
- $iter$ = number of iterations

$O(iter \times K \times p)$ : for each variable, at each iteration, we compare to the $K$ prototypes

$O(p^2)$ : we compute a **complete dissimilarity matrix** entre **all $p$ variables**. Greater $p$ lead to complexity explosion.

</details>  

---

## What happens if the user sets λ = 0 ?

Categorical block is ignored → reduces to local k-means on numeric part.

---
