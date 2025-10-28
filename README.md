# mdir

Bayesian model-based clustering for multi-modal data integration. MDI jointly analyzes multiple datasets measuring the same samples, borrowing strength across modalities to improve predictions. Particularly useful when one dataset has known labels (e.g., validated protein localizations) and others don't (e.g., GO term annotations).

## Installation

```r
# install.packages("devtools")
devtools::install_github("stcolema/mdir")
```

## Example: Protein Localization Prediction

This example integrates hyperLOPIT proteomic data with Gene Ontology cellular component annotations to predict protein subcellular localization in HEK293T cells. The model learns how strongly the datasets should be coupled and propagates label information from validated markers to unlabeled proteins.

### Load and prepare data

```r
suppressPackageStartupMessages({
  library(mdir)
  library(pRoloc)
  library(pRolocdata)
})

# Load hyperLOPIT data and GO cellular component annotations
data("HEK293T2011")
data("HEK293T2011goCC")

# Extract expression matrices
X1 <- exprs(HEK293T2011)      # LOPIT fractionation profiles
X2 <- exprs(HEK293T2011goCC)  # GO term presence/absence

# Work with proteins common to both datasets
common_proteins <- intersect(rownames(X1), rownames(X2))
X1 <- X1[common_proteins, ]
X2 <- X2[common_proteins, ]

data_list <- list(X1, X2)
```

### Set up semi-supervised learning

Only a subset of proteins have experimentally validated localizations ("markers"). These act as anchor points — the model must assign them to their known compartments while inferring labels for unmarked proteins.

```r
# Extract marker information from both datasets
markers1 <- fData(HEK293T2011)[common_proteins, "markers"]
markers2 <- fData(HEK293T2011goCC)[common_proteins, "markers"]

# Convert character labels to integer cluster assignments
label_prep <- prepareMDILabels(list(markers1, markers2))
initial_labels <- label_prep$labels_matrix

# Create fixed matrix: 1 = known label (observed), 0 = unknown (infer)
fixed <- matrix(0, nrow(X1), 2)
fixed[, 1] <- ifelse(markers1 != "unknown", 1, 0)
fixed[, 2] <- ifelse(markers2 != "unknown", 1, 0)

# Allow model to discover clusters beyond known organelles
all_markers <- getMarkerClasses(HEK293T2011)
K_classes <- length(all_markers)
```

### Run the model

`callMDI` fits the model using MCMC sampling. Key parameters:
- `types`: Data types for each dataset (TAGM handles continuous LOPIT profiles with outliers, C handles categorical GO terms)
- `K`: Upper bound on clusters (model selects active number)
- `R` and `thin`: 15000 iterations, keeping every 50th (300 samples total)

```r
mdi_result <- callMDI(
  X = data_list,
  R = 15000,
  thin = 50,
  types = c("TAGM", "C"),
  K = rep(K_classes + 5, 2),
  initial_labels = initial_labels,
  fixed = fixed,
  initial_labels_as_intended = FALSE
)
```

### Process results

```r
# Discard first 3000 iterations as burn-in, compute posterior similarity matrices
processed_mdi <- processMCMCChain(mdi_result, burn = 3000, construct_psm = TRUE)

# Point estimates: most probable localization for each protein
pred_dataset1 <- all_markers[processed_mdi$pred[[1]]]
pred_dataset2 <- all_markers[processed_mdi$pred[[2]]]

# Compare predictions between datasets
comparison <- data.frame(
  protein = common_proteins,
  dataset1 = pred_dataset1,
  dataset2 = pred_dataset2,
  prob1 = processed_mdi$prob[[1]],  # Confidence in prediction
  prob2 = processed_mdi$prob[[2]]
)

# Examine discordant predictions (biological outliers or poor integration)
discordant <- comparison[pred_dataset1 != pred_dataset2, ]
head(discordant)
```

### Assess integration strength

The parameter φ (phi) quantifies coupling between datasets. High values indicate strong agreement (datasets consistently cluster proteins together), while values near 0 suggest independent structure.

```r
phi_samples <- processed_mdi$phis
cat("Median φ:", median(phi_samples), "\n")
cat("95% CI:", quantile(phi_samples, c(0.025, 0.975)), "\n")

# Visualize: narrow distribution = high confidence in coupling strength
hist(phi_samples, 
     main = "Dataset Integration Strength", 
     xlab = "φ (higher = stronger coupling)",
     breaks = 30)
abline(v = median(phi_samples), col = "red", lwd = 2)
```

### Visualize co-clustering patterns

Posterior similarity matrices (PSMs) show the probability each pair of items (here, gene products) belongs to the same organelle. Proteins in the same compartment form dark blocks along the diagonal.

```r
library(pheatmap)

col_pal <- grDevices::colorRampPalette(c("white", "#146EB4"))(100)

organelle_colors <- c(
  "Chromatin associated" = "#CC79A7",
  "Nucleus"              = "#AA3377",
  "Cytosol"              = "#56B4E9",
  "Cytosol/Nucleus"      = "#88CCEE",
  "Endosome"             = "#009E73",
  "ER"                   = "#E69F00",
  "Golgi"                = "#F0E442",
  "Lysosome"             = "#0072B2",
  "Mitochondrion"        = "#D55E00",
  "PM"                   = "#999933",
  "Ribosome 40S"         = "#882255",
  "Ribosome 60S"         = "#AA4499",
  "unknown"              = "#FFFFFF00"
)

annotation_colors <- list(Marker = organelle_colors)

pheatmap(
  processed_mdi$psms[[1]],
  main = "PSM: LOPIT data",
  annotation_row = data.frame(Marker = markers1, row.names = common_proteins),
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = col_pal,
  annotation_colors = annotation_colors
)

pheatmap(
  processed_mdi$psms[[2]],
  main = "PSM: GO terms",
  annotation_row = data.frame(Marker = markers2, row.names = common_proteins),
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = col_pal,
  annotation_colors = annotation_colors
)
```

### Running multiple chains for diagnostics

Running multiple chains with different starting points helps verify the model has converged to the true posterior rather than a local optimum. Well-mixed chains should explore the same parameter space regardless of initialization.

```r
# Run 3 chains with different random seeds
n_chains <- 3

chains <- runMCMCChains(
  X = data_list, 
  n_chains = n_chains,
  R = 15000,
  thin = 50,
  types = c("TAGM", "C"),
  K = rep(K_classes + 10, 2),
  initial_labels = initial_labels,
  fixed = fixed,
  initial_labels_as_intended = FALSE
  )

# Process each chain
processed_chains <- lapply(chains, function(chain) {
  processMCMCChain(chain, burn = 300, construct_psm = TRUE)
})
```

### Compare trace plots

Trace plots show parameter evolution across MCMC iterations. Good convergence: chains overlap and explore the same range. Poor convergence: chains separate or trend upward/downward.

```r
# Extract phi samples from each chain
phi_traces <- data.frame(
  iteration = rep(1:30, n_chains),
  phi = c(processed_chains[[1]]$phis,
          processed_chains[[2]]$phis,
          processed_chains[[3]]$phis),
  chain = rep(1:n_chains, each = 30)
)

# Visualize mixing
library(ggplot2)

ggplot(phi_traces, aes(x = iteration, y = phi, color = factor(chain))) +
  geom_line(alpha = 0.7) +
  labs(title = "φ trace plot: assessing convergence",
       x = "Iteration (post burn-in)",
       y = "φ",
       color = "Chain") +
  theme_minimal()

# Compare posterior distributions across chains
ggplot(phi_traces, aes(x = phi, fill = factor(chain))) +
  geom_density(alpha = 0.5) +
  labs(title = "φ posterior distributions by chain",
       x = "φ",
       y = "Density",
       fill = "Chain") +
  theme_minimal()
```

**What to look for:**
- **Overlapping traces**: Chains agree on parameter range (good)
- **Chains stuck at different values**: Model hasn't converged, run longer or consider consensus clustering (see below)
- **Trends up/down**: Burn-in too short, discard more iterations
- **Similar posterior densities**: Robust inference across initializations

### Consensus clustering: handling local modes

When chains become trapped in different local modes (common in high-dimensional clustering), consensus clustering offers a pragmatic solution. Rather than attempting to run a single very long chain that may never reach the global mode, consensus clustering runs many short chains from different initializations and aggregates their results — an ensemble approach that explores multiple regions of the posterior (and thus with enough chains should capture the global mode and a better description of the uncertainty in the data than a single chain).

```r
# If chains are trapped in local modes, run many short chains
n_consensus_chains <- 100  # More chains compensate for shorter length
R_short <- 500

consensus_chains <- vector("list", n_consensus_chains)

consensus_chains <- runMCMCChains(
  X = data_list, 
  n_chains = n_consensus_chains,
  R = R_short,
  thin = R_short,
  types = c("TAGM", "C"), 
  K = rep(K_classes + 10, 2),
  initial_labels = initial_labels,
  fixed = fixed,
  initial_labels_as_intended = FALSE
)

# Aggregate results across chains
consensus_result <- compileConsensusClustering(consensus_chains)

# Or predict from multiple chains (alternative approach)
# predictions <- predictFromMultipleChains(consensus_chains)
```

**Trade-offs**: This approach sacrifices strict Bayesian interpretation (samples may not be correctly weighted) but often yields accurate expected values and explores the posterior more thoroughly than a single long chain. Parallelization across chains offers substantial computational speedup over sequential sampling.

**When to use**: Poor chain mixing (chains converge to different clusterings), multimodal posteriors, or when you need results faster than a single very long chain would provide.

## Interpreting Results

- **High φ**: Datasets agree in underlying structure
- **Low φ**: Datasets capture independent biological variation, consider analyzing separately
- **PSM blocks**: Well-separated compartments; diffuse PSMs suggest overlapping biology or noise
- **Poor chain mixing**: Re-run with more iterations or check for data quality issues
- **Chains in local modes**: Consider consensus clustering to explore multiple posterior regions

## Model Types

- **G**: Gaussian (diagonal covariance)
- **MVN**: Multivariate normal (full covariance)
- **C**: Categorical
- **GP**: Gaussian process (squared exponential kernel)
- **TAGM**: t-augmented Gaussian (MVN with global MVT for outliers)
- **TAGPM**: t-augmented Gaussian Process (GP with global MVT for outliers)

## References

**Consensus clustering**: Coleman, S., Kirk, P.D.W. & Wallace, C. (2022). [Consensus clustering for Bayesian mixture models](https://doi.org/10.1186/s12859-022-04830-8). *BMC Bioinformatics* 23, 290.

**Practical guide**: Coleman, S. (2024). [Circumventing poor mixing in Bayesian model-based clustering](https://stcolema.github.io/posts/consensusClustering/consensus_clustering.html).

**Semi-supervised MDI**: Coleman, S., et al. (2024). Semi-supervised integration of single-cell transcriptomic and protein spatial distribution data. *bioRxiv*. doi: [10.1101/2024.02.08.579519](https://doi.org/10.1101/2024.02.08.579519)