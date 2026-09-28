# nigeria-mics2021-disease-mapping
Spatial Bayesian modeling and disease risk mapping of pediatric fever prevalence in Nigeria using MICS 2021 survey microdata, R-INLA, BYM2, and SPDE frameworks.
# Spatial Bayesian Modeling of Pediatric Disease Risk in Nigeria (MICS 2021)

This repository contains the complete R analysis pipeline for spatial modeling and risk mapping of childhood fever prevalence across Nigeria using microdata from the 2021 Multiple Indicator Cluster Survey (MICS 6). 

The framework incorporates complex survey weighting, spatial autocorrelation diagnostics (Moran's I & LISA), Bayesian Integrated Nested Laplace Approximations (INLA), spatial GLMMs (BYM2 model), non-linear continuous covariates (RW1), and continuous spatial point-process modeling (SPDE Matérn field).

---

## 📌 Features & Methodology

- **Survey-Weighted Microdata Aggregation**: Ingests SPSS datasets (`ch.sav` and `hh.sav`), adjusts for complex sampling design (PSU, strata, and normalized child weights), and handles lonely PSU variance adjustments using `srvyr` and `survey`.
- **Exploratory Spatial Data Analysis (ESDA)**: Computes Global Moran's I and Local Indicators of Spatial Association (LISA) to categorize spatial clusters (High-High, Low-Low, etc.).
- **Spatial Bayesian Modeling via INLA**:
  - **Model 0**: Fixed Effects (Baseline Logistic GLM).
  - **Model 1**: Unstructured Independent and Identically Distributed (IID) Random Effects.
  - **Model 2**: Besag-York-Mollié (BYM2) Spatial Model with Penalized Complexity (PC) priors.
  - **Extension Models**: Non-linear covariate terms (Random Walk 1) and continuous spatial SPDE mesh modeling at the PSU/cluster level.
- **Model Evaluation & Sensitivity**: Evaluates goodness-of-fit via DIC, WAIC, CPO/PIT calibration histograms, and conducts PC hyperprior sensitivity checks.
- **Automated Outputs**: Generates multi-panel maps using `ggplot2` and `patchwork`, calculates adjusted Odds Ratios (ORs) with 95% Credible Intervals, and exports automated LaTeX dissertation tables via `kableExtra`.

---

## 📁 Repository Structure

```text
.
├── scripts/
│   ├── implementation_script.R        # Complete 3-part spatial Bayesian pipeline
├── data/                             # (Git-ignored) Raw MICS 6 survey data directory
│   ├── Nigeria MICS6 SPSS Datasets/
│   │   ├── ch.sav                    # Children under 5 dataset
│   │   └── hh.sav                    # Household dataset
│   └── shapefiles/
│       └── mics_boundaries3.shp      # Nigeria administrative spatial boundaries
├── output/                           # Output visualizations, datasets, and LaTeX code
│   ├── Thesis_Figure_Spatial_Fever_Analysis_4Panel.png
│   ├── Thesis_PIT_Calibration_Histogram.png
│   ├── Thesis_Nigeria_MICS2021_Full_Spatial_Results.csv
│   └── Thesis_Table_Model_Comparison.tex
├── README.md
└── .gitignore
