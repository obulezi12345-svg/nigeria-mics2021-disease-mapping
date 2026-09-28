# ==============================================================================
# PhD Thesis Implementation Script: Spatial Bayesian Modelling of Disease Risk
# Project: MICS 2021 Spatial Modeling & Social Determinants of Health in Nigeria
# Author: Ifenna PhD / Dr. O. J. Obulezi
# ==============================================================================

library(sf)
library(readxl)
library(haven)
library(tidyverse)
library(spdep)
library(INLA)
library(ggplot2)
library(viridis)

# Set Working Directory
setwd("C:/Users/Dr. O. J. Obulezi/Documents/R projects/Ifenna PhD")

# ==============================================================================
# --- STEP 1: Load and Process MICS Survey Microdata ---
# ==============================================================================

cat("Loading survey microdata from subfolder...\n")

data_dir <- file.path(getwd(), "Nigeria MICS6 SPSS Datasets")

# 1.1 Load Children Under 5 Microdata (ch.sav)
ch_path <- file.path(data_dir, "ch.sav")
if (!file.exists(ch_path)) stop(paste("File not found at:", ch_path))

ch_raw <- read_sav(ch_path)

ch_cleaned <- ch_raw %>%
  mutate(
    State_Name = as_factor(HH7),                        # Variable HH7 represents State
    Fever_Status = if_else(CA1 == 1, 1, 0, missing = 0) # CA1: Fever in past 2 weeks
  )

# 1.2 Load Household Microdata (hh.sav)
hh_path <- file.path(data_dir, "hh.sav")
if (!file.exists(hh_path)) stop(paste("File not found at:", hh_path))

hh_raw <- read_sav(hh_path)

hh_cleaned <- hh_raw %>%
  mutate(
    State_Name = as_factor(HH7),
    Wealth_Score = as.numeric(windex5), # Wealth quintile
    Improved_Sanitation = if_else(WS11 %in% c(11, 12, 13, 21, 22), 1, 0, missing = 0)
  )

# 1.3 Aggregate Data to State Level
survey_summary <- ch_cleaned %>%
  group_by(State_Name) %>%
  summarise(
    Total_Children = n(),
    Fever_Cases = sum(Fever_Status, na.rm = TRUE),
    Fever_Prevalence = mean(Fever_Status, na.rm = TRUE),
    .groups = "drop"
  )

hh_summary <- hh_cleaned %>%
  group_by(State_Name) %>%
  summarise(
    Mean_Wealth_Quintile = mean(Wealth_Score, na.rm = TRUE),
    Prop_Improved_Sanitation = mean(Improved_Sanitation, na.rm = TRUE),
    .groups = "drop"
  )

mics_survey_summary <- survey_summary %>%
  left_join(hh_summary, by = "State_Name") %>%
  mutate(State_Clean = toupper(trimws(as.character(State_Name))))

cat("Data aggregation complete!\n")

# ==============================================================================
# --- STEP 2: Spatial Data Processing & Map Generation ---
# ==============================================================================

states_sf <- st_read("shapefiles/mics_boundaries3.shp")

states_sf <- states_sf %>%
  mutate(State_Clean = toupper(trimws(as.character(GEONAMES))))

# Merge shapefile with aggregated survey metrics
states_analysis <- states_sf %>%
  left_join(mics_survey_summary, by = "State_Clean") %>%
  filter(!is.na(Fever_Cases)) # Exclude any administrative regions without survey data

# Plot Initial Exploratory Map
windows(width = 10, height = 8)

ggplot(data = states_analysis) +
  geom_sf(aes(fill = Fever_Prevalence * 100), color = "#2b2b2b", size = 0.25) +
  geom_sf_text(
    aes(label = GEONAMES), 
    size = 2.2, 
    color = "black",
    check_overlap = TRUE
  ) +
  scale_fill_viridis_c(
    option = "magma", 
    direction = -1, 
    name = "Fever Prevalence (%)"
  ) +
  labs(
    title = "Spatial Distribution of Pediatric Fever Prevalence in Nigeria",
    subtitle = "Source: Nigeria MICS 2021 Survey",
    caption = ""
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

# ==============================================================================
# --- STEP 3: Spatial Graph & Bayesian BYM2 Model Fitting ---
# ==============================================================================

# Construct adjacency matrix from spatial geometries
nb <- poly2nb(states_analysis, queen = TRUE)
nb2INLA("nigeria_adj.graph", nb)

# Create ID index required for BYM2 model component
states_analysis$spatial_id <- 1:nrow(states_analysis)

# Corrected BYM2 Prior specification
formula_bym2 <- Fever_Cases ~ 1 + 
  scale(Mean_Wealth_Quintile) + 
  scale(Prop_Improved_Sanitation) +
  f(spatial_id, 
    model = "bym2", 
    graph = "nigeria_adj.graph",
    hyper = list(
      phi = list(prior = "pc.cor1", param = c(0.5, 0.66)),
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    )
  )

# Fit Bayesian Spatial GLMM via INLA (with marginals requested)
bym2_model <- inla(
  formula = formula_bym2,
  family = "binomial",
  Ntrials = Total_Children,
  data = as.data.frame(states_analysis),
  control.compute = list(
    dic = TRUE, 
    waic = TRUE, 
    cpo = TRUE, 
    return.marginals.predictor = TRUE # Required for exceedance probabilities
  ),
  control.predictor = list(compute = TRUE, link = 1)
)

# Output Model Summaries
cat("\n=================== MODEL SUMMARY ===================\n")
summary(bym2_model)

# ==============================================================================
# --- STEP 4: Posterior Relative Risks & Spatial Mapping ---
# ==============================================================================

cat("\nExtracting spatial relative risks and exceedance probabilities...\n")

n_states <- nrow(states_analysis)

# 4.1 Extract Fitted Mean Probabilities / Relative Risk
fitted_summary <- bym2_model$summary.fitted.values[1:n_states, ]

states_analysis$RR_mean <- fitted_summary$mean
states_analysis$RR_q025 <- fitted_summary$`0.025quant`
states_analysis$RR_q975 <- fitted_summary$`0.975quant`

# 4.2 Compute Exceedance Probabilities safely
baseline_risk <- mean(states_analysis$Fever_Prevalence, na.rm = TRUE)

# Check if marginals exist before computing exceedance
if (!is.null(bym2_model$marginals.fitted.values)) {
  states_analysis$Exceedance_Prob <- sapply(
    bym2_model$marginals.fitted.values[1:n_states], 
    function(m) 1 - inla.pmarginal(marginal = m, q = baseline_risk)
  )
} else {
  # Fallback empirical probability approximation using normal distribution
  states_analysis$Exceedance_Prob <- pnorm(
    q = (states_analysis$RR_mean - baseline_risk) / fitted_summary$sd
  )
}

# 4.3 Extract Spatial Random Effects
states_analysis$Spatial_Effect <- bym2_model$summary.random$spatial_id$mean[1:n_states]

# ==============================================================================
# --- STEP 5: Visualizing Bayesian Spatial Outputs & Data Export ---
# ==============================================================================

cat("\nGenerating maps and exporting thesis deliverables...\n")

# Map A: Posterior Expected Disease Risk
p1 <- ggplot(data = states_analysis) +
  geom_sf(aes(fill = RR_mean * 100), color = "#1a1a1a", size = 0.2) +
  scale_fill_viridis_c(
    option = "plasma", 
    name = "Modeled Risk (%)"
  ) +
  labs(
    title = "A) Bayesian Smoothed Fever Risk (BYM2)",
    subtitle = "Adjusted for Wealth & Sanitation"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(), 
    plot.title = element_text(face = "bold", size = 12)
  )

# Map B: Exceedance Probability Map (Hotspot Detection)
p2 <- ggplot(data = states_analysis) +
  geom_sf(aes(fill = Exceedance_Prob), color = "#1a1a1a", size = 0.2) +
  scale_fill_gradient2(
    low = "#2c7bb6", 
    mid = "#ffffbf", 
    high = "#d7191c", 
    midpoint = 0.50,
    name = "Probability"
  ) +
  labs(
    title = "B) Spatial Hotspot Exceedance Probability",
    subtitle = "P(State Risk > National Mean)"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(), 
    plot.title = element_text(face = "bold", size = 12)
  )

# Combine and render plots safely without gridExtra namespace conflicts
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  windows(width = 14, height = 6)
  print(p1 + p2)
} else {
  windows(width = 14, height = 6)
  suppressMessages(print(p1))
  windows(width = 14, height = 6)
  suppressMessages(print(p2))
}

# Export Results to CSV
st_drop_geometry(states_analysis) %>%
  select(
    State_Clean, Total_Children, Fever_Cases, Fever_Prevalence,
    Mean_Wealth_Quintile, Prop_Improved_Sanitation, 
    RR_mean, RR_q025, RR_q975, Exceedance_Prob, Spatial_Effect
  ) %>%
  write_csv("Nigeria_MICS2021_BYM2_Results.csv")

cat("\nAnalysis fully executed! Maps rendered and 'Nigeria_MICS2021_BYM2_Results.csv' exported.\n")










# ==============================================================================
# PhD Thesis Implementation Pipeline: Spatial Bayesian Modelling of Disease Risk
# Project: MICS 2021 Spatial Modeling & Social Determinants of Health in Nigeria
# Author: Ifenna PhD / Dr. O. J. Obulezi
# Framework: INLA + Complex Survey Weights + Spatial Diagnostics (Moran's I/LISA)
# ==============================================================================

# Clear Environment and Graphics Devices
rm(list = ls())
graphics.off()

# Set Working Directory
setwd("C:/Users/Dr. O. J. Obulezi/Documents/R projects/Ifenna PhD")

# ==============================================================================
# --- STEP 0: Load Required Packages & Global Survey Settings ---
# ==============================================================================

required_packages <- c(
  "sf", "readxl", "haven", "tidyverse", "srvyr", "survey", 
  "spdep", "INLA", "ggplot2", "viridis", "patchwork", "gridExtra"
)

# Install missing packages if necessary
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) install.packages(new_packages, dependencies = TRUE)

# Safe package loading without namespace collision
for (pkg in required_packages) {
  if (!pkg %in% (.packages())) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# --- FIX FOR SINGLE-PSU STRATA ERROR ---
# Adjusts variance calculation for strata with only 1 PSU by centering around sample mean
options(survey.lonely.psu = "adjust")

cat("\n[0/6] All packages loaded and survey settings configured successfully.\n")

# ==============================================================================
# --- STEP 1: Load Survey Microdata & Apply Complex Survey Weights ---
# ==============================================================================

cat("\n[1/6] Ingesting survey microdata and setting up weighted design...\n")

data_dir <- file.path(getwd(), "Nigeria MICS6 SPSS Datasets")

# 1.1 Ingest Children Under 5 Microdata (ch.sav)
ch_path <- file.path(data_dir, "ch.sav")
if (!file.exists(ch_path)) stop(paste("File not found at:", ch_path))
ch_raw <- read_sav(ch_path)

ch_cleaned <- ch_raw %>%
  mutate(
    State_Name = as_factor(HH7),
    Fever_Status = if_else(CA1 == 1, 1, 0, missing = 0),
    chweight_norm = chweight / 1000000 # Normalize MICS 6-implicit-decimal weights
  )

# Define Complex Survey Design for Children
ch_design <- ch_cleaned %>%
  as_survey_design(
    ids = HH1,             # Primary Sampling Units / Clusters
    strata = HH2,          # Survey Strata
    weights = chweight_norm,
    nest = TRUE            # Prevents non-nested cluster/strata ID conflicts
  )

# Aggregate Children Microdata with Survey Weights
ch_summary <- ch_design %>%
  group_by(State_Name) %>%
  summarise(
    Total_Children = n(), # Raw sample size (N) for binomial denominator
    Fever_Prevalence_Weighted = survey_mean(Fever_Status, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Derive survey-adjusted case count for INLA binomial likelihood
    Fever_Cases_INLA = round(Total_Children * Fever_Prevalence_Weighted)
  )

# 1.2 Ingest Household Microdata (hh.sav)
hh_path <- file.path(data_dir, "hh.sav")
if (!file.exists(hh_path)) stop(paste("File not found at:", hh_path))
hh_raw <- read_sav(hh_path)

hh_cleaned <- hh_raw %>%
  mutate(
    State_Name = as_factor(HH7),
    Wealth_Score = as.numeric(windex5),
    Improved_Sanitation = if_else(WS11 %in% c(11, 12, 13, 21, 22), 1, 0, missing = 0),
    hhweight_norm = hhweight / 1000000
  )

hh_design <- hh_cleaned %>%
  as_survey_design(
    ids = HH1, 
    strata = HH2, 
    weights = hhweight_norm, 
    nest = TRUE
  )

hh_summary <- hh_design %>%
  group_by(State_Name) %>%
  summarise(
    Mean_Wealth_Quintile = survey_mean(Wealth_Score, na.rm = TRUE),
    Prop_Improved_Sanitation = survey_mean(Improved_Sanitation, na.rm = TRUE),
    .groups = "drop"
  )

# Merge Aggregated Metrics
mics_survey_summary <- ch_summary %>%
  left_join(hh_summary, by = "State_Name") %>%
  mutate(State_Clean = toupper(trimws(as.character(State_Name))))

# ==============================================================================
# --- STEP 2: Spatial Data Processing & Pre-Modeling Spatial Diagnostics ---
# ==============================================================================

cat("\n[2/6] Spatial join and calculating Global Moran's I & LISA Clusters...\n")

states_sf <- st_read("shapefiles/mics_boundaries3.shp", quiet = TRUE) %>%
  mutate(State_Clean = toupper(trimws(as.character(GEONAMES))))

states_analysis <- states_sf %>%
  left_join(mics_survey_summary, by = "State_Clean") %>%
  filter(!is.na(Fever_Cases_INLA))

# Construct Spatial Adjacency Matrix
nb <- poly2nb(states_analysis, queen = TRUE)
weights_list <- nb2listw(nb, style = "W", zero.policy = TRUE)

# --- DIAGNOSTIC A: Global Moran's I ---
moran_res <- moran.test(states_analysis$Fever_Prevalence_Weighted, weights_list)
cat("\n=== GLOBAL MORAN'S I TEST FOR SPATIAL AUTOCORRELATION ===\n")
print(moran_res)

# --- DIAGNOSTIC B: Local Indicators of Spatial Association (LISA) ---
local_m <- localmoran(states_analysis$Fever_Prevalence_Weighted, weights_list)

# Standardize variables for LISA quadrant calculation
quadrant <- vector(mode = "character", length = nrow(states_analysis))
scaled_prev <- scale(states_analysis$Fever_Prevalence_Weighted)
lag_scaled_prev <- scale(lag.listw(weights_list, states_analysis$Fever_Prevalence_Weighted))

pval <- local_m[, 5]
quadrant[scaled_prev > 0 & lag_scaled_prev > 0] <- "High-High"
quadrant[scaled_prev < 0 & lag_scaled_prev < 0] <- "Low-Low"
quadrant[scaled_prev < 0 & lag_scaled_prev > 0] <- "Low-High"
quadrant[scaled_prev > 0 & lag_scaled_prev < 0] <- "High-Low"
quadrant[pval > 0.05] <- "Not Significant"

states_analysis$LISA_Cluster <- factor(
  quadrant, 
  levels = c("High-High", "Low-Low", "Low-High", "High-Low", "Not Significant")
)

# Save INLA Adjacency Graph
nb2INLA("nigeria_adj.graph", nb)
states_analysis$spatial_id <- 1:nrow(states_analysis)

# ==============================================================================
# --- STEP 3: Model Comparison Framework (M0 vs M1 vs M2 BYM2) ---
# ==============================================================================

cat("\n[3/6] Fitting candidate models (Fixed Effects, IID, BYM2) via INLA...\n")

# Model 0: Fixed Effects (No Random Effects)
f_m0 <- Fever_Cases_INLA ~ 1 + scale(Mean_Wealth_Quintile) + scale(Prop_Improved_Sanitation)
m0 <- inla(
  f_m0, family = "binomial", Ntrials = Total_Children,
  data = as.data.frame(states_analysis),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# Model 1: Unstructured Random Effects (Binomial IID / Normal random intercept)
f_m1 <- Fever_Cases_INLA ~ 1 + scale(Mean_Wealth_Quintile) + scale(Prop_Improved_Sanitation) +
  f(spatial_id, model = "iid")
m1 <- inla(
  f_m1, family = "binomial", Ntrials = Total_Children,
  data = as.data.frame(states_analysis),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# Model 2: Full BYM2 Spatial Model (PC Priors)
f_m2 <- Fever_Cases_INLA ~ 1 + scale(Mean_Wealth_Quintile) + scale(Prop_Improved_Sanitation) +
  f(spatial_id, model = "bym2", graph = "nigeria_adj.graph",
    hyper = list(
      phi = list(prior = "pc.cor1", param = c(0.5, 0.66)),
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    )
  )

m2_bym2 <- inla(
  f_m2, family = "binomial", Ntrials = Total_Children,
  data = as.data.frame(states_analysis),
  control.compute = list(
    dic = TRUE, waic = TRUE, cpo = TRUE, 
    return.marginals.predictor = TRUE
  ),
  control.predictor = list(compute = TRUE, link = 1)
)

# --- THESIS TABLE: Model Comparison Summary ---
model_comparison <- data.frame(
  Model = c("Model 0: Fixed Effects Only", "Model 1: Unstructured IID", "Model 2: BYM2 Spatial Structured"),
  DIC = c(m0$dic$dic, m1$dic$dic, m2_bym2$dic$dic),
  Effective_pD = c(m0$dic$p.eff, m1$dic$p.eff, m2_bym2$dic$p.eff),
  WAIC = c(m0$waic$waic, m1$waic$waic, m2_bym2$waic$waic),
  Marginal_LogLik = c(m0$mlik[1,1], m1$mlik[1,1], m2_bym2$mlik[1,1])
)

cat("\n=================== MODEL SELECTION & COMPARISON TABLE ===================\n")
print(model_comparison)

# ==============================================================================
# --- STEP 4: Extraction of Posteriors, Relative Risk & Exceedance ---
# ==============================================================================

cat("\n[4/6] Extracting posterior statistics from winning BYM2 model...\n")

n_states <- nrow(states_analysis)
fitted_summary <- m2_bym2$summary.fitted.values[1:n_states, ]

states_analysis$BYM2_Risk_Mean <- fitted_summary$mean
states_analysis$BYM2_Risk_Q025 <- fitted_summary$`0.025quant`
states_analysis$BYM2_Risk_Q975 <- fitted_summary$`0.975quant`

# Compute Exceedance Probabilities P(State Risk > National Weighted Mean)
national_weighted_mean <- weighted.mean(
  states_analysis$Fever_Prevalence_Weighted, 
  w = states_analysis$Total_Children
)

if (!is.null(m2_bym2$marginals.fitted.values)) {
  states_analysis$Exceedance_Prob <- sapply(
    m2_bym2$marginals.fitted.values[1:n_states], 
    function(m) 1 - inla.pmarginal(marginal = m, q = national_weighted_mean)
  )
} else {
  states_analysis$Exceedance_Prob <- pnorm(
    q = (states_analysis$BYM2_Risk_Mean - national_weighted_mean) / fitted_summary$sd
  )
}

# ==============================================================================
# --- STEP 5: Model Calibration Diagnostics (PIT & CPO Analysis) ---
# ==============================================================================

cat("\n[5/6] Generating CPO/PIT Model Calibration Diagnostic Plots...\n")

# Check for CPO failure flags
cpo_failures <- sum(m2_bym2$cpo$failure > 0, na.rm = TRUE)
cat(sprintf("CPO Calculation Warnings/Failures: %d\n", cpo_failures))

# Plot PIT Histogram
png("Thesis_PIT_Calibration_Histogram.png", width = 800, height = 600, res = 120)
hist(
  m2_bym2$cpo$pit[1:n_states], 
  breaks = seq(0, 1, by = 0.1),
  probability = TRUE,
  main = "Probability Integral Transform (PIT) Calibration Histogram",
  xlab = "PIT Values",
  ylab = "Density",
  col = "#4575b4", border = "white"
)
abline(h = 1, lty = 2, col = "red", lwd = 2)
dev.off()

# ==============================================================================
# --- STEP 6: Thesis Visualizations & Data Export ---
# ==============================================================================

cat("\n[6/6] Generating multi-panel thesis figures and exporting final datasets...\n")

# Map 1: Raw Survey-Weighted Prevalence
p_raw <- ggplot(data = states_analysis) +
  geom_sf(aes(fill = Fever_Prevalence_Weighted * 100), color = "black", size = 0.2) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Prevalence (%)") +
  labs(title = "(A) Weighted Empirical Prevalence", subtitle = "Direct Survey Estimates") +
  theme_minimal() + theme(axis.title = element_blank())

# Map 2: BYM2 Smoothed Risk
p_bym2 <- ggplot(data = states_analysis) +
  geom_sf(aes(fill = BYM2_Risk_Mean * 100), color = "black", size = 0.2) +
  scale_fill_viridis_c(option = "plasma", name = "Modeled Risk (%)") +
  labs(title = "(B) BYM2 Bayesian Posterior Risk", subtitle = "Adjusted for Covariates & Spatial Structure") +
  theme_minimal() + theme(axis.title = element_blank())

# Map 3: Hotspot Exceedance Probabilities
p_exceed <- ggplot(data = states_analysis) +
  geom_sf(aes(fill = Exceedance_Prob), color = "black", size = 0.2) +
  scale_fill_gradient2(
    low = "#2c7bb6", mid = "#ffffbf", high = "#d7191c", 
    midpoint = 0.50, name = "Probability"
  ) +
  labs(
    title = "(C) Hotspot Exceedance Probability", 
    subtitle = sprintf("P(State Risk > National Mean [%.1f%%])", national_weighted_mean * 100)
  ) +
  theme_minimal() + theme(axis.title = element_blank())

# Map 4: LISA Cluster Map
lisa_colors <- c(
  "High-High" = "#d7191c", "Low-Low" = "#2c7bb6", 
  "Low-High" = "#abd9e9", "High-Low" = "#fdae61", 
  "Not Significant" = "#eeeeee"
)

p_lisa <- ggplot(data = states_analysis) +
  geom_sf(aes(fill = LISA_Cluster), color = "black", size = 0.2) +
  scale_fill_manual(values = lisa_colors, name = "LISA Cluster") +
  labs(title = "(D) Local Moran's I Clusters", subtitle = "Unadjusted Spatial Autocorrelation") +
  theme_minimal() + theme(axis.title = element_blank())

# Combine and save 4-panel dissertation figure
thesis_4panel_map <- (p_raw + p_bym2) / (p_exceed + p_lisa)

ggsave(
  filename = "Thesis_Figure_Spatial_Fever_Analysis_4Panel.png", 
  plot = thesis_4panel_map, 
  width = 14, height = 12, dpi = 300
)

# Export Full Diagnostic CSV Dataset
st_drop_geometry(states_analysis) %>%
  select(
    State_Clean, Total_Children, Fever_Cases_INLA, 
    Fever_Prevalence_Weighted, Mean_Wealth_Quintile, Prop_Improved_Sanitation,
    BYM2_Risk_Mean, BYM2_Risk_Q025, BYM2_Risk_Q975, Exceedance_Prob, LISA_Cluster
  ) %>%
  write_csv("Thesis_Nigeria_MICS2021_Full_Spatial_Results.csv")

# Export Model Selection Table for LaTeX/Doc Inclusion
write_csv(model_comparison, "Thesis_Model_Comparison_Table.csv")

cat("\n=================== PIPELINE COMPLETED SUCCESSFULLY ===================\n")
cat("1. 4-Panel Map saved to: 'Thesis_Figure_Spatial_Fever_Analysis_4Panel.png'\n")
cat("2. PIT Histogram saved to: 'Thesis_PIT_Calibration_Histogram.png'\n")
cat("3. Full State Dataset exported to: 'Thesis_Nigeria_MICS2021_Full_Spatial_Results.csv'\n")
cat("4. Model Selection Metrics exported to: 'Thesis_Model_Comparison_Table.csv'\n")



# ==============================================================================
# PhD Thesis Extension Script: Advanced Spatial Bayesian Modeling Diagnostics
# Project: MICS 2021 Spatial Modeling & Social Determinants of Health in Nigeria
# Author: Ifenna PhD / Dr. O. J. Obulezi
# Framework: Sensitivity Analysis, Non-linear Covariates, ORs & Automated LaTeX Tables
# ==============================================================================

library(INLA)
library(tidyverse)
library(knitr)
library(kableExtra)
library(texreg)

# Assumes 'states_analysis' sf object and 'nigeria_adj.graph' are loaded from previous step

# ==============================================================================
# --- EXTENSION 1: Posterior Odds Ratios (OR) & Credible Intervals ---
# ==============================================================================

cat("\n[EXT 1/4] Extracting Fixed Effects and Calculating Adjusted Odds Ratios...\n")

# Function to extract and exponentiate fixed effects from an INLA model
get_fixed_or <- function(inla_model, model_name = "Model") {
  fixed_tab <- inla_model$summary.fixed %>%
    as.data.frame() %>%
    rownames_to_column(var = "Parameter") %>%
    mutate(
      Model = model_name,
      OR_Mean  = exp(mean),
      OR_Q025  = exp(`0.025quant`),
      OR_Q975  = exp(`0.975quant`),
      `OR (95% CrI)` = sprintf("%.2f (%.2f, %.2f)", OR_Mean, OR_Q025, OR_Q975)
    ) %>%
    select(Model, Parameter, mean, sd, OR_Mean, OR_Q025, OR_Q975, `OR (95% CrI)`)
  
  return(fixed_tab)
}

# Example usage with fitted BYM2 model (m2_bym2)
fixed_effects_or <- get_fixed_or(m2_bym2, "BYM2 Main Model")
print(fixed_effects_or %>% select(Parameter, `OR (95% CrI)`))

# Export to CSV
write_csv(fixed_effects_or, "Thesis_Fixed_Effects_Odds_Ratios.csv")


# ==============================================================================
# --- EXTENSION 2: Non-Linear Covariate Effects (Random Walk / RW1 / RW2) ---
# ==============================================================================

cat("\n[EXT 2/4] Testing Non-linear Relationships via Random Walk (RW1) Priors...\n")

# Discretize continuous wealth and sanitation metrics into fine grid bins for INLA f()
states_analysis <- states_analysis %>%
  mutate(
    wealth_bin = as.numeric(cut(Mean_Wealth_Quintile, breaks = 20, labels = FALSE)),
    sanitation_bin = as.numeric(cut(Prop_Improved_Sanitation, breaks = 20, labels = FALSE))
  )

# Non-linear Formula: RW1 on Wealth + Linear Sanitation + BYM2 Spatial
f_nonlinear <- Fever_Cases_INLA ~ 1 + 
  f(wealth_bin, model = "rw1") +
  scale(Prop_Improved_Sanitation) +
  f(spatial_id, model = "bym2", graph = "nigeria_adj.graph",
    hyper = list(
      phi  = list(prior = "pc.cor1", param = c(0.5, 0.66)),
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    )
  )

m_nonlinear <- inla(
  f_nonlinear, family = "binomial", Ntrials = Total_Children,
  data = as.data.frame(states_analysis),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

cat(sprintf("Linear BYM2 DIC: %.2f | Non-Linear Wealth RW1 DIC: %.2f\n", 
            m2_bym2$dic$dic, m_nonlinear$dic$dic))

# Plot non-linear smooth effect of Wealth
wealth_rw1_effect <- m_nonlinear$summary.random$wealth_bin

png("Thesis_Nonlinear_Wealth_Effect_RW1.png", width = 800, height = 600, res = 120)
plot(wealth_rw1_effect$mean, type = "l", lwd = 2, col = "darkblue",      ylim = range(c(wealth_rw1_effect$`0.025quant`, wealth_rw1_effect$`0.975quant`)),
     main = "Non-linear Effect of Wealth Quintile on Pediatric Fever Risk",
     xlab = "Binned Wealth Index", ylab = "Log-Odds Contribution")
lines(wealth_rw1_effect$`0.025quant`, lty = 2, col = "red")
lines(wealth_rw1_effect$`0.975quant`, lty = 2, col = "red")
legend("topright", legend = c("Posterior Mean", "95% Credible Interval"),
       col = c("darkblue", "red"), lty = c(1, 2), bty = "n")
dev.off()


# ==============================================================================
# --- EXTENSION 3: Prior Sensitivity Analysis Across Hyperparameters ---
# ==============================================================================

cat("\n[EXT 3/4] Performing Prior Sensitivity Analysis on Spatial Precision & PC Cor...\n")

# Define grid of PC hyperprior variations for phi and precision
prior_specs <- list(
  "Base PC Prior"         = list(phi = c(0.5, 0.66), prec = c(1, 0.01)),
  "Informative Spatial"   = list(phi = c(0.8, 0.50), prec = c(0.5, 0.01)),
  "Diffuse Spatial"       = list(phi = c(0.2, 0.80), prec = c(2, 0.01)),
  "High Precision Variance"= list(phi = c(0.5, 0.66), prec = c(0.1, 0.01))
)

sensitivity_results <- data.frame()

for (p_name in names(prior_specs)) {
  p_val <- prior_specs[[p_name]]
  
  f_sens <- Fever_Cases_INLA ~ 1 + scale(Mean_Wealth_Quintile) + scale(Prop_Improved_Sanitation) +
    f(spatial_id, model = "bym2", graph = "nigeria_adj.graph",
      hyper = list(
        phi  = list(prior = "pc.cor1", param = p_val$phi),
        prec = list(prior = "pc.prec", param = p_val$prec)
      )
    )
  
  mod_sens <- inla(
    f_sens, family = "binomial", Ntrials = Total_Children,
    data = as.data.frame(states_analysis),
    control.compute = list(dic = TRUE, waic = TRUE)
  )
  
  # Extract posterior estimates for spatial fraction phi and fixed effects
  sens_row <- data.frame(
    Prior_Specification = p_name,
    DIC = mod_sens$dic$dic,
    WAIC = mod_sens$waic$waic,
    Intercept_Mean = mod_sens$summary.fixed["(Intercept)", "mean"],
    Wealth_Beta = mod_sens$summary.fixed["scale(Mean_Wealth_Quintile)", "mean"],
    Sanitation_Beta = mod_sens$summary.fixed["scale(Prop_Improved_Sanitation)", "mean"]
  )
  
  sensitivity_results <- rbind(sensitivity_results, sens_row)
}

cat("\n=================== SENSITIVITY ANALYSIS SUMMARY ===================\n")
print(sensitivity_results)

write_csv(sensitivity_results, "Thesis_Sensitivity_Analysis_Results.csv")


# ==============================================================================
# --- EXTENSION 4: Automated LaTeX / HTML Table Generator for Dissertation ---
# ==============================================================================

cat("\n[EXT 4/4] Formatting LaTeX Tables for Direct Copy-Paste into LaTeX Thesis...\n")

# 4.1 Export Model Selection Comparison Table to LaTeX
latex_comparison_table <- model_comparison %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    caption = "Model Comparison and Goodness-of-Fit Diagnostics for Pediatric Fever Risk",
    label = "tab:model_comparison",
    col.names = c("Model Specification", "DIC", "Effective Parameters ($p_d$)", "WAIC", "Marginal Log-Likelihood"),
    digits = 2
  ) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

# Save LaTeX Code to file
writeLines(latex_comparison_table, "Thesis_Table_Model_Comparison.tex")

# 4.2 Export Sensitivity Analysis Table to LaTeX
latex_sensitivity_table <- sensitivity_results %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    caption = "Prior Sensitivity Analysis Across Various PC Prior Hyperparameters",
    label = "tab:sensitivity_analysis",
    digits = 3
  ) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))

writeLines(latex_sensitivity_table, "Thesis_Table_Sensitivity_Analysis.tex")

cat("\n=================== EXTENSION PIPELINE COMPLETED ===================\n")
cat("1. Adjusted Odds Ratios exported to: 'Thesis_Fixed_Effects_Odds_Ratios.csv'\n")
cat("2. Non-linear RW1 plot saved to: 'Thesis_Nonlinear_Wealth_Effect_RW1.png'\n")
cat("3. Sensitivity Analysis exported to: 'Thesis_Sensitivity_Analysis_Results.csv'\n")
cat("4. LaTeX Dissertation Tables written to '.tex' files in working directory.\n")






# ==============================================================================
# PhD Thesis Script: Residual Diagnostics & Cluster-Level Spatial INLA Framework
# Project: MICS 2021 Spatial Modeling & Social Determinants of Health in Nigeria
# Author: Ifenna PhD / Dr. O. J. Obulezi
# Framework: Cluster-Level (EA/PSU) INLA SPDE (Matern Field) + Residual Analysis
# ==============================================================================

library(sf)
library(readxl)
library(haven)
library(tidyverse)
library(srvyr)
library(survey)
library(spdep)
library(INLA)
library(ggplot2)
library(viridis)

setwd("C:/Users/Dr. O. J. Obulezi/Documents/R projects/Ifenna PhD")

# ==============================================================================
# --- PART 1: RESIDUAL DIAGNOSTICS FOR EXISTING STATE-LEVEL BYM2 MODEL ---
# ==============================================================================

cat("\n[PART 1] Running residual diagnostics and spatial autocorrelation tests on BYM2...\n")

# 1.1 Compute Deviance and Pearson Residuals for State-Level Model (m2_bym2)
# Observed binomial events (y) and denominator (N)
y_obs <- states_analysis$Fever_Cases_INLA
n_obs <- states_analysis$Total_Children
p_hat <- states_analysis$BYM2_Risk_Mean

# Fitted expected cases
mu_hat <- n_obs * p_hat

# Pearson Residuals: (y - mu) / sqrt(Var(y))
var_y <- n_obs * p_hat * (1 - p_hat)
states_analysis$Pearson_Residuals <- (y_obs - mu_hat) / sqrt(var_y)

# Deviance Residuals for Binomial Likelihood
sign_diff <- sign(y_obs - mu_hat)
term1 <- ifelse(y_obs == 0, 0, y_obs * log(y_obs / mu_hat))
term2 <- ifelse(y_obs == n_obs, 0, (n_obs - y_obs) * log((n_obs - y_obs) / (n_obs - mu_hat)))
dev_sq <- 2 * (term1 + term2)
states_analysis$Deviance_Residuals <- sign_diff * sqrt(pmax(0, dev_sq))

# 1.2 Test Residuals for Remaining Spatial Autocorrelation (Moran's I)
cat("\n=== MORAN'S I TEST ON BYM2 DEVIANCE RESIDUALS ===\n")
res_moran <- moran.test(states_analysis$Deviance_Residuals, weights_list, zero.policy = TRUE)
print(res_moran)

# 1.3 Export Residual Diagnostic Plots
png("Thesis_Residual_Diagnostics_BYM2.png", width = 1200, height = 400, res = 120)
par(mfrow = c(1, 3))

# Plot 1: Fitted vs Deviance Residuals
plot(states_analysis$BYM2_Risk_Mean, states_analysis$Deviance_Residuals,
     xlab = "Fitted Prevalence (Risk)", ylab = "Deviance Residuals",
     main = "Residuals vs Fitted Risk", pch = 19, col = "#2b5c8f")
abline(h = 0, lty = 2, col = "red", lwd = 1.5)

# Plot 2: Normal Q-Q Plot of Residuals
qqnorm(states_analysis$Deviance_Residuals, main = "Normal Q-Q Plot (Deviance Residuals)", pch = 19, col = "#2b5c8f")
qqline(states_analysis$Deviance_Residuals, col = "red", lwd = 1.5)

# Plot 3: Residual Distribution Histogram
hist(states_analysis$Deviance_Residuals, breaks = 12, probability = TRUE,
     main = "Residual Distribution", xlab = "Deviance Residuals", col = "#a6bddb", border = "white")
lines(density(states_analysis$Deviance_Residuals), col = "darkblue", lwd = 2)

dev.off()
par(mfrow = c(1, 1))

cat("State-level residual diagnostics completed. Diagnostic plot saved to 'Thesis_Residual_Diagnostics_BYM2.png'.\n")

# ==============================================================================
# --- PART 2: EXPANDING TO CLUSTER-LEVEL (PSU/EA) SPATIAL RESOLUTION ---
# ==============================================================================

cat("\n[PART 2] Processing microdata at Cluster Level (HH1 Primary Sampling Units)...\n")

# 2.1 Aggregate Microdata to Primary Sampling Unit / Cluster Level (HH1)
data_dir <- file.path(getwd(), "Nigeria MICS6 SPSS Datasets")

ch_raw <- read_sav(file.path(data_dir, "ch.sav"))
hh_raw <- read_sav(file.path(data_dir, "hh.sav"))

# Process Child Microdata by Cluster
ch_cluster <- ch_raw %>%
  mutate(
    Cluster_ID = as.character(HH1),
    State_Name = as_factor(HH7),
    Fever_Status = if_else(CA1 == 1, 1, 0, missing = 0),
    chweight_norm = chweight / 1000000
  ) %>%
  group_by(Cluster_ID, State_Name) %>%
  summarise(
    Cluster_Total_Children = n(),
    Cluster_Fever_Cases = sum(Fever_Status, na.rm = TRUE),
    Cluster_Weighted_Prevalence = weighted.mean(Fever_Status, w = chweight_norm, na.rm = TRUE),
    .groups = "drop"
  )

# Process Household Microdata by Cluster
hh_cluster <- hh_raw %>%
  mutate(
    Cluster_ID = as.character(HH1),
    Wealth_Score = as.numeric(windex5),
    Improved_Sanitation = if_else(WS11 %in% c(11, 12, 13, 21, 22), 1, 0, missing = 0),
    hhweight_norm = hhweight / 1000000
  ) %>%
  group_by(Cluster_ID) %>%
  summarise(
    Cluster_Mean_Wealth = weighted.mean(Wealth_Score, w = hhweight_norm, na.rm = TRUE),
    Cluster_Prop_Sanitation = weighted.mean(Improved_Sanitation, w = hhweight_norm, na.rm = TRUE),
    .groups = "drop"
  )

# Combine Cluster-Level Data
cluster_dataset <- ch_cluster %>%
  left_join(hh_cluster, by = "Cluster_ID") %>%
  filter(!is.na(Cluster_Fever_Cases) & Cluster_Total_Children > 0)

cat(sprintf("Aggregated dataset to %d Primary Sampling Unit (PSU) clusters.\n", nrow(cluster_dataset)))

# ==============================================================================
# --- PART 3: SPATIAL MESH & SPDE MODELING FOR CLUSTER DATA ---
# ==============================================================================

cat("\n[PART 3] Constructing Spatial Mesh & SPDE Model Framework...\n")

# Note: In standard MICS workflows, cluster GPS points (lat/long) are joined here.
# If cluster shapefile/coordinates are available, load cluster_coords sf object.
# Here we demonstrate setting up the INLA SPDE mesh for Nigeria using state boundaries.

nigeria_domain <- st_union(states_sf)
nigeria_sp <- as(nigeria_domain, "Spatial")

# Extract State Centroids as cluster coordinate proxies if point GPS are unjittered/aggregated
centroids <- st_centroid(states_analysis)
coords_matrix <- st_coordinates(centroids)

# 3.1 Construct Constrained Delaunay Triangulation Mesh
mesh <- inla.mesh.2d(
  loc = coords_matrix,
  max.edge = c(0.5, 2.0), # Inner and outer max edge lengths (in degrees)
  cutoff = 0.1,           # Minimum distance between mesh vertices
  crs = st_crs(states_sf)
)

# 3.2 Define SPDE Model using Matérn Covariance and PC Priors
# Prior: P(spatial range < 1.0 degree) = 0.01; P(spatial sigma > 1.0) = 0.01
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(1.0, 0.01),
  prior.sigma = c(1.0, 0.01)
)

# 3.3 Create Projector Matrix A and Spatial Index
A_matrix <- inla.spde.make.A(mesh = mesh, loc = coords_matrix)
s_index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

# 3.4 Build INLA Data Stack for Cluster Model
# Map cluster covariates to state-level spatial coordinates (if GPS cluster level unavailable)
cluster_spatial_data <- cluster_dataset %>%
  mutate(State_Clean = toupper(trimws(as.character(State_Name)))) %>%
  left_join(
    states_analysis %>% st_drop_geometry() %>% select(State_Clean, spatial_id),
    by = "State_Clean"
  ) %>%
  filter(!is.na(spatial_id))

# Project cluster coordinates matching spatial index
coords_cluster <- coords_matrix[cluster_spatial_data$spatial_id, ]
A_cluster <- inla.spde.make.A(mesh = mesh, loc = coords_cluster)

stack_cluster <- inla.stack(
  data = list(y = cluster_spatial_data$Cluster_Fever_Cases, Ntrials = cluster_spatial_data$Cluster_Total_Children),
  A = list(A_cluster, 1),
  effects = list(
    s_index,
    list(
      Intercept = 1,
      Wealth = scale(cluster_spatial_data$Cluster_Mean_Wealth),
      Sanitation = scale(cluster_spatial_data$Cluster_Prop_Sanitation)
    )
  ),
  tag = "cluster_spde_fit"
)

# 3.5 Fit Spatial Continuous SPDE GLMM Model via INLA
cat("\nFitting Cluster-Level Continuous SPDE Model via INLA...\n")

formula_spde <- y ~ 0 + Intercept + Wealth + Sanitation + f(spatial.field, model = spde)

m_spde_cluster <- inla(
  formula_spde,
  family = "binomial",
  Ntrials = Ntrials,
  data = inla.stack.data(stack_cluster),
  control.predictor = list(A = inla.stack.A(stack_cluster), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

cat("\n=================== CLUSTER SPDE MODEL SUMMARY ===================\n")
summary(m_spde_cluster)

# ==============================================================================
# --- PART 4: RESIDUAL DIAGNOSTICS FOR CLUSTER-LEVEL SPDE MODEL ---
# ==============================================================================

cat("\n[PART 4] Running residual diagnostics on Cluster-Level SPDE Model...\n")

# Extract predicted probabilities for cluster data
fitted_idx <- inla.stack.index(stack_cluster, "cluster_spde_fit")$data
cluster_p_hat <- m_spde_cluster$summary.fitted.values$mean[fitted_idx]

y_cl <- cluster_spatial_data$Cluster_Fever_Cases
n_cl <- cluster_spatial_data$Cluster_Total_Children
mu_cl <- n_cl * cluster_p_hat

# Cluster Pearson Residuals
cluster_spatial_data$Pearson_Residuals <- (y_cl - mu_cl) / sqrt(mu_cl * (1 - cluster_p_hat))

# Cluster Deviance Residuals
sign_cl <- sign(y_cl - mu_cl)
t1_cl <- ifelse(y_cl == 0, 0, y_cl * log(y_cl / mu_cl))
t2_cl <- ifelse(y_cl == n_cl, 0, (n_cl - y_cl) * log((n_cl - y_cl) / (n_cl - mu_cl)))
cluster_spatial_data$Deviance_Residuals <- sign_cl * sqrt(pmax(0, 2 * (t1_cl + t2_cl)))

# Export Cluster Results
write_csv(cluster_spatial_data, "Thesis_Cluster_Level_SPDE_Results.csv")

cat("\n=================== EXECUTION COMPLETED ===================\n")
cat("1. State-level BYM2 residual diagnostics plotted & tested with Moran's I.\n")
cat("2. Cluster dataset aggregated across primary sampling units (PSUs).\n")
cat("3. Spatial SPDE mesh & continuous Matérn field fitted to cluster microdata.\n")
cat("4. Cluster SPDE diagnostic outputs exported to 'Thesis_Cluster_Level_SPDE_Results.csv'.\n")