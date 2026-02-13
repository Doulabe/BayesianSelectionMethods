# ==============================================================================
# File: Mice_Data_Application.R
# Author: N. K. Doulabe
# Description: Application of the Vectorized QR Variable Selection method to 
#              real-world genetic data.
#              - Dataset: Heterogeneous Stock (HS) Mice (Valdar et al., 2006)
#              - Trait: Obesity/BMI
#              - Goal: Demonstrate that integrating biological priors (known QTLs)
#                improves detection power compared to uniform priors.
#              - Output: Generates data for Table 6 and performance metrics.
# ==============================================================================

# ==============================================================================
# 1. SETUP AND LIBRARIES
# ==============================================================================

# setwd("path/to/your/project") 
#Ensure_packages(required_libs)

library(doParallel) # For parallel processing if implemented in helper functions
library(MASS)       # Standard statistical functions
library(BGLR)       # Contains the 'mice' dataset (Valdar et al., 2006)
library(dplyr)      # Data manipulation

# --- Load Project Functions ---
source("Test_statistics_functions.R")          # Core algorithms (NaW, ExW, etc.)
source("Test_statistics_functions_Criteria.R") # Model selection (BIC/AIC/R2)
source("Data_Analysis_Functions.R")            # Helper functions specific to real data (e.g., Get_Informed_Prior)

# ==============================================================================
# 2. DATA LOADING AND PRE-PROCESSING
# ==============================================================================

message("Loading Heterogeneous Stock Mice Data...")
data(mice) 
# 'mice.X' is the genotype matrix (n=1814 individuals, p=10346 SNPs)
# 'mice.pheno' contains various phenotypes (BMI, Insulin, etc.)
# 'mice.map' contains chromosome and position info for every SNP

# Define Target Trait: Obesity.BMI
# We use the corrected residuals or raw values as provided in the package
Y_Obesity.BMI <- mice.pheno$Obesity.BMI
X_full <- mice.X

# Check dimensions
n_mice <- nrow(X_full)
p_snps <- ncol(X_full)
message(paste("Data loaded: n =", n_mice, ", p =", p_snps))

# ==============================================================================
# 3. PRIOR CONSTRUCTION (The Biological Input)
# ==============================================================================

# We use literature-based knowledge to construct an "Informative Prior".
# For HS Mice BMI, QTLs are known to exist on Chromosomes 1, 2, 6, 7, 11, 15, 17.
map.data <- mice.map
known_targets <- c(1, 2, 6, 7, 11, 15, 17) 

# Construct the probability vector:
# - SNPs on known target chromosomes get High Probability (0.8)
# - Other SNPs get Low Probability (0.2)
informative_prior <- Get_Informed_Prior(
  map_data = map.data, 
  target_chroms = known_targets, 
  p_high = 0.8, 
  p_low = 0.2
)

# Construct the baseline: Uniform Prior (0.5 for everyone)
uniform_prior <- rep(0.5, p_snps)

# ==============================================================================
# 4. COMPARATIVE ANALYSIS (Fixed K)
# ==============================================================================
# We compare the selection performance when fixing the number of selected variants (K=20).

message("Running Analysis with Uniform Prior (Baseline)...")
results_uniform <- Real_Data_Comparison(
  X = X_full, 
  Y = Y_Obesity.BMI, 
  a.priori = uniform_prior, 
  K = 20, 
  tau = 0.2, 
  S = 1000 # MC samples for multiple tests
)

message("Running Analysis with Informative Prior (Proposed)...")
results_informative <- Real_Data_Comparison(
  X = X_full, 
  Y = Y_Obesity.BMI, 
  a.priori = informative_prior, 
  K = 20, 
  tau = 0.2, 
  S = 1000
)

# ==============================================================================
# 5. PERFORMANCE EVALUATION (Generates Table 6)
# ==============================================================================
# We evaluate success by checking if selected SNPs fall within the known target chromosomes.

message("Calculating Performance Metrics...")

# Performance for Uniform Prior
perf_table_uniform <- Calculate_Performance_Table(
  results = results_uniform, 
  map_data = mice.map, 
  known_targets = known_targets
)

# Performance for Informative Prior
perf_table_informative <- Calculate_Performance_Table(
  results = results_informative, 
  map_data = mice.map, 
  known_targets = known_targets
)

# Combine and View (This represents Table 6 in the manuscript)
final_comparison_table <- rbind(
  cbind(Prior = "Uniform", perf_table_uniform),
  cbind(Prior = "Informative", perf_table_informative)
)

print(final_comparison_table)

# ==============================================================================
# 6. SELECTION CRITERIA ANALYSIS (Adjusted R2 / BIC)
# ==============================================================================
# Here we allow the algorithm to determine K automatically using Adjusted R-squared.
# This tests the method's ability to stop overfitting.

message("Running Analysis with Stopping Criteria (Adjusted R2)...")

results_R2_inf <- Real_Data_Comparison_criteria(
  X = X_full, 
  Y = Y_Obesity.BMI, 
  a.priori = informative_prior, 
  tau = 0.2,
  criteria = "R2" # Can be switched to "BIC"
)

# Evaluate the automated selection
perf_table_R2_inf <- Calculate_Performance_Table(
  results= results_R2_inf, 
  map_data = mice.map, 
  known_targets = known_targets
)

print(perf_table_R2_inf)

# Save results for paper
# saveRDS(final_comparison_table, "Table6_RealData_Results.rds")