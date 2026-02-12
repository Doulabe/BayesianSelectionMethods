# ==============================================================================
# File: main.R
# Project: Bayesian Variable Selection with Vectorized QR Decomposition
# Author: N.K. Doulabe
# Repository: https://github.com/Doulabe/BayesianSelectionMethods
# 
# Description: 
#   This is the master script to reproduce the results presented in the manuscript.
#   It serves three purposes:
#   1. Installs/Loads necessary dependencies.
#   2. "Quick Start": A demo example to demonstrate the method mechanics.
#   3. "Reproduction": Blocks to run the full simulations (Tables 2, 3, 4) 
#      and Real Data Analysis (Table 6).
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP & DEPENDENCIES
# ------------------------------------------------------------------------------

# Load Project Functions
# Ensure these files are in your working directory
source("GenerateXY.R")                         # Data generation
source("Test_statistics_functions.R")          # Core Algorithms (NaW, AsW, ExW, etc.)
source("Test_statistics_functions_Criteria.R") # Algorithms with Stopping Criteria (BIC/R2)
source("Simulation_Main.R")                    # Functions for Simulation Tables
source("Real_Data_Analysis.R")                 # Functions for Mice Data Analysis

# List of required libraries
required_libs <- c("base", "MASS", "stats", "Matrix", "data.table", "doParallel", 
                   "mvtnorm", "BGLR", "dplyr") # BGLR needed for Mice Data
Ensure_packages(required_libs)
# ------------------------------------------------------------------------------
# 2. DEMO EXAMPLE (Quick Start)
# ------------------------------------------------------------------------------
message("--- Running Demo Example ---")

# A. Generate Synthetic Data
n <- 200; m <- 500; r <- 0.2; f <- 0.2
sigma <- 1
set.seed(123)

# Generate Genotypes (X) with LD and Phenotypes (Y)
X_demo <- Generate.X(n, r, m, f)
# Create 5 causal variants with effect size 0.5
beta_demo <- c(1, rep(0.5, 5), rep(0, m - 5)) 
Y_demo <- Generate.Y(n, X_demo, beta_demo, sigma)

# B. Define Priors (Example: Informative Prior)
# High probability (0.8) for first 10 SNPs, Low (0.2) for rest
prior_demo <- c(rep(0.8, 10), rep(0.2, m - 10))

# C. Run Exact Wald Test (ExW)
# We ask it to select top K=5 variables
results_demo <- Exact_Wald(X_demo, Y_demo, a.priori = prior_demo, K = 5, tau = 0.2)

print("Selected Indices (demo Example):")
print(results_demo$indices)
print("Posterior Probabilities:")
print(results_demo$post.probs)

# ------------------------------------------------------------------------------
# 3. REPRODUCING SIMULATION TABLES (Tables 2, 3, 4)
# ------------------------------------------------------------------------------
# NOTE: These steps are computationally intensive. 
# Adjust 'B' (repetitions) and 'num_cores' in the source files for testing.

# --- Table 2 & 3: Fixed K Comparison ---
# Parameters from manuscript
p_vals <- c(0.2, 0.9); tau_vals <- c(0.2, 0.5)
# Run Simulation
# results_fixed_K <- Simulation(X_demo, p_vals, tau_vals, prior_demo, alpha=0.05, K=5, B=100, S=100)
# print(results_fixed_K$main_results)

# --- Table 4: Stopping Criteria (BIC) ---
# results_BIC <- Criteria_Simulation_Table(B=100, X_demo, sigma, prior_demo, p_vals, tau_vals, criteria="BIC")
# print(results_BIC)

# ------------------------------------------------------------------------------
# 4. REPRODUCING REAL DATA ANALYSIS (Table 6)
# ------------------------------------------------------------------------------
# Requires 'BGLR' package and the 'mice' dataset.

# if(require(BGLR)) {
#   data(mice)
#   Y_bmi <- mice.pheno$Obesity.BMI
#   X_mice <- mice.X
#   
#   # Run Analysis (Commented out due to runtime)
#   # results_mice <- Real_Data_Comparison(X_mice, Y_bmi, rep(0.5, ncol(X_mice)), K=10, tau=0.2, S=500)
#   # print(results_mice)
# }