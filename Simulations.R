# ==============================================================================
# File: Simulation.R
# Author: N. K. Doulabe
# Description: Main execution script for the simulation study.
#              - Generates synthetic genotype/phenotype data with LD.
#              - Runs Naive, Asymptotic, and Exact tests (Wald & Score).
#              - Compares fixed K selection vs. Information Criteria (BIC).
#              - Reproduces Tables 1, 2, 3, and 4 from the manuscript.
# ==============================================================================

# ==============================================================================
# 1. SETUP AND DEPENDENCIES
# ==============================================================================

# Set your working directory to the folder containing the source files
# setwd("path/to/your/project") 

# Load helper functions
source("GenerateXY.R")                  # Data generation (Genotypes X and Phenotypes Y)
source("Test_statistics_functions.R")   # Core algorithms (NaW, AsW, ExW, AsS, NaM, ExM)
source("Test_statistics_functions_Criteria.R") # Model selection criteria (BIC/AIC) wrapper

# ==============================================================================
# 2. SIMULATION PARAMETERS (Corresponds to Table 1)
# ==============================================================================

# --- Experimental Grid ---
p_values   <- c(0.2, 0.6, 0.9)  # Target Power levels to calibrate effect sizes
tau_values <- c(0.2, 0.5, 1)    # Prior Standard Deviations (signal strength variability)

# --- Dataset Dimensions & Genetics ---
n <- 400                        # Sample size (Number of individuals)
m <- 1000                       # High-dimensional setting (Number of SNPs)
r <- 0.2                        # Pairwise correlation (Simulating Linkage Disequilibrium)
f <- 0.2                        # Minor Allele Frequency (MAF)

# --- Algorithm Settings ---
alpha <- 0.05                   # Significance level
K <- 15                         # True number of causal SNPs (Sparse setting)
B <- 1000                       # Number of simulation repetitions (for robust error estimation)
S <- 1000                       # Monte Carlo samples for Exact Test integral approximation

# --- Random Seed ---
# Ensures reproducibility of the X matrix and simulation results
sigma <- 1
seed <- 160323
set.seed(seed)

# ==============================================================================
# 3. PRIOR ARCHITECTURES (Scenarios I - IV)
# ==============================================================================

# Scenario I: Informative Prior (Ideal Case)
# High probability assigned to the true causal variants (first 15).
a.priori1 <- c(rep(0.8, 15), rep(0.2, m - 15))

# Scenario II: Uniform Prior (Uninformative)
# All variants have equal probability; tests the method's baseline performance.
a.priori2 <- rep(0.5, m)

# Scenario III: Misspecified Prior (Adversarial Case)
# Low probability assigned to true causal variants; tests robustness.
a.priori3 <- c(rep(0.2, 15), rep(0.8, 15), rep(0.5, m - 30))

# Scenario IV: Mixed Prior (Block Structure)
# Alternating blocks of probabilities; simulates realistic genomic uncertainty.
a.priori4 <- c(rep(c(rep(0.2, 5), rep(0.5, 5), rep(0.8, 5)), 2), rep(0.5, m - 30))

# ==============================================================================
# 4. DATA GENERATION
# ==============================================================================

message("Generating Genotype Matrix X with LD structure...")
# Generates centered X matrix with AR(1) correlation structure to mimic LD
X <- Generate.X(n, r, m, f) 


# ==============================================================================
# 5. EXECUTION: SINGLE RUN CHECKS (Unit Testing)
# ==============================================================================
# Use this section to verify code is working before running the full simulation loop.

## 5.1. Test with Fixed K (Top 15 variables)
# ------------------------------------------------------------------------------
# Example using Naive Wald in Scenario I
# Note: 'tau' must be defined (e.g., tau <- 0.2) before running this single line
tau_test <- 0.2
# Generate a temporary Y for testing
beta_test <- c(1, rep(0.5, 15), rep(0, m - 15))
Y_test <- Generate.Y(n, X, beta_test, sigma)

res_Naive_Wald_fixed_K <- Naive_Wald(X, Y_test, a.priori = a.priori1, K = 15, tau = tau_test)

# Inspect outputs:
# res_Naive_Wald_fixed_K$post.probs   # Posterior probabilities of selected models
# res_Naive_Wald_fixed_K$indices      # Indices of the 15 selected SNPs

## 5.2. Test with Stopping Criteria (BIC)
# ------------------------------------------------------------------------------
res_Naive_Wald_det_K <- Naive_Wald_Criteria(
  X, Y_test, 
  a.priori = a.priori1, 
  tau = tau_test, 
  threshold = 16, 
  target_prob = 0.95, 
  criteria = "BIC" # Switch to "AIC" or "R2_adj" if needed
)

# Inspect outputs:
# res_Naive_Wald_det_K$indices        # Selected indices (stopped automatically)
# res_Naive_Wald_det_K$bic_values     # Trajectory of BIC values

# ==============================================================================
# 6. FULL SIMULATION STUDY (Reproducing Paper Results)
# ==============================================================================
# WARNING: These functions use parallel processing and may take time.

## 6.1. Fixed K Analysis (Generates Tables 2 & 3)
# ------------------------------------------------------------------------------
message("Running full simulation for Scenario I (Fixed K)...")

# Runs B=1000 repetitions across the p/tau grid
results_p_tau_fixed_K_I <- Simulation(
  X = X, 
  p_values = p_values, 
  tau_values = tau_values,
  a.priori = a.priori1, 
  alpha = alpha, 
  K = K, 
  B = B, 
  S = S, 
  threshold = 16
)

# Save computationally expensive results
saveRDS(results_p_tau_fixed_K_I, "results_p_tau_fixed_K_I_per_p_tau_B1000.rds")

# --- Analyze Results ---

# 1. Type I Error and True Discovery Rates
# This output corresponds to the columns in Table 2 of the manuscript
print(results_p_tau_fixed_K_I$main_results)

# 2. Probability of Superiority (Pairwise Comparison)
# This output corresponds to Table 3 of the manuscript
superior_probabilty_results <- Compare_Methods_Matrix(results_p_tau_fixed_K_I$selection_summary)

# Example: View comparison matrix for specific parameters
# print(superior_probabilty_results$`p=0.2_tau=0.2`)


## 6.2. Criteria-Based Analysis (Generates Table 4)
# ------------------------------------------------------------------------------
message("Running simulation with BIC stopping criteria...")

# Runs simulation where K is determined by BIC minimization
results_p_tau_det_K_BIC <- Criteria_Simulation_Table(
  B = B, 
  X = X, 
  sigma = sigma, 
  a.priori = a.priori1, 
  p_values = p_values, 
  tau_values = tau_values, 
  S = S,  
  threshold = 16, 
  critere = "BIC"
)

# This output corresponds to **Table 4** (Precision/Recall/F1 scores with BIC)
print(results_p_tau_det_K_BIC)

saveRDS(results_p_tau_det_K_BIC, "results_p_tau_det_K_BIC.rds")