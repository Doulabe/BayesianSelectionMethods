##### -------- REAL DATA ANALYSIS ------------
#setwd("") set directory 

# --- 1. Load Libraries & Data ---
library(doParallel)
library(MASS)
library(BGLR)      # For the mice dataset
library(ggplot2)   # For plotting
library(dplyr)     # For data manipulation

#--- required functions
source("Test_statistics_functions.R")
source("Test_statistics_functions_Criteria.R")
source("Data_Analysis_Functions.R")

# Load the mice dataset
data(mice) # 'mice.X' is the genotype matrix (n=1814, p=10346)
            # 'mice.pheno' contains the phenotypes
            # We focus on 'Obesity.BMI'
Y_Obesity.BMI <- mice.pheno$Obesity.BMI
X_full <- mice.X

# --- 2. Pre-processing ---
map.data <- mice.map
known_targets <- c(1, 2, 6, 7, 11, 15, 17) # from literature
informative_prior <- Get_Informed_Prior(map.data, known_targets, p_high = 0.8, p_low = 0.2)

results_uniform <- Real_Data_Comparison(
  X = X_full, 
  Y = Y_Obesity.BMI, 
  a.priori = rep(0.5, ncol(X_full)), 
  K = 20, 
  tau = 0.2, 
  S = 1000
)

results_informative <- Real_Data_Comparison(
  X = X_full, 
  Y = Y_Obesity.BMI, 
  a.priori = informative_prior, 
  K = 20, 
  tau = 0.2, 
  S = 1000
)

perf_table_uniform <- Calculate_Performance_Table(results_uniform, mice.map, known_targets) # produce Table 6 
perf_table_informative <- Calculate_Performance_Table(results_informative, mice.map, known_targets)

#---- Selection criteria
# Example with R2 and informed prior

results_R2_inf <- Real_Data_Comparison_criteria(
  X = X_full, 
  Y = Y_Obesity.BMI, 
  a.priori = informative_prior, 
  tau = 0.2,
  critere = "R2"
)
# View Summary Table
known_targets <- c(1, 2, 6, 7, 11, 15, 17)
perf_table_R2_inf <- Calculate_Performance_Table(results_R2_inf, mice.map, known_targets)


