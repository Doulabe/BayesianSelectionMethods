##### -------- SIMULATION SET-UP ------------
#setwd("") set directory 
source("GenerateXY.R")
source("Test_statistics_functions.R")
source("Test_statistics_functions_Criteria.R")

#--- Simulations parameters (Table 1)

p_values <- c(0.2, 0.6, 0.9)  # power values
tau_values <- c(0.2, 0.5, 1)  # prior SD tau values
n <- 400                      # number of individuals
r <- 0.2                      # pairwise correlation   
m <- 1000                     # number of SNPs
f <- 0.2                      # MAF
alpha <- 0.05                 # alpha value
K <- 15                       # specified number of causal SNPs 
B <- 1000                     # number of simulations repetitions
S <- 1000                     # Monte Carlo for Multiple tests integrals calculation (Fisher and Chisq)

##---- Prior architecture
a.priori1 <- c(rep(0.8, 15), rep(0.2, m - 15))  # Scenario I   : Informative prior architecture
a.priori2 <- rep(0.5, m)                        # Scenario  II : Uniform prior architecture
a.priori3 <- c(rep(0.2, 15), rep(0.8, 15), rep(0.5, m-30)) # Scenario III: Misspecified prior architecture
a.priori4 <- c(rep(c(rep(0.2, 5), rep(0.5, 5), rep(0.8, 5)), 2), rep(0.5, m-30)) # Scenario IV : Mixed blocks prior architecture

sigma <- 1
seed <- 160323
set.seed(seed)
##--- Generate centered X
X <- Generate.X(n, r, m, f)
##
library(doParallel)
library(MASS) 
library(mvtnorm)

#1--- Check results per method
## Example with  Naive Wald in scenario I
###1.1- Fixed K
res_Naive_Wald_fixed_K <- Naive_Wald(X, Y, a.priori = a.priori1, K = 15, tau)
res_Naive_Wald_fixed_K$post.probs   # for posterior probablities 
res_Naive_Wald_fixed_K$indices.     # for selected indices 

###1.2- determination of K with criteria (BIC)
res_Naive_Wald_det_K <- Naive_Wald_Criteria(X, Y, a.priori = a.priori1, tau, threshold = 16, target_prob = 0.95, criteria = "BIC")
res_Naive_Wald_det_K$post.probs           # for posterior probablities 
res_Naive_Wald_det_K$indices.             # for selected indices 
res_Naive_Wald_det_K$r_squared_adj        # for R2_adj values     
res_Naive_Wald_det_K$bic_values           # for BIC values


#2--- all results 
## Example of Scenario I
###2.1- Fixed K
results_p_tau_fixed_K_I <- Simulation(X, p_values, tau_values,a.priori = a.priori1, alpha, K, B, S, threshold = 16)
saveRDS(results_p_tau_fixed_K_I, "results_p_tau_fixed_K_I_per_p_tau_B1000")
results_p_tau_fixed_K_I$main_results  # to get results in Table 2
superior_probabilty_results <- Compare_Methods_Matrix(results_p_tau_fixed_K_I$selection_summary) # to get results in Table 3
#ask for example: superior_probabilty_results$`p=0.2_tau=0.2`

###2.2- determination of K with criteria (BIC) 
results_p_tau_det_K_BIC <- Criteria_Simulation_Table(B, X, sigma, a.priori = a.priori1, p_values , tau_values, S,  threshold = 16, critere = "BIC")
results_p_tau_det_K_BIC. # For Table 4
