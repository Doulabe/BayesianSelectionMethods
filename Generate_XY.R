# ==============================================================================
# File: GenerateXY.R
# Author: N. K. Doulabe
# Description: Functions to generate synthetic Genotype (X) and Phenotype (Y) 
#              data for GWAS simulation studies. Includes power calculation 
#              for effect size calibration.
# ==============================================================================

# Required libraries
library(mvtnorm) # For rmvnorm
library(stats)   # For pnorm, qnorm, rnorm

# ==============================================================================
# 1. Genotype Generation
# ==============================================================================

#' Generate Genotype Matrix (X)
#'
#' @description Generates a centered and scaled genotype matrix with a specified 
#' correlation structure and Minor Allele Frequency (MAF).
#'
#' @param n Integer. Number of individuals (sample size).
#' @param r Numeric. Pairwise correlation coefficient between SNPs. 
#'          Note: This implementation creates a Compound Symmetry structure 
#'          (all SNPs correlated with all others by r).
#' @param m Integer. Number of SNPs (variables).
#' @param f Numeric. Minor Allele Frequency (MAF), e.g., 0.2.
#'
#' @return A centered and scaled matrix X of dimensions (n x m).
#' @export 
Generate.X <- function(n, r, m, f) {
  
  # 1. Create Correlation Matrix V (Compound Symmetry)
  # Diagonal = 1, Off-diagonal = r
  un <- rep(1, m)
  V <- r * un %*% t(un) + (1 - r) * diag(un) 
  
  # 2. Generate Latent Normal Variables
  # We generate two sets (Z1, Z2) to simulate diploidy (two alleles per SNP)
  Z1 <- rmvnorm(n, sigma = V)
  Z2 <- rmvnorm(n, sigma = V)
  
  # 3. Transform to Uniform Distribution
  U1 <- pnorm(Z1)
  U2 <- pnorm(Z2)
  
  # 4. Threshold to Discrete Genotypes {0, 1} per allele
  # If U < f, the allele is the minor allele (coded as 1)
  G1 <- matrix(as.numeric(U1 < f), ncol = m)
  G2 <- matrix(as.numeric(U2 < f), ncol = m)
  
  # 5. Combine to form Genotype {0, 1, 2}
  G <- G1 + G2
  
  # 6. Scale and Center the matrix
  # This is standard for LASSO/GWAS regression inputs
  X <- scale(G)
  return(X)
}

# ==============================================================================
# 2. Phenotype Generation
# ==============================================================================

#' Generate Phenotype Vector (Y)
#'
#' @description Generates continuous phenotypes based on a linear model:
#'              Y = Beta_0 + X*Beta + Epsilon
#'
#' @param n Integer. Sample size.
#' @param X Matrix. Genotype matrix (n x m).
#' @param beta Vector. Effect sizes. MUST have length m+1 (Index 1 is Intercept).
#' @param sigma Numeric. Standard deviation of the error term (noise).
#'
#' @return A vector Y of length n.
#' @export
Generate.Y <- function(n, X, beta, sigma) {
  # Add a column of 1s for the Intercept
  Design_Matrix <- cbind(rep(1, n), X)
  # Linear combination + Gaussian Noise
  Y <- Design_Matrix %*% beta + rnorm(n, mean = 0, sd = sigma)
  # Return as vector (drop matrix dimensions)
  return(Y[, 1])
}

# ==============================================================================
# 3. Power Calculation
# ==============================================================================

#' Power Function for Effect Size Calibration
#'
#' @description Calculates the difference between the theoretical power of a 
#'              Z-test and a target power 'p'. Used within uniroot() to find 
#'              the required effect size (beta) for a simulation scenario.
#'
#' @param beta Numeric. The effect size (signal strength).
#' @param n Integer. Sample size.
#' @param alpha Numeric. Significance level (Type I error).
#' @param p Numeric. Target probability (Power). Default is 0 (returns raw power).
#'
#' @return Numeric value. If p is provided, returns (Calculated_Power - Target_Power).
#' @export
Power_Fun <- function(beta, n, alpha, p = 0) {
  # Critical value for two-sided test
  z <- qnorm(1 - alpha / 2)
  # Power calculation using Normal distribution approximation
  # Power = P(Z > z - beta*sqrt(n)) + P(Z < -z - beta*sqrt(n))
  calc_power <- 1 - pnorm(z, mean = sqrt(n) * beta) + 
    pnorm(-z, mean = sqrt(n) * beta)
  # Return difference for root finding
  return(calc_power - p)
}