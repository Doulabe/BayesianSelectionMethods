# ==============================================================================
# File: Test_statistics_functions.R
# Author: N. K. Doulabe
# Description: Implements Bayesian Variable Selection algorithms for GWAS using
#              Vectorized QR decomposition. Includes Exact Wald, Asymptotic Score,
#              and Multiple-SNP joint tests.
# ==============================================================================

# Load necessary libraries (ensure these are installed via Ensure_packages function)
#Ensure_packages(required_libs)
library(MASS)
library(stats)
library(Matrix)
library(doParallel)
library(mvtnorm)
library(dplyr)

# ==============================================================================
# Section 1: Integrand Functions for Marginal Likelihoods
# ==============================================================================

#' Naive Wald Integrand
#' @description Computes the density for the marginal likelihood integration (Naive approach).
#' @param lambda Effect size parameter to integrate over.
#' @param Z Observed Z-statistic.
#' @param n Sample size.
#' @param rho Correlation between candidate and included variables.
#' @param tau Prior standard deviation for effect size.
Naive_Wald_Integrand <- function(lambda, Z, n, rho, tau){
  # Assumes standard normal likelihood for Z without variance correction
  dnorm(Z, mean = lambda * sqrt(n * (1 - rho))) * dnorm(lambda, sd = tau)
} 

#' Asymptotic Wald Integrand
#' @description Integrand for the Asymptotic Wald test with variance correction.
Asymptotic_Wald_Integrand <- function(lambda, Z, n, rho, tau){
  mean.Z <- lambda * sqrt(n * (1 - rho))
  # Variance inflation due to correlation (rho)
  var.Z <- 1 + 0.5 * (1 - rho) * lambda^2
  dnorm(Z, mean = mean.Z, sd = sqrt(var.Z)) * dnorm(lambda, sd = tau)
} 

#' Exact Wald Integrand
#' @description Integrand for the Exact Wald test using the non-central t-distribution.
#' @param k Number of currently included variables (degrees of freedom adjustment).
Exact_Wald_Integrand <- function(lambda, Z, n, rho, tau, k){
  # Uses non-central t-distribution (dt) for exact finite-sample inference
  dt(Z, df = n - (k + 1), ncp = lambda * sqrt(n * (1 - rho))) * dnorm(lambda, sd = tau)
}

#' Asymptotic Score Integrand
#' @description Integrand for the Score test statistic.
Asymptotic_Score_Integrand <- function(lambda, Q, n, rho, tau){
  mean.Q <- lambda * sqrt(n) * (1 - rho) / sqrt(1 + (1 - rho) * lambda^2)
  var.Q <- (2 * (1 - rho) + ((1 - rho)^2) * lambda^2) / (2 * (1 + (1 - rho) * lambda^2)^3)
  dnorm(Q, mean = mean.Q, sd = sqrt(var.Q)) * dnorm(lambda, sd = tau)
} 


# ==============================================================================
# Section 2: Monte Carlo Approximations for Multiple-SNP Tests
# ==============================================================================

#' Fisher MC Approximation (Exact Multiple Test)
#' @description Approximates the posterior probability for Joint Tests using 
#'              Monte Carlo simulation of the non-central Fisher distribution.
#' @param Q_W The observed F-statistic (or transformed statistic).
#' @param XtX_tilde Gram matrix of the candidate model.
#' @param XtX_k Gram matrix of the null model (already included variables).
#' @param S Number of Monte Carlo samples.
Fisher_MC_approximation <- function(Q_W, XtX_tilde, XtX_k, n, p, tau2, S) {
  stopifnot(is.matrix(XtX_tilde))
  rA  <- ncol(XtX_tilde)
  df1 <- rA
  df2 <- n - (rA + 1)
  
  # 1. Simulate prior effects (Lambda) to estimate non-centrality parameter (delta)
  Lambda_A <- matrix(rnorm(rA * S, sd = sqrt(tau2)), nrow = rA)
  delta_A  <- 0.5 * colSums(Lambda_A * (XtX_tilde %*% Lambda_A))
  delta_A  <- pmax(delta_A, 0) # Ensure non-negative
  
  # 2. Compute density under H1 (Model A)
  densA <- df(Q_W, df1 = df1, df2 = df2, ncp = delta_A)
  densA[!is.finite(densA)] <- 0
  A_hat <- mean(densA)
  
  # 3. Compute density under H0 (Model B)
  if (is.null(XtX_k) || ncol(XtX_k) == 0) {
    # If no previous variables, H0 is the central F-distribution
    densB <- df(Q_W, df1 = df1, df2 = df2, ncp = 0)
    densB <- if (is.finite(densB)) densB else 0
    B_hat <- densB
  } else {
    # If variables exist, simulate H0 non-centrality
    rB       <- ncol(XtX_k)
    Lambda_B <- matrix(rnorm(rB * S, sd = sqrt(tau2)), nrow = rB)
    delta_B  <- 0.5 * colSums(Lambda_B * (XtX_k %*% Lambda_B))
    delta_B  <- pmax(delta_B, 0)
    densB <- df(Q_W, df1 = df1, df2 = df2, ncp = delta_B)
    densB[!is.finite(densB)] <- 0
    B_hat <- mean(densB)
  }
  
  # 4. Compute Posterior Probability
  num <- p * A_hat
  den <- num +  (1 - p)* B_hat
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

#' Posterior Wrapper for Exact Multiple Test
#' @description Prepares matrices and calls the Fisher MC approximation.
ExM_posterior <- function(Q_W, X_tilde_jk, X_k, n, p, tau2, S) {
  XtX_tilde <- crossprod(X_tilde_jk)                 
  XtX_k     <- if (is.null(X_k)) NULL else crossprod(X_k) 
  Fisher_MC_approximation(Q_W = Q_W,
                          XtX_tilde = XtX_tilde, XtX_k = XtX_k,
                          n = n, p = p, tau2 = tau2, S = S)
}

#' Chi-Square MC Approximation (Naive Multiple Test)
#' @description Approximates posterior probability using Non-central Chi-square.
#'              Used for the Naive Multiple (NaM) method.
Chisq_MC_approximation <- function(Q_chi, XtX_tilde, XtX_k, n, p, tau2, S) {
  stopifnot(is.matrix(XtX_tilde))
  rA  <- ncol(XtX_tilde)                 
  df  <- rA
  
  # A-hat: Approximate density under Alternative Hypothesis
  Lambda_A <- matrix(rnorm(rA * S, sd = sqrt(tau2)), nrow = rA)
  delta_A  <- 0.5 * colSums(Lambda_A * (XtX_tilde %*% Lambda_A))
  delta_A  <- pmax(delta_A, 0)
  
  densA <- dchisq(Q_chi, df = df, ncp = delta_A)
  densA[!is.finite(densA)] <- 0
  A_hat <- mean(densA)
  
  # B-hat: Approximate density under Null Hypothesis
  if (is.null(XtX_k) || ncol(XtX_k) == 0L) {
    B_hat <- dchisq(Q_chi, df = df, ncp = 0)
    B_hat <- if (is.finite(B_hat)) B_hat else 0
  } else {
    rB       <- ncol(XtX_k)
    Lambda_B <- matrix(rnorm(rB * S, sd = sqrt(tau2)), nrow = rB)
    delta_B  <- 0.5 * colSums(Lambda_B * (XtX_k %*% Lambda_B))
    delta_B  <- pmax(delta_B, 0)
    densB <- dchisq(Q_chi, df = df, ncp = delta_B)
    densB[!is.finite(densB)] <- 0
    B_hat <- mean(densB)
  }
  
  # Posterior calculation
  num <- p * A_hat
  den <- num + (1 - p) * B_hat
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

NaM_posterior <- function(Q_chi, X_tilde_jk, X_k, n, p, tau2, S) {
  XtX_tilde <- crossprod(X_tilde_jk)                 
  XtX_k     <- if (is.null(X_k)) NULL else crossprod(X_k)
  Chisq_MC_approximation(Q_chi = Q_chi, XtX_tilde = XtX_tilde, XtX_k = XtX_k,
                         n = n, p = p, tau2 = tau2, S = S)
}


# ==============================================================================
# Section 3: Main Selection Algorithms (Single SNP)
# ==============================================================================

#' Naive Wald Test (NaW)
#' @description Performs stepwise variable selection using the Naive Wald statistic.
#' @param X Genotype matrix (n x m).
#' @param Y Phenotype vector (n x 1).
#' @param a.priori Vector of prior probabilities for each SNP.
#' @param K Maximum number of steps/variables to select.
#' @param tau Prior standard deviation.
Naive_Wald <- function(X, Y, a.priori, K, tau) {
  m <- ncol(X)
  n <- nrow(X)
  
  # Center data to simplify intercept handling
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  indices.inclus <- integer(0)
  indices_remain_to_test <- 1:m
  post.probs <- numeric(K)
  
  for (k in 1:K) {
    if (length(indices_remain_to_test) == 0) break
    
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    # --- Step 1: Compute Statistics ---
    if (k == 1) {
      # Base case: Correlation
      r <- cor(Y, X_cand)[1, ]
      r[is.na(r)] <- 0
      Z_vector <- ifelse(abs(r) == 1, sign(r) * Inf, r * sqrt(n - 2) / sqrt(1 - r^2))
      rho_vector <- numeric(length(indices_remain_to_test))
    } else {
      # Iterative case: Use QR decomposition to project out included variables
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Q <- qr.Q(qr_inclus)
      
      # Compute residuals (Vectorized for all candidates)
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      X_cand_res <- X_cand - Q %*% (t(Q) %*% X_cand)
      
      ss_X_cand_res <- colSums(X_cand_res^2)
      ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      
      beta_cand <- as.vector((t(Y_res) %*% X_cand_res) / ss_X_cand_res)
      
      # Variance estimation
      RSS_current <- sum(Y_res^2)
      RSS_new <- RSS_current - beta_cand^2 * ss_X_cand_res
      sigma2_new <- RSS_new / (n - k - 1)
      se_beta_cand <- sqrt(sigma2_new / ss_X_cand_res)
      Z_vector <- beta_cand / se_beta_cand
      
      # Calculate Rho (Correlation correction factor)
      XtX_inclus_inv <- chol2inv(qr.R(qr_inclus))
      V_inv <- n * XtX_inclus_inv
      R_matrix <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_matrix
      rho_vector <- rowSums(t(R_matrix) * t(M)) 
      rho_vector <- pmax(0, pmin(rho_vector, 1 - 1e-9))
    }
    
    Z_vector[!is.finite(Z_vector)] <- 0
    
    # --- Step 2: Integrate Marginal Likelihoods ---
    integrals <- mapply(function(z, r) {
      tryCatch(
        integrate(Naive_Wald_Integrand, lower = -Inf, upper = Inf, Z = z, n = n, rho = r, tau = tau)$value,
        error = function(e) NA
      )
    }, Z_vector, rho_vector)
    
    integrals[is.na(integrals)] <- 0
    
    # --- Step 3: Compute Posterior Probabilities ---
    densities <- dnorm(Z_vector)
    p_vector <- a.priori[indices_remain_to_test]
    numerator <- p_vector * integrals
    denominator <- numerator + (1 - p_vector) * densities
    a.posteriori <- numerator / (denominator + 1e-20) 
    
    if(all(!is.finite(a.posteriori))) break
    
    # Select best candidate
    max_idx_local <- which.max(a.posteriori)
    post.probs[k] <- a.posteriori[max_idx_local]
    a.inclure <- indices_remain_to_test[max_idx_local]
    indices.inclus <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }
  
  return(list(post.probs = post.probs, indices = indices.inclus))
}

#' Asymptotic Wald Test (AsW)
#' @description Performs selection using the Asymptotic Wald statistic.
#' @details Uses `Asymptotic_Wald_Integrand` which accounts for variance inflation.
Asymptotic_Wald <- function(X, Y, a.priori, K, tau) {
  m <- ncol(X)
  n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  indices.inclus <- integer(0)
  indices_remain_to_test <- 1:m
  post.probs <- numeric(K)
  
  for (k in 1:K) {
    if (length(indices_remain_to_test) == 0) break
    
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      r <- cor(Y, X_cand)[1, ]
      r[is.na(r)] <- 0
      Z_vector <- ifelse(abs(r) == 1, sign(r) * Inf, r * sqrt(n - 2) / sqrt(1 - r^2))
      rho_vector <- numeric(length(indices_remain_to_test))
    } else {
      # Standard Vectorized QR approach (see Naive_Wald for details)
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Q <- qr.Q(qr_inclus)
      
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      X_cand_res <- X_cand - Q %*% (t(Q) %*% X_cand)
      
      ss_X_cand_res <- colSums(X_cand_res^2)
      ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector((t(Y_res) %*% X_cand_res) / ss_X_cand_res)
      RSS_current <- sum(Y_res^2)
      RSS_new <- RSS_current - beta_cand^2 * ss_X_cand_res
      sigma2_new <- RSS_new / (n - k - 1)
      se_beta_cand <- sqrt(sigma2_new / ss_X_cand_res)
      Z_vector <- beta_cand / se_beta_cand
      
      # Rho Calculation
      XtX_inclus_inv <- chol2inv(qr.R(qr_inclus))
      V_inv <- n * XtX_inclus_inv
      R_matrix <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_matrix
      rho_vector <- rowSums(t(R_matrix) * t(M)) 
      rho_vector <- pmax(0, pmin(rho_vector, 1 - 1e-9))
    }
    
    Z_vector[!is.finite(Z_vector)] <- 0
    
    integrals <- mapply(function(z, r) {
      tryCatch(
        integrate(Asymptotic_Wald_Integrand, lower = -Inf, upper = Inf, Z = z, n = n, rho = r, tau = tau)$value,
        error = function(e) NA
      )
    }, Z_vector, rho_vector)
    
    integrals[is.na(integrals)] <- 0
    densities <- dnorm(Z_vector)
    p_vector <- a.priori[indices_remain_to_test]
    
    numerator <- p_vector * integrals
    denominator <- numerator + (1 - p_vector) * densities
    a.posteriori <- numerator / (denominator + 1e-20) 
    
    if(all(!is.finite(a.posteriori))) break
    max_idx_local <- which.max(a.posteriori)
    post.probs[k] <- a.posteriori[max_idx_local]
    a.inclure <- indices_remain_to_test[max_idx_local]
    indices.inclus <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }
  return(list(post.probs = post.probs, indices = indices.inclus))
}

#' Exact Wald Test (ExW)
#' @description The proposed method: Exact Wald test using Vectorized QR decomposition.
#' @details Uses `dt` (Student's t) for exact finite-sample inference, robust to high LD.
Exact_Wald <- function(X, Y, a.priori, K, tau) {
  m <- ncol(X)
  n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  indices.inclus <- integer(0)
  indices_remain_to_test <- 1:m
  post.probs <- numeric(K)
  
  for (k in 1:K) {
    if (length(indices_remain_to_test) == 0) break
    
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      r <- cor(Y, X_cand)[1, ]
      r[is.na(r)] <- 0
      Z_hat_vector <- ifelse(abs(r) == 1, sign(r) * Inf, r * sqrt(n - 2) / sqrt(1 - r^2))
      rho_vector <- numeric(length(indices_remain_to_test)) 
    } else {
      # QR Projection
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Q <- qr.Q(qr_inclus)
      
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      X_cand_res <- X_cand - Q %*% (t(Q) %*% X_cand)
      
      ss_X_cand_res <- colSums(X_cand_res^2)
      ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector((t(Y_res) %*% X_cand_res) / ss_X_cand_res)
      RSS_current <- sum(Y_res^2)
      RSS_new <- RSS_current - beta_cand^2 * ss_X_cand_res
      sigma2_new <- RSS_new / (n - k - 1)
      se_beta_cand <- sqrt(sigma2_new / ss_X_cand_res)
      Z_vector <- beta_cand / se_beta_cand
      
      # Rho Calculation
      XtX_inclus_inv <- chol2inv(qr.R(qr_inclus))
      V_inv <- n * XtX_inclus_inv
      R_matrix <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_matrix
      rho_vector <- rowSums(t(R_matrix) * t(M))
      rho_vector <- pmax(0, pmin(rho_vector, 1 - 1e-9))
      
      # --- Exact Correction ---
      # Z_hat is scaled by (1-rho) to align with t-distribution assumptions
      Z_hat_vector <- Z_vector * sqrt(1 - rho_vector)
    }
    
    Z_hat_vector[!is.finite(Z_hat_vector)] <- 0
    
    # Exact Wald Integration (using n-k-1 degrees of freedom)
    integrals <- mapply(function(z, r) {
      tryCatch(
        integrate(Exact_Wald_Integrand, lower = -Inf, upper = Inf, Z = z, n = n, rho = r, tau = tau, k = k)$value,
        error = function(e) NA
      )
    }, Z_hat_vector, rho_vector)
    
    integrals[is.na(integrals)] <- 0
    
    # Null distribution is t-distribution
    densities <- dt(Z_hat_vector, df = n - (k + 1))
    p_vector <- a.priori[indices_remain_to_test]
    
    numerator <- p_vector * integrals
    denominator <- numerator + (1 - p_vector) * densities
    a.posteriori <- numerator / (denominator + 1e-20) 
    
    if(all(!is.finite(a.posteriori))) break
    max_idx_local <- which.max(a.posteriori)
    post.probs[k] <- a.posteriori[max_idx_local]
    a.inclure <- indices_remain_to_test[max_idx_local]
    indices.inclus <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }
  return(list(post.probs = post.probs, indices = indices.inclus))
}

#' Asymptotic Score Test (AsS)
#' @description Performs variable selection using the efficient Score statistic.
#' @details Does not require recalculating model coefficients for candidates, making it faster.
Asymptotic_Score <- function(X, Y, a.priori, K, tau) {
  m <- ncol(X)
  n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  indices.inclus <- integer(0)
  indices_remain_to_test <- 1:m
  post.probs <- numeric(K)
  
  for (k in 1:K) {
    if (length(indices_remain_to_test) == 0) break
    
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      numerator <- as.vector(t(Y) %*% X_cand)
      beta <- numerator / colSums(X_cand^2)
      RSS <- sum(Y^2) - beta * numerator
      sigma_vec <- sqrt(RSS / (n - 2))
      Q_hat_vector <- numerator / (sqrt(n) * sigma_vec)
      rho_vector <- numeric(length(indices_remain_to_test))
    } else {
      # Use QR to work on residuals
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Q <- qr.Q(qr_inclus)
      
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      X_cand_res <- X_cand - Q %*% (t(Q) %*% X_cand)
      
      numerator <- as.vector(t(Y) %*% X_cand_res) # Score numerator
      
      # Calculate sigma for every candidate model
      ss_X_cand_res <- colSums(X_cand_res^2)
      ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector((t(Y) %*% X_cand_res) / ss_X_cand_res)
      RSS_current <- sum(Y_res^2)
      RSS_new <- RSS_current - beta_cand^2 * ss_X_cand_res
      sigma_new_vec <- sqrt(RSS_new / (n - k - 1))
      
      Q_hat_vector <- numerator / (sqrt(n) * sigma_new_vec)
      
      # Rho Calculation
      XtX_inclus_inv <- chol2inv(qr.R(qr_inclus))
      V_inv <- n * XtX_inclus_inv
      R_matrix <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_matrix
      rho_vector <- rowSums(t(R_matrix) * t(M))
      rho_vector <- pmax(0, pmin(rho_vector, 1 - 1e-9))
    }
    
    Q_hat_vector[!is.finite(Q_hat_vector)] <- 0
    
    integrals <- mapply(function(q, r) {
      tryCatch(
        integrate(Asymptotic_Score_Integrand, lower = -Inf, upper = Inf, Q = q, n = n, rho = r, tau = tau)$value,
        error = function(e) NA
      )
    }, Q_hat_vector, rho_vector)
    
    integrals[is.na(integrals)] <- 0
    densities <- dnorm(Q_hat_vector)
    p_vector <- a.priori[indices_remain_to_test]
    
    numerator <- p_vector * integrals
    denominator <- numerator + (1 - p_vector) * densities
    a.posteriori <- numerator / (denominator + 1e-20)
    
    if(all(!is.finite(a.posteriori))) break
    max_idx_local <- which.max(a.posteriori)
    post.probs[k] <- a.posteriori[max_idx_local]
    a.inclure <- indices_remain_to_test[max_idx_local]
    indices.inclus <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }
  return(list(post.probs = post.probs, indices = indices.inclus))
}


# ==============================================================================
# Section 4: Multiple-SNP (Joint) Tests
# ==============================================================================

#' Naive Multiple Test (NaM)
#' @description Joint test using Chi-squared approximations.
Naive_Multiple <- function(X, Y, a.priori, K, tau2, S) {
  m <- ncol(X)
  n <- nrow(X)
  
  # Center data to make no-intercept model equivalent
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  indices.inclus <- integer(0)
  indices_remain_to_test <- seq_len(m)
  post.probs <- rep(NA_real_, K)
  TSS <- sum(Y^2)
  
  for (k in 1:K) {
    if (length(indices_remain_to_test) == 0) break
    
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    rA <- k 
    df1 <- rA
    df2 <- n - rA - 1L
    if (df2 <= 0) break
    
    # --- Vectorized F-statistic calculation ---
    if (k == 1) {
      r_sq <- (as.vector(t(Y) %*% X_cand) / sqrt(TSS * colSums(X_cand^2)))^2
      r_sq[is.na(r_sq)] <- 0
      F_block_vector <- r_sq * (n - 2) / (1 - r_sq)
      X_inclus <- NULL
    } else {
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      
      Y_res <- qr.resid(qr_inclus, Y)
      X_cand_res <- qr.resid(qr_inclus, X_cand)
      
      ss_X_cand_res <- colSums(X_cand_res^2)
      ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector(t(Y_res) %*% X_cand_res) / ss_X_cand_res
      RSS_current <- sum(Y_res^2)
      RSS_new_vector <- RSS_current - beta_cand^2 * ss_X_cand_res
      SSR_new_vector <- TSS - RSS_new_vector
      
      MSE_new_vector <- RSS_new_vector / df2
      MSR_new_vector <- SSR_new_vector / df1
      F_block_vector <- MSR_new_vector / MSE_new_vector
    }
    
    # Convert F-stats to Chi-squared stats for approximation
    Q_chi_vector <- rA * F_block_vector
    Q_chi_vector[!is.finite(Q_chi_vector) | Q_chi_vector < 0] <- -1
    
    # Monte Carlo posterior calculation
    p_vector <- a.priori[indices_remain_to_test]
    a.posteriori <- mapply(
      function(q_chi, j, p) {
        if (q_chi < 0) return(NA_real_)
        X.modele <- if (is.null(X_inclus)) X[, j, drop = FALSE] else cbind(X_inclus, X[, j, drop = FALSE])
        NaM_posterior(Q_chi = q_chi, X_tilde_jk = X.modele, X_k = X_inclus, n = n, p = p, tau2 = tau2, S = S)
      },
      Q_chi_vector, indices_remain_to_test, p_vector
    )
    
    if (all(!is.finite(a.posteriori))) break
    best_idx_local <- which.max(a.posteriori)
    post.probs[k] <- a.posteriori[best_idx_local]
    a.inclure <- indices_remain_to_test[best_idx_local]
    indices.inclus <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }
  list(post.probs = post.probs, indices = indices.inclus)
}

#' Exact Multiple Test (ExM)
#' @description Joint test using Exact Fisher distributions (via MC approximation).
Exact_Multiple <- function(X, Y, a.priori, K, tau2, S) {
  m <- ncol(X)
  n <- nrow(X)
  
  # Center data to make no-intercept model equivalent to intercept model
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  indices.inclus <- integer(0)
  indices_remain_to_test <- seq_len(m)
  post.probs <- rep(NA_real_, K)
  TSS <- sum(Y^2)
  
  for (k in 1:K) {
    if (length(indices_remain_to_test) == 0) break
    
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    rA <- k
    df1 <- rA
    df2 <- n - rA - 1L
    if (df2 <= 0) break
    
    if (k == 1) {
      r_sq <- (as.vector(t(Y) %*% X_cand) / sqrt(TSS * colSums(X_cand^2)))^2
      r_sq[is.na(r_sq)] <- 0
      Z_vector <- r_sq * (n - 2) / (1 - r_sq)
      X_inclus <- NULL
    } else {
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      
      Y_res <- qr.resid(qr_inclus, Y)
      X_cand_res <- qr.resid(qr_inclus, X_cand)
      
      ss_X_cand_res <- colSums(X_cand_res^2)
      ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector(t(Y_res) %*% X_cand_res) / ss_X_cand_res
      RSS_current <- sum(Y_res^2)
      RSS_new_vector <- RSS_current - beta_cand^2 * ss_X_cand_res
      SSR_new_vector <- TSS - RSS_new_vector
      MSE_new_vector <- RSS_new_vector / df2
      MSR_new_vector <- SSR_new_vector / df1
      Z_vector <- MSR_new_vector / MSE_new_vector
    }
    
    Z_vector[!is.finite(Z_vector) | Z_vector < 0] <- -1
    p_vector <- a.priori[indices_remain_to_test]
    
    # Exact Posterior using Fisher MC
    a.posteriori <- mapply(
      function(z, j, p) {
        if (z < 0) return(NA_real_)
        X.modele <- if (is.null(X_inclus)) X[, j, drop = FALSE] else cbind(X_inclus, X[, j, drop = FALSE])
        ExM_posterior(Q_W = z, X_tilde_jk = X.modele, X_k = X_inclus, n = n, p = p, tau2 = tau2, S = S)
      },
      Z_vector, indices_remain_to_test, p_vector
    )
    
    if (all(!is.finite(a.posteriori))) break
    best_idx_local <- which.max(a.posteriori)
    post.probs[k] <- a.posteriori[best_idx_local]
    a.inclure <- indices_remain_to_test[best_idx_local]
    indices.inclus <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }
  list(post.probs = post.probs, indices = indices.inclus)
}


# ==============================================================================
# Section 5: Simulation and Evaluation Framework
# ==============================================================================

#' Simulation Wrapper
#' @description Runs parallel simulations to compare Type I error and Power across methods.
#' @param B Number of simulation iterations per setting.
#' @param S Number of Monte Carlo samples for multiple tests.
#' @param threshold Cutoff for selected variables (default 16).
Simulation <- function(X, p_values, tau_values, a.priori, alpha, K, B, S, threshold = 16) {
  m <- ncol(X)
  n <- nrow(X)
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Load dependencies on worker nodes
  clusterEvalQ(cl, {
    library(doParallel)
    library(MASS)
    library(mvtnorm)
    library(stats)
  })
  
  # Export necessary functions to worker nodes
  clusterExport(cl, c(
    "Generate.Y", "Naive_Wald", "Asymptotic_Wald", "Exact_Wald", "Asymptotic_Score", "Exact_Multiple", 
    "Naive_Multiple", "Power_Fun", "Naive_Wald_Integrand", "Asymptotic_Wald_Integrand",
    "Exact_Wald_Integrand", "Asymptotic_Score_Integrand", "Fisher_MC_approximation",
    "Chisq_MC_approximation", "ExM_posterior", "NaM_posterior"
  ))
  
  param_grid <- expand.grid(p = p_values, tau = tau_values, stringsAsFactors = FALSE)
  
  final_main_results <- list()
  final_selection_summaries <- list()
  method_names <- c("NaW", "AsW", "ExW", "AsS", "NaM", "ExM")
  
  for (i in seq_len(nrow(param_grid))) {
    p <- param_grid$p[i]
    tau <- param_grid$tau[i]
    
    # Determine effect size (b_val) for target power
    b_val <- uniroot(Power_Fun, lower = 0, upper = 1, n = n, alpha = alpha, p = p)$root
    beta <- c(1, rep(b_val, 15), rep(0, m - 15))
    sigma <- 1
    
    # Parallel Loop
    mc_results <- foreach(b = 1:B, .combine = rbind, .packages = c("stats")) %dopar% {
      Y <- Generate.Y(n, X, beta, sigma)
      
      methods_results <- list(
        NaW   = Naive_Wald(X, Y, a.priori, K, tau),
        AsW   = Asymptotic_Wald(X, Y, a.priori, K, tau), 
        ExW   = Exact_Wald(X, Y, a.priori, K, tau),
        AsS   = Asymptotic_Score(X, Y, a.priori, K, tau), 
        NaM   = Naive_Multiple(X, Y, a.priori, K, tau^2, S),
        ExM   = Exact_Multiple(X, Y, a.priori, K, tau^2, S)
      )
      
      c(
        sapply(methods_results, function(res) sum(res$indices < threshold)),
        sapply(methods_results, function(res) sum(res$indices %in% 1:15))
      )
    }
    
    # Aggregating Results
    mean_counts <- colMeans(mc_results[, 1:6])
    mean_results <- data.frame(p = p, tau = tau)
    for (j in seq_along(method_names)) {
      mean_results[[paste0("mean_", method_names[j])]] <- mean_counts[j]
    }
    final_main_results[[i]] <- mean_results
    
    # Aggregating Selection Summaries
    for (j in seq_along(method_names)) {
      method_col <- mc_results[, 6 + j] 
      summary <- tabulate(method_col + 1, nbins = 16) / B
      summary_df <- as.data.frame(t(summary))
      colnames(summary_df) <- paste0("choix.", 0:15)
      summary_df$p <- p
      summary_df$tau <- tau
      summary_df$method <- method_names[j]
      final_selection_summaries[[length(final_selection_summaries) + 1]] <- summary_df
    }
  }
  
  stopCluster(cl)
  
  main_results <- do.call(rbind, final_main_results)
  selection_summary <- do.call(rbind, final_selection_summaries)
  row.names(main_results) <- NULL
  row.names(selection_summary) <- NULL
  
  return(list(main_results = main_results, selection_summary = selection_summary))
}

#' Calculate Probability of Superiority
#' @description Computes the probability that Method X selects more true variables than Method Y.
Get_Superior_Probability <- function(x, y) {
  sum(sapply(0:15, function(i) {
    xi <- x[i + 1] 
    sum(sapply(i:15, function(j) {
      yj <- y[j + 1]
      weight <- if (i == j) 0.5 else 1
      xi * yj * weight
    }))
  }))
}

#' Comparison Matrix
#' @description Generates pairwise comparison matrices for all methods.
Compare_Methods_Matrix <- function(selection_summary) {
  methods <- unique(selection_summary$method)
  param_groups <- unique(selection_summary[, c("p", "tau")])
  all_matrices <- list()
  
  for (i in seq_len(nrow(param_groups))) {
    p_val <- param_groups$p[i]
    tau_val <- param_groups$tau[i]
    
    subset_data <- selection_summary %>%
      filter(p == p_val, tau == tau_val)
    
    matrix_result <- matrix(NA, nrow = length(methods), ncol = length(methods))
    rownames(matrix_result) <- methods
    colnames(matrix_result) <- methods
    
    for (m1 in methods) {
      for (m2 in methods) {
        if (m1 != m2) {
          x <- as.numeric(subset_data %>% filter(method == m1) %>% select(starts_with("choix.")))
          y <- as.numeric(subset_data %>% filter(method == m2) %>% select(starts_with("choix.")))
          prob <- Get_Superior_Probability(x, y)
          matrix_result[m1, m2] <- round(prob, 3)
        }
      }
    }
    all_matrices[[paste0("p=", p_val, "_tau=", tau_val)]] <- matrix_result
  }
  return(all_matrices)
}