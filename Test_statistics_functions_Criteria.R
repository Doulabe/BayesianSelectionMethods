# ==============================================================================
# File: Test_statistics_functions_Criteria.R
# Author: N. K. Doulabe
# Description: Implements Variable Selection algorithms with automatic stopping 
#              criteria (BIC, AIC, Adj-R2, Posterior Probability).
#              Contains two sections:
#              A. Functions optimized for Simulation Studies (returning metrics).
#              B. Functions optimized for Real Data Analysis (handling thresholds).
# ==============================================================================

source("Test_statistics_functions.R") # Loads integrands and helper functions

# ==============================================================================
# SECTION A: Criteria Functions for Simulation Study (Table 4)
# ==============================================================================

#' Naive Wald with Stopping Criteria (Simulation)
#' @description Performs stepwise selection, calculating Information Criteria at each step.
#' @param criteria Stopping metric: "BIC", "AIC", "R2" (Adjusted R-squared), or "Posterior".
#' @param threshold Cutoff index to count "True Positives" (e.g., indices < 16).
#' @return List containing selected indices, metric history, and count of true positives.
Naive_Wald_Criteria <- function(X, Y, a.priori, tau, threshold = 16, target_prob = 0.95, criteria = "BIC") {
  
  # --- 1. Validation & Setup ---
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria invalid. Must be one of: 'R2', 'AIC', 'BIC', 'Posterior'.")
  }
  
  m <- ncol(X)
  n <- nrow(X)
  
  # Center Data (Standardize scale)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X) 
  
  # Initialize containers
  indices.inclus <- integer(0)
  indices_remain_to_test <- 1:m
  
  r_squared_adj_hist <- numeric(0)
  aic_hist <- numeric(0)
  bic_hist <- numeric(0)
  post.probs <- numeric(0)
  
  # Total Sum of Squares (TSS) for R2 calculation
  TSS <- sum(Y^2)
  RSS_current <- TSS # Start with Null model RSS
  
  # --- 2. Main Selection Loop ---
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    # --- A. Compute Statistics (Vectorized) ---
    if (k == 1) {
      # Base case: Simple Correlation
      nums <- as.vector(t(Y) %*% X_cand)
      dens <- sqrt(TSS * colSums(X_cand^2))
      r_vec <- nums / (dens + 1e-12)
      
      Z_vector <- r_vec * sqrt(n - 2) / sqrt(1 - r_vec^2)
      rho_vector <- numeric(length(indices_remain_to_test)) 
      
      # Store sum of squares for RSS update
      ss_X_cand <- colSums(X_cand^2)
      beta_cand <- nums / ss_X_cand
    } else {
      # Iterative case: Vectorized QR
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Q <- qr.Q(qr_inclus)
      
      # Residuals of Y and Candidates
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      XtX_cand <- t(Q) %*% X_cand
      X_cand_res <- X_cand - Q %*% XtX_cand
      
      ss_X_cand_res <- colSums(X_cand_res^2)
      ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      
      beta_cand <- as.vector((t(Y_res) %*% X_cand_res) / ss_X_cand_res)
      
      # Provisional RSS for Z-score variance
      RSS_provisional <- sum(Y_res^2) - beta_cand^2 * ss_X_cand_res
      sigma2_new <- RSS_provisional / (n - k - 1)
      se_beta <- sqrt(sigma2_new / ss_X_cand_res)
      Z_vector <- beta_cand / se_beta
      
      # Rho Calculation
      XtX_inclus_inv <- chol2inv(qr.R(qr_inclus))
      V_inv <- n * XtX_inclus_inv
      R_matrix <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_matrix
      rho_vector <- rowSums(t(R_matrix) * t(M))
      rho_vector <- pmax(0, pmin(rho_vector, 1 - 1e-9))
    }
    
    Z_vector[!is.finite(Z_vector)] <- 0
    
    # --- B. Posterior Probability ---
    sd_eff <- sqrt(1 + n * (1 - rho_vector) * tau^2)
    integrals <- dnorm(Z_vector, mean = 0, sd = sd_eff)
    densities <- dnorm(Z_vector, mean = 0, sd = 1)
    p_vec <- a.priori[indices_remain_to_test]
    numerator <- p_vec * integrals
    denominator <- numerator + (1 - p_vec) * densities
    a.posteriori <- numerator / (denominator + 1e-20)
    
    # Identify best candidate
    max_idx_local <- which.max(a.posteriori)
    current_max_prob <- a.posteriori[max_idx_local]
    best_var_global <- indices_remain_to_test[max_idx_local]
    
    # --- C. Calculate Metrics for the Potential New Model ---
    best_beta <- beta_cand[max_idx_local]
    
    if (k == 1) {
      best_ss_res <- colSums(X_cand^2)[max_idx_local]
    } else {
      best_ss_res <- colSums(X_cand_res^2)[max_idx_local]
    }
    
    RSS_new <- RSS_current - (best_beta^2 * best_ss_res)
    if(RSS_new < 0) RSS_new <- 1e-9 
    
    df_model <- k + 1 
    
    # Compute Criteria
    r2_adj_val <- 1 - (RSS_new / (n - df_model)) / (TSS / (n - 1))
    aic_val <- n * log(RSS_new / n) + 2 * df_model
    bic_val <- n * log(RSS_new / n) + log(n) * df_model
    
    # --- D. Stopping Criteria Check ---
    stop_flag <- FALSE
    
    if (k > 1) {
      if (criteria == "R2") {
        # Stop if Adjusted R2 decreases (worse fit)
        if (r2_adj_val < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      } else if (criteria == "AIC") {
        # Stop if AIC increases
        if (aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      } else if (criteria == "BIC") {
        # Stop if BIC increases (penalty outweighs fit)
        if (bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      } else if (criteria == "Posterior") {
        # Stop if max probability is below confidence threshold
        if (current_max_prob < target_prob) stop_flag <- TRUE
      }
    }
    
    if (stop_flag) {
      break # Stop BEFORE adding the variable
    }
    
    # --- E. Update Lists ---
    indices.inclus <- c(indices.inclus, best_var_global)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var_global)
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj_val)
    aic_hist <- c(aic_hist, aic_val)
    bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, current_max_prob)
    RSS_current <- RSS_new 
  }
  
  # For Simulation: Count how many selected variables are 'True' (indices < threshold)
  indices_under_16 <- sum(indices.inclus < threshold) 
  
  return(list(
    indices = indices.inclus,
    post.probs = post.probs,
    r_squared_adj = r_squared_adj_hist,
    aic_values = aic_hist,
    bic_values = bic_hist,
    indices_under_16 = indices_under_16,
    K = length(indices.inclus)
  ))
}

#' Asymptotic Wald with Stopping Criteria (Simulation)
Asymptotic_Wald_Criteria <- function(X, Y, a.priori, tau, threshold = 16, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) stop("Criteria invalid.")
  
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  indices.inclus <- integer(0)
  indices_remain_to_test <- 1:m
  r_squared_adj_hist <- numeric(0)
  aic_hist <- numeric(0)
  bic_hist <- numeric(0)
  post.probs <- numeric(0)
  TSS <- sum(Y^2)
  RSS_current <- TSS
  
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      nums <- as.vector(t(Y) %*% X_cand)
      ss_X_cand <- colSums(X_cand^2)
      beta_cand <- nums / ss_X_cand
      RSS_prov <- TSS - beta_cand^2 * ss_X_cand
      sigma2 <- RSS_prov / (n - 2)
      Z_vec <- beta_cand / sqrt(sigma2 / ss_X_cand)
      rho_vec <- numeric(length(indices_remain_to_test))
    } else {
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Q <- qr.Q(qr_inclus)
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      XtX_cand <- t(Q) %*% X_cand
      X_cand_res <- X_cand - Q %*% XtX_cand
      ss_X_cand_res <- colSums(X_cand_res^2)
      ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector((t(Y_res) %*% X_cand_res) / ss_X_cand_res)
      RSS_prov <- sum(Y_res^2) - beta_cand^2 * ss_X_cand_res
      sigma2 <- RSS_prov / (n - k - 1)
      Z_vec <- beta_cand / sqrt(sigma2 / ss_X_cand_res)
      V_inv <- n * chol2inv(qr.R(qr_inclus))
      R_mat <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_mat
      rho_vec <- colSums(R_mat * M) 
      rho_vec <- pmax(0, pmin(rho_vec, 1 - 1e-9))
    }
    
    Z_vec[!is.finite(Z_vec)] <- 0
    integrals <- mapply(function(z, r) {
      tryCatch(integrate(Asymptotic_Wald_Integrand, lower=-Inf, upper=Inf, Z=z, n=n, rho=r, tau=tau)$value, error=function(e) 0)
    }, Z_vec, rho_vec)
    
    densities <- dnorm(Z_vec)
    p_vec <- a.priori[indices_remain_to_test]
    numerator <- p_vec * integrals
    denominator <- numerator + (1 - p_vec) * densities
    post <- numerator / (denominator + 1e-20)
    best_idx_local <- which.max(post)
    best_prob <- post[best_idx_local]
    best_var <- indices_remain_to_test[best_idx_local]
    
    # Criteria Calculation
    if(k==1) ss_res_best <- colSums(X_cand^2)[best_idx_local]
    else ss_res_best <- colSums(X_cand_res^2)[best_idx_local]
    
    RSS_new <- RSS_current - (beta_cand[best_idx_local]^2 * ss_res_best)
    if(RSS_new < 0) RSS_new <- 1e-9
    df <- k + 1
    r2_adj <- 1 - (RSS_new/(n-df)) / (TSS/(n-1))
    aic_val <- n * log(RSS_new/n) + 2 * df
    bic_val <- n * log(RSS_new/n) + log(n) * df
    
    stop_flag <- FALSE
    if (k > 1) {
      if (criteria == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (criteria == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (criteria == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (criteria == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    indices.inclus <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj)
    aic_hist <- c(aic_hist, aic_val)
    bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }
  
  indices_under_16 <- sum(indices.inclus < threshold)
  return(list(indices = indices.inclus,
              post.probs = post.probs,
              r_squared_adj = r_squared_adj_hist,
              aic_values = aic_hist,
              bic_values = bic_hist,
              indices_under_16 = indices_under_16,
              K = length(indices.inclus)))
}

#' Exact Wald with Stopping Criteria (Simulation)
Exact_Wald_Criteria <- function(X, Y, a.priori, tau, threshold = 16, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) stop("Criteria invalid.")
  
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  indices.inclus <- integer(0); indices_remain_to_test <- 1:m
  r_squared_adj_hist <- c(); aic_hist <- c(); bic_hist <- c(); post.probs <- c()
  TSS <- sum(Y^2); RSS_current <- TSS
  
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      nums <- as.vector(t(Y) %*% X_cand)
      ss_X_cand <- colSums(X_cand^2)
      beta_cand <- nums / ss_X_cand
      RSS_prov <- TSS - beta_cand^2 * ss_X_cand
      sigma2 <- RSS_prov / (n - 2)
      Z_vec <- beta_cand / sqrt(sigma2 / ss_X_cand)
      rho_vec <- numeric(length(indices_remain_to_test))
      Z_hat_vec <- Z_vec
    } else {
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus); Q <- qr.Q(qr_inclus)
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      X_cand_res <- X_cand - Q %*% (t(Q) %*% X_cand)
      ss_X_cand_res <- colSums(X_cand_res^2); ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector((t(Y_res) %*% X_cand_res) / ss_X_cand_res)
      RSS_prov <- sum(Y_res^2) - beta_cand^2 * ss_X_cand_res
      sigma2 <- RSS_prov / (n - k - 1)
      Z_vec <- beta_cand / sqrt(sigma2 / ss_X_cand_res)
      
      V_inv <- n * chol2inv(qr.R(qr_inclus))
      R_mat <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_mat
      rho_vec <- colSums(R_mat * M)
      rho_vec <- pmax(0, pmin(rho_vec, 1 - 1e-9))
      Z_hat_vec <- Z_vec * sqrt(1 - rho_vec)
    }
    
    Z_hat_vec[!is.finite(Z_hat_vec)] <- 0
    integrals <- mapply(function(z, r) {
      tryCatch(integrate(Exact_Wald_Integrand, lower=-Inf, upper=Inf, Z=z, n=n, rho=r, tau=tau, k=k)$value, error=function(e) 0)
    }, Z_hat_vec, rho_vec)
    
    densities <- dt(Z_hat_vec, df = n - (k + 1))
    p_vec <- a.priori[indices_remain_to_test]
    numerator <- p_vec * integrals
    denominator <- numerator + (1 - p_vec) * densities
    post <- numerator / (denominator + 1e-20)
    best_idx_local <- which.max(post); best_prob <- post[best_idx_local]
    
    if(k==1) ss_res_best <- colSums(X_cand^2)[best_idx_local]
    else ss_res_best <- colSums(X_cand_res^2)[best_idx_local]
    RSS_new <- RSS_current - (beta_cand[best_idx_local]^2 * ss_res_best)
    if(RSS_new < 0) RSS_new <- 1e-9
    df <- k + 1
    r2_adj <- 1 - (RSS_new/(n-df)) / (TSS/(n-1))
    aic_val <- n * log(RSS_new/n) + 2 * df
    bic_val <- n * log(RSS_new/n) + log(n) * df
    
    stop_flag <- FALSE
    if (k > 1) {
      if (criteria == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (criteria == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (criteria == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (criteria == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    indices.inclus <- c(indices.inclus, indices_remain_to_test[best_idx_local])
    indices_remain_to_test <- setdiff(indices_remain_to_test, indices_remain_to_test[best_idx_local])
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj); aic_hist <- c(aic_hist, aic_val); bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }
  indices_under_16 <- sum(indices.inclus < threshold)
  return(list(indices = indices.inclus, post.probs = post.probs, r_squared_adj = r_squared_adj_hist, aic_values = aic_hist, bic_values = bic_hist, indices_under_16 = indices_under_16, K = length(indices.inclus)))
}

#' Asymptotic Score with Stopping Criteria (Simulation)
Asymptotic_Score_Criteria <- function(X, Y, a.priori, tau, threshold = 16, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) stop("Criteria invalid.")
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  indices.inclus <- integer(0); indices_remain_to_test <- 1:m
  r_squared_adj_hist <- c(); aic_hist <- c(); bic_hist <- c(); post.probs <- c()
  TSS <- sum(Y^2); RSS_current <- TSS
  
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      numerator <- as.vector(t(Y) %*% X_cand)
      ss_X <- colSums(X_cand^2)
      beta <- numerator / ss_X
      RSS <- TSS - beta * numerator
      sigma_vec <- sqrt(RSS / (n - 2))
      Q_hat <- numerator / (sqrt(n) * sigma_vec)
      rho_vec <- numeric(length(indices_remain_to_test))
      beta_cand <- beta
    } else {
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus); Q <- qr.Q(qr_inclus)
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      X_cand_res <- X_cand - Q %*% (t(Q) %*% X_cand)
      numerator <- as.vector(t(Y) %*% X_cand_res)
      ss_X_cand_res <- colSums(X_cand_res^2); ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector((t(Y) %*% X_cand_res) / ss_X_cand_res)
      RSS_prov <- sum(Y_res^2) - beta_cand^2 * ss_X_cand_res
      sigma_new_vec <- sqrt(RSS_prov / (n - k - 1))
      Q_hat <- numerator / (sqrt(n) * sigma_new_vec)
      
      V_inv <- n * chol2inv(qr.R(qr_inclus))
      R_mat <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_mat
      rho_vec <- colSums(R_mat * M)
      rho_vec <- pmax(0, pmin(rho_vec, 1 - 1e-9))
    }
    
    Q_hat[!is.finite(Q_hat)] <- 0
    integrals <- mapply(function(q, r) {
      tryCatch(integrate(Asymptotic_Score_Integrand, lower=-Inf, upper=Inf, Q=q, n=n, rho=r, tau=tau)$value, error=function(e) 0)
    }, Q_hat, rho_vec)
    
    densities <- dnorm(Q_hat)
    p_vec <- a.priori[indices_remain_to_test]
    numerator <- p_vec * integrals
    denominator <- numerator + (1 - p_vec) * densities
    post <- numerator / (denominator + 1e-20)
    best_idx_local <- which.max(post); best_prob <- post[best_idx_local]
    
    if(k==1) ss_res_best <- colSums(X_cand^2)[best_idx_local]
    else ss_res_best <- colSums(X_cand_res^2)[best_idx_local]
    RSS_new <- RSS_current - (beta_cand[best_idx_local]^2 * ss_res_best)
    if(RSS_new < 0) RSS_new <- 1e-9
    df <- k + 1
    r2_adj <- 1 - (RSS_new/(n-df)) / (TSS/(n-1))
    aic_val <- n * log(RSS_new/n) + 2 * df
    bic_val <- n * log(RSS_new/n) + log(n) * df
    
    stop_flag <- FALSE
    if (k > 1) {
      if (criteria == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (criteria == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (criteria == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (criteria == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    indices.inclus <- c(indices.inclus, indices_remain_to_test[best_idx_local])
    indices_remain_to_test <- setdiff(indices_remain_to_test, indices_remain_to_test[best_idx_local])
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj); aic_hist <- c(aic_hist, aic_val); bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }
  indices_under_16 <- sum(indices.inclus < threshold)
  return(list(indices = indices.inclus, post.probs = post.probs, r_squared_adj = r_squared_adj_hist, aic_values = aic_hist, bic_values = bic_hist, indices_under_16 = indices_under_16, K = length(indices.inclus)))
}

#' Exact Multiple Test (Fisher) with Stopping Criteria
Exact_Multiple_Criteria <- function(X, Y, a.priori, tau2, S, threshold = 16, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) stop("Criteria invalid.")
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  indices.inclus <- integer(0); indices_remain_to_test <- 1:m
  r_squared_adj_hist <- c(); aic_hist <- c(); bic_hist <- c(); post.probs <- c()
  TSS <- sum(Y^2); RSS_current <- TSS
  X_inclus <- NULL
  
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    rA <- k; df1 <- rA; df2 <- n - rA - 1
    if (df2 <= 0) break
    
    if (k == 1) {
      r_sq <- (as.vector(t(Y) %*% X_cand) / sqrt(TSS * colSums(X_cand^2)))^2
      r_sq[is.na(r_sq)] <- 0
      Z_vector <- r_sq * (n - 2) / (1 - r_sq)
      nums <- as.vector(t(Y) %*% X_cand)
      beta_cand <- nums / colSums(X_cand^2)
    } else {
      if(is.null(X_inclus)) X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Y_res <- qr.resid(qr_inclus, Y)
      X_cand_res <- qr.resid(qr_inclus, X_cand)
      ss_X_cand_res <- colSums(X_cand_res^2); ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector(t(Y_res) %*% X_cand_res) / ss_X_cand_res
      RSS_prov <- sum(Y_res^2)
      RSS_new_vec <- RSS_prov - beta_cand^2 * ss_X_cand_res
      SSR_new_vec <- TSS - RSS_new_vec
      Z_vector <- (SSR_new_vec / df1) / (RSS_new_vec / df2)
    }
    
    Z_vector[!is.finite(Z_vector) | Z_vector < 0] <- -1
    p_vec <- a.priori[indices_remain_to_test]
    
    post <- mapply(function(z, j, p) {
      if (z < 0) return(NA)
      X_mod <- if(is.null(X_inclus)) X[, j, drop=FALSE] else cbind(X_inclus, X[, j, drop=FALSE])
      ExM_posterior(Q_W=z, X_tilde_jk=X_mod, X_k=X_inclus, n=n, p=p, tau2=tau2, S=S)
    }, Z_vector, indices_remain_to_test, p_vec)
    
    if (all(!is.finite(post))) break
    best_idx_local <- which.max(post); best_prob <- post[best_idx_local]
    
    if(k==1) ss_res_best <- colSums(X_cand^2)[best_idx_local]
    else ss_res_best <- colSums(X_cand_res^2)[best_idx_local]
    RSS_new <- RSS_current - (beta_cand[best_idx_local]^2 * ss_res_best)
    if(RSS_new < 0) RSS_new <- 1e-9
    df <- k + 1
    r2_adj <- 1 - (RSS_new/(n-df)) / (TSS/(n-1))
    aic_val <- n * log(RSS_new/n) + 2 * df
    bic_val <- n * log(RSS_new/n) + log(n) * df
    
    stop_flag <- FALSE
    if (k > 1) {
      if (criteria == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (criteria == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (criteria == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (criteria == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    best_var <- indices_remain_to_test[best_idx_local]
    indices.inclus <- c(indices.inclus, best_var)
    if(is.null(X_inclus)) X_inclus <- X[, best_var, drop=FALSE] else X_inclus <- cbind(X_inclus, X[, best_var, drop=FALSE])
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj); aic_hist <- c(aic_hist, aic_val); bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }
  indices_under_16 <- sum(indices.inclus < threshold)
  return(list(indices = indices.inclus, post.probs = post.probs, r_squared_adj = r_squared_adj_hist, aic_values = aic_hist, bic_values = bic_hist, indices_under_16 = indices_under_16, K = length(indices.inclus)))
}

#' Naive Multiple Test (Chi-Sq) with Stopping Criteria
Naive_Multiple_Criteria <- function(X, Y, a.priori, tau2, S, threshold = 16, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) stop("Criteria invalid.")
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  indices.inclus <- integer(0); indices_remain_to_test <- 1:m
  r_squared_adj_hist <- c(); aic_hist <- c(); bic_hist <- c(); post.probs <- c()
  TSS <- sum(Y^2); RSS_current <- TSS
  X_inclus <- NULL
  
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    rA <- k; df1 <- rA; df2 <- n - rA - 1
    if (df2 <= 0) break
    
    if (k == 1) {
      r_sq <- (as.vector(t(Y) %*% X_cand) / sqrt(TSS * colSums(X_cand^2)))^2
      r_sq[is.na(r_sq)] <- 0
      F_stat <- r_sq * (n - 2) / (1 - r_sq)
      nums <- as.vector(t(Y) %*% X_cand)
      beta_cand <- nums / colSums(X_cand^2)
    } else {
      if(is.null(X_inclus)) X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Y_res <- qr.resid(qr_inclus, Y)
      X_cand_res <- qr.resid(qr_inclus, X_cand)
      ss_X_cand_res <- colSums(X_cand_res^2); ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector(t(Y_res) %*% X_cand_res) / ss_X_cand_res
      RSS_prov <- sum(Y_res^2)
      RSS_new_vec <- RSS_prov - beta_cand^2 * ss_X_cand_res
      SSR_new_vec <- TSS - RSS_new_vec
      F_stat <- (SSR_new_vec / df1) / (RSS_new_vec / df2)
    }
    
    Q_chi <- rA * F_stat
    Q_chi[!is.finite(Q_chi) | Q_chi < 0] <- -1
    p_vec <- a.priori[indices_remain_to_test]
    post <- mapply(function(q, j, p) {
      if (q < 0) return(NA)
      X_mod <- if(is.null(X_inclus)) X[, j, drop=FALSE] else cbind(X_inclus, X[, j, drop=FALSE])
      NaM_posterior(Q_chi=q, X_tilde_jk=X_mod, X_k=X_inclus, n=n, p=p, tau2=tau2, S=S)
    }, Q_chi, indices_remain_to_test, p_vec)
    
    if (all(!is.finite(post))) break
    best_idx_local <- which.max(post); best_prob <- post[best_idx_local]
    
    if(k==1) ss_res_best <- colSums(X_cand^2)[best_idx_local]
    else ss_res_best <- colSums(X_cand_res^2)[best_idx_local]
    RSS_new <- RSS_current - (beta_cand[best_idx_local]^2 * ss_res_best)
    if(RSS_new < 0) RSS_new <- 1e-9
    
    df <- k + 1
    r2_adj <- 1 - (RSS_new/(n-df)) / (TSS/(n-1))
    aic_val <- n * log(RSS_new/n) + 2 * df
    bic_val <- n * log(RSS_new/n) + log(n) * df
    
    stop_flag <- FALSE
    if (k > 1) {
      if (criteria == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (criteria == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (criteria == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (criteria == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    best_var <- indices_remain_to_test[best_idx_local]
    indices.inclus <- c(indices.inclus, best_var)
    if(is.null(X_inclus)) X_inclus <- X[, best_var, drop=FALSE] else X_inclus <- cbind(X_inclus, X[, best_var, drop=FALSE])
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj); aic_hist <- c(aic_hist, aic_val); bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }
  indices_under_16 <- sum(indices.inclus < threshold)
  return(list(indices = indices.inclus, post.probs = post.probs, r_squared_adj = r_squared_adj_hist, aic_values = aic_hist, bic_values = bic_hist, indices_under_16 = indices_under_16, K = length(indices.inclus)))
}


# ==============================================================================
# Simulation Wrapper (Table 4 Generator)
# ==============================================================================

#' Simulation Wrapper for Stopping Criteria
#' @description Runs parallel simulations using criteria-based selection.
#' @return A summary table with Average K, CV, and Performance metrics.
Criteria_Simulation_Table <- function(B, X, sigma, a.priori, p_values, tau_values, S = 1000, target_prob = 0.95, threshold = 16, criteria = "BIC") {
  
  m <- ncol(X)
  n <- nrow(X)
  
  # --- 1. Parallel Setup ---
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  clusterEvalQ(cl, {
    library(MASS)
    library(mvtnorm)
    library(stats)
  })
  
  clusterExport(cl, c(
    "Generate.Y", "Naive_Wald_Criteria", "Asymptotic_Wald_Criteria", "Exact_Wald_Criteria", "Asymptotic_Score_Criteria", 
    "Exact_Multiple_Criteria", "Naive_Multiple_Criteria", "Power_Fun", "Naive_Wald_Integrand", 
    "Asymptotic_Wald_Integrand", "Exact_Wald_Integrand", "Asymptotic_Score_Integrand", "Fisher_MC_approximation",
    "Chisq_MC_approximation", "ExM_posterior", "NaM_posterior"
  ))
  
  final_summary_table <- data.frame()
  
  # --- 2. Main Loops ---
  for (p in p_values) {
    for (tau in tau_values) {
      
      # Calculate effect size
      b_root <- try(uniroot(Power_Fun, lower = 0, upper = 5, n = n, alpha = 0.05, p = p), silent = TRUE)
      if (inherits(b_root, "try-error")) b.val <- 0.5 else b.val <- b_root$root
      beta <- c(1, rep(b.val, 15), rep(0, m - 15))
      
      cat(sprintf("Running: p=%.2f, tau=%.2f\n", p, tau))
      
      # --- 3. Parallel Execution ---
      raw_results <- foreach(b = 1:B, .combine = rbind) %dopar% {
        
        safe_get <- function(vec) { if (length(vec) > 0) return(tail(vec, 1)) else return(NA) }
        
        Y <- Generate.Y(n, X, beta, sigma)
        
        iter_mat <- matrix(NA, nrow = 6, ncol = 7)
        colnames(iter_mat) <- c("Method", "K", "Posterior", "R2_adj", "AIC", "BIC", "indices_under_16")
        
        # Define the 6 runs (Passing 'criteria' correctly to all functions)
        runs <- list(
          list("NaW",  function() Naive_Wald_Criteria(X, Y, a.priori, tau, threshold, target_prob, criteria)),
          list("AsW",  function() Asymptotic_Wald_Criteria(X, Y, a.priori, tau, threshold, target_prob, criteria)),
          list("ExW",  function() Exact_Wald_Criteria(X, Y, a.priori, tau, threshold, target_prob, criteria)),
          list("AsS",  function() Asymptotic_Score_Criteria(X, Y, a.priori, tau, threshold, target_prob, criteria)),
          list("NaM",  function() Naive_Multiple_Criteria(X, Y, a.priori, tau^2, S, threshold, target_prob, criteria)),
          list("ExM",  function() Exact_Multiple_Criteria(X, Y, a.priori, tau^2, S, threshold, target_prob, criteria))
        )
        
        for(i in 1:6) {
          Method <- runs[[i]][[1]]
          func <- runs[[i]][[2]]
          
          # Try/Catch block prevents simulation crash
          res <- tryCatch(func(), error = function(e) return(NULL))
          
          iter_mat[i, 1] <- Method
          if (!is.null(res)) {
            iter_mat[i, 2] <- length(res$indices)
            iter_mat[i, 3] <- safe_get(res$post.probs)
            iter_mat[i, 4] <- safe_get(res$r_squared_adj)
            iter_mat[i, 5] <- safe_get(res$aic_values)
            iter_mat[i, 6] <- safe_get(res$bic_values)
            iter_mat[i, 7] <- res$indices_under_16
          }
        }
        iter_mat
      }
      
      # --- 4. Summarization ---
      raw_df <- as.data.frame(raw_results, stringsAsFactors = FALSE)
      # Ensure numeric columns are actually numeric
      for(j in 2:7) raw_df[, j] <- suppressWarnings(as.numeric(raw_df[, j]))
      method_Methods <- c("NaW", "AsW", "ExW", "AsS", "NaM", "ExM")
      
      for(m_Method in method_Methods){
        sub_dat <- raw_df[raw_df[, 1] == m_Method, ]
        
        avg_K <- mean(sub_dat[, 2], na.rm = TRUE)
        median_K <- median(sub_dat[, 2], na.rm = TRUE)
        
        if (is.na(avg_K) || is.nan(avg_K) || avg_K == 0) {
          cv_K <- 0
        } else {
          cv_K <- sd(sub_dat[, 2], na.rm = TRUE) / avg_K
        }
        
        final_summary_table <- rbind(final_summary_table, data.frame(
          Method = m_Method,
          p = p,
          tau = tau,
          Average_K = avg_K,
          Median_K = median_K,
          CV_K = cv_K,
          indices_under_16 = mean(sub_dat[, 7], na.rm = TRUE), # Precision Metric
          Average_prob = mean(sub_dat[, 3], na.rm = TRUE),
          Average_R2.adj = mean(sub_dat[, 4], na.rm = TRUE),
          Average_Aic = mean(sub_dat[, 5], na.rm = TRUE),
          Average_Bic = mean(sub_dat[, 6], na.rm = TRUE)
        ))
      }
    }
  }
  stopCluster(cl)
  return(final_summary_table)
}


# ==============================================================================
# SECTION B: Criteria functions for Real Data Analysis (Table 5/6)
# ==============================================================================
# Note: These are identical to Section A but allow a hard stop 'threshold_K'
#       to prevent infinite loops in real data if criteria fail to converge.

Naive_Wald_Criteria_For_Data <- function(X, Y, a.priori, tau, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) stop("Criteria invalid.")
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  indices.inclus <- integer(0); indices_remain_to_test <- 1:m
  r_squared_adj_hist <- numeric(0); aic_hist <- numeric(0); bic_hist <- numeric(0); post.probs <- numeric(0)
  TSS <- sum(Y^2); RSS_current <- TSS 
  
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      nums <- as.vector(t(Y) %*% X_cand)
      dens <- sqrt(TSS * colSums(X_cand^2))
      r_vec <- nums / (dens + 1e-12)
      Z_vector <- r_vec * sqrt(n - 2) / sqrt(1 - r_vec^2)
      rho_vector <- numeric(length(indices_remain_to_test))
      ss_X_cand <- colSums(X_cand^2)
      beta_cand <- nums / ss_X_cand
    } else {
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus); Q <- qr.Q(qr_inclus)
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      XtX_cand <- t(Q) %*% X_cand
      X_cand_res <- X_cand - Q %*% XtX_cand
      ss_X_cand_res <- colSums(X_cand_res^2); ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector((t(Y_res) %*% X_cand_res) / ss_X_cand_res)
      RSS_provisional <- sum(Y_res^2) - beta_cand^2 * ss_X_cand_res
      sigma2_new <- RSS_provisional / (n - k - 1)
      se_beta <- sqrt(sigma2_new / ss_X_cand_res)
      Z_vector <- beta_cand / se_beta
      
      V_inv <- n * chol2inv(qr.R(qr_inclus))
      R_mat <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_mat
      rho_vector <- rowSums(t(R_matrix) * t(M))
      rho_vector <- pmax(0, pmin(rho_vector, 1 - 1e-9))
    }
    
    Z_vector[!is.finite(Z_vector)] <- 0
    sd_eff <- sqrt(1 + n * (1 - rho_vector) * tau^2)
    integrals <- dnorm(Z_vector, mean = 0, sd = sd_eff)
    densities <- dnorm(Z_vector, mean = 0, sd = 1)
    p_vec <- a.priori[indices_remain_to_test]
    numerator <- p_vec * integrals
    denominator <- numerator + (1 - p_vec) * densities
    a.posteriori <- numerator / (denominator + 1e-20)
    
    max_idx_local <- which.max(a.posteriori)
    current_max_prob <- a.posteriori[max_idx_local]
    best_var_global <- indices_remain_to_test[max_idx_local]
    
    best_beta <- beta_cand[max_idx_local]
    if (k == 1) best_ss_res <- colSums(X_cand^2)[max_idx_local]
    else best_ss_res <- colSums(X_cand_res^2)[max_idx_local]
    
    RSS_new <- RSS_current - (best_beta^2 * best_ss_res)
    if(RSS_new < 0) RSS_new <- 1e-9 
    df_model <- k + 1 
    
    r2_adj_val <- 1 - (RSS_new / (n - df_model)) / (TSS / (n - 1))
    aic_val <- n * log(RSS_new / n) + 2 * df_model
    bic_val <- n * log(RSS_new / n) + log(n) * df_model
    
    stop_flag <- FALSE
    if (k > 1) {
      if (criteria == "R2" && r2_adj_val < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (criteria == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (criteria == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (criteria == "Posterior" && current_max_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    indices.inclus <- c(indices.inclus, best_var_global)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var_global)
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj_val)
    aic_hist <- c(aic_hist, aic_val)
    bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, current_max_prob)
    RSS_current <- RSS_new 
    
    if (!is.null(threshold_K) && (length(indices.inclus) >= threshold_K)) break 
  }
  
  return(list(indices = indices.inclus, post.probs = post.probs, r_squared_adj = r_squared_adj_hist, aic_values = aic_hist, bic_values = bic_hist, K = length(indices.inclus)))
}

Asymptotic_Wald_Criteria_For_Data <- function(X, Y, a.priori, tau, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) stop("Criteria invalid.")
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  indices.inclus <- integer(0); indices_remain_to_test <- 1:m
  r_squared_adj_hist <- numeric(0); aic_hist <- numeric(0); bic_hist <- numeric(0); post.probs <- numeric(0)
  TSS <- sum(Y^2); RSS_current <- TSS
  
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      nums <- as.vector(t(Y) %*% X_cand)
      ss_X_cand <- colSums(X_cand^2)
      beta_cand <- nums / ss_X_cand
      RSS_prov <- TSS - beta_cand^2 * ss_X_cand
      sigma2 <- RSS_prov / (n - 2)
      Z_vec <- beta_cand / sqrt(sigma2 / ss_X_cand)
      rho_vec <- numeric(length(indices_remain_to_test))
    } else {
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus); Q <- qr.Q(qr_inclus)
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      XtX_cand <- t(Q) %*% X_cand
      X_cand_res <- X_cand - Q %*% XtX_cand
      ss_X_cand_res <- colSums(X_cand_res^2); ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector((t(Y_res) %*% X_cand_res) / ss_X_cand_res)
      RSS_prov <- sum(Y_res^2) - beta_cand^2 * ss_X_cand_res
      sigma2 <- RSS_prov / (n - k - 1)
      Z_vec <- beta_cand / sqrt(sigma2 / ss_X_cand_res)
      V_inv <- n * chol2inv(qr.R(qr_inclus))
      R_mat <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_mat
      rho_vec <- colSums(R_mat * M) 
      rho_vec <- pmax(0, pmin(rho_vec, 1 - 1e-9))
    }
    
    Z_vec[!is.finite(Z_vec)] <- 0
    integrals <- mapply(function(z, r) {
      tryCatch(integrate(Asymptotic_Wald_Integrand, lower=-Inf, upper=Inf, Z=z, n=n, rho=r, tau=tau)$value, error=function(e) 0)
    }, Z_vec, rho_vec)
    
    densities <- dnorm(Z_vec)
    p_vec <- a.priori[indices_remain_to_test]
    numerator <- p_vec * integrals
    denominator <- numerator + (1 - p_vec) * densities
    post <- numerator / (denominator + 1e-20)
    best_idx_local <- which.max(post); best_prob <- post[best_idx_local]
    best_var <- indices_remain_to_test[best_idx_local]
    
    if(k==1) ss_res_best <- colSums(X_cand^2)[best_idx_local]
    else ss_res_best <- colSums(X_cand_res^2)[best_idx_local]
    RSS_new <- RSS_current - (beta_cand[best_idx_local]^2 * ss_res_best)
    if(RSS_new < 0) RSS_new <- 1e-9
    df <- k + 1
    r2_adj <- 1 - (RSS_new/(n-df)) / (TSS/(n-1))
    aic_val <- n * log(RSS_new/n) + 2 * df
    bic_val <- n * log(RSS_new/n) + log(n) * df
    
    stop_flag <- FALSE
    if (k > 1) {
      if (criteria == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (criteria == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (criteria == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (criteria == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    indices.inclus <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj)
    aic_hist <- c(aic_hist, aic_val)
    bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
    
    if (!is.null(threshold_K) && (length(indices.inclus) >= threshold_K)) break 
  }
  return(list(indices = indices.inclus, post.probs = post.probs, r_squared_adj = r_squared_adj_hist, aic_values = aic_hist, bic_values = bic_hist, K = length(indices.inclus)))
}

Exact_Wald_Criteria_For_Data <- function(X, Y, a.priori, tau, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) stop("Criteria invalid.")
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  indices.inclus <- integer(0); indices_remain_to_test <- 1:m
  r_squared_adj_hist <- c(); aic_hist <- c(); bic_hist <- c(); post.probs <- c()
  TSS <- sum(Y^2); RSS_current <- TSS
  
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      nums <- as.vector(t(Y) %*% X_cand)
      ss_X_cand <- colSums(X_cand^2)
      beta_cand <- nums / ss_X_cand
      RSS_prov <- TSS - beta_cand^2 * ss_X_cand
      sigma2 <- RSS_prov / (n - 2)
      Z_vec <- beta_cand / sqrt(sigma2 / ss_X_cand)
      rho_vec <- numeric(length(indices_remain_to_test))
      Z_hat_vec <- Z_vec
    } else {
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus); Q <- qr.Q(qr_inclus)
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      X_cand_res <- X_cand - Q %*% (t(Q) %*% X_cand)
      ss_X_cand_res <- colSums(X_cand_res^2); ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector((t(Y_res) %*% X_cand_res) / ss_X_cand_res)
      RSS_prov <- sum(Y_res^2) - beta_cand^2 * ss_X_cand_res
      sigma2 <- RSS_prov / (n - k - 1)
      Z_vec <- beta_cand / sqrt(sigma2 / ss_X_cand_res)
      V_inv <- n * chol2inv(qr.R(qr_inclus))
      R_mat <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_mat
      rho_vec <- colSums(R_mat * M)
      rho_vec <- pmax(0, pmin(rho_vec, 1 - 1e-9))
      Z_hat_vec <- Z_vec * sqrt(1 - rho_vec)
    }
    
    Z_hat_vec[!is.finite(Z_hat_vec)] <- 0
    integrals <- mapply(function(z, r) {
      tryCatch(integrate(Exact_Wald_Integrand, lower=-Inf, upper=Inf, Z=z, n=n, rho=r, tau=tau, k=k)$value, error=function(e) 0)
    }, Z_hat_vec, rho_vec)
    
    densities <- dt(Z_hat_vec, df = n - (k + 1))
    p_vec <- a.priori[indices_remain_to_test]
    numerator <- p_vec * integrals
    denominator <- numerator + (1 - p_vec) * densities
    post <- numerator / (denominator + 1e-20)
    best_idx_local <- which.max(post); best_prob <- post[best_idx_local]
    
    if(k==1) ss_res_best <- colSums(X_cand^2)[best_idx_local]
    else ss_res_best <- colSums(X_cand_res^2)[best_idx_local]
    RSS_new <- RSS_current - (beta_cand[best_idx_local]^2 * ss_res_best)
    if(RSS_new < 0) RSS_new <- 1e-9
    df <- k + 1
    r2_adj <- 1 - (RSS_new/(n-df)) / (TSS/(n-1))
    aic_val <- n * log(RSS_new/n) + 2 * df
    bic_val <- n * log(RSS_new/n) + log(n) * df
    
    stop_flag <- FALSE
    if (k > 1) {
      if (criteria == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (criteria == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (criteria == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (criteria == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    indices.inclus <- c(indices.inclus, indices_remain_to_test[best_idx_local])
    indices_remain_to_test <- setdiff(indices_remain_to_test, indices_remain_to_test[best_idx_local])
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj); aic_hist <- c(aic_hist, aic_val); bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
    if (!is.null(threshold_K) && (length(indices.inclus) >= threshold_K)) break 
  }
  return(list(indices = indices.inclus, post.probs = post.probs, r_squared_adj = r_squared_adj_hist, aic_values = aic_hist, bic_values = bic_hist, K = length(indices.inclus)))
}

Asymptotic_Score_Criteria_For_Data <- function(X, Y, a.priori, tau, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) stop("Criteria invalid.")
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  indices.inclus <- integer(0); indices_remain_to_test <- 1:m
  r_squared_adj_hist <- c(); aic_hist <- c(); bic_hist <- c(); post.probs <- c()
  TSS <- sum(Y^2); RSS_current <- TSS
  
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      numerator <- as.vector(t(Y) %*% X_cand)
      ss_X <- colSums(X_cand^2)
      beta <- numerator / ss_X
      RSS <- TSS - beta * numerator
      sigma_vec <- sqrt(RSS / (n - 2))
      Q_hat <- numerator / (sqrt(n) * sigma_vec)
      rho_vec <- numeric(length(indices_remain_to_test))
      beta_cand <- beta
    } else {
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus); Q <- qr.Q(qr_inclus)
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      X_cand_res <- X_cand - Q %*% (t(Q) %*% X_cand)
      numerator <- as.vector(t(Y) %*% X_cand_res)
      ss_X_cand_res <- colSums(X_cand_res^2); ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector((t(Y) %*% X_cand_res) / ss_X_cand_res)
      RSS_prov <- sum(Y_res^2) - beta_cand^2 * ss_X_cand_res
      sigma_new_vec <- sqrt(RSS_prov / (n - k - 1))
      Q_hat <- numerator / (sqrt(n) * sigma_new_vec)
      V_inv <- n * chol2inv(qr.R(qr_inclus))
      R_mat <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_mat
      rho_vec <- colSums(R_mat * M)
      rho_vec <- pmax(0, pmin(rho_vec, 1 - 1e-9))
    }
    
    Q_hat[!is.finite(Q_hat)] <- 0
    integrals <- mapply(function(q, r) {
      tryCatch(integrate(Asymptotic_Score_Integrand, lower=-Inf, upper=Inf, Q=q, n=n, rho=r, tau=tau)$value, error=function(e) 0)
    }, Q_hat, rho_vec)
    
    densities <- dnorm(Q_hat)
    p_vec <- a.priori[indices_remain_to_test]
    numerator <- p_vec * integrals
    denominator <- numerator + (1 - p_vec) * densities
    post <- numerator / (denominator + 1e-20)
    best_idx_local <- which.max(post); best_prob <- post[best_idx_local]
    
    if(k==1) ss_res_best <- colSums(X_cand^2)[best_idx_local]
    else ss_res_best <- colSums(X_cand_res^2)[best_idx_local]
    RSS_new <- RSS_current - (beta_cand[best_idx_local]^2 * ss_res_best)
    if(RSS_new < 0) RSS_new <- 1e-9
    df <- k + 1
    r2_adj <- 1 - (RSS_new/(n-df)) / (TSS/(n-1))
    aic_val <- n * log(RSS_new/n) + 2 * df
    bic_val <- n * log(RSS_new/n) + log(n) * df
    
    stop_flag <- FALSE
    if (k > 1) {
      if (criteria == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (criteria == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (criteria == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (criteria == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    indices.inclus <- c(indices.inclus, indices_remain_to_test[best_idx_local])
    indices_remain_to_test <- setdiff(indices_remain_to_test, indices_remain_to_test[best_idx_local])
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj); aic_hist <- c(aic_hist, aic_val); bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
    if (!is.null(threshold_K) && (length(indices.inclus) >= threshold_K)) break 
  }
  return(list(indices = indices.inclus, post.probs = post.probs, r_squared_adj = r_squared_adj_hist, aic_values = aic_hist, bic_values = bic_hist, K = length(indices.inclus)))
}

Exact_Multiple_Criteria_For_Data <- function(X, Y, a.priori, tau2, S, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) stop("Criteria invalid.")
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  indices.inclus <- integer(0); indices_remain_to_test <- 1:m
  r_squared_adj_hist <- c(); aic_hist <- c(); bic_hist <- c(); post.probs <- c()
  TSS <- sum(Y^2); RSS_current <- TSS
  X_inclus <- NULL
  
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    rA <- k; df1 <- rA; df2 <- n - rA - 1
    if (df2 <= 0) break
    
    if (k == 1) {
      r_sq <- (as.vector(t(Y) %*% X_cand) / sqrt(TSS * colSums(X_cand^2)))^2
      r_sq[is.na(r_sq)] <- 0
      Z_vector <- r_sq * (n - 2) / (1 - r_sq)
      nums <- as.vector(t(Y) %*% X_cand)
      beta_cand <- nums / colSums(X_cand^2)
    } else {
      if(is.null(X_inclus)) X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Y_res <- qr.resid(qr_inclus, Y)
      X_cand_res <- qr.resid(qr_inclus, X_cand)
      ss_X_cand_res <- colSums(X_cand_res^2); ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector(t(Y_res) %*% X_cand_res) / ss_X_cand_res
      RSS_prov <- sum(Y_res^2)
      RSS_new_vec <- RSS_prov - beta_cand^2 * ss_X_cand_res
      SSR_new_vec <- TSS - RSS_new_vec
      Z_vector <- (SSR_new_vec / df1) / (RSS_new_vec / df2)
    }
    
    Z_vector[!is.finite(Z_vector) | Z_vector < 0] <- -1
    p_vec <- a.priori[indices_remain_to_test]
    
    post <- mapply(function(z, j, p) {
      if (z < 0) return(NA)
      X_mod <- if(is.null(X_inclus)) X[, j, drop=FALSE] else cbind(X_inclus, X[, j, drop=FALSE])
      ExM_posterior(Q_W=z, X_tilde_jk=X_mod, X_k=X_inclus, n=n, p=p, tau2=tau2, S=S)
    }, Z_vector, indices_remain_to_test, p_vec)
    
    if (all(!is.finite(post))) break
    best_idx_local <- which.max(post); best_prob <- post[best_idx_local]
    
    if(k==1) ss_res_best <- colSums(X_cand^2)[best_idx_local]
    else ss_res_best <- colSums(X_cand_res^2)[best_idx_local]
    RSS_new <- RSS_current - (beta_cand[best_idx_local]^2 * ss_res_best)
    if(RSS_new < 0) RSS_new <- 1e-9
    df <- k + 1
    r2_adj <- 1 - (RSS_new/(n-df)) / (TSS/(n-1))
    aic_val <- n * log(RSS_new/n) + 2 * df
    bic_val <- n * log(RSS_new/n) + log(n) * df
    
    stop_flag <- FALSE
    if (k > 1) {
      if (criteria == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (criteria == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (criteria == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (criteria == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    best_var <- indices_remain_to_test[best_idx_local]
    indices.inclus <- c(indices.inclus, best_var)
    if(is.null(X_inclus)) X_inclus <- X[, best_var, drop=FALSE] else X_inclus <- cbind(X_inclus, X[, best_var, drop=FALSE])
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj); aic_hist <- c(aic_hist, aic_val); bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
    if (!is.null(threshold_K) && (length(indices.inclus) >= threshold_K)) break
  }
  return(list(indices = indices.inclus, post.probs = post.probs, r_squared_adj = r_squared_adj_hist, aic_values = aic_hist, bic_values = bic_hist, K = length(indices.inclus)))
}

Naive_Multiple_Criteria_For_Data <- function(X, Y, a.priori, tau2, S, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) stop("Criteria invalid.")
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  indices.inclus <- integer(0); indices_remain_to_test <- 1:m
  r_squared_adj_hist <- c(); aic_hist <- c(); bic_hist <- c(); post.probs <- c()
  TSS <- sum(Y^2); RSS_current <- TSS
  X_inclus <- NULL
  
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    rA <- k; df1 <- rA; df2 <- n - rA - 1
    if (df2 <= 0) break
    
    if (k == 1) {
      r_sq <- (as.vector(t(Y) %*% X_cand) / sqrt(TSS * colSums(X_cand^2)))^2
      r_sq[is.na(r_sq)] <- 0
      F_stat <- r_sq * (n - 2) / (1 - r_sq)
      nums <- as.vector(t(Y) %*% X_cand)
      beta_cand <- nums / colSums(X_cand^2)
    } else {
      if(is.null(X_inclus)) X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Y_res <- qr.resid(qr_inclus, Y)
      X_cand_res <- qr.resid(qr_inclus, X_cand)
      ss_X_cand_res <- colSums(X_cand_res^2); ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector(t(Y_res) %*% X_cand_res) / ss_X_cand_res
      RSS_prov <- sum(Y_res^2)
      RSS_new_vec <- RSS_prov - beta_cand^2 * ss_X_cand_res
      SSR_new_vec <- TSS - RSS_new_vec
      F_stat <- (SSR_new_vec / df1) / (RSS_new_vec / df2)
    }
    
    Q_chi <- rA * F_stat
    Q_chi[!is.finite(Q_chi) | Q_chi < 0] <- -1
    p_vec <- a.priori[indices_remain_to_test]
    post <- mapply(function(q, j, p) {
      if (q < 0) return(NA)
      X_mod <- if(is.null(X_inclus)) X[, j, drop=FALSE] else cbind(X_inclus, X[, j, drop=FALSE])
      NaM_posterior(Q_chi=q, X_tilde_jk=X_mod, X_k=X_inclus, n=n, p=p, tau2=tau2, S=S)
    }, Q_chi, indices_remain_to_test, p_vec)
    
    if (all(!is.finite(post))) break
    best_idx_local <- which.max(post); best_prob <- post[best_idx_local]
    
    if(k==1) ss_res_best <- colSums(X_cand^2)[best_idx_local]
    else ss_res_best <- colSums(X_cand_res^2)[best_idx_local]
    RSS_new <- RSS_current - (beta_cand[best_idx_local]^2 * ss_res_best)
    if(RSS_new < 0) RSS_new <- 1e-9
    
    df <- k + 1
    r2_adj <- 1 - (RSS_new/(n-df)) / (TSS/(n-1))
    aic_val <- n * log(RSS_new/n) + 2 * df
    bic_val <- n * log(RSS_new/n) + log(n) * df
    
    stop_flag <- FALSE
    if (k > 1) {
      if (criteria == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (criteria == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (criteria == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (criteria == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    best_var <- indices_remain_to_test[best_idx_local]
    indices.inclus <- c(indices.inclus, best_var)
    if(is.null(X_inclus)) X_inclus <- X[, best_var, drop=FALSE] else X_inclus <- cbind(X_inclus, X[, best_var, drop=FALSE])
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj); aic_hist <- c(aic_hist, aic_val); bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
    if (!is.null(threshold_K) && (length(indices.inclus) >= threshold_K)) break
  }
  return(list(indices = indices.inclus, post.probs = post.probs, r_squared_adj = r_squared_adj_hist, aic_values = aic_hist, bic_values = bic_hist, K = length(indices.inclus)))
}