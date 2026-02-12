source("Test_statistics_functions.R")

#---- Functions with stopping criteria 
###NOTE: Keep threshold ONLY for simulation 


###-------------  A. Criteria functions for simulation study ---------------- 
### creating Table 4
## Naive Wald
Naive_Wald_Criteria <- function(X, Y, a.priori, tau, threshold = 16, target_prob = 0.95, criteria = "BIC") {
  
  # --- 1. Validation & Setup ---
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria ivalid, shoud be one of 'R2', 'AIC', 'BIC' ou 'Posterior'.")
  }
  
  m <- ncol(X)
  n <- nrow(X)
  
  # Center Data 
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X) # to be removed if Simulated X is scaled 
  
  # Initialize containers
  indices.inclus <- integer(0)
  indices_remain_to_test <- 1:m
  
  r_squared_adj_hist <- numeric(0)
  aic_hist <- numeric(0)
  bic_hist <- numeric(0)
  post.probs <- numeric(0)
  
  # Calculate Total Sum of Squares (TSS) for R2 calculation
  TSS <- sum(Y^2)
  RSS_current <- TSS # Start with Null model RSS
  
  # --- 2. Main Selection Loop ---
  # We loop up to 'm' (or stop early), but usually bounded by threshold in practice
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    # --- A. Compute Statistics (Vectorized) ---
    if (k == 1) {
      # Base case: Correlations
      # Vectorized calc: (y'x) / sqrt(sum(x^2))
      nums <- as.vector(t(Y) %*% X_cand)
      dens <- sqrt(TSS * colSums(X_cand^2))
      r_vec <- nums / (dens + 1e-12)
      
      # Z-scores
      Z_vector <- r_vec * sqrt(n - 2) / sqrt(1 - r_vec^2)
      rho_vector <- numeric(length(indices_remain_to_test)) # rho=0 at start
      # Store sum of squares of candidates for later RSS update
      ss_X_cand <- colSums(X_cand^2)
      beta_cand <- nums / ss_X_cand
    } else {
      # General case: QR Update
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Q <- qr.Q(qr_inclus)
      # Residuals of Y
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      # Orthogonalize candidates (Gram-Schmidt step)
      XtX_cand <- t(Q) %*% X_cand
      X_cand_res <- X_cand - Q %*% XtX_cand
      ss_X_cand_res <- colSums(X_cand_res^2)
      ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      # Beta candidates on the residual space
      beta_cand <- as.vector((t(Y_res) %*% X_cand_res) / ss_X_cand_res)
      # Variance estimation for Z-score
      # Note: We use provisional RSS for each candidate to get t-stat
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
    
    # Handle Infs
    Z_vector[!is.finite(Z_vector)] <- 0
    
    # --- B. Posterior Probability (Analytic Integration) ---
    # Analytic solution for Naive_Wald_Integrand:
    # Marginal distribution is N(0, 1 + n*(1-rho)*tau^2)
    sd_eff <- sqrt(1 + n * (1 - rho_vector) * tau^2)
    integrals <- dnorm(Z_vector, mean = 0, sd = sd_eff)
    densities <- dnorm(Z_vector, mean = 0, sd = 1)
    p_vec <- a.priori[indices_remain_to_test]
    numerator <- p_vec * integrals
    denominator <- numerator + (1 - p_vec) * densities
    a.posteriori <- numerator / (denominator + 1e-20)
    # Pick Best
    max_idx_local <- which.max(a.posteriori)
    current_max_prob <- a.posteriori[max_idx_local]
    best_var_global <- indices_remain_to_test[max_idx_local]
    
    # --- C. Calculate Metrics for the potential New Model ---
    # To get exact metrics, we update RSS based on the chosen variable
    # RSS_new = RSS_old - (beta_cand^2 * ss_X_cand_res) [of the best var]
    
    # We need the specific stats of the best variable from the vector calculation above
    best_beta <- beta_cand[max_idx_local]
    # Get Sum of Squares of the residual of the candidate
    if (k == 1) {
      best_ss_res <- colSums(X_cand^2)[max_idx_local]
    } else {
      best_ss_res <- colSums(X_cand_res^2)[max_idx_local]
    }
    
    RSS_new <- RSS_current - (best_beta^2 * best_ss_res)
    if(RSS_new < 0) RSS_new <- 1e-9 # Safety
    # Degrees of freedom: k variables + 1 intercept
    df_model <- k + 1 
    
    # Calculate Metrics
    # R2 Adj: 1 - (RSS/(n-p)) / (TSS/(n-1))
    r2_adj_val <- 1 - (RSS_new / (n - df_model)) / (TSS / (n - 1))
    # AIC: n * log(RSS/n) + 2*p
    aic_val <- n * log(RSS_new / n) + 2 * df_model
    # BIC: n * log(RSS/n) + log(n)*p
    bic_val <- n * log(RSS_new / n) + log(n) * df_model
    
    # --- D. Stopping Criteria Check ---
    stop_flag <- FALSE
    
    if (k > 1) {
      if (critere == "R2") {
        # Stop if Adjusted R2 decreases
        if (r2_adj_val < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      } else if (critere == "AIC") {
        # Stop if AIC increases
        if (aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      } else if (critere == "BIC") {
        # Stop if BIC increases
        if (bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      } else if (critere == "Posterior") {
        # Stop if Probability drops below target
        if (current_max_prob < target_prob) stop_flag <- TRUE
      }
    }
    
    if (stop_flag) {
      break # Stop BEFORE adding the variable (or after, depending on preference)
      # Usually stepwise stops avoiding the 'bad' step. 
      # We break here, so the variable is NOT added to indices.inclus
    }
    
    # --- E. Update Lists ---
    indices.inclus <- c(indices.inclus, best_var_global)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var_global)
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj_val)
    aic_hist <- c(aic_hist, aic_val)
    bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, current_max_prob)
    RSS_current <- RSS_new # Update for next iteration
  }
  # Count selected variables < threshold (e.g. 16)
  indices_under_16 <- sum(indices.inclus < threshold) # to be removed for data analysis
  return(list(
    indices = indices.inclus,
    post.probs = post.probs,
    r_squared_adj = r_squared_adj_hist,
    aic_values = aic_hist,
    bic_values = bic_hist,
    K = length(indices.inclus)
  ))
}

## Asymptotic Wald
Asymptotic_Wald_Criteria <- function(X, Y, a.priori, tau, threshold = 16, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria ivalid, shoud be one of 'R2', 'AIC', 'BIC' ou 'Posterior'.")
  }
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
      # M is (k x p). R_mat is (k x p). We multiply element-wise and sum columns.
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
      if (critere == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (critere == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (critere == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (critere == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
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
              K = length(indices.inclus)))
}
## Exact Wald
Exact_Wald_Criteria <- function(X, Y, a.priori, tau, threshold = 16, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria ivalid, shoud be one of 'R2', 'AIC', 'BIC' ou 'Posterior'.")
  }
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
      if (critere == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (critere == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (critere == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (critere == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    indices.inclus <- c(indices.inclus, indices_remain_to_test[best_idx_local])
    indices_remain_to_test <- setdiff(indices_remain_to_test, indices_remain_to_test[best_idx_local])
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj); aic_hist <- c(aic_hist, aic_val); bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }
  indices_under_16 <- sum(indices.inclus < threshold)
  return(list(indices = indices.inclus,
              post.probs = post.probs,
              r_squared_adj = r_squared_adj_hist,
              aic_values = aic_hist,
              bic_values = bic_hist,
              K = length(indices.inclus)))
}

## Asymptotic Score
Asymptotic_Score_Criteria <- function(X, Y, a.priori, tau, threshold = 16, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria ivalid, shoud be one of 'R2', 'AIC', 'BIC' ou 'Posterior'.")
  }
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
      if (critere == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (critere == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (critere == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (critere == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    
    indices.inclus <- c(indices.inclus, indices_remain_to_test[best_idx_local])
    indices_remain_to_test <- setdiff(indices_remain_to_test, indices_remain_to_test[best_idx_local])
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj); aic_hist <- c(aic_hist, aic_val); bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }
  indices_under_16 <- sum(indices.inclus < threshold)
  return(list(indices = indices.inclus,
              post.probs = post.probs,
              r_squared_adj = r_squared_adj_hist,
              aic_values = aic_hist,
              bic_values = bic_hist,
              K = length(indices.inclus)))
}
## Multi wald

Exact_Multiple_Criteria <- function(X, Y, a.priori, tau2, S, threshold = 16, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria ivalid, shoud be one of 'R2', 'AIC', 'BIC' ou 'Posterior'.")
  }
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
      # For criteria calc
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
    # Vectorized Posterior (using mapply wrapper for approximation function)
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
      if (critere == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (critere == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (critere == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (critere == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
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
  return(list(indices = indices.inclus,
              post.probs = post.probs,
              r_squared_adj = r_squared_adj_hist,
              aic_values = aic_hist,
              bic_values = bic_hist,
              K = length(indices.inclus)))
}

## Naive multi wald

Naive_Multiple_Criteria <- function(X, Y, a.priori, tau2, S, threshold = 16, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria ivalid, shoud be one of 'R2', 'AIC', 'BIC' ou 'Posterior'.")
  }
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
      if (critere == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (critere == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (critere == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (critere == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
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
  return(list(indices = indices.inclus,
              post.probs = post.probs,
              r_squared_adj = r_squared_adj_hist,
              aic_values = aic_hist,
              bic_values = bic_hist,
              K = length(indices.inclus)))
}

#---- Simulation function

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
  
  # Export functions
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
      # Calculate b.val 
      b_root <- try(uniroot(Power_Fun, lower = 0, upper = 5, n = n, alpha = 0.05, p = p), silent = TRUE)
      if (inherits(b_root, "try-error")) b.val <- 0.5 else b.val <- b_root$root
      beta <- c(1, rep(b.val, 15), rep(0, m - 15))
      cat(sprintf("Running: p=%.2f, tau=%.2f\n", p, tau))
      # --- 3. Parallel Execution ---
      raw_results <- foreach(b = 1:B, .combine = rbind) %dopar% {
        # fun to safely extract scalars
        safe_get <- function(vec) { if (length(vec) > 0) return(tail(vec, 1)) else return(NA) }
  
        Y <- Generate.Y(n, X, beta, sigma)
        
        # We define a matrix to hold this iteration's results (6 methods, 7 columns)
        # Cols: Method, K, Posterior, R2_adj, AIC, BIC, indices_under_16
        iter_mat <- matrix(NA, nrow = 6, ncol = 7)
        colMethods(iter_mat) <- c("Method", "K", "Posterior", "R2_adj", "AIC", "BIC", "indices_under_16")
        # Define the 6 runs
        runs <- list(
          list("NaW",  function() Naive_Wald_Criteria(X, Y, a.priori, tau, threshold, target_prob, critere)),
          list("AsW",  function() Asymptotic_Wald_Criteria(X, Y, a.priori, tau, threshold, target_prob, critere)),
          list("ExW",  function() Exact_Wald_Criteria(X, Y, a.priori, tau, threshold, target_prob, critere)),
          list("AsS",  function() Asymptotic_Score_Criteria(X, Y, a.priori, tau, threshold, target_prob, critere)),
          list("NaM",  function() Naive_Multiple_Criteria(X, Y, a.priori, tau^2, S, threshold, target_prob, critere)),
          list("ExM",  function() Exact_Multiple_Criteria(X, Y, a.priori, tau^2, S, threshold, target_prob, critere))
        )
        # Execute each method
        for(i in 1:6) {
          Method <- runs[[i]][[1]]
          func <- runs[[i]][[2]]
          # Try/Catch to prevent one method crashing the whole batch
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
      # Convert result columns to numeric (Columns 2 to 7)
      # Suppress warnings about NAs introduced by coercion
      for(j in 2:7) raw_df[, j] <- suppressWarnings(as.numeric(raw_df[, j]))
      method_Methods <- c("NaW", "AsW", "ExW", "AsS", "NaM", "ExM")
      
      for(m_Method in method_Methods){
        # Filter data for this method
        # Column 1 is "Method" (V1 in raw_df usually)
        # Note: foreach .combine=rbind on matrices creates a big matrix.
        # as.data.frame makes columns V1, V2, etc.
        sub_dat <- raw_df[raw_df[, 1] == m_Method, ]
        # Calculate Average K
        avg_K <- mean(sub_dat[, 2], na.rm = TRUE)
        median_K <- median(sub_dat[, 2], na.rm = TRUE)
        # If avg_K is NaN (e.g. all runs failed or resulted in NA), we set CV to 0
        if (is.na(avg_K) || is.nan(avg_K) || avg_K == 0) {
          cv_K <- 0
        } else {
          cv_K <- sd(sub_dat[, 2], na.rm = TRUE) / avg_K
        }
        # -----------------------------------------
        final_summary_table <- rbind(final_summary_table, data.frame(
          Method = m_Method,
          p = p,
          tau = tau,
          Average_K = avg_K,
          Median_K = median_K,
          CV_K = cv_K,
          indices_under_16 = mean(sub_dat[, 7], na.rm = TRUE),
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

###-------------  A. Criteria functions for Real data analysis ---------------- 
### creating Table 5

## Naive Wald
Naive_Wald_Criteria_For_Data <- function(X, Y, a.priori, tau, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  # --- 1. Validation & Setup ---
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria ivalid, shoud be one of 'R2', 'AIC', 'BIC' ou 'Posterior'.")
  }
  m <- ncol(X)
  n <- nrow(X)
  # Center Data
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  # Initialize containers
  indices.inclus <- integer(0)
  indices_remain_to_test <- 1:m
  r_squared_adj_hist <- numeric(0)
  aic_hist <- numeric(0)
  bic_hist <- numeric(0)
  post.probs <- numeric(0)
  TSS <- sum(Y^2)
  RSS_current <- TSS 
  # --- 2. Main Selection Loop ---
  for (k in 1:m) {
    if (length(indices_remain_to_test) == 0) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    if (k == 1) {
      nums <- as.vector(t(Y) %*% X_cand)
      dens <- sqrt(TSS * colSums(X_cand^2))
      r_vec <- nums / (dens + 1e-12)
      Z_vector <- r_vec * sqrt(n - 2) / sqrt(1 - r_vec^2)
      rho_vector <- numeric(length(indices_remain_to_test)) # rho=0 at start
      ss_X_cand <- colSums(X_cand^2)
      beta_cand <- nums / ss_X_cand
    } else {
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Q <- qr.Q(qr_inclus)
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      # Orthogonalize candidates (Gram-Schmidt step)
      XtX_cand <- t(Q) %*% X_cand
      X_cand_res <- X_cand - Q %*% XtX_cand
      ss_X_cand_res <- colSums(X_cand_res^2)
      ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector((t(Y_res) %*% X_cand_res) / ss_X_cand_res)
      # Variance estimation for Z-score
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
    # --- B. Posterior Probability (Analytic Integration) ---
    sd_eff <- sqrt(1 + n * (1 - rho_vector) * tau^2)
    integrals <- dnorm(Z_vector, mean = 0, sd = sd_eff)
    densities <- dnorm(Z_vector, mean = 0, sd = 1)
    p_vec <- a.priori[indices_remain_to_test]
    numerator <- p_vec * integrals
    denominator <- numerator + (1 - p_vec) * densities
    a.posteriori <- numerator / (denominator + 1e-20)
    # Pick Best
    max_idx_local <- which.max(a.posteriori)
    current_max_prob <- a.posteriori[max_idx_local]
    best_var_global <- indices_remain_to_test[max_idx_local]
    # --- C. Calculate Metrics for the POTENTIAL New Model ---
    # We need the specific stats of the best variable from the vector calculation above
    best_beta <- beta_cand[max_idx_local]
    if (k == 1) {
      best_ss_res <- colSums(X_cand^2)[max_idx_local]
    } else {
      best_ss_res <- colSums(X_cand_res^2)[max_idx_local]
    }
    RSS_new <- RSS_current - (best_beta^2 * best_ss_res)
    if(RSS_new < 0) RSS_new <- 1e-9 
    # Degrees of freedom: k variables + 1 intercept
    df_model <- k + 1 
    # Calculate Metrics
    # R2 Adj:
    r2_adj_val <- 1 - (RSS_new / (n - df_model)) / (TSS / (n - 1))
    # AIC: 
    aic_val <- n * log(RSS_new / n) + 2 * df_model
    # BIC: 
    bic_val <- n * log(RSS_new / n) + log(n) * df_model
    # --- D. Stopping Criteria Check ---
    stop_flag <- FALSE
    if (k > 1) {
      if (critere == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (critere == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (critere == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (critere == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    # --- E. Update Lists ---
    indices.inclus <- c(indices.inclus, best_var_global)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var_global)
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj_val)
    aic_hist <- c(aic_hist, aic_val)
    bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, current_max_prob)
    RSS_current <- RSS_new # Update for next iteration
    # Hard stop threshold if needed
    if (!is.null(threshold_K) && (length(indices.inclus) >= threshold_K)) break 
  }
  return(list(
    indices = indices.inclus,
    post.probs = post.probs,
    r_squared_adj = r_squared_adj_hist,
    aic_values = aic_hist,
    bic_values = bic_hist,
    K = length(indices.inclus)
  ))
}

## Asymptotic Wald
Asymptotic_Wald_Criteria_For_Data <- function(X, Y, a.priori, tau, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria ivalid, shoud be one of 'R2', 'AIC', 'BIC' ou 'Posterior'.")
  }
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
      if (critere == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (critere == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (critere == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (critere == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
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
  return(list(indices = indices.inclus,
              post.probs = post.probs,
              r_squared_adj = r_squared_adj_hist,
              aic_values = aic_hist,
              bic_values = bic_hist,
              K = length(indices.inclus)))
}

## Exact Wald
Exact_Wald_Criteria_For_Data <- function(X, Y, a.priori, tau, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria ivalid, shoud be one of 'R2', 'AIC', 'BIC' ou 'Posterior'.")
  }
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
      if (critere == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (critere == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (critere == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (critere == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    indices.inclus <- c(indices.inclus, indices_remain_to_test[best_idx_local])
    indices_remain_to_test <- setdiff(indices_remain_to_test, indices_remain_to_test[best_idx_local])
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj); aic_hist <- c(aic_hist, aic_val); bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
    
    if (!is.null(threshold_K) && (length(indices.inclus) >= threshold_K)) break 
  }
  return(list(indices = indices.inclus,
              post.probs = post.probs,
              r_squared_adj = r_squared_adj_hist,
              aic_values = aic_hist,
              bic_values = bic_hist,
              K = length(indices.inclus)))
}

## Asymptotic score
Asymptotic_Score_Criteria_For_Data <- function(X, Y, a.priori, tau, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria ivalid, shoud be one of 'R2', 'AIC', 'BIC' ou 'Posterior'.")
  }
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
      if (critere == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (critere == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (critere == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (critere == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
    }
    if (stop_flag) break
    indices.inclus <- c(indices.inclus, indices_remain_to_test[best_idx_local])
    indices_remain_to_test <- setdiff(indices_remain_to_test, indices_remain_to_test[best_idx_local])
    r_squared_adj_hist <- c(r_squared_adj_hist, r2_adj); aic_hist <- c(aic_hist, aic_val); bic_hist <- c(bic_hist, bic_val)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
    if (!is.null(threshold_K) && (length(indices.inclus) >= threshold_K)) break 
  }
  return(list(indices = indices.inclus,
              post.probs = post.probs,
              r_squared_adj = r_squared_adj_hist,
              aic_values = aic_hist,
              bic_values = bic_hist,
              K = length(indices.inclus)))
}
## Multi wald

Exact_Multiple_Criteria_For_Data <- function(X, Y, a.priori, tau2, S, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria ivalid, shoud be one of 'R2', 'AIC', 'BIC' ou 'Posterior'.")
  }
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
      if (critere == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (critere == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (critere == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (critere == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
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
  return(list(indices = indices.inclus,
              post.probs = post.probs,
              r_squared_adj = r_squared_adj_hist,
              aic_values = aic_hist,
              bic_values = bic_hist,
              K = length(indices.inclus)))
}

Naive_Multiple_Criteria_For_Data <- function(X, Y, a.priori, tau2, S, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  if (!criteria %in% c("R2", "AIC", "BIC", "Posterior")) {
    stop("Criteria ivalid, shoud be one of 'R2', 'AIC', 'BIC' ou 'Posterior'.")
  }
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
    # Criteria
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
      if (critere == "R2" && r2_adj < tail(r_squared_adj_hist, 1)) stop_flag <- TRUE
      if (critere == "AIC" && aic_val > tail(aic_hist, 1)) stop_flag <- TRUE
      if (critere == "BIC" && bic_val > tail(bic_hist, 1)) stop_flag <- TRUE
      if (critere == "Posterior" && best_prob < target_prob) stop_flag <- TRUE
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
  return(list(indices = indices.inclus,
              post.probs = post.probs,
              r_squared_adj = r_squared_adj_hist,
              aic_values = aic_hist,
              bic_values = bic_hist,
              K = length(indices.inclus)))
}

