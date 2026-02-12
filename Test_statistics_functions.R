######---- Mise a jour : Extension to complex null hypothesis (reprogramme)
Naive_Wald_Integrand <- function(lambda,Z,n,rho,tau){
  dnorm(Z,mean=lambda*sqrt(n*(1-rho)))*dnorm(lambda,sd=tau)
} 

Asymptotic_Wald_Integrand <- function(lambda,Z,n,rho,tau){
  mean.Z <- lambda*sqrt(n*(1-rho))
  var.Z <- 1 + 0.5*(1 - rho)*lambda^2
  dnorm(Z,mean = mean.Z, sd = sqrt(var.Z))*dnorm(lambda, sd = tau)
} 

Exact_Wald_Integrand <- function(lambda,Z,n,rho,tau,k){
  dt(Z, df=n-(k+1),ncp=lambda*sqrt(n*(1-rho)))*dnorm(lambda,sd=tau)
}

Asymptotic_Score_Integrand <- function(lambda,Q,n,rho,tau){
  mean.Q <- lambda*sqrt(n)*(1 - rho)/sqrt(1 + (1-rho)*lambda^2)
  var.Q <- (2*(1 - rho) + ((1-rho)^2)*lambda^2)/(2*(1 + (1-rho)*lambda^2)^3)
  dnorm(Q, mean=mean.Q, sd = sqrt(var.Q))*dnorm(lambda, sd=tau)
} 

# Exact Multiple test (ExM): Non central Fisher distribution
Fisher_MC_approximation <- function(Q_W, XtX_tilde, XtX_k, n, p, tau2, S) {
  stopifnot(is.matrix(XtX_tilde))
  rA  <- ncol(XtX_tilde)
  df1 <- rA
  df2 <- n - (rA + 1)
  
  # draw Lambda in matching dims
  Lambda_A <- matrix(rnorm(rA * S, sd = sqrt(tau2)), nrow = rA)
  delta_A  <- 0.5 * colSums(Lambda_A * (XtX_tilde %*% Lambda_A))
  delta_A  <- pmax(delta_A, 0)
  
  densA <- df(Q_W, df1 = df1, df2 = df2, ncp = delta_A)
  densA[!is.finite(densA)] <- 0
  A_hat <- mean(densA)
  
  # B term: central F if no previously included preds
  if (is.null(XtX_k) || ncol(XtX_k) == 0) {
    densB <- df(Q_W, df1 = df1, df2 = df2, ncp = 0)
    densB <- if (is.finite(densB)) densB else 0
    B_hat <- densB
  } else {
    rB       <- ncol(XtX_k)
    Lambda_B <- matrix(rnorm(rB * S, sd = sqrt(tau2)), nrow = rB)
    delta_B  <- 0.5 * colSums(Lambda_B * (XtX_k %*% Lambda_B))
    delta_B  <- pmax(delta_B, 0)
    densB <- df(Q_W, df1 = df1, df2 = df2, ncp = delta_B)
    densB[!is.finite(densB)] <- 0
    B_hat <- mean(densB)
  }
  
  # posterior: q = (1-p)A / ((1-p)A + pB)
  num <- p * A_hat
  den <- num +  (1 - p)* B_hat
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

ExM_posterior <- function(Q_W, X_tilde_jk, X_k, n, p, tau2, S) {
  XtX_tilde <- crossprod(X_tilde_jk)                 # rA x rA
  XtX_k     <- if (is.null(X_k)) NULL else crossprod(X_k)  # (rA-1) x (rA-1)
  Fisher_MC_approximation(Q_W = Q_W,
                         XtX_tilde = XtX_tilde, XtX_k = XtX_k,
                         n = n, p = p, tau2 = tau2, S = S)
}

##### Naive Multiple Test (NaM): Non central chi-square distribution
Chisq_MC_approximation <- function(Q_chi, XtX_tilde, XtX_k, n, p, tau2, S) {
  stopifnot(is.matrix(XtX_tilde))
  rA  <- ncol(XtX_tilde)                 # block size in the *full* model
  df  <- rA
  df2 <- n - (rA + 1)
  
  # A-hat: Λ_A ~ N(0, τ^2 I_{rA}), delta_A = 1/2 Λ' XtX_tilde Λ
  Lambda_A <- matrix(rnorm(rA * S, sd = sqrt(tau2)), nrow = rA)
  delta_A  <- 0.5 * colSums(Lambda_A * (XtX_tilde %*% Lambda_A))
  delta_A  <- pmax(delta_A, 0)
  
  densA <- dchisq(Q_chi, df = df, ncp = delta_A)
  densA[!is.finite(densA)] <- 0
  A_hat <- mean(densA)
  
  # B-hat: same chi-square df (df = rA), but noncentrality from previous block
  if (is.null(XtX_k) || ncol(XtX_k) == 0L) {
    # No previously included predictors -> central chi-square under B
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
  
  # Prior clamp and posterior combine: q = (1-p)A / ((1-p)A + pB)
  num <- p * A_hat
  den <- num + (1 - p) * B_hat
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

# Wrapper: accepts design matrices (X_k can be NULL)
NaM_posterior <- function(Q_chi, X_tilde_jk, X_k, n, p, tau2, S) {
  XtX_tilde <- crossprod(X_tilde_jk)                 # rA x rA
  XtX_k     <- if (is.null(X_k)) NULL else crossprod(X_k)  # (rA-1) x (rA-1)
  Chisq_MC_approximation(Q_chi = Q_chi, XtX_tilde = XtX_tilde, XtX_k = XtX_k,
                       n = n, p = p, tau2 = tau2, S = S)
}

##---- distributional and selection functions

# NaW
Naive_Wald <- function(X, Y, a.priori, K, tau) {
  m <- ncol(X)
  n <- nrow(X)
  
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  indices.inclus <- integer(0)
  indices_remain_to_test <- 1:m
  post.probs <- numeric(K)
  
  for (k in 1:K) {
    if (length(indices_remain_to_test) == 0) {
      if (k <= K) post.probs[k:K] <- NA
      break
    }
    
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      # Base case is simple correlation
      r <- cor(Y, X_cand)[1, ]
      r[is.na(r)] <- 0
      # Handle perfect correlation case where sqrt denominator is 0
      Z_vector <- ifelse(abs(r) == 1, sign(r) * Inf, r * sqrt(n - 2) / sqrt(1 - r^2))
      rho_vector <- numeric(length(indices_remain_to_test))
    } else {
      # General case using QR decomposition for speed and stability
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
      
      # --- CORRECTED RHO CALCULATION ---
      # This is the corrected implementation to ensure rho_vector has the right length.
      XtX_inclus_inv <- chol2inv(qr.R(qr_inclus))
      V_inv <- n * XtX_inclus_inv
      R_matrix <- t(X_inclus) %*% X_cand / n
      
      # Correctly compute diag(t(R) %*% V_inv %*% R) in a vectorized way
      M <- V_inv %*% R_matrix
      rho_vector <- rowSums(t(R_matrix) * t(M)) # The result is a vector of length p_cand
      
      # Clamp rho to prevent NaN errors from numerical imprecision
      rho_vector <- pmax(0, pmin(rho_vector, 1 - 1e-9))
    }
    
    # Handle potential Inf/-Inf in Z_vector (e.g., from perfect correlation)
    Z_vector[!is.finite(Z_vector)] <- 0
    
    integrals <- mapply(function(z, r) {
      tryCatch(
        integrate(Naive_Wald_Integrand, lower = -Inf, upper = Inf, Z = z, n = n, rho = r, tau = tau)$value,
        error = function(e) NA
      )
    }, Z_vector, rho_vector)
    
    integrals[is.na(integrals)] <- 0
    
    densities <- dnorm(Z_vector)
    p_vector <- a.priori[indices_remain_to_test]
    
    numerator <- p_vector * integrals
    denominator <- numerator + (1 - p_vector) * densities
    a.posteriori <- numerator / (denominator + 1e-20) 
    
    # Handle cases where all posterior probabilities might be NaN or NA
    if(all(!is.finite(a.posteriori))) {
      # If all calculations fail, stop and return what we have so far
      if (k <= K) post.probs[k:K] <- NA
      break
    }
    max_idx_local <- which.max(a.posteriori)
    
    post.probs[k] <- a.posteriori[max_idx_local]
    
    a.inclure <- indices_remain_to_test[max_idx_local]
    indices.inclus <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }
  
  return(list(post.probs = post.probs, indices = indices.inclus))
}

# AsW
Asymptotic_Wald <- function(X, Y, a.priori, K, tau) {
  m <- ncol(X)
  n <- nrow(X)
  
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  indices.inclus <- integer(0)
  indices_remain_to_test <- 1:m
  post.probs <- numeric(K)
  
  for (k in 1:K) {
    if (length(indices_remain_to_test) == 0) {
      if (k <= K) post.probs[k:K] <- NA
      break
    }
    
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      # Base case is simple correlation
      r <- cor(Y, X_cand)[1, ]
      r[is.na(r)] <- 0
      # Handle perfect correlation case where sqrt denominator is 0
      Z_vector <- ifelse(abs(r) == 1, sign(r) * Inf, r * sqrt(n - 2) / sqrt(1 - r^2))
      rho_vector <- numeric(length(indices_remain_to_test))
    } else {
      # General case using QR decomposition for speed and stability
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
      
      
      # corrected implementation to ensure rho_vector has the right length.
      XtX_inclus_inv <- chol2inv(qr.R(qr_inclus))
      V_inv <- n * XtX_inclus_inv
      R_matrix <- t(X_inclus) %*% X_cand / n
      
      # Correctly compute diag(t(R) %*% V_inv %*% R) in a vectorized way
      M <- V_inv %*% R_matrix
      rho_vector <- rowSums(t(R_matrix) * t(M)) # The result is a vector of length p_cand
      # --- END OF CORRECTION ---
      
      # Clamp rho to prevent NaN errors from numerical imprecision
      rho_vector <- pmax(0, pmin(rho_vector, 1 - 1e-9))
    }
    
    # Handle potential Inf/-Inf in Z_vector (e.g., from perfect correlation)
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
    
    # Handle cases where all posterior probabilities might be NaN or NA
    if(all(!is.finite(a.posteriori))) {
      # If all calculations fail, stop and return what we have so far
      if (k <= K) post.probs[k:K] <- NA
      break
    }
    max_idx_local <- which.max(a.posteriori)
    
    post.probs[k] <- a.posteriori[max_idx_local]
    
    a.inclure <- indices_remain_to_test[max_idx_local]
    indices.inclus <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }
  
  return(list(post.probs = post.probs, indices = indices.inclus))
}

#ExW
Exact_Wald <- function(X, Y, a.priori, K, tau) {
  m <- ncol(X)
  n <- nrow(X)
  
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  indices.inclus <- integer(0)
  indices_remain_to_test <- 1:m
  post.probs <- numeric(K)
  
  for (k in 1:K) {
    if (length(indices_remain_to_test) == 0) {
      if (k <= K) post.probs[k:K] <- NA
      break
    }
    
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    # --- Vectorized calculations for all candidates ---
    if (k == 1) {
      # Base case: rho is 0, Z.hat is the standard t-statistic
      r <- cor(Y, X_cand)[1, ]
      r[is.na(r)] <- 0
      # Handle perfect correlation case where sqrt denominator is 0
      Z_hat_vector <- ifelse(abs(r) == 1, sign(r) * Inf, r * sqrt(n - 2) / sqrt(1 - r^2))
      rho_vector <- numeric(length(indices_remain_to_test)) # rho is 0
    } else {
      # General case for k > 1 using QR decomposition
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
      
      # Vectorized rho calculation (same as before)
      XtX_inclus_inv <- chol2inv(qr.R(qr_inclus))
      V_inv <- n * XtX_inclus_inv
      R_matrix <- t(X_inclus) %*% X_cand / n
      M <- V_inv %*% R_matrix
      rho_vector <- rowSums(t(R_matrix) * t(M))
      
      # Clamp rho to prevent numerical errors
      rho_vector <- pmax(0, pmin(rho_vector, 1 - 1e-9))
      
      # --- Key difference for exact.wald ---
      # Scale the t-statistic to get the final Z.hat
      Z_hat_vector <- Z_vector * sqrt(1 - rho_vector)
    }
    
    # Clean up any non-finite values before integration
    Z_hat_vector[!is.finite(Z_hat_vector)] <- 0
    
    # --- Posterior probability calculation ---
    integrals <- mapply(function(z, r) {
      tryCatch(
        integrate(Exact_Wald_Integrand, lower = -Inf, upper = Inf, Z = z, n = n, rho = r, tau = tau, k = k)$value,
        error = function(e) NA
      )
    }, Z_hat_vector, rho_vector)
    
    integrals[is.na(integrals)] <- 0
    
    # The null distribution is dt() with appropriate degrees of freedom
    densities <- dt(Z_hat_vector, df = n - (k + 1))
    p_vector <- a.priori[indices_remain_to_test]
    
    numerator <- p_vector * integrals
    denominator <- numerator + (1 - p_vector) * densities
    a.posteriori <- numerator / (denominator + 1e-20) 
    
    if(all(!is.finite(a.posteriori))) {
      if (k <= K) post.probs[k:K] <- NA
      break
    }
    max_idx_local <- which.max(a.posteriori)
    
    post.probs[k] <- a.posteriori[max_idx_local]
    
    a.inclure <- indices_remain_to_test[max_idx_local]
    indices.inclus <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }
  
  return(list(post.probs = post.probs, indices = indices.inclus))
}

#AsS
Asymptotic_Score <- function(X, Y, a.priori, K, tau) {
  m <- ncol(X)
  n <- nrow(X)
  
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  indices.inclus <- integer(0)
  indices_remain_to_test <- 1:m
  post.probs <- numeric(K)
  
  for (k in 1:K) {
    if (length(indices_remain_to_test) == 0) {
      if (k <= K) post.probs[k:K] <- NA
      break
    }
    
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    if (k == 1) {
      # Vectorized calculation for the first variable
      numerator <- as.vector(t(Y) %*% X_cand)
      
      # Vectorized sigma calculation
      beta <- numerator / colSums(X_cand^2)
      
      # --- FIX: Use sum() for a vector, not colSums() ---
      RSS <- sum(Y^2) - beta * numerator
      
      sigma_vec <- sqrt(RSS / (n - 2))
      
      Q_hat_vector <- numerator / (sqrt(n) * sigma_vec)
      rho_vector <- numeric(length(indices_remain_to_test)) # rho is 0
    } else {
      # General case using QR decomposition
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      Q <- qr.Q(qr_inclus)
      
      Y_res <- Y - Q %*% (t(Q) %*% Y)
      X_cand_res <- X_cand - Q %*% (t(Q) %*% X_cand)
      
      # Vectorized Q.hat calculation
      numerator <- as.vector(t(Y) %*% X_cand_res)
      
      # The sigma for each potential model is needed for the denominator
      ss_X_cand_res <- colSums(X_cand_res^2)
      ss_X_cand_res[ss_X_cand_res < 1e-9] <- 1e-9
      beta_cand <- as.vector((t(Y) %*% X_cand_res) / ss_X_cand_res)
      RSS_current <- sum(Y_res^2)
      RSS_new <- RSS_current - beta_cand^2 * ss_X_cand_res
      sigma_new_vec <- sqrt(RSS_new / (n - k - 1))
      
      Q_hat_vector <- numerator / (sqrt(n) * sigma_new_vec)
      
      # Vectorized rho calculation
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
    
    if(all(!is.finite(a.posteriori))) {
      if (k <= K) post.probs[k:K] <- NA
      break
    }
    max_idx_local <- which.max(a.posteriori)
    
    post.probs[k] <- a.posteriori[max_idx_local]
    
    a.inclure <- indices_remain_to_test[max_idx_local]
    indices.inclus <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }
  
  return(list(post.probs = post.probs, indices = indices.inclus))
}

#NaM
Naive_Multiple <- function(X, Y, a.priori, K, tau2, S) {
  m <- ncol(X)
  n <- nrow(X)
  
  # --- Crucial Step: Center data to align model assumptions ---
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  
  indices.inclus <- integer(0)
  indices_remain_to_test <- seq_len(m)
  post.probs <- rep(NA_real_, K)
  
  # Total Sum of Squares for the centered Y
  TSS <- sum(Y^2)
  
  for (k in 1:K) {
    if (length(indices_remain_to_test) == 0) {
      if (k <= K) post.probs[k:K] <- NA
      break
    }
    
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    rA <- k # Number of predictors in the candidate models
    df1 <- rA
    df2 <- n - rA - 1L # Degrees of freedom for error (assumes intercept)
    if (df2 <= 0) break
    
    # --- Vectorized F-statistic calculation (same as fast_multi.wald) ---
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
    
    # --- Key difference for naive.multi.wald ---
    # Convert the vector of F-stats to a vector of Chi-squared stats
    Q_chi_vector <- rA * F_block_vector
    
    Q_chi_vector[!is.finite(Q_chi_vector) | Q_chi_vector < 0] <- -1
    
    # --- Loop only for the Monte Carlo approximation part ---
    p_vector <- a.priori[indices_remain_to_test]
    
    a.posteriori <- mapply(
      function(q_chi, j, p) {
        if (q_chi < 0) return(NA_real_)
        
        X.modele <- if (is.null(X_inclus)) X[, j, drop = FALSE] else cbind(X_inclus, X[, j, drop = FALSE])
        
        # Call the naive version of the posterior function
        NaM_posterior(
          Q_chi = q_chi,
          X_tilde_jk = X.modele,
          X_k = X_inclus,
          n = n, p = p, tau2 = tau2, S = S
        )
      },
      Q_chi_vector, indices_remain_to_test, p_vector
    )
    
    # --- Select best variable and update for the next iteration ---
    if (all(!is.finite(a.posteriori))) { break }
    
    best_idx_local <- which.max(a.posteriori)
    post.probs[k] <- a.posteriori[best_idx_local]
    
    a.inclure <- indices_remain_to_test[best_idx_local]
    indices.inclus <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }
  
  list(post.probs = post.probs, indices = indices.inclus)
}

#ExM
Exact_Multiple <- function(X, Y, a.priori, K, tau2, S) {
  m <- ncol(X)
  n <- nrow(X)
  
  # --- FIX: Center data to make no-intercept model equivalent to intercept model ---
  # This makes the vectorized F-statistic calculation mathematically sound.
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)
  # --- END OF FIX ---
  
  indices.inclus <- integer(0)
  indices_remain_to_test <- seq_len(m)
  post.probs <- rep(NA_real_, K)
  
  # Total Sum of Squares is now simply sum of squares of the centered Y
  TSS <- sum(Y^2)
  
  for (k in 1:K) {
    if (length(indices_remain_to_test) == 0) {
      if (k <= K) post.probs[k:K] <- NA
      break
    }
    
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    
    rA <- k
    df1 <- rA
    df2 <- n - rA - 1L # Degrees of freedom for error in a model with an intercept
    if (df2 <= 0) break
    
    if (k == 1) {
      # For centered data, cor^2 * (n-2) / (1-r^2) is equivalent to the F-stat
      r_sq <- (as.vector(t(Y) %*% X_cand) / sqrt(TSS * colSums(X_cand^2)))^2
      r_sq[is.na(r_sq)] <- 0
      Z_vector <- r_sq * (n - 2) / (1 - r_sq)
      X_inclus <- NULL
    } else {
      X_inclus <- X[, indices.inclus, drop = FALSE]
      qr_inclus <- qr(X_inclus)
      
      # Residuals after regressing on already included variables
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
    a.posteriori <- mapply(
      function(z, j, p) {
        if (z < 0) return(NA_real_)
        X.modele <- if (is.null(X_inclus)) X[, j, drop = FALSE] else cbind(X_inclus, X[, j, drop = FALSE])
        ExM_posterior(
          Q_W = z,
          X_tilde_jk = X.modele,
          X_k = X_inclus,
          n = n, p = p, tau2 = tau2, S = S
        )
      },
      Z_vector, indices_remain_to_test, p_vector
    )
    
    if (all(!is.finite(a.posteriori))) { break }
    
    best_idx_local <- which.max(a.posteriori)
    post.probs[k] <- a.posteriori[best_idx_local]
    
    a.inclure <- indices_remain_to_test[best_idx_local]
    indices.inclus <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }
  
  list(post.probs = post.probs, indices = indices.inclus)
}


#### Parallel simulation avec nouvelles focnctions

Simulation <- function(X, p_values, tau_values, a.priori, alpha, K, B, S, threshold = 16) {
  m <- ncol(X)
  n <- nrow(X)
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  clusterEvalQ(cl, {
             library(doParallel)
             library(MASS)
             library(mvtnorm)
             library(stats)
         })
       clusterExport(cl, c(
             "Generate.Y", "Naive_Wald", "Asymptotic_Wald", "Exact_Wald", "Asymptotic_Score", "Exact_Multiple", 
             "Naive_Multiple", "Power_Fun", "Naive_Wald_Integrand", "Asymptotic_Wald_Integrand",
             "Exact_Wald_Integrand", "Asymptotic_Score_Integrand", "Fisher_MC_approximation",
             "Chisq_MC_approximation", "ExM_posterior", "NaM_posterior"
         ))
  
  param_grid <- expand.grid(p = p_values, tau = tau_values, stringsAsFactors = FALSE)
  
  final_main_results <- list()
  final_selection_summaries <- list()
  
  method_names <- c("NaW", "AsW", "ExW", 
                    "AsS", "NaM", "ExM")
  
  for (i in seq_len(nrow(param_grid))) {
    p <- param_grid$p[i]
    tau <- param_grid$tau[i]
    
    b_val <- uniroot(Power_Fun, lower = 0, upper = 1, n = n, alpha = alpha, p = p)$root
    beta <- c(1, rep(b_val, 15), rep(0, m - 15))
    sigma <- 1
    
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
      
      # Return both: count under threshold and count true variables selected (1:15)
      c(
        sapply(methods_results, function(res) sum(res$indices < threshold)),
        sapply(methods_results, function(res) sum(res$indices %in% 1:15))
      )
    }
    
    # Separate columns
    mean_counts <- colMeans(mc_results[, 1:6])
    mean_results <- data.frame(p = p, tau = tau)
    for (j in seq_along(method_names)) {
      mean_results[[paste0("mean_", method_names[j])]] <- mean_counts[j]
    }
    
    final_main_results[[i]] <- mean_results
    
    # Create selection summary for each method (0 to 15 selected true vars)
    for (j in seq_along(method_names)) {
      method_col <- mc_results[, 6 + j]  # Because the second block starts at column 7
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

Get_Superior_Probability <- function(x, y) {
  sum(sapply(0:15, function(i) {
    xi <- x[i + 1]  # x0 is at index 1
    sum(sapply(i:15, function(j) {
      yj <- y[j + 1]
      weight <- if (i == j) 0.5 else 1
      xi * yj * weight
    }))
  }))
}

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
#comp_res$`p=0.2_tau=0.2`


