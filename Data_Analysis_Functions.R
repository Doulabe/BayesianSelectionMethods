source("Test_statistics_functions.R")

#---- For fixed K (Table 6)
Real_Data_Comparison <- function(X, Y, a.priori, K, tau, S = 1000) {
  
  # --- Helper: Uniform Metric Calculator ---
  # Fits a linear model on selected indices to get comparable RMSE/AIC/BIC
  get_metrics <- function(y, x, indices, method_name) {
    if (length(indices) == 0) {
      # Null model
      fit <- lm(y ~ 1)
      k <- 0
    } else {
      # Subset X
      x_sub <- x[, indices, drop = FALSE]
      fit <- lm(y ~ x_sub)
      k <- length(indices)
    }
    
    summ <- summary(fit)
    y_hat <- fitted(fit)
    resids <- residuals(fit)
    n <- length(y)
    
    # Calculate Metrics
    r2_adj   <- summ$adj.r.squared
    
    # AIC/BIC (Calculated manually to ensure consistency)
    log_lik <- as.numeric(logLik(fit))
    aic_val <- -2 * log_lik + 2 * (k + 1) # +1 for intercept
    bic_val <- -2 * log_lik + log(n) * (k + 1)
    
    return(list(
      Method = method_name,
      K = k,
      R2_adj = r2_adj,
      AIC = aic_val,
      BIC = bic_val
    ))
  }
  
  # --- Setup ---
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  
  summary_list <- list()
  indices_list <- list()
  pip_list     <- list()
  
  # --- 1. Run Methods ---
  runs <- list(
    list(Name = "NaW",  Func = function() Naive_Wald(X, Y, a.priori, K, tau)),
    list(Name = "AsW",  Func = function() Asymptotic_Wald(X, Y, a.priori, K, tau)),
    list(Name = "ExW",  Func = function() Exact_Wald(X, Y, a.priori, K, tau)),
    list(Name = "AsS",  Func = function() Asymptotic_Score(X, Y, a.priori, K, tau)),
    list(Name = "NaM",  Func = function() Naive_Multiple(X, Y, a.priori, K, tau^2, S)),
    list(Name = "ExM",  Func = function() Exact_Multiple(X, Y, a.priori, K, tau^2, S))
  )
  
  cat("Processing Custom Methods...\n")
  for(run in runs) {
    cat(paste0("  Running ", run$Name, "...\n"))
    res <- run$Func()
    
    # Check if result structure is valid
    if(is.null(res) || is.null(res$indices)) {
      warning(paste(run$Name, "returned NULL or no indices."))
    } else {
      metrics <- get_metrics(Y, X, res$indices, run$Name)
      summary_list[[run$Name]] <- metrics
      indices_list[[run$Name]] <- res$indices
      pip_list[[run$Name]]     <- res$post.probs
    }
  }
  
  # --- 2. Final Compilation ---
  final_df <- do.call(rbind, lapply(summary_list, as.data.frame))
  return(list(Summary = final_df, Indices = indices_list, PIP = pip_list))
}

#---- For K determination criteria (Table 5)
Real_Data_Comparison_criteria <- function(X, Y, a.priori, tau, S = 1000, threshold_K = NULL, target_prob = 0.95, criteria = "BIC") {
  
  # --- 1. Basic Checks ---
  if (nrow(X) != length(Y)) stop("Dimensions of X and Y do not match.")
  # Ensure inputs are matrix/vector
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  
  # --- 2. Define the Methods ---
  # Note: Multi Wald methods use tau^2 based on previous logic
  runs <- list(
    list(Method = "NaW",  Func = function() Naive_Wald_Criteria_For_Data(X, Y, a.priori, tau, threshold_K , target_prob, criteria)),
    list(Method = "AsW",  Func = function() Asymptotic_Wald_Criteria_For_Data(X, Y, a.priori, tau, threshold_K , target_prob, criteria)),
    list(Method = "ExW",  Func = function() Exact_Wald_Criteria_For_Data(X, Y, a.priori, tau, threshold_K , target_prob, criteria)),
    list(Method = "AsS",  Func = function() Asymptotic_Score_Criteria_For_Data(X, Y, a.priori, tau, threshold_K , target_prob, criteria)),
    list(Method = "NaM",  Func = function() Naive_Multiple_Criteria_For_Data(X, Y, a.priori, tau^2, S, threshold_K , target_prob, criteria)),
    list(Method = "ExM",  Func = function() Exact_Multiple_Criteria_For_Data(X, Y, a.priori, tau^2, S, threshold_K , target_prob, criteria))
  )
  
  # Initialize storage
  summary_df <- data.frame(
    Method = character(),
    K = integer(),
    PIP = numeric(),
    R2_adj = numeric(),
    AIC = numeric(),
    BIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  indices_list <- list()
  # Helper to safely extract the last element of a vector (final model stats)
  safe_get <- function(vec) { 
    if (length(vec) > 0) return(tail(vec, 1)) else return(NA) 
  }
  cat("Processing Real Data...\n")
  # --- 3. Execution Loop ---
  for(run in runs) {
    method_Method <- run$Method
    method_func <- run$Func
    cat(sprintf("  Running %s...\n", method_Method))
    # Run the method safely
    res <- tryCatch({
      method_func()
    }, error = function(e) {
      warning(paste("Error in", method_Method, ":", e$message))
      return(NULL)
    })
    if (!is.null(res)) {
      # Extract scalars
      k_val   <- length(res$indices)
      pip_val <- safe_get(res$post.probs)
      r2_val  <- safe_get(res$r_squared_adj)
      aic_val <- safe_get(res$aic_values)
      bic_val <- safe_get(res$bic_values)
      
      # Add to Summary Table
      summary_df <- rbind(summary_df, data.frame(
        Method = method_Method,
        K = k_val,
        PIP = pip_val,
        R2_adj = r2_val,
        AIC = aic_val,
        BIC = bic_val
      ))
      # Store Indices separately
      indices_list[[method_Method]] <- res$indices
    } else {
      # Handle Failure
      summary_df <- rbind(summary_df, data.frame(
        Method = method_Method, K = NA, PIP = NA, R2_adj = NA, AIC = NA, BIC = NA
      ))
      indices_list[[method_Method]] <- NA
    }
  }
  return(list(Summary = summary_df, Selected_Indices = indices_list))
}

##-- Calculate Pbio
Calculate_Performance_Table <- function(results, map_data, known_targets) {
  
  # --- 1. Get Statistical Metrics ---
  # Create a copy to avoid modifying the original object
  stat_df <- results$Summary[, c("Method", "K")]
  
  # --- 2. Calculate Biological coverage (Pbio) ---
  known_targets <- as.character(known_targets)
  
  # Ensure map has 'Chr' column
  if(!all(c("chr", "pos") %in% colnames(map_data))) {
    colnames(map_data)[1:2] <- c("Chr", "Pos") 
  }
  
  bio_stats <- data.frame()
  # We loop through the names in Indices 
  methods_in_indices <- names(results$Indices)
  
  for(m in methods_in_indices) {
    idx <- results$Indices[[m]]
    
    if(length(idx) == 0 || is.na(idx[1])) {
      hits <- 0
      prec <- 0
    } else {
      # Get Chromosomes for selected SNPs
      # Ensure we look up the correct rows
      selected_chrs <- as.character(map_data[idx, "Chr"])
      
      # Count hits on known targets
      hits <- sum(selected_chrs %in% known_targets)
      
      # Calculate coverage %
      prec <- (hits / length(idx)) * 100
    }
    
    bio_stats <- rbind(bio_stats, data.frame(
      Method = m,
      Target_Hits = hits,
      Pbio_Pct = prec
    ))
  }
  
  # --- 3. Merge ---
  # all.x = TRUE ensures that even if a match is imperfect, we keep the stats
  final_df <- merge(stat_df, bio_stats, by = "Method", all = TRUE)
  
  # Clean up: Remove rows where Pbio are NA (if any)
  final_df <- final_df[!is.na(final_df$Pbio_Pct), ]
  
  # Rounding
  final_df$Pbio_Pct <- round(final_df$Pbio_Pct, 2)
  
  return(final_df)
}

##----  calculate informed priors
Get_Informed_Prior <- function(map_data, target_chroms, p_high = 0.8, p_low = 0.2) {
  
  # Ensure map has proper column names
  if(!all(c("chr", "pos") %in% colnames(map_data))) {
    colnames(map_data)[1:2] <- c("Chr", "Pos") 
  }
  
  # Initialize the prior vector with the low probability
  prior_vec <- rep(p_low, nrow(map_data))
  
  # Identify indices of target chromosomes
  # We convert to character to ensure "6" matches "6"
  target_indices <- which(as.character(map_data$Chr) %in% as.character(target_chroms))
  
  # Update those specific indices to the high probability
  prior_vec[target_indices] <- p_high
  
  return(prior_vec)
}

