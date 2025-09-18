library(dplyr)
library(ivmodel)
source("estimator_helpers.R")
# implement unadjusted estimator from Chen and Li

unadjusted_ANCOVA <- function(data, Ij = NA, estimand = "ind",
                              alpha = 0.05, variance = "DB") {
  # number of clusters
  clusters <- as.integer(unique(data$cluster))
  I <- length(unique(data$cluster))
  # number of periods
  # assume periods numbered 0 thru J+1
  J <- max(data$period) - 1
  periods <- 1:J
  
  # initialize Ij
  if(is.na(Ij)) {
    result <- data %>%
      filter(treatment == 1 & period != 0 & period != (J+1)) %>%
      group_by(period) %>%
      summarise(unique_clusters = n_distinct(cluster))
    Ij <- result$unique_clusters
    
  }
  
  I_tr_seq <- Ij
  I_ct_seq <- I-I_tr_seq
  I_a <- I/(J+1)
  
  if (estimand == "ind") {
    data$weights <- rep(1, nrow(data))
  } else if (estimand == "cell") {
    data <- data %>%
      group_by(cluster, period) %>%
      mutate(weights = 1 / n()) %>%
      ungroup()
    
  } else if (estimand == "cell") {
    data <- data %>%
      group_by(period) %>%
      mutate(weights = 1 / n()) %>%
      ungroup()
    
  }

  # discard pre and post rollout periods
  data_middle <- data %>% filter(period != 0 & period != (J+1))
  
  # make schedule column as needed for variance estimation
  data_middle <- data_middle %>%
    group_by(cluster) %>%
    mutate(schedule = ifelse(any(treatment == 1), min(period[treatment == 1]), J + 1)) %>%
    ungroup()
  
  # Convert period to a factor if it's not already
  data_middle$period <- as.factor(data_middle$period)
  
  # Create dummy variables for period
  dummy_period <- model.matrix(~ period - 1, data = data_middle)
  
  # Get the names of the period dummy variables
  period_vars <- colnames(dummy_period)
  
  # Create the interaction terms
  interaction_terms <- paste(period_vars, "treatment", sep = ":")
  
  # Combine the period dummy variables and interaction terms
  all_terms <- c(period_vars, interaction_terms)
  
  # Add the dummy variables to the original data frame
  data_middle <- cbind(data_middle, dummy_period)
  
  # fit linear model
  # Dynamically create the formula
  formula <- as.formula(paste("outcome ~", paste(all_terms, collapse = " + "), "- 1"))
  
  # Run the regression model with the dynamically created formula
  model <- lm(formula, weights = weights, data = data_middle)
  
  # reconvert period to numeric
  data_middle$period <- as.numeric(as.character(data_middle$period))
  
  # compute the point estimate
  varpi <- tapply(data_middle$weights, data_middle$period, sum) / sum(data_middle$weights)
  point_estimate <- coef(model)[(J+1):(2*J)]%*%varpi
  
  # variance estimation
  ind_mat <- cbind(data_middle[,period_vars],
                      data_middle$treatment*data_middle[,period_vars])
  
  res_x <- residuals(model) + as.matrix(ind_mat) %*% coef(model)
  
  covDB <- plug_in_cov(data = data_middle, residual = res_x, 
                       weights = data_middle$weights, I = I, J = J,
                       I_tr_seq = I_tr_seq, I_ct_seq = I_ct_seq, I_a = I_a) 
  DB_sd_estimate <- sqrt(c(varpi %*% covDB %*% varpi))
  
  covCR0 <- vcovCR(model, cluster = data_middle$cluster, type = "CR0")[(J+1):(2*J),(J+1):(2*J)]
  CR0_sd_estimate <- sqrt(c(varpi %*% covCR0 %*% varpi))
  
  covCR3 <- try({
    vcovCR(model, cluster = data_middle$cluster, type = "CR3")[(J+1):(2*J),(J+1):(2*J)]
  }, silent = FALSE)
  
  if (inherits(covCR3, "try-error")) {
    print("CR3 fail")
    covCR3 <- try({
      vcovCR.pseudo.lm(model, cluster = data_middle$cluster, type = "CR3")[(J+1):(2*J),(J+1):(2*J)]
    }, silent = FALSE)
    if (inherits(covCR3, "try-error")) {
      print("CR3 psuedo fail")
      covCR3 <- NA
      CR3_sd_estimate <- NA
    } else {
      CR3_sd_estimate <- sqrt(c(varpi %*% covCR3 %*% varpi))
    }
  } else {
    CR3_sd_estimate <- sqrt(c(varpi %*% covCR3 %*% varpi))
  }

  return(list(est = point_estimate, DB_sd_est = DB_sd_estimate, 
              DB_deviate = point_estimate / DB_sd_estimate,
              CR0_sd_est = CR0_sd_estimate, CR0_deviate = point_estimate / CR0_sd_estimate,
              CR3_sd_est = CR3_sd_estimate, CR3_deviate = point_estimate / CR3_sd_estimate,
              DB_CI = c(point_estimate - qnorm(1 - alpha / 2) * DB_sd_estimate, 
                        point_estimate + qnorm(1 - alpha / 2) * DB_sd_estimate),
              CR0_CI = c(point_estimate - qnorm(1 - alpha / 2) * CR0_sd_estimate, 
                         point_estimate + qnorm(1 - alpha / 2) * CR0_sd_estimate),
              CR3_CI = c(point_estimate - qnorm(1 - alpha / 2) * CR3_sd_estimate, 
                         point_estimate + qnorm(1 - alpha / 2) * CR3_sd_estimate)))
}


# implement ANCOVAI from Chen and Li
# compute weighted means of individual covariates

ANCOVAI <- function(data, Ij = NA, estimand = "ind",
                              alpha = 0.05) {
  # number of clusters
  clusters <- as.integer(unique(data$cluster))
  I <- length(unique(data$cluster))
  # number of periods
  # assume periods numbered 0 thru J+1
  J <- max(data$period) - 1
  periods <- 1:J
  
  # initialize Ij
  if(is.na(Ij)) {
    result <- data %>%
      filter(treatment == 1 & period != 0 & period != (J+1)) %>%
      group_by(period) %>%
      summarise(unique_clusters = n_distinct(cluster))
    Ij <- result$unique_clusters
    
  }
  
  I_tr_seq <- Ij
  I_ct_seq <- I-I_tr_seq
  I_a <- I/(J+1)
  
  if (estimand == "ind") {
    data$weights <- rep(1, nrow(data))
  } else if (estimand == "cell") {
    data <- data %>%
      group_by(cluster, period) %>%
      mutate(weights = 1 / n()) %>%
      ungroup()
    
  } else if (estimand == "cell") {
    data <- data %>%
      group_by(period) %>%
      mutate(weights = 1 / n()) %>%
      ungroup()
    
  }
  
  # discard pre and post rollout periods
  data_middle <- data %>% filter(period != 0 & period != (J+1))
  
  # make schedule column as needed for variance estimation
  data_middle <- data_middle %>%
    group_by(cluster) %>%
    mutate(schedule = ifelse(any(treatment == 1), min(period[treatment == 1]), J + 1)) %>%
    ungroup()
  
  # Get the names of the original covariate columns
  original_cols <- colnames(data_middle)[!(colnames(data_middle) %in% 
                                             c("cluster", "period", "weights", 
                                               "treatment", "outcome", "schedule"))]
  
  # Create a list of centered column names
  centered_cols <- paste0("centered_", original_cols)
  
  # make mean covariate columns and the centered covariate columns
  data_middle <- data_middle %>%
    group_by(period) %>%
    mutate(across(all_of(original_cols), 
                  ~ weighted.mean(.x, weights, na.rm = TRUE), 
                  .names = "weighted_mean_{col}")) %>%
    mutate(across(all_of(original_cols), 
                  ~ .x - get(paste0("weighted_mean_", cur_column())), 
                  .names = "centered_{col}")) %>%
    dplyr::select(all_of(centered_cols), cluster, period, weights, treatment, outcome, schedule) %>%
    ungroup()
  
  # Convert period to a factor if it's not already
  data_middle$period <- as.factor(data_middle$period)
  
  # Create dummy variables for period
  dummy_period <- model.matrix(~ period - 1, data = data_middle)
  
  # Get the names of the period dummy variables
  period_vars <- colnames(dummy_period)
  
  # Create the interaction terms
  interaction_terms <- paste(period_vars, "treatment", sep = ":")
  
  # Get the names of the centered covariate terms
  centered_column_names <- data_middle %>%
    dplyr::select(starts_with("centered_")) %>%
    names()
  
  # Combine the period dummy variables and interaction terms
  all_terms <- c(centered_column_names, period_vars, interaction_terms)
  
  # Add the dummy variables to the original data frame
  data_middle <- cbind(data_middle, dummy_period)
  
  # fit linear model
  # Dynamically create the formula
  formula <- as.formula(paste("outcome ~", paste(all_terms, collapse = " + "), "- 1"))
  
  # Run the regression model with the dynamically created formula
  model <- lm(formula, weights = weights, data = data_middle)
  
  # reconvert period to numeric
  data_middle$period <- as.numeric(as.character(data_middle$period))
  
  # compute the point estimate
  x_size <- length(centered_column_names)
  varpi <- tapply(data_middle$weights, data_middle$period, sum) / sum(data_middle$weights)
  point_estimate <- coef(model)[(x_size + J+1):(x_size + 2*J)]%*%varpi
  
  # variance estimation
  ind_mat <- cbind(data_middle[,period_vars],
                   data_middle$treatment*data_middle[,period_vars])
  
  res_x <- residuals(model) + as.matrix(ind_mat) %*% coef(model)[-(1:x_size)]
  
  covDB <- plug_in_cov(data = data_middle, residual = res_x, 
                       weights = data_middle$weights, I = I, J = J,
                       I_tr_seq = I_tr_seq, I_ct_seq = I_ct_seq, I_a = I_a) 
  DB_sd_estimate <- sqrt(c(varpi %*% covDB %*% varpi))
  
  covCR0 <- vcovCR(model, cluster = data_middle$cluster, type = "CR0")[(x_size + J+1):(x_size + 2*J),
                                                                      (x_size + J+1):(x_size + 2*J)]
  CR0_sd_estimate <- sqrt(c(varpi %*% covCR0 %*% varpi))
  
  covCR3 <- try({
    vcovCR(model, cluster = data_middle$cluster, type = "CR3")[
      (x_size + J + 1):(x_size + 2 * J),
      (x_size + J + 1):(x_size + 2 * J)
    ]
  }, silent = FALSE)
  
  if (inherits(covCR3, "try-error")) {
    print("CR3 fail")
    covCR3 <- try({
      vcovCR.pseudo.lm(model, cluster = data_middle$cluster, type = "CR3")[
        (x_size + J + 1):(x_size + 2 * J),
        (x_size + J + 1):(x_size + 2 * J)
      ]
    }, silent = FALSE)
    if (inherits(covCR3, "try-error")) {
      print("CR3 psuedo fail")
      covCR3 <- NA
      CR3_sd_estimate <- NA
    } else {
      CR3_sd_estimate <- sqrt(c(varpi %*% covCR3 %*% varpi))
    }
  } else {
    CR3_sd_estimate <- sqrt(c(varpi %*% covCR3 %*% varpi))
  }
  
  return(list(est = point_estimate, DB_sd_est = DB_sd_estimate, 
              DB_deviate = point_estimate / DB_sd_estimate,
              CR0_sd_est = CR0_sd_estimate, CR0_deviate = point_estimate / CR0_sd_estimate,
              CR3_sd_est = CR3_sd_estimate, CR3_deviate = point_estimate / CR3_sd_estimate,
              DB_CI = c(point_estimate - qnorm(1 - alpha / 2) * DB_sd_estimate, 
                        point_estimate + qnorm(1 - alpha / 2) * DB_sd_estimate),
              CR0_CI = c(point_estimate - qnorm(1 - alpha / 2) * CR0_sd_estimate, 
                        point_estimate + qnorm(1 - alpha / 2) * CR0_sd_estimate),
              CR3_CI = c(point_estimate - qnorm(1 - alpha / 2) * CR3_sd_estimate, 
                         point_estimate + qnorm(1 - alpha / 2) * CR3_sd_estimate)))
}



# implement ANCOVAIII from Chen and Li
ANCOVAIII <- function(data, Ij = NA, estimand = "ind",
                    alpha = 0.05) {
  # number of clusters
  clusters <- as.integer(unique(data$cluster))
  I <- length(unique(data$cluster))
  # number of periods
  # assume periods numbered 0 thru J+1
  J <- max(data$period) - 1
  periods <- 1:J
  
  # initialize Ij
  if(is.na(Ij)) {
    result <- data %>%
      filter(treatment == 1 & period != 0 & period != (J+1)) %>%
      group_by(period) %>%
      summarise(unique_clusters = n_distinct(cluster))
    Ij <- result$unique_clusters
    
  }
  
  I_tr_seq <- Ij
  I_ct_seq <- I-I_tr_seq
  I_a <- I/(J+1)
  
  if (estimand == "ind") {
    data$weights <- rep(1, nrow(data))
  } else if (estimand == "cell") {
    data <- data %>%
      group_by(cluster, period) %>%
      mutate(weights = 1 / n()) %>%
      ungroup()
    
  } else if (estimand == "cell") {
    data <- data %>%
      group_by(period) %>%
      mutate(weights = 1 / n()) %>%
      ungroup()
    
  }
  
  # discard pre and post rollout periods
  data_middle <- data %>% filter(period != 0 & period != (J+1))
  
  # make schedule column as needed for variance estimation
  data_middle <- data_middle %>%
    group_by(cluster) %>%
    mutate(schedule = ifelse(any(treatment == 1), min(period[treatment == 1]), J + 1)) %>%
    ungroup()
  
  # Get the names of the original covariate columns
  original_cols <- colnames(data_middle)[!(colnames(data_middle) %in% 
                                  c("cluster", "period", "weights", 
                                    "treatment", "outcome", "schedule"))]
  
  # Create a list of centered column names
  centered_cols <- paste0("centered_", original_cols)
  
  # make mean covariate columns and the centered covariate columns
  data_middle <- data_middle %>%
    group_by(period) %>%
    mutate(across(all_of(original_cols), 
                  ~ weighted.mean(.x, weights, na.rm = TRUE), 
                  .names = "weighted_mean_{col}")) %>%
    mutate(across(all_of(original_cols), 
                  ~ .x - get(paste0("weighted_mean_", cur_column())), 
                  .names = "centered_{col}")) %>%
    dplyr::select(all_of(centered_cols), cluster, period, weights, treatment, outcome, schedule) %>%
    ungroup()
  
  # Convert period to a factor if it's not already
  data_middle$period <- as.factor(data_middle$period)
  
  # Create dummy variables for period
  dummy_period <- model.matrix(~ period - 1, data = data_middle)
  
  # Get the names of the period dummy variables
  period_vars <- colnames(dummy_period)
  
  # Create the interaction terms
  interaction_terms <- paste(period_vars, "treatment", sep = ":")
  
  # Get the names of the centered covariate terms
  centered_column_names <- data_middle %>%
    dplyr::select(starts_with("centered_")) %>%
    names()
  
  # Create covariate interaction terms
  covariate_interaction_terms <- paste(centered_column_names, "treatment", sep = ":")
  
  # Combine the period dummy variables and interaction terms
  all_terms <- c(centered_column_names, covariate_interaction_terms, period_vars, interaction_terms)
  
  # Add the dummy variables to the original data frame
  data_middle <- cbind(data_middle, dummy_period)
  
  # fit linear model
  # Dynamically create the formula
  formula <- as.formula(paste("outcome ~", paste(all_terms, collapse = " + "), "- 1"))
  
  # Run the regression model with the dynamically created formula
  model <- lm(formula, weights = weights, data = data_middle)
  
  # reconvert period to numeric
  data_middle$period <- as.numeric(as.character(data_middle$period))
  
  # compute the point estimate
  x_size <- length(centered_column_names)
  varpi <- tapply(data_middle$weights, data_middle$period, sum) / sum(data_middle$weights)
  point_estimate <- coef(model)[(2 * x_size + J+1):(2 * x_size + 2*J)]%*%varpi
  
  # variance estimation
  ind_mat <- cbind(data_middle[,period_vars],
                   data_middle$treatment*data_middle[,period_vars])
  
  res_x <- residuals(model) + as.matrix(ind_mat) %*% coef(model)[-(1:(2 * x_size))]
  
  covDB <- plug_in_cov(data = data_middle, residual = res_x, 
                       weights = data_middle$weights, I = I, J = J,
                       I_tr_seq = I_tr_seq, I_ct_seq = I_ct_seq, I_a = I_a) 
  DB_sd_estimate <- sqrt(c(varpi %*% covDB %*% varpi))
  
  covCR0 <- vcovCR(model, cluster = data_middle$cluster, type = "CR0")[(2 * x_size + J+1):(2 * x_size + 2*J),
                                                                      (2 * x_size + J+1):(2 * x_size + 2*J)]
  CR0_sd_estimate <- sqrt(c(varpi %*% covCR0 %*% varpi))
  
  covCR3 <- try({
    vcovCR(model, cluster = data_middle$cluster, type = "CR3")[
      (2 * x_size + J + 1):(2 * x_size + 2 * J),
      (2 * x_size + J + 1):(2 * x_size + 2 * J)
    ]
  }, silent = FALSE)
  
  if (inherits(covCR3, "try-error")) {
    print("CR3 fail")
    covCR3 <- try({
      vcovCR.pseudo.lm(model, cluster = data_middle$cluster, type = "CR3")[
        (2 * x_size + J + 1):(2 * x_size + 2 * J),
        (2 * x_size + J + 1):(2 * x_size + 2 * J)
      ]
    }, silent = FALSE)
    if (inherits(covCR3, "try-error")) {
      print("CR3 psuedo fail")
      covCR3 <- NA
      CR3_sd_estimate <- NA
    } else {
      CR3_sd_estimate <- sqrt(c(varpi %*% covCR3 %*% varpi))
    }
  } else {
    CR3_sd_estimate <- sqrt(c(varpi %*% covCR3 %*% varpi))
  }
  
  return(list(est = point_estimate, DB_sd_est = DB_sd_estimate, 
              DB_deviate = point_estimate / DB_sd_estimate,
              CR0_sd_est = CR0_sd_estimate, CR0_deviate = point_estimate / CR0_sd_estimate,
              CR3_sd_est = CR3_sd_estimate, CR3_deviate = point_estimate / CR3_sd_estimate,
              DB_CI = c(point_estimate - qnorm(1 - alpha / 2) * DB_sd_estimate, 
                        point_estimate + qnorm(1 - alpha / 2) * DB_sd_estimate),
              CR0_CI = c(point_estimate - qnorm(1 - alpha / 2) * CR0_sd_estimate, 
                         point_estimate + qnorm(1 - alpha / 2) * CR0_sd_estimate),
              CR3_CI = c(point_estimate - qnorm(1 - alpha / 2) * CR3_sd_estimate, 
                         point_estimate + qnorm(1 - alpha / 2) * CR3_sd_estimate)))
}



# Function for Horvitz-Thompson estimation of estimands.
# data contains outcome, treatment, cluster, period, covariates
# Ij is a vector of length J that indicates how many clusters are treated at
# each time point
stepped_wedge_HT <- function(data, Ij = NA, weights = NA, 
                             estimand = "ind", alpha = 0.05,
                             g1 = NA, g0 = NA, reg_adjust = TRUE, 
                             mixed_model = FALSE,
                             full_data = FALSE) {
  # cannot have reg_adjust and mixed_model
  if (mixed_model & reg_adjust) {
    stop("Pick one of ordinary regression adjustment and mixed model adjustment.")
  }

  # number of clusters
  clusters <- as.integer(unique(data$cluster))
  I <- length(unique(data$cluster))
  # number of periods
  # assume periods numbered 0 thru J+1
  J <- max(data$period) - 1
  periods <- 1:J

  # initialize Ij
  if(is.na(Ij)) {
    result <- data %>%
      filter(treatment == 1 & period != 0 & period != (J+1)) %>%
      group_by(period) %>%
      summarise(unique_clusters = n_distinct(cluster))

    Ij <- result$unique_clusters

  }

  # propensity scores
  p_scores <- Ij / I
  joint_p_scores_11 <- array(NA, dim = c(I, J, I, J), dimnames =
                               list(paste0("cluster_i = ", clusters),
                                    paste0("period_j = ", periods),
                                    paste0("cluster_i' = ", clusters),
                                    paste0("period_j' = ", periods)))
  joint_p_scores_10 <- joint_p_scores_11
  joint_p_scores_01 <- joint_p_scores_11
  joint_p_scores_00 <- joint_p_scores_11
  # number of subjects
  N <- nrow(data)

  # initialize cluster sizes and  avg outcomes and treatments.
  N_ij <- matrix(NA, nrow = I, ncol = J)
  N_j <- matrix(NA, nrow = I, ncol = J)
  bar_Y_ij <- matrix(NA, nrow = I, ncol = J)
  g0_ij <- matrix(NA, nrow = I, ncol = J)
  g1_ij <- matrix(NA, nrow = I, ncol = J)
  Z_ij <- matrix(NA, nrow = I, ncol = J)

  # If NA, no regression adjustment
  if (is.na(g1)) {
    g1 <- function(X) {
      return(0)
    }
  }
  if (is.na(g0)) {
    g0 <- function(X) {
      return(0)
    }
  }
  # if regression adjust, call helper function
  if (reg_adjust) {
    if (full_data) {
      reg_adjust_object <- reg_adjust_HT_full(data)
      g1 <- reg_adjust_object$g1
      g0 <- reg_adjust_object$g0
    } else {
      reg_adjust_object <- reg_adjust_HT(data)
      g1 <- reg_adjust_object$g1
      g0 <- reg_adjust_object$g0
    }
  }
  
  if (mixed_model) {
    lmm_or <- lmm_HT(data)
    g1 <- lmm_or
    g0 <- lmm_or
  }

  for (i in 1:I) {
    for (j in 1:J) {
      rel <- data %>% filter(cluster == clusters[i] & period == periods[j])
      if (reg_adjust & full_data) {
        rel_covariates <- rel %>% dplyr::select(-c(treatment, outcome)) %>% 
          mutate(period = as.factor(period))
      } else if (reg_adjust & !full_data) {
        rel_covariates <- rel %>% dplyr::select(-c(period, treatment, outcome))
      }
      
      if (mixed_model) {
        rel_covariates <- rel %>% dplyr::select(-c(outcome)) %>% 
          mutate(period = as.factor(period)) %>% 
          mutate(cluster = as.factor(cluster))
      }
      
      N_ij[i,j] <- nrow(rel)
      bar_Y_ij[i, j] <- mean(rel$outcome)
      g0_ij[i, j] <- mean(g0(rel_covariates))
      g1_ij[i, j] <- mean(g1(rel_covariates))
      Z_ij[i, j] <- any(rel$treatment == 1)

    }
  }

  for (i in 1:I) {
    for (j in 1:J) {
      N_j[i,j] <- sum(N_ij[,j])

    }
  }

  if (estimand == "ind") {
    weights <- N_ij
  } else if (estimand == "period") {
    weights <- N_ij / N_j
  } else if (estimand == "cell") {
    weights <- matrix(1, nrow = I, ncol = J)
  }

  # initialize joint p_scores
  for (i in 1:I) {
    for (j in 1:J) {
      for (i_prime in 1:I) {
        for (j_prime in 1:J) {
          # 11 case
          if (i == i_prime) {
            if (j <= j_prime) {
              joint_p_scores_11[i, j, i_prime, j_prime] <- p_scores[j]
            } else {
              joint_p_scores_11[i, j, i_prime, j_prime] <- p_scores[j_prime]
            }

          } else {
            if (j <= j_prime) {
              joint_p_scores_11[i, j, i_prime, j_prime] <- Ij[j] / I *
                (Ij[j_prime] - 1) / (I - 1)
            } else {
              joint_p_scores_11[i, j, i_prime, j_prime] <- Ij[j_prime] / I *
                (Ij[j] - 1) / (I - 1)
            }
          }
          # 00 case
          if (i == i_prime) {
            if (j <= j_prime) {
              joint_p_scores_00[i, j, i_prime, j_prime] <- 1 - p_scores[j_prime]
            } else {
              joint_p_scores_00[i, j, i_prime, j_prime] <- 1 - p_scores[j]
            }

          } else {
            if (j <= j_prime) {
              joint_p_scores_00[i, j, i_prime, j_prime] <- (1 - Ij[j_prime] / I) *
                (I - Ij[j] - 1) / (I - 1)
            } else {
              joint_p_scores_00[i, j, i_prime, j_prime] <- (1 - Ij[j] / I) *
                (I - Ij[j_prime] - 1) / (I - 1)
            }
          }
          # 10 case
          if (i == i_prime) {
            if (j <= j_prime) {
              joint_p_scores_10[i, j, i_prime, j_prime] <- 0
            } else {
              joint_p_scores_10[i, j, i_prime, j_prime] <- (Ij[j] - Ij[j_prime]) / I
            }
          } else {
            if (j == j_prime) {
              joint_p_scores_10[i, j, i_prime, j_prime] <- (1 - Ij[j] / I) *
                (Ij[j]) / (I - 1)
            } else if (j > j_prime) {
              joint_p_scores_10[i, j, i_prime, j_prime] <- (Ij[j] - Ij[j_prime]) / I *
                (Ij[j] - 1) / (I - 1) + (1 - Ij[j] / I) * (Ij[j]) / (I - 1)
            } else if (j < j_prime) {
              joint_p_scores_10[i, j, i_prime, j_prime] <- (Ij[j] / I) *
                (I - Ij[j_prime]) / (I - 1)
            }
          }

        }
      }

    }
  }

  # compute the HT estimate
  ht_hat <- 0
  for (i in 1:I) {
    for (j in 1:J) {
      if (Z_ij[i,j] == 1) {
        ht_hat =  ht_hat + weights[i,j] * (bar_Y_ij[i,j] - g1_ij[i,j]) / p_scores[j]
      } else {
        ht_hat =  ht_hat - weights[i,j] * (bar_Y_ij[i,j] - g0_ij[i,j]) / (1-p_scores[j])
      }
      ht_hat =  ht_hat + weights[i,j] * g1_ij[i,j] - weights[i,j] * g0_ij[i,j]
    }
  }
  ht_hat <- ht_hat / N

  # variance estimate
  var_treat <- 0
  var_control <- 0
  covar <- 0
  # conservative variance estimate
  cons_var_treat <- 0
  cons_var_control <- 0
  cons_covar <- 0

  for (i in 1:I) {
    for (j in 1:J) {
      for (i_prime in 1:I) {
        for (j_prime in 1:J) {
          # var_treat
          if (joint_p_scores_11[i,j,i_prime,j_prime] > 0) {
            term = Z_ij[i,j] * Z_ij[i_prime,j_prime] *
              (joint_p_scores_11[i,j,i_prime,j_prime] - p_scores[j] * p_scores[j_prime]) /
              (joint_p_scores_11[i,j,i_prime,j_prime] * p_scores[j] * p_scores[j_prime]) *
              weights[i,j] * weights[i_prime, j_prime] * (bar_Y_ij[i, j] - g1_ij[i, j]) *
              (bar_Y_ij[i_prime, j_prime] - g1_ij[i_prime, j_prime])
            var_treat = var_treat + term
            cons_var_treat = cons_var_treat + term
          } else {
            # Add conservative correction by Young's inequality
              cons_var_treat = cons_var_treat + Z_ij[i,j] * weights[i,j]^2 * (bar_Y_ij[i, j] - g1_ij[i, j])^2 /
                (2 * p_scores[j]) + Z_ij[i_prime,j_prime] * weights[i_prime,j_prime]^2 *
                (bar_Y_ij[i_prime, j_prime] - g1_ij[i_prime, j_prime])^2 / (2 * p_scores[j_prime])

            
          }

          # var_control
          if (joint_p_scores_00[i,j,i_prime,j_prime] > 0) {
            term = (1 - Z_ij[i,j]) * (1 - Z_ij[i_prime,j_prime]) *
              (joint_p_scores_00[i,j,i_prime,j_prime] - (1 - p_scores[j]) * (1 - p_scores[j_prime])) /
              (joint_p_scores_00[i,j,i_prime,j_prime] * (1 - p_scores[j]) * (1 - p_scores[j_prime])) *
              weights[i,j] * weights[i_prime, j_prime] * (bar_Y_ij[i, j] - g0_ij[i, j]) *
              (bar_Y_ij[i_prime, j_prime] - g0_ij[i_prime, j_prime])
            var_control = var_control + term
            cons_var_control = cons_var_control + term
          } else {
            # fill in conservative object by Young's inequality
              cons_var_control = cons_var_control + (1 - Z_ij[i,j]) * weights[i,j]^2 * (bar_Y_ij[i, j] - g0_ij[i, j])^2 /
                (2 * (1 - p_scores[j])) + (1 - Z_ij[i_prime,j_prime]) * weights[i_prime,j_prime]^2 *
                (bar_Y_ij[i_prime, j_prime] - g0_ij[i_prime, j_prime])^2 / (2 * (1 - p_scores[j_prime]))

            

          }

          # covar
          if (joint_p_scores_10[i,j,i_prime,j_prime] > 0) {
            term = (Z_ij[i,j]) * (1 - Z_ij[i_prime,j_prime]) *
              (joint_p_scores_10[i,j,i_prime,j_prime] - (p_scores[j]) * (1 - p_scores[j_prime])) /
              (joint_p_scores_10[i,j,i_prime,j_prime] * (p_scores[j]) * (1 - p_scores[j_prime])) *
              weights[i,j] * weights[i_prime, j_prime] * (bar_Y_ij[i, j] - g1_ij[i, j]) *
              (bar_Y_ij[i_prime, j_prime] - g0_ij[i_prime, j_prime])
            covar = covar + term
            cons_covar = cons_covar + term
          } else {
            # fill in conservative object by Young's inequality
              cons_covar = cons_covar - (Z_ij[i,j]) * weights[i,j]^2 * (bar_Y_ij[i, j] - g1_ij[i, j])^2 /
                (2 * p_scores[j]) - (1 - Z_ij[i_prime,j_prime]) * weights[i_prime,j_prime]^2 *
                (bar_Y_ij[i_prime, j_prime] - g0_ij[i_prime, j_prime])^2 / (2 * (1 - p_scores[j_prime]))
          
          }

        }
      }

    }
  }

  # normalize
  var_treat <- var_treat / N^2
  cons_var_treat <- cons_var_treat / N^2
  var_control <- var_control / N^2
  cons_var_control <- cons_var_control / N^2
  covar <- covar / N^2
  cons_covar <- cons_covar / N^2

  # put everything together
  var_hat <- var_treat + var_control - 2 * covar
  cons_var_hat <- cons_var_treat + cons_var_control - 2 * cons_covar
  # return the deviate
  deviate <- ht_hat / sqrt(var_hat)
  cons_deviate <- ht_hat / sqrt(cons_var_hat)

  return(list(est = ht_hat,  sd_est = sqrt(var_hat), deviate = deviate,
              cons_sd_est = sqrt(cons_var_hat), cons_deviate = cons_deviate,
              CI = c(ht_hat - qnorm(1 - alpha / 2) * sqrt(var_hat), ht_hat + 
                       qnorm(1 - alpha / 2) * sqrt(var_hat)),
              cons_CI = c(ht_hat - qnorm(1 - alpha / 2) * sqrt(cons_var_hat), ht_hat + 
                       qnorm(1 - alpha / 2) * sqrt(cons_var_hat))))

}

# IV model estimator
ivmodel_estimator <- function(data, adjust = TRUE) {
  data <- data %>% mutate(period = as.factor(period)) 
  if (adjust) {
    X <- model.matrix(~ . - 1, data = data %>% dplyr::select(-outcome, -D, -treatment, -cluster))
    model <- ivmodel(Y = data$outcome, D = data$D, Z = data$treatment,
                     X = X,
                     clusterID = data$cluster, intercept = FALSE)
  } else {
    model <- ivmodel(Y = data$outcome, D = data$D, Z = data$treatment,
                     clusterID = data$cluster)
  }
  return(model)

}

aer_ivreg_estimator <- function(data, adjust = TRUE) {

}

