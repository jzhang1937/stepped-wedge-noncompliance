# Simulation helpers for stepped wedge simulations
library(dplyr)

# I - number of clusters
# J - number of periods, not including pre and post rollout
# Ij - a vector controlling the number of clusters crossing over at each
# period, including pre and post, so the first element must be zero, and the last
# must be I.
# compliance_intercept - the intercept in the multinomial logistic regression
# compliance_factor - how much to shrink by to get more equal probabilities
# 0 - complier, 1 - always taker, 2 - never taker
# for compliance. smaller means more compliance
generate_data <- function(I, J, Ij, compliance_intercept = 0, 
                          compliance_factor = 1.5,
                          inf_cluster_size = FALSE) {
  
  # cluster size and binomial covariate
  N_ij <- array(NA, dim = c(I, J + 2), dimnames = 
          list(paste0("cluster ", 1:I), 
               paste0("period ", 0:(J+1))))
  X_ij1 <- array(NA, dim = c(I, J + 2), dimnames = 
                  list(paste0("cluster ", 1:I), 
                       paste0("period ", 0:(J+1))))
  for (i in 1:I) {
    for (j in 0:(J+1)) {
      N_ij[paste0("cluster ", i),paste0("period ", j)] <- 
        ceiling(runif(1, min = 10, max = 90) + 2 * (j + 1)^(1.5))
      
      X_ij1[paste0("cluster ", i),paste0("period ", j)] <- 
        rbinom(1, 1, prob = 0.5)
    }
  }
  
  # average period size
  avg_period_size <- sum(N_ij) / (J + 2)
  
  # cluster effect
  c_i <- rnorm(I, mean = 0, sd = 0.2)
  
  # put everything into dataframe form
  # Create vectors for cluster and period covariates
  clusters <- rep(1:I, each = J+2)
  periods <- rep(0:(J+1), times = I)
  
  # Expand each cluster-period by the number of individuals in each
  expanded_clusters <- rep(clusters, times = as.vector(t(N_ij)))
  expanded_periods <- rep(periods, times = as.vector(t(N_ij)))
  
  # Create the dataframe
  df <- data.frame(cluster = expanded_clusters, period = expanded_periods)
  
  cluster_effect_df <- data.frame(cluster = 1:I,
                           c_i = c_i)
  
  cluster_period_covariate_df <- expand.grid(cluster = 1:I,
                              period = 0:(J+1))
  cluster_period_covariate_df$X_ij1 <- as.vector(X_ij1)
  
  df <- df %>%
    left_join(cluster_period_covariate_df, by = c("cluster", "period")) %>% 
    left_join(cluster_effect_df, by = "cluster")
  
  
  # individual covariate and noise
  df <- df %>% group_by(cluster, period) %>%
    mutate(X_ijk2 = cluster/I + runif(n(), min = -1, max = 1)) %>%
    mutate(e_ijk = rnorm(n(), mean = 0, sd = 0.9)) %>%
    mutate(N_ij = n()) %>%
    ungroup()
  df <- df %>% group_by(period) %>% mutate(X_j2 = mean(X_ijk2)) %>% ungroup()
  
  # set up compliance predictor
  if (inf_cluster_size) {
    df <- df %>% mutate(log_1_over_0 = (compliance_intercept + (period + 1) / (J + 2) +
                                          2 * N_ij * I / avg_period_size +
                                          0.7 * X_ij1 +
                                          0.5 * (X_ijk2 - X_j2)^3 + c_i + e_ijk) / compliance_factor) %>% 
      mutate(log_2_over_0 = (compliance_intercept - (period + 1) / (J + 2) - 
                               2 * N_ij * I / avg_period_size -
                               0.4 * X_ij1 +
                               (X_ijk2 - X_j2)^2 - c_i) / compliance_factor) %>% 
      mutate(prob_0 = 1 / (1 + exp(log_2_over_0) + exp(log_1_over_0))) %>% 
      mutate(prob_1 = exp(log_1_over_0) / (1 + exp(log_2_over_0) + exp(log_1_over_0))) %>% 
      mutate(prob_2 = exp(log_2_over_0) / (1 + exp(log_2_over_0) + exp(log_1_over_0))) 
  } else{
    df <- df %>% mutate(log_1_over_0 = (compliance_intercept + (period + 1) / (J + 2) + 0.7 * X_ij1 +
                                          0.5 * (X_ijk2 - X_j2)^3 + c_i + e_ijk) / compliance_factor) %>% 
      mutate(log_2_over_0 = (compliance_intercept - (period + 1) / (J + 2) - 0.4 * X_ij1 +
                               (X_ijk2 - X_j2)^2 - c_i) / compliance_factor) %>% 
      mutate(prob_0 = 1 / (1 + exp(log_2_over_0) + exp(log_1_over_0))) %>% 
      mutate(prob_1 = exp(log_1_over_0) / (1 + exp(log_2_over_0) + exp(log_1_over_0))) %>% 
      mutate(prob_2 = exp(log_2_over_0) / (1 + exp(log_2_over_0) + exp(log_1_over_0))) 
  }
  
  
  # draw compliance class
  df <- df %>%
    rowwise() %>%
    mutate(compliance = sample(0:2, size = 1, prob = c(prob_0, prob_1, prob_2))) %>%
    ungroup()
  table(df$compliance)
  
  # potential outcomes
  if (!inf_cluster_size) {
    df <- df %>% mutate(Y_0 = (period + 1) / (J + 2) + X_ij1 +
                          (X_ijk2 - X_j2)^2 + c_i + e_ijk) %>%
      mutate(Y_1 = Y_0 + 1.25 * (0.5 * X_ij1 + (X_ijk2 - X_j2)^3))
  } else {
  df <- df %>% mutate(Y_0 = (period + 1) / (J + 2) + X_ij1 +
                        (X_ijk2 - X_j2)^2 + c_i + e_ijk) %>%
    mutate(Y_1 = Y_0 + 0.45 * (N_ij * I / avg_period_size + 0.5 * X_ij1 + (X_ijk2 - X_j2)^3))
  }
  
  # draw treatment assignment
  # random permutation of 1 thru I
  random_permutation <- sample(1:I)
  Z_ij <- array(0, dim = c(I, J + 2), dimnames = 
                  list(paste0("cluster ", 1:I), 
                       paste0("period ", 0:(J+1))))
  for (i in 1:I) {
    for (j in 0:(J+1)) {
      if (random_permutation[i] <= Ij[j + 1]) {
        Z_ij[paste0("cluster ", i),paste0("period ", j)] <- 1
      } 
      
    }
  }
  
  cluster_period_treatment_df <- expand.grid(cluster = 1:I,
                                             period = 0:(J+1))
  cluster_period_treatment_df$treatment <- as.vector(Z_ij)
  
  df <- df %>%
    left_join(cluster_period_treatment_df, by = c("cluster", "period")) 
  
  # observed outcome and received treatment D
  df <- df %>% 
    mutate(D = (compliance == 0) * treatment + (compliance == 1) * 1 + (compliance == 2) * 0) %>%
    mutate(outcome = D * Y_1 + (1 - D) * Y_0) 
  
  
  return(data = df)
  
  
}

