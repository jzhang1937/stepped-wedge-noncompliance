# general code to run a simulation with splitting into multiple fold 
# (only used for ivmodel simulation)
args <- commandArgs(trailingOnly = TRUE)
cluster_period <- as.character(args[1])
method <- as.character(args[2])
informative <- as.logical(args[3])
fold <- as.integer(args[4])

# Stepped-wedge simulation
library(dplyr)
source("simulation_helpers.R")
source("estimator_helpers.R")
source("estimators.R") 

# name
sim_name <- paste(cluster_period, method, informative, fold, sep = "_")

# simulation parameters
# Simulation hyperparameters
compliance_intercept = -0.5 
inf_cluster_size = informative
if (inf_cluster_size) {
  compliance_factor = 6
} else {
  compliance_factor = 2.5
}
alpha = 0.05
n_iter = 1000
tol = 0.0001
folds = 10

# Cluster-period setup
IJ <- cluster_period
# Use gregexpr to find the positions of the numbers after I and J
I_position <- gregexpr("I\\d+", IJ)[[1]]
J_position <- gregexpr("J\\d+", IJ)[[1]]

# Extract the matched substrings
I_match <- regmatches(IJ, I_position)
J_match <- regmatches(IJ, J_position)

# Remove the "I" and "J" and convert to numeric
I_number <- as.numeric(sub("I", "", I_match))
J_number <- as.numeric(sub("J", "", J_match))

I = I_number
J = J_number
Ij = seq(0,I,I/(J+1))


# HT function factory
stepped_wedge_HT_function <- function(reg_adjust = TRUE, mixed_model = FALSE, full_data = FALSE) {
  if (reg_adjust & !full_data) {
    placeholder <- function(data) {
      return(stepped_wedge_HT(data = data, reg_adjust = TRUE, mixed_model = FALSE,
                              full_data = FALSE))
    }
  } else if(reg_adjust & full_data) {
    placeholder <- function(data) {
      return(stepped_wedge_HT(data = data, reg_adjust = TRUE, mixed_model = FALSE,
                              full_data = TRUE))
    }
  } else if(!reg_adjust & mixed_model) {
    placeholder <- function(data) {
      return(stepped_wedge_HT(data = data, reg_adjust = FALSE, mixed_model = TRUE,
                              full_data = FALSE))
    }
  } else if(!reg_adjust & !mixed_model) {
    placeholder <- function(data) {
      return(stepped_wedge_HT(data = data, reg_adjust = FALSE, mixed_model = FALSE,
                              full_data = FALSE))
    }
  }
  return(placeholder)
}
# adjust and no adjust functions
stepped_wedge_HT_no_adjust <- stepped_wedge_HT_function(reg_adjust = FALSE, mixed_model = FALSE)
stepped_wedge_HT_reg_adjust <- stepped_wedge_HT_function(reg_adjust = TRUE, mixed_model = FALSE)
stepped_wedge_HT_lmm_adjust <- stepped_wedge_HT_function(reg_adjust = FALSE, mixed_model = TRUE)
stepped_wedge_HT_reg_full_adjust <- stepped_wedge_HT_function(reg_adjust = TRUE, mixed_model = FALSE,
                                                              full_data = TRUE)

# HT methods
HT_methods <- c("unadj_HT", "reg_adj_HT", "lmm_adj_HT", "reg_adj_full_HT")
# ANCOVA methods
ANCOVA_methods <- c("unadj_ANCOVA", "ANCOVAI", "ANCOVAIII")
# ivmodel methods
ivmodel_methods <- c("adj_ivmodel", "unadj_ivmodel")


# Run sim if an HT method
if (method %in% HT_methods) {
  metrics <- c("true.lambda", "point.est", "zero.deviate", "cons.zero.deviate", 
               "deviate", "cons.deviate")
  results <- array(data = NA, dim = c(length(cluster_period), length(metrics), 
                                      n_iter/folds))
  dimnames(results) <- list(
    Dimension1 = cluster_period,            
    Dimension2 = metrics,
    Dimension3 = seq(1:(n_iter/folds))  
  )
  
  if (method == "unadj_HT") {
    method_function <- stepped_wedge_HT_no_adjust
  } else if (method == "reg_adj_HT") {
    method_function <- stepped_wedge_HT_reg_adjust
  } else if (method == "lmm_adj_HT") {
    method_function <- stepped_wedge_HT_lmm_adjust
  } else if (method == "reg_adj_full_HT") {
    method_function <- stepped_wedge_HT_reg_full_adjust
  }
  load(paste0("simulations/reproducible_seeds/",IJ,"_",inf_cluster_size,"_",fold,".RData")) 
  for (i in 1:(n_iter/folds)){
    print(i)
    sim_data <- generate_data(I = I, J = J, Ij = Ij, 
                              compliance_intercept = compliance_intercept, 
                              compliance_factor = compliance_factor,
                              inf_cluster_size = inf_cluster_size)
    if (inf_cluster_size) {
      min_sim_data <- sim_data  |> 
        dplyr::select(cluster, period, X_ij1, X_ijk2, N_ij, treatment, outcome, D)
    } else {
      min_sim_data <- sim_data  |> 
        dplyr::select(cluster, period, X_ij1, X_ijk2, treatment, outcome, D)
    }
    
    compliers <- sim_data |> filter(compliance == 0 & period != 0 & period != (J+1))
    # true value of effect ratio
    true_lambda <- mean(compliers$Y_1 - compliers$Y_0)
    results[IJ,"true.lambda",i] <- true_lambda
    # test lambda = 0
    ER_zero <- effect_ratio_test(data = min_sim_data, 
                                 test_function = method_function, lambda = 0)
    # deviates under zero
    results[IJ,"zero.deviate",i] <- ER_zero$deviate
    results[IJ,"cons.zero.deviate",i] <- ER_zero$cons_deviate
    
    # test lambda = true lambda
    ER_true <- effect_ratio_test(data = min_sim_data, 
                                 test_function = method_function, lambda = true_lambda)
    
    # deviates under zero
    results[IJ,"deviate",i] <- ER_true$deviate
    results[IJ,"cons.deviate",i] <- ER_true$cons_deviate
    
    # point estimates
    results[IJ,"point.est",i] <- 
      effect_ratio_point_estimate(data = min_sim_data, 
                                  test_function = method_function, 
                                  tol = tol) 
  }
} else if (method %in% ANCOVA_methods) {
  metrics <- c("true.lambda", "point.est", "DB.zero.deviate", "CR0.zero.deviate",
               "CR3.zero.deviate", "DB.deviate", "CR0.deviate", "CR3.deviate")
  results <- array(data = NA, dim = c(length(cluster_period), length(metrics), 
                                      n_iter/folds))
  dimnames(results) <- list(
    Dimension1 = cluster_period,            
    Dimension2 = metrics,
    Dimension3 = seq(1:(n_iter/folds))  
  )
  
  if (method == "unadj_ANCOVA") {
    method_function <- unadjusted_ANCOVA
  } else if (method == "ANCOVAI") {
    method_function <- ANCOVAI
  } else if (method == "ANCOVAIII") {
    method_function <- ANCOVAIII
  }
  # set reproducibility
  load(paste0("simulations/reproducible_seeds/",IJ,"_",inf_cluster_size,"_",fold,".RData"))
  for (i in 1:(n_iter/folds)){
    print(i)
    sim_data <- generate_data(I = I, J = J, Ij = Ij, 
                              compliance_intercept = compliance_intercept, 
                              compliance_factor = compliance_factor,
                              inf_cluster_size = inf_cluster_size)
    if (inf_cluster_size) {
      min_sim_data <- sim_data  |> 
        dplyr::select(cluster, period, X_ij1, X_ijk2, N_ij, treatment, outcome, D)
      compliers <- sim_data |> filter(compliance == 0 & period != 0 & period != (J+1))
    } else {
      min_sim_data <- sim_data  |> 
        dplyr::select(cluster, period, X_ij1, X_ijk2,  treatment, outcome, D)
      compliers <- sim_data |> filter(compliance == 0 & period != 0 & period != (J+1))
    }
    
    # true value of effect ratio
    true_lambda <- mean(compliers$Y_1 - compliers$Y_0)
    results[IJ,"true.lambda",i] <- true_lambda
    # test lambda = 0
    ER_zero <- effect_ratio_test(data = min_sim_data, 
                                 test_function = method_function, lambda = 0)
    # deviates under zero
    results[IJ,"DB.zero.deviate",i] <- ER_zero$DB_deviate
    results[IJ,"CR0.zero.deviate",i] <- ER_zero$CR0_deviate
    results[IJ,"CR3.zero.deviate",i] <- ER_zero$CR3_deviate
    
    # test lambda = true lambda
    ER_true <- effect_ratio_test(data = min_sim_data, 
                                 test_function = method_function, lambda = true_lambda)
    
    # deviates under zero
    results[IJ,"DB.deviate",i] <- ER_true$DB_deviate
    results[IJ,"CR0.deviate",i] <- ER_true$CR0_deviate
    results[IJ,"CR3.deviate",i] <- ER_true$CR3_deviate
    
    # point estimates
    results[IJ,"point.est",i] <- 
      effect_ratio_point_estimate(data = min_sim_data, 
                                  test_function = method_function, 
                                  tol = tol) 
  }
} else if (method %in% ivmodel_methods) {
  metrics <- c("true.lambda", "lower.ci.CLR", "upper.ci.CLR", "lower.ci.AR",
               "upper.ci.AR", "lower.ci.Fuller", "upper.ci.Fuller",
               "point.est.Fuller", "lower.ci.LIML", "upper.ci.LIML",
               "point.est.LIML")
  results <- array(data = NA, dim = c(length(cluster_period), length(metrics), 
                                      n_iter/folds))
  dimnames(results) <- list(
    Dimension1 = cluster_period,            
    Dimension2 = metrics,
    Dimension3 = seq(1:(n_iter/folds))  
  )
  
  # set reproducibility
  load(paste0("simulations/reproducible_seeds/",IJ,"_",inf_cluster_size,"_",fold,".RData"))
       for (i in 1:(n_iter/folds)){
         print(i)
         sim_data <- generate_data(I = I, J = J, Ij = Ij, 
                                   compliance_intercept = compliance_intercept, 
                                   compliance_factor = compliance_factor,
                                   inf_cluster_size = inf_cluster_size)
         if (inf_cluster_size) {
           min_sim_data <- sim_data  |> 
             dplyr::select(cluster, period, X_ij1, X_ijk2, N_ij, treatment, outcome, D)
           compliers <- sim_data |> filter(compliance == 0 & period != 0 & period != (J+1))
         } else {
           min_sim_data <- sim_data  |> 
             dplyr::select(cluster, period, X_ij1, X_ijk2,  treatment, outcome, D)
           compliers <- sim_data |> filter(compliance == 0 & period != 0 & period != (J+1))
         }
         
         # true value of effect ratio
         true_lambda <- mean(compliers$Y_1 - compliers$Y_0)
         results[IJ,"true.lambda",i] <- true_lambda
         
         if (method == "adj_ivmodel") {
           fit_model <- ivmodel_estimator(data = min_sim_data, adjust = TRUE)
         } else if (method == "unadj_ivmodel") {
           fit_model <- ivmodel_estimator(data = min_sim_data, adjust = FALSE)
         }
         
         results[IJ,"lower.ci.CLR",i] <- fit_model$CLR$ci[1]
         results[IJ,"upper.ci.CLR",i] <- fit_model$CLR$ci[2]
         results[IJ,"lower.ci.AR",i] <- fit_model$AR$ci[1]
         results[IJ,"upper.ci.AR",i] <- fit_model$AR$ci[2]
         results[IJ,"lower.ci.Fuller",i] <- fit_model$Fuller$ci[1]
         results[IJ,"upper.ci.Fuller",i] <- fit_model$Fuller$ci[2]
         results[IJ,"point.est.Fuller",i] <- fit_model$Fuller$point.est
         results[IJ,"lower.ci.LIML",i] <- fit_model$LIML$ci[1]
         results[IJ,"upper.ci.LIML",i] <- fit_model$LIML$ci[2]
         results[IJ,"point.est.LIML",i] <- fit_model$LIML$point.est
         
       }
}


results_dir <- paste0("simulations/results/", sim_name)
if (!dir.exists(results_dir)) dir.create(results_dir)
saveRDS(results, file = paste0(results_dir, "/", sim_name, "_results.rds"))



