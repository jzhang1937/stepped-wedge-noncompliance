# code to get the reproducible RNG state checkpoints
# when splitting the simulation into multiple scripts
cluster_period <- c("I11J10", "I12J5", "I30J5", "I60J5", "I90J5")
informative <- c(TRUE, FALSE)
folds <- 10
n_iter <- 1000

# Stepped-wedge simulation
library(dplyr)
source("simulation_helpers.R")
source("estimator_helpers.R")
source("estimators.R") 

# simulation parameters
compliance_intercept = -0.5 

for (cp in cluster_period) {
  # Cluster-period setup
  IJ <- cp
  # find the positions of the numbers after I and J
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
  for (inf in informative) {
    
    inf_cluster_size = inf
    if (inf_cluster_size) {
      compliance_factor = 6
    } else {
      compliance_factor = 2.5
    }
    # name of cluster-period and informative or not
    cp_inf <- paste(cp, inf, sep = "_")
    # reproducibility
    set.seed(123)
    for (i in 1:folds){
      checkpoint_name <- paste0("simulations/reproducible_seeds/", cp_inf,"_", i, ".RData")
      save(.Random.seed, file = checkpoint_name)
      for (j in 1:(n_iter/folds)) {
        sim_data <- generate_data(I = I, J = J, Ij = Ij, 
                                  compliance_intercept = compliance_intercept, 
                                  compliance_factor = compliance_factor,
                                  inf_cluster_size = inf_cluster_size)
      }
    }
  }
}




