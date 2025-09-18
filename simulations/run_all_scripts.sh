#!/bin/bash

# Check if the R script is provided as the first argument
if [ -z "$1" ]; then
  echo "Usage: $0 <full_or_relative_path_to_R_script>"
  exit 1
fi

# Get the R script from the first argument
R_SCRIPT="$1"

# Define the possible choices for each argument
ARG1_CHOICES=("I11J10" "I12J5" "I30J5" "I60J5" "I90J5")
ARG2_CHOICES=("unadj_HT" "reg_adj_HT" "lmm_adj_HT" "reg_adj_full_HT" "unadj_ANCOVA" "ANCOVAI" "ANCOVAIII")
ARG3_CHOICES=("TRUE" "FALSE")

# Loop through all combinations of arg1, arg2, and arg3
for arg1 in "${ARG1_CHOICES[@]}"; do
  for arg2 in "${ARG2_CHOICES[@]}"; do
    for arg3 in "${ARG3_CHOICES[@]}"; do
      # Run each combination in the background
      qsub -l m_mem_free=8G simulations/run_r_script_args.sh "$R_SCRIPT" "$arg1" "$arg2" "$arg3"
    done
  done
done

# Wait for all background jobs to finish
wait

echo "All jobs submitted."
