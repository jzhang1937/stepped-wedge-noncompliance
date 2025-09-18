#!/bin/bash

# Check if an argument is passed
if [ -z "$1" ]; then
  echo "Usage: $0 <full_or_relative_path_to_R_script> [args_for_R_script...]"
  exit 1
fi

# Get the full or relative path to the R script from the first argument
R_SCRIPT="$1"

# Shift the argument list to exclude the first argument (R script filename)
shift

# Run the R script using Rscript command with any additional arguments
Rscript "$R_SCRIPT" "$@"

# Check if Rscript executed successfully
if [ $? -eq 0 ]; then
  echo "R script $R_SCRIPT executed successfully."
else
  echo "Error: Failed to execute $R_SCRIPT."
fi

# example usage
# ./simulations/run_r_script_args.sh ./simulations/run_simulation.R arg1 arg2 arg3
