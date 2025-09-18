Reproducing the simulations and real data analysis reported in
\`Randomization inference for stepped-wedge designs with noncompliance
with application to a palliative care pragmatic trial’
================

This repository contains code for the simulations and analyses in the
paper “Randomization inference for stepped-wedge designs with
noncompliance with application to a palliative care pragmatic trial”
(arXiv, 2025).

# Dependencies

Ensure that your system has the required software dependencies.

- R version 4.2.2 or higher

# Get started

First, clone the `stepped-wedge-noncompliance` repository onto your
machine.

    git@github.com:jzhang1937/stepped-wedge-noncompliance.git

# Simulations

## Simulation results and figures

Simulation results and figures are in directories `simulation/results`
and `figures`, respectively.

If you would like to rerun the simulations from scratch, follow the
steps in the next section.

## Run simulation scripts

Navigate to the stepped-wedge-noncompliance directory. All scripts below
must be executed from this directory.

These commands were run the Wharton HPC Cluster. For the commands below,
depending on the limits of your cluster, you may need to set memory and
time limit parameters differently within
`simulations/run_all_scripts.sh` and
`simulations/run_all_iv_model_scripts.sh`.

    # ANCOVA and Horvitz-Thompson simulations:
    qsub simulations/run_all_scripts.sh

    # ivmodel simulations:
    Rscript get_seeds.R
    qsub simulations/run_all_iv_model_scripts.sh

## Create the simulation figures

Before creating the figures, please ensure that the working directory is
`stepped-wedge-noncompliance`. The figures are placed in the `figures`
directory.

    # Figures 3, 5, 6, 7
    Rscript plot_figures.R

## Access raw tables

The tables were generated using the code from the
`process_full_results.R` file, which prints latex to the console. They
were then edited manually for style and rounding for the manuscript.

# Real data analysis

The data is not available due to privacy constraints, but we provide the
scripts used to generate results pertaining to the real data analysis in
the `analysis` folder. The code for generating Figure 1, Table 2, 3, and
the left half of Table 4 is in `analysis/primary_analyses.R`. The code
for generating the right half of Table4 is in
`analysis/readm30_analysis.R`.
