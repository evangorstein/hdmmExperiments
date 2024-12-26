This repository contains all scripts and notebooks used to produce tables, figures, and other analyses in "HighDimMixedModels.jl: Robust High Dimensional Mixed Models across Omics Data".


# Simulated data

## Simulating data

GWAS and OTU design matrices for the simulation study were generated in R by `R/simulating_gwas_matrices/gwas_sim.R` and `R/simulating_otu_matrices/spring.R`, respectively. Simulation of the response for the gene expression, gwas, and microbiome experiments (and of the design matrices for the former) are done by scripts `gen_data_gene.jl`, `gen_data_gwas.jl` and `gen_data_OTU.jl`, respectively, all found within `scripts`.  

## Fitting models to simualted data
The fitting of the high dimensional models was done by running a script on CHTC. The script for the GWAS experiment is provided at `scripts/experiment_gwas.jl` as an example. Essentially, this script goes through each simulated dataset from each setting and fits the model to the dataset over a grid of lambdas. The fitted models for all simulated data sets were then processed with, e.g. `scripts/process_results_gwas.jl`, which extracts the top two fits by BIC from each dataset and stores them in the CSV files `sim_results/gwas_new/best_results.csv` for downstream visualization

## Visualizing performance of fitted models

The visualizations of simulation results found in the paper are all generated in `scripts/plot_results.Rmd`

# Real data analyses

The real gene expression data set analyzed in the paper is saved in `data/real/gene_expressions/`. The real OTU data set is saved in `data/real/OTU`, with pre-processing of the data done by `R/real_otu_data/real_OTU.r`.  Finally, the real GWAS data set is the `mice` dataset from the R package `BGLR`. All analyses and model fitting of the real data sets are performed by the scripts `scripts/real_*.jl`, for * = `gene`, `gwas` or `OTU`, respectively.