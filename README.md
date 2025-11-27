# growth-precip-thailand
Sensitivity of tree growth to drought using dendrometer band records from Huai Kha Khaeng
2009 - 2022
[![DOI](https://zenodo.org/badge/816825930.svg)](https://doi.org/10.5281/zenodo.15777974)

**To reproduce analyses and produce figures, run scripts from the *scripts/analysis* folder, sequentially.**

Scripts under the *processing* are provided for transparency in the data processing steps that output the shared datasets, but cannot be run without necessary dependencies.

# Folder structure:

- data - contains climate data and tree data csvs for HKK
- doc - contains figures and documents
- results - contains climate data plots
- scripts - R scripts to analyse the results

# data/climate/
- CHIRPS_Daily_HKK.csv - contains file with daily precipitation values for HKK from CHIRPS 
- CHIRPS_HKK.csv - monthly precipitation summaries for HKK from CHIRPS
- ENSO_index_meiv2.txt - ENSO index for HKK
- ERA5Land_Daily_HKK.csv - daily temperature and relative humidity values for HKK from ERA5Land
- ERA5Land_HKK.csv - monthly summaries of ERA5Land variables
- SPEI_HKK_from_GEE.csv - Standardised Precipitation Evapotranpiration Index values for 1-, 3-, 6- and 12-month windows

# data/dendro/
- sensitivity_dataset.csv - main dataset containing tree by tree information
- sensitivity_metadata.csv - metadata describing columns of the main dataset
- summaries_dataset.csv - annual summaries of increment
- summaries_metadata.csv - metadata describing columns of summary dataset
- sp_iucn.csv - dataset with IUCN categories of species used in the study. Accessed in September 2025

# scripts/
- session_info.txt - text file describing software and package versions used and hardware details. Output from sessionInfo() in R

## scripts/analysis
Contains all scripts to fully reproduce analyses, results and figures
- 01_species.R - script to run Bayesian level models with only species predictors (intercept only models and deciduous/TWI models)
- 02_orderedcii.R - script to run Bayesian causal models with CII, DBH, TWI and species random effects
- 03_orderedcii_tpi.R - script to run models Bayesian causal models with CII, DBH, TPI and species random effects
- 04_causal.R - code for alternate causal inquiry models
- 05_figures.R - script to make main figures from data and model outputs.

## scripts/processing
Contains all scripts for processing raw data. Outputs datasets available in data/dendro
- 000_twi.R - script to calculate TWI and TPI across the plot
- 01_prepdata.R - script to prep data for analysis starting from dendrometer band measurements

# doc/
- manuscript_hkk_drought_sensitivity.Rmd - Main manuscript document in Rmd format
- manuscript_hkk_drought_sensitivity.docx - Main manuscript in docx format (knit version)
- manuscript_hkk_drought_sensitivity_v1.docx - Main manuscript v1 (first submission) in docx format (knit version)
- manuscript_hkk_drought_sensitivity_v2_track.docx - Main manuscript v2 in docx format with track changes from v1
- supplementary_information.Rmd - combined Supplementary Appendix
- supplementary_information.pdf - pdf version Supplmentary Appendix
- growth-precip-thailand.bib - BibTex file with references
- ecology-letters.csl - Citation Style Language file
- word-styles-reference - Style Reference to knit manuscript
- cover letter - in Rmd and pdf

## doc/display/
- Contains all figures in png and tables in csv




