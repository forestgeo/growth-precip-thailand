# growth-precip-thailand
Sensitivity of tree growth to drought using dendrometer band records from Huai Kha Khaeng
2009 - 2022
[![DOI](https://zenodo.org/badge/816825930.svg)](https://doi.org/10.5281/zenodo.15777974)

# Folder structure:

- data - contains climate data csvs for HKK
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
- sensitivity_dataset.csv - main dataset 
- sensitivity_metadata.csv - metadata describing columns of the main dataset
- summaries_dataset.csv - annual summaries of increment
- summaries_metadata.csv - metadata describing columns of summary dataset

# scripts
- 000_twi.R - script to calculate TWI and TPI across the plot
- 01_prepdata.R - script to prep data for analysis starting from dendrometer band measurements
- 02_species.R - [START HERE] script to run Bayesian level models with only species predictors (intercept only models and deciduous/TWI models)
- 03_orderedcii.R - script to run Bayesian causal models with CII, DBH, TWI and species random effects
- 04_orderedcii_tpi.R - script to run models Bayesian causal models with CII, DBH, TPI and species random effects
- 05_figures.R - script to make main figures from data and model outputs.
- 06_causal.R - code for alternate causal inquiry models
- 07_otherfigs.R - script for some SI figures

# doc/
- manuscript_hkk_drought_sensitivity.Rmd - Main manuscript document in Rmd format
- manuscript_hkk_drought_sensitivity.docx - Main manuscript in docx format (knit version)
- supplementary_information.Rmd - combined Supplementary Appendix
- supplementary_information.pdf - pdf version Supplmentary Appendix
- growth-precip-thailand.bib - BibTex file with references
- ecology-letters.csl - Citation Style Language file
- word-styles-reference - Style Reference to knit manuscript
- cover letter - in Rmd and pdf

# doc/display/
- Contains all figures in png and tables in csv




