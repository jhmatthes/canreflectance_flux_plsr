canreflectance_flux_plsr
========================
This repository contains code used within: Matthes, J. H., S. H. Knox, C. Sturtevant, O. Sonnentag, J. Verfaillie, D. Baldocchi, Predicting landscape-scale CO2 flux at a pasture and rice paddy with long-term hyperspectral canopy reflectance measurements. Biogeosciences Discussions, submitted.

This code takes long-term repeated hyperspectral canopy reflectance measurements, eddy covariance flux measurements, and conducts PLSR modeling to examine the ability of the canopy reflectance measurements to predict landscape-scale CO2 flux.

hyperspectral_flux_plsr.R: This is the main set of code used to run the analysis.

pslr_bootstrap.R: Conducts bootstrapped PLSR fitting using randomized sets of calibration/evaluation data.

plsr_validate.R: Used to evaluate the mean bootstrapped PLSR model on an independent validation dataset.

plot_plsr_validate.R: Used to create a plot of the plsr validation results.

VIP.R: Function to calculate the variable importance of projection metric.


