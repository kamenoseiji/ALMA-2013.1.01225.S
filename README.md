# ALMA-2013.1.01225.S
This repository stores reduction scripts for ALMA 2013.1.01225.S and scripts to produce figures used in Kameno+2020, ApJ, 895, 73 https://ui.adsabs.harvard.edu/link_gateway/2020ApJ...895...73K/doi:10.3847/1538-4357/ab8bd6

# CASA Reduction Scripts (prepared for CASA 5.6.0)
B6/Xa76868_X2423.py             : Band-6 reduction for ExecBlock uid://A002/Xa76868/X2423
B7CO/Xa7a216.py                 : Band-7 reduction, spectral setup#1 for ExecBlocks uid://A002/Xa7a216/X1c4c and uid://A002/Xa7c533/X1dc7
B7HCN/NGC1052.py                : Band-7 reduction, spectral setup#2 for ExecBlocks uid://A002/Xa830fc/X16e2 and uid://A002/Xa830fc/X1a2f

# Text-based spectral data, produced by CASA image viewer
SpecData/                       : text-based spectral data (frequency and velocity basis)

# R scripts 
Scripts/AngularMomentum.R       : for figure 11
       /ColumnDensity.R         : for figure 10
       /CoveringFactor.R        : to calculate covering factors
       /FileBaseline.R          : set baselines in figures 4, 5, and 6
       /GaussFit.R              : Gaussian fit for figure 8
       /LineDecompose.R         : Spectral decomposition for mixed profiles in figures 4, 5, and 6
       /LineList.R              : Line species from the splatalogue
       /MCspec.R                : identificaiton of molecular clouds in figure 1
       /plotMS2Spec.R           : for figures 4, 5, and 6
       /plotSpec.R              : prototype of plotMS2Spec.R
       /PVregression.R          : for figure 7
       /VelocityGradient.R      : determine velocity gradient from figure 7

# figures
figures/
