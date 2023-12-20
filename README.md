# ALMA-2013.1.01225.S
This repository stores reduction scripts for ALMA 2013.1.01225.S and scripts to produce figures used in Kameno+2020, ApJ, 895, 73 https://iopscience.iop.org/article/10.3847/1538-4357/ab8bd6

# CASA Reduction Scripts (prepared for CASA 5.6.0) and image files
- B6/Xa76868_X2423.py             : Band-6 reduction for ExecBlock uid://A002/Xa76868/X2423
- B6/NGC_1052_ContSubforChanMap*  : CO (J=2-1) continuum-subtracted images 
- B7CO/Xa7a216.py                 : Band-7 reduction, spectral setup#1 for ExecBlocks uid://A002/Xa7a216/X1c4c and uid://A002/Xa7c533/X1dc7
- B7CO/NGC_1052_comb_CO_Tap.*     : CO (J=3-2) continuum-subtracted and tapered images 
- B7HCN/NGC1052.py                : Band-7 reduction, spectral setup#2 for ExecBlocks uid://A002/Xa830fc/X16e2 and uid://A002/Xa830fc/X1a2f

Other large image data are stored in https://kamenoseiji.sakura.ne.jp/ALMA-2013.1.01225.S/

# Text-based spectral data, produced by CASA image viewer
- SpecData/                       : text-based spectral data (frequency and velocity basis)

# R scripts 
- Scripts/AngularMomentum.R       : for figure 11
- Scripts/ColumnDensity.R         : for figure 10
- Scripts/CoveringFactor.R        : to calculate covering factors
- Scripts/FileBaseline.R          : set baselines in figures 4, 5, and 6
- Scripts/GaussFit.R              : Gaussian fit for figure 8
- Scripts/LineDecompose.R         : Spectral decomposition for mixed profiles in figures 4, 5, and 6
- Scripts/LineList.R              : Line species from the splatalogue
- Scripts/MCspec.R                : identificaiton of molecular clouds in figure 1
- Scripts/plotMS2Spec.R           : for figures 4, 5, and 6
- Scripts/plotSpec.R              : prototype of plotMS2Spec.R
- Scripts/PVregression.R          : for figure 7
- Scripts/VelocityGradient.R      : determine velocity gradient from figure 7

# figures
- figures/
