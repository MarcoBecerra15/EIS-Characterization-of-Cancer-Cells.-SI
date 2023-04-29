# EIS-Characterization-of-Cancer-Cells.-SI
This is the supplementary information for Electrochemical Impedance Spectroscopy Characterization and Identification of Cancer Cells using Machine Learning Methods.

__________________________________________________________Code Files Definitions_______________________________________________________________________________________
MyPlotCells = Plot Zmod/Zphz vs Frequency 
My_CodeOutsCompatible = 1. Arrange GAMRY Zphz/Zmod into .mat table files // 2. Arrange filtered data into .csv files for PCA and Random Forests 
Cell_Outlier_ID_Code = Filter out outliers and remove baseline if desired
LoopPCA = Perform Principal Component Analysis (PCA)
RF_Manual = Generate Random Forests (RF) Models
LOOCV Code = Perform Leave-One-Out Cross Validation on the RF models 
Paper Plotter = Compile PCA Plots 
Table Paper = Arrange Model Accuracies into Table 
