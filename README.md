# SAR-R-Group-Decomposition

jupyternotebook named R-group-decomposition has been made to analyse the tox assays for given set of analogues and extract insights for SAR analysis. 
This notebook performs cluster analysis using RDKit functions to identify a common scaffold. In next step R-group decomposition is performed on the clusters of the analogs to extract the substitutents (R1,R2, R3 and R4) on the common scaffold.

2D- SAR is then performed to correlate the combination of R groups and identify the structures that impact the toxicity and potency.
The input files are attached along.

Further this analysis includes calculation of distribution of different property and safety assay data.

distribution-safety-data.py is used to see the distribution of LD20 safety assays.

correlation-matrix-LD20.py is used to calculate and plot the correlations between different safety assays.

property-distribution.py is used to calculate the distribution of properties. The input data used is SMILES file in the folder.

distribution-shape-smiles.py used to calculate shape descriptors for the SMILES of given candidates.

pca-shape-tox.py is used to calculate the shape descriptors for the compounds and perform PCA analysis on this data to get the heatmap.


