# sims-ISCB2023
**Simulations of two-level data sets that are analyzed by G-computation and IPTW estimators, combined with random effects models.**

*Steps to generate the simulated data sets are the following:*

# Step 1: 
Import required libraries (file: libraries.R)

# Step 2: 
Import function expit() (file: expit.R)

# Step 3: 
Import function to simulate the independent covariates (file: cancerDataDesign.R)

# Step 4: 
Import functions to simulate the treatment and the potential outcomes (files SimulTreat.R \& SimulOut.R)

# Step 5: 
Generate a .csv file with the Data Generating Mechnaism (DGM) indices (file dgm_index.R generates the respective dgm_index.csv document, which provides a detailed description for each of the 288 assumed DGMs).

# Step 6: 
Combine all functions into one and run via an HPC (file simdata_hpc.R) 

# Step 7: 
Generate the simulated data sets via the HPC (file test_script.txt) 


*Steps to analyze the simulated data sets are the following:*

# Step 1: 
Import function to exclude clusters with only one treatment level being assigned (file: mydata_excl.R)

# Step 2: 
Import functions for marginal and clustered IPTW, and DR (AIPTW) (files mrgest.R; clest.R; drest.R, respectively)

# Step 3: 
Analyze the simulated data sets (file analysis_hpc.R)
















