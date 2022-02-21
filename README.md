# Nearest Neighbour Distance Matching Leave-One-Out Cross-Validation for map validation

This repository contains the analysis code for the paper "Nearest Neighbour Distance Matching Leave-One-Out Cross-Validation for map validation" by Carles Milà, Jorge Mateu, Edzer Pebesma, and Hanna Meyer. The manuscript has been submitted to the journal *Methods in Ecology and Evolution*.

The contents of the repository are as follows:

* [sim1_analysis.R](code/sim1_analysis.R): Contains the code to run simulation 1 (random fields), which uses functions defined in [sim1_functions.R](code/sim1_functions.R) and [sim_utils.R](code/sim_utils.R).
* [sim2_analysis.R](code/sim2_analysis.R): Contains the code to run simulation 2 (virtual species), which uses a landscape dataset and simulates a species outcome created with the code included in [sim2_landscape.R](code/sim2_landscape.R), and uses functions defined in [sim2_functions.R](code/sim2_functions.R) and [sim_utils.R](code/sim_utils.R).
* [Data](data/): It contains the landscape data and simulated outcome used in simulation 2.
* [Results](results/): It contains csv files with the output of simulations 1 and 2.
* [methods_example.Rmd](reports/methods_example.Rmd): It contains the analysis and visualization code to generate the figures included in section 1 and 2.1 of the manuscript that serve as an illustration of the proposed methods. A pdf version of the file can be found in [methods_example.pdf](reports/methods_example.pdf).
* [sims_workflow.Rmd](reports/sims_workflow.Rmd): It contains the analysis and visualization code to generate simulation workflow figures included in sections 2.2 and 2.3 of the manuscript as well as supporting figure S1. A pdf version of the file can be found in [sims_workflow.pdf](reports/sims_workflow.pdf).
* [sims1_fields.Rmd](reports/sims1_fields.Rmd): It contains the code to generate the results figures of simulation 1 included in section 3.1 of the manuscript and the appendix. A pdf version of the file can be found in [sims1_fields.pdf](reports/sims1_fields.pdf).
* [sims2_species.Rmd](reports/sims2_species.Rmd): It contains the code to generate the results figures of simulation 2 included in section 3.2 of the manuscript and the appendix. A pdf version of the file can be found in [sims2_species.pdf](reports/sims2_species.pdf).
