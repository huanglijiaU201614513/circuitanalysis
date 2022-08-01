# Circuit analysis GitHub repository
This README file is a description script for the GitHub repository of "What makes a functional gene regulatory network? A circuit motif analysis" written by Lijia Huang,  Benjamin Clauss and Mingyang Lu.
## HH_FF_scores.numbers
This is a numbers file containing a 60212*2 dimentional matrix in which the first colomn is H score (multiplicity score), and the second one is F score (flexibility score). There are 60212 rows in total, corresponding to 60212 4-node gene circuits.
## Rdata.zip
This is a zip file of all 60212 Rdata files. Each Rdata file contains a information matrix of one gene circuit.
## code_for_HH_FF_calculation.R
This is a R file(3.6.2) containing one template code for calculating the H and F scores for all 60212 4-node gene circuits. Detailed numerical methods used in this code can be found in the paper.

## code_for_enrichment_scores_visualization.R
This is a R file(3.6.2) containing the code for caculating and visualizing(including bar plots and heatmaps) the enrichment scores for all 2-node motifs. Detailed numerical methods used in this code can be found in the paper.

## code_for_large_scale_HH_FF.R
This is a R file(3.6.2) containing one template code for caculating the H score (multiplicity score) and F score (flexibility score) for all large-scale gene regulatory networks. Detailed numerical methods used in this code can be found in the paper.
