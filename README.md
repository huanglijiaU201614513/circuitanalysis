# Circuit analysis GitHub repository
This README file describes the contents of the GitHub repository for the paper "What makes a functional gene regulatory network? A circuit motif analysis" written by Lijia Huang, Benjamin Clauss and Mingyang Lu.
## HH_FF_scores.numbers
A file containing a 60212*2 dimensional matrix, in which the first column shows the multiplicity scores (H), and the second one shows the flexibility scores (F). There are 60212 rows in total, corresponding to all the 60212 non-redundant 4-node gene circuits.
## Rdata.zip
A compressed data file containing all 60212 Rdata files, each of which specifies the network topology of a 4-node gene circuit.
## code_for_HH_FF_calculation.R
This is a R file (3.6.2) containing a template code for calculating the H and F scores for all 60212 4-node gene circuits. Detailed numerical methods used in this code can be found in the paper.

## code_for_enrichment_scores_visualization.R
This is a R file (3.6.2) containing the code for calculating and visualizing (including bar plots and heatmaps) the enrichment scores for all 2-node circuit motifs. Detailed numerical methods used in this code can be found in the paper.

## code_for_large_scale_HH_FF.R
This is a R file(3.6.2) containing one template code for caculating the H score (multiplicity score) and F score (flexibility score) for all large-scale gene regulatory networks. Detailed numerical methods used in this code can be found in the paper.
