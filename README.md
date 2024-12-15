# TGFb1_HTM_interactome
miRNA-gene interactome of normal human trabecular meshwork cells treated with TGFb1

This repository contains all necessary files to create a miRNA-gene interactome of normal human trabecular meshwork cells treated with TGFb1 and show the results either in R (Rstudio) or in cytoscape.

There are two steps :

#1 create object iact_sel : This is achieved by running script TGFb_miRNAGene_corr_github.r
The options of this script need to be set in section 'dashboard' of this script.

#2 plot networks and matrices : This is done by running script iact2network_github.r
Once again the options can be set in section 'dashboard'. If one desires to plot a network in cytoscape then first cytoscape needs to be started.
