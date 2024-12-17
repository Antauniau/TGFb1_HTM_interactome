# TGFb1_HTM_interactome
miRNA-gene interactome of normal human trabecular meshwork cells treated with TGFb1

# Basic flow
This repository contains all necessary files to create a miRNA-gene interactome of normal human trabecular meshwork cells treated with TGFb1 and show the results either in R (Rstudio) or in cytoscape.

There are two steps :

#1 create object iact_sel : This is achieved by running script TGFb_miRNAGene_corr_github.r
The options of this script need to be set in section 'dashboard' of this script.

#2 plot networks and matrices : This is done by running script iact2network_github.r
Once again the options can be set in section 'dashboard'. If one desires to plot a network in cytoscape then cytoscape needs to be started first.

# more detailed description of script TGFb_miRNAGene_corr.r
	• purpose : determine miR-gene interactions using both inverse foldchange and correlation
	• inputs :
		○ gene information such as log2FC, log2CPM, FDR (csv)
		○ genes.counts (needed for correlation calculation)
		○ mirs information such as log2FC, log2CPM, FDR (csv)
		○ mirs.counts
		○ iact.db : interaction database obtained from online databases like miRTarBase etc
		○ cor.lim : correlation limit
		○ iact.db.pval.lim : limit of interaction p-value
		○ mirdip results : *.Rdata file
	• outputs :
		○ iact.sel
			§ files created with this script : 
				□ iact_sel_github.csv
    				□ iact_sel_github.Rdata    
			§ fields :
				□ gene info
				□ mir info including FC etc
				□ correlation info
				□ iact info : which databases

# more detailed description of script iact2network.r
	• purposes : 
		○ automatically draw networks in cytoscape !
		○ also draw heatmap of the same interaction submatrix
	• inputs : 
		○ iact_sel
			§ from : iact_sel_github
		○ enrichr.gem
			§ from : "enrichr.gem.Rdata"
		○ clusters.raw
			§ from : "AutoAnnotate - Summary Network default node.csv"
	• method : 
		○ select submatrix of total interaction matrix using enrichment results
		○ using RCy3 control to cytoscape and automatically draw networks
	• outputs : 
		○ create cytoscape network
		○ heatmap per 
			§ miR-gene per theme
			§ miR-cluster
			§ miR-term
		○ interaction table 
			§ name : iact_sel_<filter$name>.txt such as : 
				□ iact_sel_table_gene
			§ fields :
				□ miRNA, geneSymb
				□ correlation
				□ q (=1-p_interaction)
				□ db binary vector : interaction is in what databases ? ; db.sum
				□ gene.and.mir.de : are gene and miR both FDR<0.05 ?
				□ evidence.exper  : is there experimental evidence ?
				□ score : interaction score
				□ miRNA log2cpm
				□ iact pval and padj : interaction pval and FDR
				□ invFC y/n : is there inversely correlated expression for this gene-miRNA pair ?
				□ q per database (=1-p_interaction)
				□ gene info : log2FC, log2CPM, FDR etc
				□ miR info  :  ,,

