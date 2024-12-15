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
				□ /TGFb_interactome/tgfb1_n4_run2_deseq2/iact_sel_noCorrFilter
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
			§ from : /TGFb_interactome/tgfb1_n4_run2_deseq2/iact_sel_noCorrFilter.Rdata
		○ enrichr.gem
			§ from : "/PXFG_paper/n4_deseq2_enrichr/results/enrichr.gem.Rdata"
		○ clusters.raw
			§ from : "PXFG_paper/n4_deseq2_enrichr/EnrichmentMap3_new_gmt/results/AutoAnnotate - Summary Network default node.csv"
	• method : 
		○ select submatrix using enrichment results
		○ using RCy3 control to cytoscape
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
				□ q
				□ db binary vector ; db.sum
				□ gene.and.mir.de
				□ evidence.exper
				□ score
				□ miRNA log2cpm
				□ iact pval and padj
				□ invFC y/n
				□ q per database
				□ gene info
				□ miR info

