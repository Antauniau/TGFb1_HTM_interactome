# iact2network_github.R
# Anton Roodnat, 2024, Ulster University
#
# • purposes of this script : 
#   ○ automatically draw networks in cytoscape !
#   ○ also draw heatmap of the same interaction submatrix
# • inputs : 
#   ○ iact_sel
#     § from : /TGFb_interactome/tgfb1_n4_run2_deseq2/iact_sel_noCorrFilter.Rdata
#   ○ enrichr.gem
#     § from : "/PXFG_paper/n4_deseq2_enrichr/results/enrichr.gem.Rdata"
#   ○ clusters.raw
#     § from : "PXFG_paper/n4_deseq2_enrichr/EnrichmentMap3_new_gmt/results/AutoAnnotate - Summary Network default node.csv"
# • method : 
#   ○ select submatrix using enrichment results
#   ○ using RCy3 control to cytoscape
# • outputs : 
#   ○ create cytoscape network
#   ○ heatmap per 
#     § miR-gene per theme
#     § miR-cluster
#     § miR-term
#   ○ interaction table 
#     § name : iact_sel_<filter$name>.txt such as : 
#       □ iact_sel_table_gene
#     § fields :
#       □ miRNA, geneSymb
#       □ correlation
#       □ q
#       □ db binary vector ; db.sum
#       □ gene.and.mir.de
#       □ evidence.exper
#       □ score
#       □ miRNA log2cpm
#       □ iact pval and padj
#       □ invFC y/n
#       □ q per database
#       □ gene info
#       □ miR info


rm(list = ls())
library(RCy3)
library(pheatmap)
library(ggplot2)
#library(magrittr)
#library(dplyr)


#*******************************************************************************
#---- dashboard ----
#*******************************************************************************

# output file from previous script TGFb_miRNAGene_corr_github.r
iact.sel.file = "iact_sel_github.Rdata"

# determine what analysis to perform and how
network.analysis     = T # determine connectiveness, nodes, edges, hubs etc

create.mir.clus      = F # create miR - enrichment cluster network
create.mir.term      = F # create miR - enrichment term network
show.heatmap         = T 
  heatmap.show.corr  = F # show correlation (T) or interaction score (F)
  heatmap.save       = F
send.to.cytoscape    = F # draw a network using cytoscape
use.iact.all.DE      = T # if true then only genes & mirs that are DE (FDR<0.05) will be used

do.tf.loops          = F # check for miR->TF->miR loops
do.tf2mirs           = F # what TFs may be responsible for DE miRs ? and are they DE as well?

# determine what to save
write.gene.list              = F # write iact.sel genes to file iact_sel_genelist.txt
write.gene.id.list           = F # write iact.sel geneids to file iact_sel_genelist.txt
write.iact.sel               = F # write interactions to file iact_sel.txt in WD
write.iact.table.mir.first   = F # create interaction table for report either miR or gene sorting key
write.iact.table.gene.first  = F # create interaction table for report either miR or gene sorting key
write.cluster.gene.mir       = F # create cluster-gene-mir table

# uncomment one line to select a certain biological theme
#filter.sel = "all"
filter.sel = c("XFM")           # exfolation material genes
#filter.sel = c("OX")            # oxidative stress
#filter.sel = c("ESM")           # extracellular signal molecules
#filter.sel  = c("G1/S cluster") # cell growth & arrest
#filter.sel  = c("UPR")          # unfolded protein response
#filter.sel  = c("ECM")          # extra-cellular matrix
#filter.sel  = c("EMT")          # epithelial to mesenchymal transition
#filter.sel  = c("TF")           # transcription factors

# enrichment data
terms.file    = "enrichr.gem.Rdata"
clusters.file = "AutoAnnotate - Summary Network default node.csv"
# read gmt file containing all genes per enrichment term
load("gmt_df_cor.Rdata")

# read TF->miR interactions
tf2mir <- read.csv( "miRNet-mir-tf-hsa.csv",header=T,sep=",")


#*******************************************************************************
#---- filter definitions ----
#*******************************************************************************

filters = list()

# exfoliation material
filters[["all"]]   = list(name="all",type="all",plot.DE.only=F,list=c())

# exfoliation material
filters[["XFM"]]   = list(name="XFM",type="gene.sel",plot.DE.only=F,list=c()) # list will be added below

# extracellular signal molecules
filters[["ESM"]]   = list(name="extracellular signal molecules",type="term.sel",plot.DE.only=T,list=c("GO:0008083","GO:0048018","KEGG:hsa04390"))

# extracellular matrix
filters[["ECM"]]   = list(name="extracelullar matrix (ECM)",type="term.sel",plot.DE.only=T,list=c("GO:0030198","REAC:R-HSA-1474244"))

# extracellular matrix
filters[["EMT"]]   = list(name="EMT",type="term.sel",plot.DE.only=T,list=c("GO:0051240","GO:0090189","GO:0010717","GO:1905332","GO:0007399","GO:0051241"))

# unfolded protein response
filters[["UPR"]]   = list(name="Unfolded protein response",type="term.sel",plot.DE.only=F,list=c("REAC:R-HSA-381038","REAC:R-HSA-381070","REAC:R-HSA-381119","GO:0036498"))

# oxidative stress (last entry is antiox activity)
filters[["OX"]]    = list(name="Oxidative stress", type="term.sel",plot.DE.only=F,list=c("GO:0062197","GO:0034599","GO:0016209"))

# all miR29 targets
filters[["mir29"]] = list(name="mir29",type="mir.sel",plot.DE.only=F,list=c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"))

# cluster : fibroblast proliferation ; G1/S
filters[["G1/S cluster"]] = list(name="G1S",type="cluster.sel",plot.DE.only=T,list=c("fibroblast proliferation ; G1/S"))

# transcription factors
filters[["TF"]] = list(name="TF",type="gene.sel",plot.DE.only=F,list=c()) # list will be added below


#*******************************************************************************
#---- create custom gene lists ----
#*******************************************************************************

#**** create XFM genelist ******************************************************

XFM.genes.duplicates <- c("EMILIN1","APCS","ELN","LAMA1","LAMB1","LAMC1","APCS",
                          "FN1","LAMA1","LAMB1","LAMC1","NID1","SDC2","EMILIN1",
                          "FBN1","LTBP2","FBN1","LTBP1","LTBP2","MFAP2","VTN",
                          "ADAM19","ADAM21","ADAMTS8","APCS","C1QA","C1QB","C1QC",
                          "C3","C4A","C4B","CLU","DSC2","DSC3","FBLN2","FBN1","FN1",
                          "LAMA1","LAMB1","LAMC1","NID1","SDC3","TIMP3","VCN","VTN",
                          "LOXL1","LTBP2","C3","CLU","APOE","LOXL1","ALDH3A1","ANXA1",
                          "ANXA7","APOA1","APOA2","APOA4","APOC3","APOE","ARHGAP42",
                          "BFSP1","BFSP2","C2CD4A","C3","CLU","COL18A1","CRYAA","CRYAB",
                          "CRYBB2","EMILIN1","ENO1","FBN1","FGA","FGB","FGG","FN1","FTH1",
                          "HBA1","HBB","HBD","HBE1","HBG1","HEATR1","HIST1H2BO","HIST1H4A",
                          "HIST2H2AC","HSPB1","IGHG1","IGKC","ITIH2","ITIH4","KERA","LDHA",
                          "LOXL1","LTBP2","MRVI1","MYH7B","MYH9","MYL6","PKM","PRDX2",
                          "PTPN5","RASD1","S100A6","SASH1","SERPINA1","SLC4A1","SPTA1",
                          "SPTB","TIMP3","TKT","TRIM31","TRRAP","VIM","VTN","YWHAB")
# add corrections for gene symbols that are outdated or plain wrong
XFM.genes.duplicates <- c(XFM.genes.duplicates,"H2BC17","H4C1","H2AC20","IRAG1","VCAN")
XFM.genes.pub <- sort(unique(XFM.genes.duplicates))
filters[["XFM"]]$list = XFM.genes.pub

#**** create TF genelist from GO:MF ********************************************

GO.MF.file   <- "GO_Molecular_Function_2021.txt"
GO.MF        <- read.csv( GO.MF.file,header=F,sep=")") # not so handy separators in this file
GO.MF.genestring <- GO.MF$V2[grep("GO:0000978",GO.MF$V1)]
TF.genes     <- strsplit(GO.MF.genestring, split = "\t")[[1]]
TF.genes     <- TF.genes[!TF.genes==""]
filters[["TF"]]$list = TF.genes




#*******************************************************************************
#---- load data and further filtering ----
#*******************************************************************************

# gene-mir interactions
load(file=iact.sel.file)
iact.sel.old <- iact.sel
# xx klopt dit wel :
iact.sel$gene.log2cpm = log2(iact.sel$gene.log2cpm)

# add if there is experimental evidence xx could also be added to the *corr.r script
iact.sel$evidence.exper= (iact.sel$q.mirtarbase>0 | iact.sel$q.tarbase>0)

# rather arbitrary attempt to determine strongest miR
iact.sel$score = iact.sel$evidence.exper + 0.5*iact.sel$gene.and.mir.de + 2*abs(iact.sel$cor.sam) + 0.2*iact.sel$mir.log2cpm

# for each gene, remove interactions for miRs with log2CPM lower than 8 orders wrt max(log2cpm) so 256x
idx.skip = rep(F,nrow(iact.sel))
for (gene in unique(iact.sel$gene.Symbol)){
  idx = which(iact.sel$gene.Symbol == gene)
  mir.log2cpm.max = max(iact.sel$mir.log2cpm[idx])
  idx.low = which(iact.sel$mir.log2cpm[idx] < mir.log2cpm.max - 10)
  idx.skip[idx[idx.low]] = T
}
cat("#interactions before removal of low-cpm miRs=",nrow(iact.sel),"\n")
iact.sel <- iact.sel[!idx.skip,]
cat("#interactions after  removal of low-cpm miRs=",nrow(iact.sel),"\n")


# reshuffle columns a bit for easier excel browsing
iact.sel.temp = iact.sel[,c("miRNA","gene.Symbol","cor.sam","q","db","db.sum","gene.and.mir.de","evidence.exper","score","mir.log2cpm",
                            "iact.padj","iact.pval","invFC","q.mirtarbase","q.targetscan.default","q.targetscan.noncons","q.tarbase","q.mirwalk","q.mirdip",
                            "gene.ID","gene.log2FC","gene.log2cpm","gene.pval","gene.padj",
                            "mir.log2FC","mir.pval","mir.padj")]
iact.sel <- iact.sel.temp ; rm(iact.sel.temp)

if(use.iact.all.DE) {
  # use conservative method : genes and mirs with FDR<0.05
  iact.sel <- iact.sel[which(iact.sel$gene.and.mir.de==T),]
} else {
  # use FDR for interaction rather than genes and mirs
  iact.sel <- iact.sel[which(iact.sel$iact.padj<0.05),]
}
iact.sel.all <- iact.sel

genes.all = unique(iact.sel.all$gene.Symbol)
mirs.all  = unique(iact.sel.all$miRNA)


# TGFb1 enrichment terms
load(terms.file)
terms <- enrichr.gem
rm(enrichr.gem)

# enrichment clusters
clusters.raw <- read.csv( clusters.file,header=T,sep=",")
clusters <- data.frame(nr=clusters.raw$X__mclCluster,
                       name=clusters.raw$shared.name,
                       term.list=clusters.raw$EnrichmentMap..Name,
                       genes.map=clusters.raw$EnrichmentMap..Genes,
                       genes="")

for (k in 1:nrow(clusters)) {
  # what terms are in this cluster ?
  cluster = clusters[k,]
  cluster.terms = strsplit(cluster$term.list,split = ",")[[1]]
    # what genes are in each term ?
  cluster.genes = c()  
  for (cluster.term in cluster.terms){
      term.idx = grep(cluster.term,terms$GO.ID,ignore.case = T)
      term.genes = strsplit(terms$Genes[term.idx],split = "\t")[[1]]
      cluster.genes <- c(cluster.genes,term.genes)
  }
  cluster.genes <- unique(cluster.genes)
  clusters$genes[k] <- paste(unlist(cluster.genes),collapse=",")
  # does my list agree with the list from cytoscape ?
}


#*******************************************************************************
#---- network analysis of whole interactome ----
#*******************************************************************************

if (network.analysis) {
  
  cat("number of genes =",length(genes.all),"\n")
  cat("number of miRs  =",length(mirs.all),"\n")
  cat("number of edges =",nrow(iact.sel.all),"\n")
  Emax = length(genes.all)*length(mirs.all)
  Emin = max(length(genes.all),length(mirs.all))
  cat("E=",nrow(iact.sel.all)," Emin=",Emin," Emax=",Emax,"\n")
  density = (nrow(iact.sel.all)-Emin)/(Emax-Emin)
  cat("density=",100*density,"%\n")
  # degree per gene
  
  M.iact = matrix(F,nrow=length(genes.all),ncol=length(mirs.all))
  rownames(M.iact)=genes.all ; colnames(M.iact)=mirs.all
  for (k in 1:nrow(iact.sel.all)) {
    r = which(iact.sel.all$gene.Symbol[k]==genes.all)
    c = which(iact.sel.all$miRNA[k]==mirs.all)
    M.iact[r,c]=T
  }
  genes.degree = sort(rowSums(M.iact+0),decreasing = T)
  mirs.degree  = sort(colSums(M.iact+0),decreasing = T)
  #barplot(height = mirs.degree, names.arg = names(mirs.degree),horiz = T,las=2)
  print(genes.degree[1:20])
  print(mirs.degree[1:20])
  cat("average gene degree =",mean(genes.degree),"=",
      nrow(iact.sel.all)/length(genes.all),"\n")
  cat("average miR degree  =",mean(mirs.degree),"\n")
  
}

#*******************************************************************************
#---- miR -> TF-> miR loops ----
#*******************************************************************************


if(do.tf.loops) {
  genes.subset <- filters[["TF"]]$list
  iact.sel     <- iact.sel[iact.sel$gene.Symbol %in% genes.subset,]
  iact.tf      <- unique(iact.sel$gene.Symbol)
  cat("miR -| TF interactions : ",nrow(iact.sel),"for ",length(iact.tf),"TFs")
  
  # create dataframe to store loop
  tf.loop.df <- data.frame(miRNA="",gene="",miRNet_id="",action_type="")
  
  for(k in 1:nrow(iact.sel)) {
    # for each interaction : is the reverse in tf2mir ?
    
    # try to find stemloop symbol (perhaps there is a database to do this)
    miRNA_stemloop <- tolower(strsplit(iact.sel$miRNA[k],split = "-5p|-3p")[[1]][1])
    tf2mir.idx = which(tf2mir$mir_id==miRNA_stemloop & tf2mir$symbol==iact.sel$gene.Symbol[k])
    if(length(tf2mir.idx)>0) {
      print(tf2mir[tf2mir.idx,])
      tf.loop.df.x <- data.frame(miRNA=iact.sel$miRNA[k],gene=iact.sel$gene.Symbol[k],
                                 miRNet_id=tf2mir$mirnet[tf2mir.idx],
                                 action_type=tf2mir$action_type[tf2mir.idx])
      tf.loop.df <- rbind(tf.loop.df,tf.loop.df.x)
    }
  }
}


#*******************************************************************************
#---- create miR-gene subnetwork based on cluster/term/gene/miR ----
#*******************************************************************************

for (filter.name in filter.sel) {
  
  filter = filters[[filter.name]]

  if(filter$type == "cluster.sel") {
    #== use cluster ==
    # find idx
    cluster.idx <- grep(filter$list,clusters$name,ignore.case = T)
    genes.subset <- unique(strsplit(clusters$genes[cluster.idx],split = ",")[[1]])

    # alternative way of obtaining all genes in this cluster
    terms.this.cluster <- strsplit(clusters$term.list[cluster.idx],split = ",")[[1]]
    # find genes in gmt.df
    genes.this.cluster = c()
    for (term in terms.this.cluster) {
      gmt.df.idx = which(term==gmt.df$id)
      genes.this.term = strsplit(gmt.df$genes[gmt.df.idx],split = "\t")[[1]]
      genes.this.cluster = c(genes.this.cluster,genes.this.term)
    }
    genes.subset <- unique(genes.this.cluster)
    
    #== create subset ==
    iact.sel <- iact.sel.all[iact.sel.all$gene.Symbol %in% genes.subset,]
  }

  if(filter$type == "term.sel") {
    GO.ID.subset.idx = which(terms$GO.ID %in% filter$list)
    # create one string out of several strings :
    genes.subset.s <- paste(unlist(terms$Genes[GO.ID.subset.idx]),collapse = ",")
    genes.subset   <- unique(strsplit(genes.subset.s,split=",")[[1]])
    
    # alternative way of obtaining all genes in this cluster
    # find genes in gmt.df
    genes.this.cluster = c()
    for (term in filter$list) {
      gmt.df.idx = which(term==gmt.df$id)
      genes.this.term = strsplit(gmt.df$genes[gmt.df.idx],split = "\t")[[1]]
      genes.this.cluster = c(genes.this.cluster,genes.this.term)
    }
    genes.subset <- unique(genes.this.cluster)
  
    iact.sel <- iact.sel[iact.sel$gene.Symbol %in% genes.subset,]
  }

  if(filter$type == "mir.sel") {
    mir.v = c()
    for (mir in filter$list) {
      idx = grep(mir,mirs.all)
      if (length(idx)>0) {
        mir.v = c(mir.v,mirs.all[idx])
      }
    }
    iact.sel <- iact.sel.all[iact.sel.all$miRNA %in% mir.v,]
  }

  if(filter$type == "gene.sel") {
    iact.sel.idx = which(iact.sel.all$gene.Symbol %in% filter$list)
    iact.sel <- iact.sel.all[iact.sel.idx,]
  }
}


#---->> send subnetwork to cytoscape ----

if(send.to.cytoscape) {
  deleteAllNetworks()
  
  iact.sel.idx <- sort(iact.sel$miRNA,index.return=T)
  iact.sel <- iact.sel[iact.sel.idx$ix,]
  
  
  node.mirs   <- iact.sel[,c("miRNA","mir.log2FC","mir.log2cpm")]
  #node.mirs   <- unique(node.mirs)
  node.genes  <- iact.sel[,c("gene.Symbol","gene.log2FC","gene.log2cpm")]
  #node.genes  <- unique(node.genes)
  node.groups <- c(rep("mir",nrow(node.mirs)),rep("gene",nrow(node.genes)))
  nodes.id    <- c(node.mirs$miRNA,node.genes$gene.Symbol)
  nodes <- data.frame(id=nodes.id,
                      group=node.groups,
                      log2FC=c(node.mirs$mir.log2FC,node.genes$gene.log2FC),
                      log2cpm=c(node.mirs$mir.log2cpm,node.genes$gene.log2cpm),
                      #score=as.integer(c(20,10,15,5)), # integers
                      stringsAsFactors=FALSE)
  
  
  
  
  # construct edges
  edge.sources <- iact.sel$miRNA
  edge.targets <- iact.sel$gene.Symbol
  
  edges <- data.frame(source=edge.sources,
                      target=edge.targets,
                      interaction=rep("inhibits",nrow(iact.sel)),
                      pval.iact=1-iact.sel$q,
                      cor.iact=iact.sel$cor.sam,
                      #weight=c(5.1,3.0,5.2,9.9), # numeric
                      stringsAsFactors=FALSE)
  
  createNetworkFromDataFrames(nodes=nodes, edges=edges[,1:2], title="Anton", collection="miR-gene interactome")
  
  
  #createNetworkFromDataFrames(edges=edges[,1:2])
  #edges <- edges %>%
  key.list=paste(edge.sources, "(interacts with)",edge.targets)
  #mutate(key=paste(edge.sources, "(interacts with)",edge.targets))
  
  edges$key.list = key.list
  edges.col.idx  = !colnames(edges) %in% c("source","target")
  edges.attr = edges[,edges.col.idx]
  
  loadTableData(edges.attr,data.key.column = "key.list", table = "edge")
  #pipo123
  
  #---->> cytoscape style and mapping ----
  # mapping :
  # node colour -> log2FC
  # node size   -> log2cpm?
  # node shape  -> circles for gene & mir, something else for enrichment
  # edge colour -> correlation ?
  # edge size   -> pval of interaction ? number of databases ?
  
  style.name = "mirGene"
  if (style.name %in% getVisualStyleNames()){
    deleteVisualStyle(style.name)  
  }
  
  getVisualPropertyNames()
  
  defaults      <- list(NODE_SHAPE="ELLIPSE",
                   NODE_SIZE=20,
                   EDGE_TRANSPARENCY=120,
                   NODE_BORDER_WIDTH=0,
                   NODE_TRANSPARENCY=100,
                   EDGE_WIDTH=1.1,
                   NETWORK_FORCE_HIGH_DETAIL=T,
                   EDGE_STROKE_UNSELECTED_PAINT=rgb(0.2,0.7,0.2))
                   #NODE_LABEL_POSITION="N,S,c,0.00,0.00")
  nodeLabels    <- mapVisualProperty('node label','id','p')
  
  col.down      <- colorRampPalette(c(rgb(0.0,1.0,1.0), rgb(0.8,1.0,1.0)), alpha = F)(5)
  col.up        <- colorRampPalette(c(rgb(1.0,0.9,0.8), rgb(1.0,0.5,0.0)), alpha = F)(5)
  
  # colours to indicate strength
  col.strength  <- colorRampPalette(c(rgb(1.0,0.7,1.0), rgb(0.9,0.9,0.9)), alpha = F)(2)
  
  
  #nodeFills.col <- colorRampPalette(c(rgb(0,1,1), rgb(1,0.5,0)), alpha = F)(10)
  nodeFills.col <- c(col.down,col.up)
  nodeFills.val <- c(seq(-3, 3, length=10))
  nodeFills     <- mapVisualProperty('node fill color','log2FC','c',nodeFills.val,nodeFills.col)
  
  nodeSize.size <- c(seq(30,50,length=10))
  nodeSize.val  <- c(seq(1,11,length=10))
  nodeSize      <- mapVisualProperty('node size','log2cpm','c',nodeSize.val,nodeSize.size)
  arrowShapes   <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',
                                     c("activates","inhibits","interacts"),
                                     c("None","None","None"))
                                     #c("Arrow","T","None")) # or T=inhibit
  
  edgeWidth.val <- c(-1.0,-0.3)
  edgeWidth.w   <- c( 2.5,0.5)
  edgeWidth     <- mapVisualProperty('Edge Width','cor.iact','c',edgeWidth.val,edgeWidth.w)
  
  edgePaint.val <- c(seq(0.0,0.1,length=2))
  edgePaint     <- mapVisualProperty('Edge Target Arrow Unselected Paint','pval.iact','c',edgePaint.val,col.strength)
  
  
  createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,nodeSize,arrowShapes,edgeWidth,edgePaint))
  #createVisualStyle(style.name, defaults)
  setVisualStyle(style.name)
  
  #---->> layout
  layoutNetwork('force-directed defaultSpringCoefficient=.001 defaultSpringLength=40')
  # remove overlaps
  # organic
}

#*******************************************************************************
#---- >> show miR-gene heatmap ----
#*******************************************************************************
if(show.heatmap) {
  
  # heatmap colours

  # define colours
  breaks       = c(seq(-2, -0.1, length=10), 0, seq(0.1,2,length=10))
  breaks_ext   = c(seq(-4, -0.1, length=10), 0, seq(0.1,4,length=10))
  
  #set different color vectors for each interval
  col_down    <- colorRampPalette(c(rgb(0,0,1,0.1), rgb(0,0,1,0.8)), alpha = TRUE)(9)
  col_none    <- rep("white", 2)
  col_up      <- colorRampPalette(c(rgb(1,0,0,0.1), rgb(1,0,0,0.8)), alpha = TRUE)(10)
  col_up_ext  <- colorRampPalette(c(rgb(1,0,0,0.1), rgb(1,0,0,0.8)), alpha = TRUE)(20)
  colours_bin <- colorRampPalette(c(rgb(0,1,0,0.1), rgb(0,0,1,0.1)), alpha = TRUE)(10)
  colours_de  <- c(col_down,col_none,col_up)
  colours_ext <- c(col_down,col_none,col_up_ext)
  # map correlation to green 
  #colours_cor_neg <- colorRampPalette(c(rgb(0,0.5,0,1.0), "white"), alpha = TRUE)(10)
  colours_cor_neg2 <- colorRampPalette(c(rgb(1.0,0.3,0.0,0.8), rgb(1.0, 0.3, 0.0, 0.3) ), alpha = TRUE)(20)
  colours_cor_neg <- colorRampPalette(c(rgb(0,0.5,0,0.8), rgb(0.0, 0.5, 0.0, 0.3) ), alpha = TRUE)(20)
  #colours_cor_pos <- colorRampPalette(c("white",rgb(0.7,0.7,1.0,1.0)) , alpha = TRUE)(9)
  colours_cor_pos <- colorRampPalette(c(rgb(0.7,0.7,1.0,0.3),rgb(0.7,0.7,1.0,0.8)) , alpha = TRUE)(20)
  colours_cor_pos2 <- colorRampPalette(c(rgb(1.0,0.7,1.0,0.3),rgb(1.0,0.7,1.0,0.8)) , alpha = TRUE)(20) 
  
  colours_cor     <- c(colours_cor_neg2,"white",colours_cor_neg,"white",colours_cor_pos,"white",colours_cor_pos2)
  # number of breaks is number of colours + 1 because outer breaks are in fact limits:
  breaks_cor  <- c(seq(-2.0,-1.05,length=21),seq(-1.0,-0.05,length=21),seq(0.05,1.0,length=21),seq(1.05,2.0,length=21))
  
  colours_score_iact <- c("white",colorRampPalette(c(rgb(0,1,0,0.3), rgb(0,0,1,0.8)), alpha = TRUE)(20))
  breaks_score_iact  <- c(seq(0,7,length=22))
  
  colours_cpm <- c("white","grey")
  
  iact.sel.mirs  <- unique(iact.sel$miRNA)
  iact.sel.genes <- unique(iact.sel$gene.Symbol)
  M.mirs.genes = matrix(0,nrow=length(iact.sel.mirs),ncol=length(iact.sel.genes))
  mir.ann=c();gen.ann=c();mir.cpm=c()
  
  for (k in 1:nrow(iact.sel)) {
    mir.idx  <- which(iact.sel$miRNA[k]==iact.sel.mirs)
    gen.idx  <- which(iact.sel$gene.Symbol[k]==iact.sel.genes)
    if(heatmap.show.corr) {
      M.mirs.genes[mir.idx,gen.idx] = (2*iact.sel$gene.and.mir.de[k]-1)*min(-0.15,iact.sel$cor.sam[k])      
      if(iact.sel$gene.and.mir.de[k]==F & iact.sel$evidence.exper[k]){
        M.mirs.genes[mir.idx,gen.idx] =  1 + M.mirs.genes[mir.idx,gen.idx]  
      }
      if(iact.sel$gene.and.mir.de[k]==T & iact.sel$evidence.exper[k]){
        M.mirs.genes[mir.idx,gen.idx] = -1 + M.mirs.genes[mir.idx,gen.idx]  
      }
    } else {
      M.mirs.genes[mir.idx,gen.idx] = iact.sel$score[k]
    }
    mir.ann[mir.idx] = iact.sel$mir.log2FC[k]
    mir.cpm[mir.idx] = iact.sel$mir.log2cpm[k]
    gen.ann[gen.idx] = iact.sel$gene.log2FC[k]
  }
  
  #mir.ann = seq(-1,10,length=length(mir.ann))
  #mir.ann[3] = 4
  #mir.ann = seq(1,length(mir.ann))
  
  rownames(M.mirs.genes) <- iact.sel.mirs
  colnames(M.mirs.genes) <- iact.sel.genes
  
  annotation.colours <- list()
  # map colour range to variable range : log2FC=4 is max/min for colours_de
  col.len     = length(colours_de)
  col.len.mid = round((col.len-1)/2)
  
  mir.ann.pk  = max(abs(mir.ann))
  col.mir.idx.min = max(1,col.len.mid + round(0.5*col.len*min(mir.ann)/mir.ann.pk))
  col.mir.idx.max = col.len.mid + round(0.5*col.len*max(mir.ann)/mir.ann.pk)
  col.mir.idx = col.mir.idx.min:col.mir.idx.max
  #col.mir.idx = round(seq(col.mir.idx.min,col.mir.idx.max,length=10))
  
  gen.ann.pk  = max(abs(gen.ann))
  col.gen.idx.min = col.len.mid + round(0.5*col.len*min(gen.ann)/gen.ann.pk)
  col.gen.idx.max = col.len.mid + round(0.5*col.len*max(gen.ann)/gen.ann.pk)
  col.gen.idx = round(seq(col.gen.idx.min,col.gen.idx.max,length=10))
  
  
  #mir.ann.col.min.idx = max(1,col.len.mid + round(((min(mir.ann)/4.0))*col.len.mid))
  #mir.ann.col.max.idx = min(col.len,col.len.mid + round(((max(mir.ann)/4.0))*col.len.mid))
  annotation.colours$mir.ann <- colours_de[col.mir.idx]
  annotation.colours$gen.ann <- colours_de[col.gen.idx] #seq(-2,2,length=length(iact.sel.mirs))
  annotation.colours$mir.cpm <- colours_cpm
  

  #annotation.row <- data.frame(mir.ann=mir.ann)
  annotation.col <- data.frame(gen.ann=gen.ann)
  #annotation.row.max =  3
  #annotation.row.min = -3
  annotation.row <- data.frame(mir.ann=mir.ann,mir.cpm=mir.cpm)
  
  rownames(annotation.row) <- rownames(M.mirs.genes)
  rownames(annotation.col) <- colnames(M.mirs.genes) 
  
  
  if(heatmap.show.corr) {
    colours.heatmap <- colours_cor
    breaks.heatmap  <- breaks_cor
    heatmap.title <- paste0(filters[[filter.sel]]$name,"\n",
                            "standard : experimental(orange), predicted(green) ;\n",
                            "extended : experimental(pink), predicted(blue)")
    
  } else {
    colours.heatmap <- colours_score_iact
    breaks.heatmap  <- breaks_score_iact
    heatmap.title <- paste0(filters[[filter.sel]]$name,"\n","interaction score")
  }
  
  if (heatmap.save) {
    heatmap.filename <- paste0(filter$name,"_heatmap.png")
    pheatmap(M.mirs.genes, 
             main = heatmap.title,
             fontsize = 8,fontsize_row = 10,fontsize_col = 10,
             cellwidth = 12, cellheight = 12,
             color = colours.heatmap,breaks = breaks.heatmap,
             #cellwidth=12,cellheight = 12,
             cluster_rows=T,cluster_cols=T, # group genes
             clustering_distance_rows="binary",
             clustering_distance_cols="binary",
             clustering_method = "complete", #"ward.D",
             annotation_row = annotation.row,
             annotation_col = annotation.col,
             annotation_colors = annotation.colours,
             treeheight_row = 10,
             treeheight_col = 10,
             filename = heatmap.filename
    )
  } else {
    pheatmap(M.mirs.genes, 
             main = heatmap.title,
             fontsize = 8,fontsize_row = 10,fontsize_col = 10,
             cellwidth = 12, cellheight = 12,
             color = colours.heatmap,breaks = breaks.heatmap,
             #cellwidth=12,cellheight = 12,
             cluster_rows=T,cluster_cols=T, # group genes
             clustering_distance_rows="binary",
             clustering_distance_cols="binary",
             clustering_method = "complete", #"ward.D",
             annotation_row = annotation.row,
             annotation_col = annotation.col,
             annotation_colors = annotation.colours,
             treeheight_row = 10,
             treeheight_col = 10
             )
  }
}

#*******************************************************************************
#---- create miR-cluster network / matrix for miR (sub)set ----
#*******************************************************************************

if(create.mir.clus) {

  cat("creating miR-cluster matrix for",nrow(clusters),"clusters and ",length(mirs.all),"miRs\n")
  M.mirs.genes.all = matrix(data=0,ncol=length(mirs.all),nrow=length(genes.all))
  for(k in 1:nrow(iact.sel.all)) {
    r = which(genes.all %in% iact.sel.all$gene.Symbol[k])
    c = which(mirs.all  %in% iact.sel.all$miRNA[k])
    M.mirs.genes.all[r,c] = M.mirs.genes.all[r,c] + 1
  }
  rownames(M.mirs.genes.all) = genes.all
  colnames(M.mirs.genes.all) = mirs.all
  
  M.mirs.clus = matrix(data=0,ncol=length(mirs.all),nrow=nrow(clusters))
  M.mirs.clus.abs = M.mirs.clus
  
  clus.gene.mir.df <- data.frame(clusters=clusters$name,genes=clusters$genes,mirs="")

  for(k in 1:nrow(clusters)) {
    
    # filter genes
    genes.this.clus = unique(strsplit(clusters$genes[k],split = ",")[[1]])
    gene.idx = genes.all %in% genes.this.clus
    
    # create submatrix from total gene-miR matrix, some genes not regulated
    M.mirs.genes.sub = matrix(M.mirs.genes.all[gene.idx,],nrow=length(which(gene.idx)))
    rownames(M.mirs.genes.sub) = genes.all[gene.idx]
    colnames(M.mirs.genes.sub) = mirs.all
    
    # number of genes per miR-cluster combination
    genes.reg.idx.v  = which(rowSums(M.mirs.genes.sub)>0)
    genes.reg.clus.v = rownames(M.mirs.genes.sub[genes.reg.idx.v,])
    mirs.this.clus.v = colSums(M.mirs.genes.sub)
    # normalise such that sum of each cluster is 1 : show dominant miRs
    #M.mirs.clus[k,] = mirs.this.clus.v/sum(mirs.this.clus.v)
    if(length(genes.reg.clus.v)>1) {
      M.mirs.clus[k,] = mirs.this.clus.v/length(genes.reg.clus.v)
    } else {
      M.mirs.clus[k,] = mirs.this.clus.v
    }
    #M.mirs.clus[k,] = mirs.this.clus.v/length(genes.reg.clus.v)
    M.mirs.clus.abs[k,] = mirs.this.clus.v
    #if(k==31){pipo} # test for ECM
    
    # create dataframe for supplemental cluster-genes-miRs
    #clus.gene.mir.df[k]$cluster=clusters$name[k]
    #clus.gene.mir.df[k]$genes=clusters$genes[k]
    clus.gene.mir.df$mirs[k]=paste(unlist(mirs.all[mirs.this.clus.v>0]),collapse=",")
  }
  
  if(write.cluster.gene.mir){
    write.table(clus.gene.mir.df,file=paste0("cluster_gene_mir_table.csv"),
                row.names = F,col.names = T,quote = F,sep = "\t")
    
  }
  
  rownames(M.mirs.clus) = clusters$name
  colnames(M.mirs.clus) = mirs.all

  rownames(M.mirs.clus.abs) = clusters$name
  colnames(M.mirs.clus.abs) = mirs.all
  
    
  # breaks & colours : |breaks| = |colours|+1
  breaks_terms  = seq(0,0.6,length=51) # for miR-cluster map use with col_up_ext
  colours_terms = c("white",colorRampPalette(c(rgb(0.4,0.0,0.5,0.1), rgb(0.4,0.0,0.5,0.8)), alpha = TRUE)(25),
                            colorRampPalette(c(rgb(0.4,0.0,0.5,0.8), rgb(0.4,0.0,0.5,1.0)), alpha = TRUE)(25))
  pheatmap(M.mirs.clus,breaks = breaks_terms,color = colours_terms,
           cluster_rows = T,cluster_cols = T,
           #file = "Figxx_miR_vs_cluster.tiff",
           #clustering_distance_rows="binary",clustering_distance_cols="binary",
           #clustering_method = "complete", #"ward.D",
  )
  
  cluster.regulation.v = colSums(M.mirs.clus)
  
  plot(cluster.regulation.v)

}

#*******************************************************************************
#---- create miR-term network / matrix for miR (sub)set ----
#*******************************************************************************
# The idea : given a list of miRs :
# - what terms may be regulated by these ?
# - do miRs cooperate ?
# - how strongly ? how many genes for instance ?
#   . map #genes to edge width (network) or to intensity (heatmap)

if(create.mir.term) {

  M.mirs.terms = matrix(data=0,ncol=length(mirs.all),nrow=nrow(terms))
  for(k in 1:nrow(terms)) {
    # filter genes
    genes.this.term = strsplit(terms$Genes[k],split = ",")[[1]]
    gene.idx = genes.all %in% genes.this.term
    M.mirs.genes.sub = matrix(M.mirs.genes.all[gene.idx,],nrow=length(which(gene.idx)))
    rownames(M.mirs.genes.sub) = genes.all[gene.idx]
    colnames(M.mirs.genes.sub) = mirs.all
    M.mirs.terms[k,] = colSums(M.mirs.genes.sub) / length(genes.this.term)
  }
  rownames(M.mirs.terms) = terms$Description
  colnames(M.mirs.terms) = mirs.all
  pheatmap(M.mirs.terms,
           cellwidth = 10,cellheight = 10,fontsize = 8,filename = "pipo.png")
}

if(write.gene.list) {
  write.table(unique(iact.sel$gene.Symbol),file=paste0("iact_sel_genelist_",filter$name,".txt"),
              row.names = F,col.names = F,quote = F)
}

if(write.gene.id.list) {
  write.table(unique(iact.sel$gene.ID),file=paste0("iact_sel_geneidlist_",filter$name,".txt"),
              row.names = F,col.names = F,quote = F)
}


if(write.iact.sel){
  iact.sel.s <- iact.sel
  iact.sel.s$cor.sam = sprintf("%1.2f",iact.sel$cor.sam)
  iact.sel.s$q = sprintf("%1.2f",iact.sel$q)
  iact.sel.s$mir.log2cpm = sprintf("%1.2f",iact.sel$mir.log2cpm)
  iact.sel.s$score = sprintf("%1.2f",iact.sel$score)
  iact.sel.s$iact.padj = sprintf("%1.2e",iact.sel$iact.padj)
  iact.sel.s$iact.pval = sprintf("%1.2e",iact.sel$iact.pval)
  iact.sel.s$db = paste0("x",iact.sel$db)
  write.table(iact.sel.s,file=paste0("iact_sel_",filter$name,".txt"),
              row.names = F,col.names = T,quote = F)
}

if(write.iact.table.mir.first){
  iact.sel.table=c()
  mirs.sel = unique(iact.sel$miRNA)
  iact.sel.df = data.frame(mir=mirs.sel,genes.basic="",genes.extended="")
  mirs.sel.nr = c()
  for (k in 1:length(mirs.sel)) {
    mir=mirs.sel[k]
    genes.this.mir.basic  = iact.sel$gene.Symbol[which(iact.sel$miRNA==mir & iact.sel$gene.and.mir.de==T)]
    genes.this.mir.extend = iact.sel$gene.Symbol[which(iact.sel$miRNA==mir & iact.sel$gene.and.mir.de==F)]
    iact.sel.df$mir[k] = mir
    iact.sel.df$genes.basic[k]    = paste(unlist(genes.this.mir.basic),collapse=",")
    iact.sel.df$genes.extended[k] = paste(unlist(genes.this.mir.extend),collapse=",")
    mirs.sel.nr[k] = strsplit(mir,"-")[[1]][3]
  }
  # nicely sort by miR number so 146 should come after 29
  mirs.sel.nr = sub("a",".1",mirs.sel.nr,fixed = T)
  mirs.sel.nr = sub("b",".2",mirs.sel.nr,fixed = T)
  mirs.sel.nr = sub("c",".3",mirs.sel.nr,fixed = T)
  mirs.sel.ix = sort.int(as.double(mirs.sel.nr),index.return=T)
  iact.sel.df = iact.sel.df[mirs.sel.ix$ix,]

  write.table(iact.sel.df,file="iact_sel_table_mirFirst.csv",row.names = F)
}

#---- create interaction table, gene first ----
if(write.iact.table.gene.first){
  iact.sel.table=c()
  genes.sel = unique(iact.sel$gene.Symbol)
  #iact.sel.df2 = data.frame(gene=genes.sel,standard_interactions="",extended_interactions="")
  iact.sel.df2 = data.frame(gene=genes.sel,interactions="")
  for (k in 1:length(genes.sel)) {
    gene=genes.sel[k]
    mirs.this.gene.standard = iact.sel$miRNA[which(iact.sel$gene.Symbol==gene & iact.sel$gene.and.mir.de==T)]
    mirs.this.gene.extended = iact.sel$miRNA[which(iact.sel$gene.Symbol==gene & iact.sel$gene.and.mir.de==F)]
    iact.sel.df2$gene[k] = gene
    
    if (length(mirs.this.gene.standard)>0) {
      iact.sel.df2$interactions[k]= paste(unlist(mirs.this.gene.standard),collapse=", ")
    }
    if (length(mirs.this.gene.standard)>0 & length(mirs.this.gene.extended)>0) {
      iact.sel.df2$interactions[k] = paste0(iact.sel.df2$interactions[k],", ")
    }
    if (length(mirs.this.gene.extended)>0) {
      iact.sel.df2$interactions[k] = paste0(iact.sel.df2$interactions[k],"(",
                                            paste(unlist(mirs.this.gene.extended),collapse=", "),")")
    }
    #iact.sel.df2$interactions[k]= paste0(paste(unlist(mirs.this.gene.standard),collapse=", "),", (",
    #                                     paste(unlist(mirs.this.gene.extended),collapse=", "),")")
    #iact.sel.df2$extended_interactions[k]  = paste(unlist(mirs.this.gene.extended),collapse=", ")
  }
  iact.sel.df2 = iact.sel.df2[sort.list(iact.sel.df2$gene),]
  iact_sel_table_geneFirst_name = paste0("iact_sel_table_geneFirst_",filter$name,".csv")
  cat("->FILE writing gene |- mirs table to :",iact_sel_table_geneFirst_name,"\n")
  write.table(iact.sel.df2,file=iact_sel_table_geneFirst_name,row.names = F,sep="\t",quote = F)
}

#---- network analysis of subnetwork ----
if(network.analysis){
  cat("-----------------------------------------------\n")
  cat("network analysis of subnetwork",filter$name,"\n")
  cat("-----------------------------------------------\n")
  cat(nrow(iact.sel)," interactions (",length(unique(iact.sel$gene.Symbol))," genes",length(unique(iact.sel$miRNA))," miRs)\n")
  cat(length(which(iact.sel$gene.and.mir.de==T)),"standard interactions (",
      length(unique(iact.sel$gene.Symbol[which(iact.sel$gene.and.mir.de==T)])),"genes",
      length(unique(iact.sel$miRNA[which(iact.sel$gene.and.mir.de==T)])),"miRs )\t",
      length(which(iact.sel$gene.and.mir.de==T & iact.sel$evidence.exper==T)),"experimental,\t",
      length(which(iact.sel$gene.and.mir.de==T & iact.sel$evidence.exper==F)),"predicted\n")
  cat(length(which(iact.sel$gene.and.mir.de==F)),"extended interactions (",
      length(unique(iact.sel$gene.Symbol[which(iact.sel$gene.and.mir.de==F)])),"genes",
      length(unique(iact.sel$miRNA[which(iact.sel$gene.and.mir.de==F)])),"miRs )\t",
      length(which(iact.sel$gene.and.mir.de==F & iact.sel$evidence.exper==T)),"experimental,\t",
      length(which(iact.sel$gene.and.mir.de==F & iact.sel$evidence.exper==F)),"predicted\n")

  # create submatrix
  M.iact.sel = M.iact[which(rownames(M.iact) %in% unique(iact.sel$gene.Symbol)),
                      which(colnames(M.iact) %in% unique(iact.sel$miRNA))]
  
    # gene with highest number of edges
  
  genes.degree = sort(rowSums(M.iact.sel+0),decreasing = T)
  mirs.degree  = sort(colSums(M.iact.sel+0),decreasing = T)
  #barplot(height = mirs.degree, names.arg = names(mirs.degree),horiz = T,las=2)
  # for paper :
  print(paste(unlist(paste0(" ",names(mirs.degree[1:20])," (",mirs.degree[1:20]," interactions)")),collapse=","))
  
  print(genes.degree[1:20])
  print(mirs.degree[1:20])
  cat("average gene degree =",mean(genes.degree),"\n")
  cat("average miR degree  =",mean(mirs.degree),"\n")

  # mir with highest number of edges
  # strongest mir per gene
  if(F){
    genes = unique(iact.sel$gene.Symbol)
    for (gene in genes) {
      idx = which(iact.sel$gene.Symbol==gene)
      cat(gene," : ",iact.sel$score[idx],"\n")
    }
  }
  
}

