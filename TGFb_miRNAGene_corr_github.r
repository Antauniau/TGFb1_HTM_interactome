# TGFb_miRNAGene_corr_github
# 2024 Anton Roodnat, Ulster University
#
# • purpose : determine miR-gene interactions using both inv FC and corr
# • inputs :
#   ○ genes (csv)
#   ○ genes.counts (needed for correlation calculation)
#   ○ mirs (csv)
#   ○ mirs.counts
#   ○ iact.db : interaction database
#   ○ iact.db.pval.lim : limit of interaction p-value
#   ○ mirdip results : *.Rdata file
#   ○ published_PEXG_miRs.xlsx : list of published miRNA-gene pairs + biological context
# • outputs :
#   ○ iact.sel
#     § files created with this script : 
#       iact_sel_github
#     § object fields :
#       □ gene info
#       □ mir info including FC etc
#       □ correlation info
#       □ iact info : which databases

rm(list = ls())
library(limma)
library(edgeR)

# correlation between mRNAs and miRs

#*******************************************************************************
#---- dashboard ----
#*******************************************************************************

# count limits
min.count.val.mir = 0
min.count.val.gen = 50

# DE limits (both miR and gene)
padj.lim          = 0.05
pval.lim          = 1.00

# interactions
iact.db.pval.lim  = 0.10
cor.lim           = 0.1

# save output
save.iact.sel     = T
iact.sel.file     = "iact_sel_github"

# compare results with publications
test.results      = T

# network analysis : how many genes, mirs, what hub nodes ?
network.analysis  = T

# inputs :
# genes : pval, padj, counts
# miRs  : pval, padj, counts

#*******************************************************************************
#---- retrieve input files ----
#*******************************************************************************

# genes
#******
genes.file   <- "TGFb1_htm_genetable_deseq2_n4.csv"
genes        <- read.csv(genes.file,header=T,sep="\t")
genes$gene_id <- substr(genes$gene_id,1,15) # remove last part
cat("read ",nrow(genes)," genes\n")

genes.counts.file <- "tgfb1_counts_genenames1.txt"
genes.counts <- read.csv(genes.counts.file,header=T,sep="\t")
genes.counts$gene_id <- substr(genes.counts$gene_id,1,15)
genes.counts.cols.discard <- c("TM_C_3","TM_T_3")
genes.counts <- genes.counts[,!(colnames(genes.counts) %in% genes.counts.cols.discard)]

# miRs
#******
mirs.file    <- "TGFb1_miRNA_DEGs_v2.csv"
mirs         <- read.csv(mirs.file,header=T,sep="\t")
mirs.counts.file  <- "miR_counts.txt"
mirs.counts  <- read.csv(mirs.counts.file,header=T,sep="\t")
# filter out columns :
mirs.counts.cols.keep <- c("Geneid","X021_control","X021_TGFb1",
                                    "X029_control","X029_TGFb1",
                                    "LGP3_control","LGP3_TGFb1",
                                    "LGP5_control","LGP5_TGFb1")
mirs.counts <- mirs.counts[,colnames(mirs.counts) %in% mirs.counts.cols.keep]
mirs

# database
#**********
# database has been split up due to github filesize limitations
database.files = c("miRtargetsCombined1.Rdata","miRtargetsCombined2.Rdata","miRtargetsCombined3.Rdata")

#*******************************************************************************
#---- normalisation ----
#*******************************************************************************

# genes
y  <- DGEList(counts=genes.counts[,3:ncol(genes.counts)], genes=genes.counts[,1:2]) 
o  <- order(rowSums(y$counts), decreasing=T)
y  <-y[o,]
d  <- duplicated(y$genes$gene_id)
y  <- y[!d,]
y$samples$lib.size <- colSums( y$counts)
donor <- factor(rep(c(1,2,4,5),2))
treat <- factor(c(rep("C",4),rep("T",4)))
design <- model.matrix(~donor+treat)
rownames(design) <- colnames(y)
keep <- filterByExpr(y,design,min.count=min.count.val.mir)
y <- y[keep, , keep.lib.sizes = F] 
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y,method="TMM") 
genes.log.M <- cpm(y ,log=T)
#rownames(mirs.log) <- y$genes$genes
colnames(genes.log.M) <- colnames(genes.counts[3:ncol(genes.counts)])
genes.log.counts <- data.frame(y$genes,genes.log.M)

# mirs
rm(y)
y  <- DGEList(counts=mirs.counts[,2:ncol(mirs.counts)], genes=mirs.counts[,1])  
o  <- order(rowSums(y$counts), decreasing=T)
y  <-y[o,]
d  <- duplicated(y$genes$genes)
y  <- y[!d,]
y$samples$lib.size <- colSums( y$counts)
donor <- factor(c("X021","X021","X029","X029","LGP3","LGP3","LGP5","LGP5"))
treat <- factor(rep(c("C","T"),4))
design <- model.matrix(~donor+treat)
rownames(design) <- colnames(y)
keep <- filterByExpr(y,design,min.count=min.count.val.mir)
y <- y[keep, , keep.lib.sizes = F] 
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y,method="TMM") 
mirs.log.M <- cpm(y ,log=T)
#rownames(mirs.log) <- y$genes$genes
colnames(mirs.log.M) <- colnames(mirs.counts[2:ncol(mirs.counts)])
mirs.log.counts <- data.frame(mir_id=y$genes$genes,mirs.log.M)

#*******************************************************************************
# ---- filter out on padj & pval ----
#*******************************************************************************

genes.de.idx <- which(genes$pval < pval.lim & genes$padj < padj.lim )
genes.de     <- genes[genes.de.idx,]
genes.log.counts.idx <- which(genes.log.counts$gene_id %in% genes$gene_id[genes.de.idx])
genes.log.counts <- genes.log.counts[genes.log.counts.idx,]

mirs.de.idx <- which(mirs$PValue < pval.lim & mirs$FDR < padj.lim )
mirs.de     <- mirs[mirs.de.idx,]
mirs.log.counts.idx <- which(mirs.log.counts$mir_id %in% mirs$genes[mirs.de.idx])
mirs.log.counts <- mirs.log.counts[mirs.log.counts.idx,]

#*******************************************************************************
#---- calc log2cpm per sample for mir,gene ----
#*******************************************************************************

genes.log2cnt.M <- genes.log.counts[,c("TM_T_1","TM_T_2","TM_T_4","TM_T_5",
                                       "TM_C_1","TM_C_2","TM_C_4","TM_C_5")]
colnames(genes.log2cnt.M) <- c("LGP3_T","LGP5_T","X021_T","X029_T",
                               "LGP3_C","LGP5_C","X021_C","X029_C")
rownames(genes.log2cnt.M) <- genes.log.counts$gene_id

genes.log2FC.M  <- genes.log.counts[,c("TM_T_1","TM_T_2","TM_T_4","TM_T_5")] -
                   genes.log.counts[,c("TM_C_1","TM_C_2","TM_C_4","TM_C_5")]
colnames(genes.log2FC.M) <- c("LGP3","LGP5","X021","X029")
rownames(genes.log2FC.M) <- genes.log.counts$gene_id
#genes.log2FC <- data.frame(genes.log.counts[,1:2],genes.log2FC.M)

mirs.log2cnt.M <- mirs.log.counts[,c("LGP3_TGFb1","LGP5_TGFb1","X021_TGFb1","X029_TGFb1",
                                     "LGP3_control","LGP5_control","X021_control","X029_control")]
colnames(mirs.log2cnt.M) <- c("LGP3_T","LGP5_T","X021_T","X029_T",
                              "LGP3_C","LGP5_C","X021_C","X029_C")
rownames(mirs.log2cnt.M)  <- mirs.log.counts$mir_id

mirs.log2FC.M <- mirs.log.counts[,c("LGP3_TGFb1","LGP5_TGFb1","X021_TGFb1","X029_TGFb1")] -
                 mirs.log.counts[,c("LGP3_control","LGP5_control","X021_control","X029_control")]
colnames(mirs.log2FC.M)  <- c("LGP3","LGP5","X021","X029")
rownames(mirs.log2FC.M)  <- mirs.log.counts$mir_id
#mirs.log2FC <- data.frame(mirs.log.counts[1],mirs.log2FC.M)

#*******************************************************************************
#---- determine correlation matrices for log2cnt and log2FC ----
#*******************************************************************************

# correlation based on FC
C  <- cor(t(genes.log2FC.M),t(mirs.log2FC.M))
# correlation based on sample expression (log2CPM)
C2 <- cor(t(genes.log2cnt.M),t(mirs.log2cnt.M))

# test correlation for specific pairs
#gene.cor = ""

#*******************************************************************************
#---- create dataframe iact : possible interactions from databases ----
#*******************************************************************************

iact.db = data.frame()
for (database.file in database.files) {
  cat("load ",database.file,"\n")
  load(database.file)
  iact.db <- rbind(iact.db,dbt)
}

#load(database.file)
#iact.db <- dbt
print(sprintf("loaded interaction database with %d interactions",nrow(iact.db)))

# select interactions containing our mirs
iact.db.sel.idx1 <- which(iact.db$miRNA %in% mirs.de$genes)
iact             <- iact.db[iact.db.sel.idx1,]
iact.db.sel.idx2 <- which(iact$gene.ID %in% genes.de$gene_id)
iact             <- iact[iact.db.sel.idx2,]

#---- add mirdip results to this database ----

cat("selected ",nrow(iact),"interactions from integrated database using DE genes,mirs\n")
load("mirdippie.Rdata")
cat("loaded ",nrow(iact.mirdip),"interactions from mirdip\n")

cat("adding mirdip to iact object. iact starts with",nrow(iact)," rows\n")
iact$db.mirdip   = rep(F,nrow(iact))
iact$qval.mirdip = rep(0,nrow(iact))
n.iact = nrow(iact)
for(k in 1:nrow(iact.mirdip)){
  iact.idx = which(iact$miRNA==iact.mirdip$mir[k] & iact$gene.Symbol==iact.mirdip$gene_symbol[k])
  if(length(iact.idx)==1) {
    iact$db.mirdip[iact.idx] = T
    iact$qval.mirdip[iact.idx] = 1-iact.mirdip$conf_score[k]
    # adjust qvalue
    iact$qval[iact.idx] = max(iact$qval[iact.idx], iact$qval.mirdip[iact.idx])
    iact$pval[iact.idx] = 1 - iact$qval[iact.idx]
  }
  if(length(iact.idx)>1) {
    print("something is really wrong")
    pipo654
  }
  if(length(iact.idx)==0 & F) { # always false
    # new entry : just skip because only mirdip is a bit risky perhaps
    iact[n.iact+1,] = iact[n.iact,]
    n.iact = n.iact + 1
    iact$miRNA[n.iact] = iact.mirdip$mir[k]
    iact$gene.ID[n.iact] = genes$gene_id[which(genes$gene_symbol==iact.mirdip$gene_symbol[k])]
    iact$gene.Symbol[n.iact] = iact.mirdip$gene_symbol[k]
    iact$db.miRTarBase_hsa_MTI[n.iact]              = F
    iact$db.TargetScan_DefaultPredictions[n.iact]   = F
    iact$db.TargetScan_Nonconserved[n.iact]         = F
    iact$db.TarBase[n.iact]                         = F
    iact$db.miRWalk[n.iact]                         = F
    iact$db.mirdip[n.iact]                          = T
    
    iact$qval.miRTarBase_hsa_MTI[n.iact]            = 0
    iact$qval.TarBase[n.iact]                       = 0
    iact$qval.TargetScan_DefaultPredictions[n.iact] = 0
    iact$qval.TargetScan_Nonconserved[n.iact]       = 0
    iact$qval.miRWalk[n.iact]                       = 0
    iact$qval.mirdip[n.iact] = 1-iact.mirdip$conf_score[k]
    
    iact$pval[n.iact] = iact.mirdip$conf_score[k]
    iact$qval[n.iact] = 1-iact$pval[n.iact]
  }
}

# remove iacts with total pval > iact.db.pval.lim
iact.db.sel.idx3 <- which(iact$pval <= iact.db.pval.lim)
iact             <- iact[iact.db.sel.idx3,]
cat("adding mirdip to iact object. iact ends with",nrow(iact)," rows\n")

# Check : are all genes in the database ?
genes.de.in.db.idx <- (genes.de$gene_id %in% iact.db$gene.ID)
print(sprintf("%d DE genes out of %d found in interactions using ENSG ID",
              length(which(genes.de.in.db.idx)),
              nrow(genes.de)))
#print(paste0("not found in database :",genes.de$gene_symbol[!genes.de.in.db.idx]))
# Check : are all mirs in the database ?
mirs.de.in.db.idx <- (mirs.de$genes %in% iact.db$miRNA)
print(sprintf("%d DE mirs out of %d found in interactions using miR id",
              length(which(mirs.de.in.db.idx)),
              nrow(mirs.de)))
#print(paste0("not found in database :",mirs.de$genes[!mirs.de.in.db.idx]))
print(sprintf("number of genes in db=%d, number of miRs in db=%d",
              length(unique(iact.db$gene.ID)),
              length(unique(iact.db$miRNA)) ))


# add gene, mir, cor info per interaction
iact$gene.log2FC  <- rep(0,nrow(iact))
iact$gene.log2cpm <- rep(0,nrow(iact))
iact$gene.pval    <- rep(0,nrow(iact))
iact$gene.padj    <- rep(0,nrow(iact))
iact$mir.log2FC   <- rep(0,nrow(iact))
iact$mir.log2cpm  <- rep(0,nrow(iact))
iact$mir.pval     <- rep(0,nrow(iact))
iact$mir.padj     <- rep(0,nrow(iact))
iact$cor.fc       <- rep(0,nrow(iact))
iact$cor.sam      <- rep(0,nrow(iact))

print("add gene, mir, correlation fields to iact object, this may take a while")

for (k in 1:nrow(iact)) {
  gene.idx <- which(iact$gene.ID[k] == genes.de$gene_id)
  iact$gene.log2FC[k]  <- genes.de$log2FC[gene.idx]
  iact$gene.log2cpm[k] <- genes.de$basemean[gene.idx]
  iact$gene.pval[k]   <- genes.de$pval[gene.idx]
  iact$gene.padj[k]   <- genes.de$padj[gene.idx]
  
  mir.idx  <- which(iact$miRNA[k]   == mirs.de$genes)
  iact$mir.log2FC[k]  <- mirs.de$logFC[mir.idx]
  iact$mir.log2cpm[k] <- mirs.de$logCPM[mir.idx]
  iact$mir.pval[k]    <- mirs.de$PValue[mir.idx]
  iact$mir.padj[k]    <- mirs.de$FDR[mir.idx]
  
  iact$cor.fc[k]  <-  C[iact$gene.ID[k],iact$miRNA[k]]
  iact$cor.sam[k] <- C2[iact$gene.ID[k],iact$miRNA[k]]  
  
}
print("done adding gene, mir, correlation fields to iact object")

#*******************************************************************************
#---- filter on correlation and inv FC : iact -> iact.sel ----
#*******************************************************************************

print(sprintf("-- starting with %d interactions",nrow(iact)))
# filter on correlation
# what pairs have neg correlation but no inv DE ?
iact$invFC <- sign(iact$gene.log2FC)*sign(iact$mir.log2FC) == -1

# size of subsets
invFC.cor.TT.idx <- which(iact$invFC==T &  iact$cor.sam<cor.lim)
invFC.cor.TF.idx <- which(iact$invFC==T & !iact$cor.sam<cor.lim)
invFC.cor.FT.idx <- which(iact$invFC==F &  iact$cor.sam<cor.lim)
invFC.cor.FF.idx <- which(iact$invFC==F & !iact$cor.sam<cor.lim)

iact.venn.m <- matrix(nrow=3,ncol=3)
iact.venn.m[1,1] = length(invFC.cor.TT.idx)
iact.venn.m[1,2] = length(invFC.cor.TF.idx)
iact.venn.m[2,1] = length(invFC.cor.FT.idx)
iact.venn.m[2,2] = length(invFC.cor.FF.idx)
iact.venn.m[1:2,3] = rowSums(iact.venn.m[1:2,1:2])
iact.venn.m[3,1:2] = colSums(iact.venn.m[1:2,1:2])
iact.venn.m[3,3]   = sum(iact.venn.m[1:2,3])
rownames(iact.venn.m) = c("invFC=T","invFC=F","sum")
colnames(iact.venn.m) = c("cor<lim","cor>lim","sum")
View(iact.venn.m)

# check known interactions : create df of known interactions + ref
iact.sel <- iact[invFC.cor.TT.idx,]

#*******************************************************************************
#---- construct binary database presence word  ----
#*******************************************************************************

# per entry determine word of the form 00100 to assess database presence

mirtarbase.b         = 0+(iact.sel$qval.miRTarBase_hsa_MTI>0)
targetscan.default.b = 0+(iact.sel$qval.TargetScan_DefaultPredictions>0)
targetscan.noncons.b = 0+(iact.sel$qval.TargetScan_Nonconserved>0)
tarbase.b            = 0+(iact.sel$qval.TarBase>0)
mirwalk.b            = 0+(iact.sel$qval.miRWalk>0)
mirdip.b             = 0+(iact.sel$qval.mirdip>0)

db.word = paste0(mirtarbase.b,targetscan.default.b,targetscan.noncons.b,
                 tarbase.b,mirwalk.b,mirdip.b)
db.sum  = mirtarbase.b+targetscan.default.b+targetscan.noncons.b+tarbase.b+
          mirwalk.b+mirdip.b


#*******************************************************************************
#---- reorder iact.sel ----
#*******************************************************************************
iact.sel.extended = iact.sel
rm(iact.sel)
iact.sel <- data.frame(  miRNA=iact.sel.extended$miRNA,
                         gene.ID=iact.sel.extended$gene.ID,
                         gene.Symbol=iact.sel.extended$gene.Symbol,
                         invFC=iact.sel.extended$invFC,
                         cor.fc=iact.sel.extended$cor.fc,
                         cor.sam=iact.sel.extended$cor.sam,
                         q=iact.sel.extended$qval,
                         db=db.word,
                         db.sum=db.sum,
                         q.mirtarbase=iact.sel.extended$qval.miRTarBase_hsa_MTI,
                         q.targetscan.default=iact.sel.extended$qval.TargetScan_DefaultPredictions,
                         q.targetscan.noncons=iact.sel.extended$qval.TargetScan_Nonconserved,
                         q.tarbase=iact.sel.extended$qval.TarBase,
                         q.mirwalk=iact.sel.extended$qval.miRWalk,
                         q.mirdip=iact.sel.extended$qval.mirdip,
                         gene.log2FC=iact.sel.extended$gene.log2FC,
                         gene.log2cpm=iact.sel.extended$gene.log2cpm,
                         gene.pval=iact.sel.extended$gene.pval,
                         gene.padj=iact.sel.extended$gene.padj,
                         mir.log2FC=iact.sel.extended$mir.log2FC,
                         mir.log2cpm=iact.sel.extended$mir.log2cpm,
                         mir.pval=iact.sel.extended$mir.pval,
                         mir.padj=iact.sel.extended$mir.padj)
View(iact.sel)

#*******************************************************************************
#---- filter out if only in 1 predictive database mirwalk, targetscan, mirdip
#*******************************************************************************
#iact.sel <- iact.sel$db==""
cat("Filter out interactions which originate from only 1 predictive database",
    " starting with",nrow(iact.sel),"interactions\n")
iact.sel <- iact.sel[!(iact.sel$q.mirwalk>0 & iact.sel$db.sum==1),]
iact.sel <- iact.sel[!(iact.sel$q.targetscan.noncons>0 & iact.sel$db.sum==1),]
iact.sel <- iact.sel[!(iact.sel$q.targetscan.default>0 & iact.sel$db.sum==1),]
iact.sel <- iact.sel[!(iact.sel$q.mirdip>0 & iact.sel$db.sum==1),]
cat("Filter out interactions which only occur in tarbase",
    " starting with",nrow(iact.sel),"interactions\n")
iact.sel <- iact.sel[!(iact.sel$q.tarbase>0 & iact.sel$db.sum==1),]
cat("Remaining :" ,nrow(iact.sel),"interactions\n")

#*******************************************************************************
#---- do network analysis ----
#*******************************************************************************
if(network.analysis) {
  print(sprintf("iact.sel : %d/%d genes and %d/%d mirs in %d interactions",
                length(unique(iact.sel$gene.ID )),nrow(genes.de),
                length(unique(iact.sel$miRNA)),nrow(mirs.de),
                nrow(iact.sel)))
  max.nr.edges = nrow(genes.de)*nrow(mirs.de)
  cat("fully connected graph has ",max.nr.edges," edges so ",
      100*nrow(iact.sel)/max.nr.edges,"% connected\n")
}

#*******************************************************************************
#---- construct Karnaugh ----
#*******************************************************************************
if(F){
library(pheatmap)
Mkar = matrix(0,8,8)
# 000 001 010 100 011 101 110 111
#   0   1   2   4   3   5   6   7
bin.order = c(0,1,2,4,3,5,6,7) + 1

for (k in 0:63) {
  msb = floor(k/8)
  lsb = k - 8*msb
  r = bin.order[msb+1]
  c = bin.order[lsb+1]
  db.dec = strtoi(iact.sel$db,2)
  Mkar[r,c] = length(which(db.dec==k))
}

# Karnaugh-type 

kar.targetscan_noncons = c(F,T,F,F,T,T,F,T)
kar.targetscan_default = c(F,F,T,F,T,F,T,T)
kar.mirtarbase         = c(F,F,F,T,F,T,T,T)
kar.mirdip             = c(F,T,F,F,T,T,F,T)
kar.mirwalk            = c(F,F,T,F,T,F,T,T)
kar.tarbase            = c(F,F,F,T,F,T,T,T)

ann.targetscan_noncons = c()
ann.targetscan_noncons[ kar.targetscan_noncons]="y"
ann.targetscan_noncons[!kar.targetscan_noncons]="n"

ann.targetscan_default = c()
ann.targetscan_default[ kar.targetscan_default]="y"
ann.targetscan_default[!kar.targetscan_default]="n"

ann.mirtarbase = c()
ann.mirtarbase[ kar.mirtarbase]="y"
ann.mirtarbase[!kar.mirtarbase]="n"

ann.mirdip = c()
ann.mirdip[ kar.mirdip]="y"
ann.mirdip[!kar.mirdip]="n"

ann.mirwalk = c()
ann.mirwalk[ kar.mirwalk]="y"
ann.mirwalk[!kar.mirwalk]="n"

ann.tarbase = c()
ann.tarbase[ kar.tarbase]="y"
ann.tarbase[!kar.tarbase]="n"

#ann.gene.fdr.txt = c()
#ann.gene.fdr.txt[ann.gene.fdr] = "yes"
#ann.gene.fdr.txt[!ann.gene.fdr] = "no"
annotation.colours <- list()
annotation.colours$targetscan_noncons <- c(n="white",y="grey")
annotation.colours$targetscan_default <- c(n="white",y="grey")
annotation.colours$mirtarbase         <- c(n="white",y="grey")
annotation.colours$mirdip             <- c(n="white",y="grey")
annotation.colours$mirwalk            <- c(n="white",y="grey")
annotation.colours$tarbase            <- c(n="white",y="grey")

annotation.row <- data.frame(targetscan_noncons = ann.targetscan_noncons,
                             targetscan_default = ann.targetscan_default,
                             mirtarbase         = ann.mirtarbase
                             )
annotation.col <- data.frame(mirdip             = ann.mirdip,
                             mirwalk            = ann.mirwalk,
                             tarbase            = ann.tarbase
)

rownames(Mkar) = c("none","targetscan.noncons","targetscan.default","mirtarbase",
                   "targetscan.all","mirtarbase+targetscan.noncons",
                   "mirtarbase+targetscan.default","mirtarbase+targetscan.all")
colnames(Mkar) = c("none","mirdip","mirwalk","tarbase","mirwalk+mirdip",
                   "tarbase+mirdip","tarbase+mirwalk","tarbase+mirwalk+mirdip")

rownames(annotation.row) <- rownames(Mkar)
rownames(annotation.col) <- colnames(Mkar)

pheatmap(mat=Mkar,cluster_rows = F,cluster_cols = F,display_numbers = T,
         annotation_row = annotation.row,
         annotation_col = annotation.col,
         annotation_colors = annotation.colours,
         show_rownames = F,show_colnames = F
)
View(Mkar)
}

#*******************************************************************************
#---- compare results with published miR-gene pairs in glaucoma ----
#*******************************************************************************

if(test.results) {
  publ.mirs.file <- "published_PEXG_miRs.xlsx"
  publ.mirs      <- as.data.frame(readxl::read_xlsx( publ.mirs.file,sheet='ref_xcheck'))
  testset.idx <- which(publ.mirs$testset=="x")
  testset <- data.frame(mir=publ.mirs$miR[testset.idx],
                        gen=publ.mirs$targ_gene[testset.idx],
                        outcome="?",
                        in.iact.sel="",                        
                        in.genes="",
                        in.mirs="",
                        in.genes.de="",
                        in.mirs.de="",                        
                        in.iact.db="",
                        invFC="",
                        corr="",cor.sam="")
  
  for (k in 1:nrow(testset)){
    # test if in iact.sel
    
    mir.ref=testset$mir[k]
    gen.ref=testset$gen[k]
  
    # interrogation flow 
    q0 = length(which(mir.ref==iact.sel$miRNA & gen.ref==iact.sel$gene.Symbol))>0
    q1 = length(which(mir.ref==mirs$genes))>0
    q2 = length(which(gen.ref==genes$gene_symbol))>0
    if(q1 & q2) {
      # if gene and mir exist, do other tests : 
      q3 = mir.ref %in% mirs.de$genes
      q4 = gen.ref %in% genes.de$gene_symbol
      q5 = length(which(iact.db$gene.Symbol==gen.ref & iact.db$miRNA==mir.ref))>0
      q6 = ( sign(genes$log2FC[which(gen.ref==genes$gene_symbol)])*
             sign(mirs$logFC[which(mir.ref==mirs$genes)]) ) < 0
      if(q3 & q4 & q5 & q6) {
        iact.idx = which(iact$miRNA==mir.ref & iact$gene.Symbol==gen.ref)
        testset$cor.sam[k] = iact$cor.sam[iact.idx]
        q7 = (iact$cor.sam[iact.idx] < cor.lim)
      } else {
        q7 = NA
      }
    } else {
      q3=NA;q4=NA;q5=NA;q6 = NA;q7=NA
    }

    if(q0) {
      testset$outcome[k] ="OK:in sel"
    } else if(!q1) {
      testset$outcome[k] ="miR not known"
    } else if(!q2) {
      testset$outcome[k] ="gene not known"
    } else if(!q3) {
      testset$outcome[k] ="miR is not DE"
    } else if(!q4) {
      testset$outcome[k] ="gene is not DE"
    } else if(!q5) {
      testset$outcome[k] ="interaction not in database"
    } else if(!q6) {
      testset$outcome[k] ="FC not inverse"
    } else if(!q7) {
      testset$outcome[k] ="insufficient correlation"
    }
  
    
    
    # assign answers to classify according to pattern

    testset$in.iact.sel[k] = q0
    testset$in.mirs[k]     = q1
    testset$in.genes[k]    = q2
    testset$in.mirs.de[k]  = q3
    testset$in.genes.de[k] = q4
    testset$in.iact.db[k]  = q5
    testset$invFC[k]       = q6
    testset$corr[k]        = q7


  }
  View(testset)
}

#---- calculate pval and padj per interaction ----
iact.sel$iact.pval = iact.sel$gene.pval + iact.sel$mir.pval
iact.sel$iact.padj = p.adjust(iact.sel$iact.pval,method = "fdr")
iact.sel$gene.and.mir.de = (iact.sel$gene.padj<0.05 & iact.sel$mir.padj<0.05)


#---- save interatome ----
if(save.iact.sel){
  print(paste0("saving iact.sel to ",iact.sel.file))
  save(iact.sel,file=paste0(iact.sel.file,".Rdata"))
  write.table(iact.sel,file=paste0(iact.sel.file,".csv"),quote=F,sep="\t")
}

