#############################
# most packages do not support the use of TPM or FPKM for differential expression testing. This means that e.g. 
# it's completely wrong to feed them to programs expecting counts (e.g. DESeq2 or EdgeR). One reason for this is that these measures are normalized.
# What I mean by this is that, for example, if you sum the TPM of all genes/transcripts in a sample, the sum will always be 1,000,000; this is a direct result of the way TPM is calculated. 
# FPKM will behave in a similar manner (though, as many have pointed out, TPM should always be preferred to FPKM which has a more arbitrary and less stable normalization term).


## ##################################################################### ##
############################# ---RNA-Seq --- ##############################
## ##################################################################### ##
#
#### **** ---RNA_Seq data. RNAi NULL & WT ------**** ####
#
################################################ ##
library(edgeR)
library(DESeq2)
library(genefilter)
library(MethylIT)
# devtools::install_git("https://github.com/genomaths/MethylIT.git")

setwd("/data5/tomato_RNAseq/QoRTs_out")

sampleFiles <- c("RRCtlP1_QC.geneCounts.formatted.for.DESeq.txt.gz",
                 "RRCtlP2_QC.geneCounts.formatted.for.DESeq.txt.gz",
                 "RRCtlP3_QC.geneCounts.formatted.for.DESeq.txt.gz",
                 "RDrP1_QC.geneCounts.formatted.for.DESeq.txt.gz",
                 "RDrP2_QC.geneCounts.formatted.for.DESeq.txt.gz",   
                 "RDrP3_QC.geneCounts.formatted.for.DESeq.txt.gz")

sampleName <- c("RRCtlP1","RRCtlP2","RRCtlP3","RDrP1","RDrP2","RDrP3")


sampleCondition <- factor( c( "RRCtl", "RRCtl", "RRCtl",
                              "RDr", "RDr", "RDr"),
                           levels = c("RRCtl", "RDr"))

sampleTable <- data.frame( sampleName = sampleName,
                           fileName = sampleFiles,
                           condition = sampleCondition )
# sampleName                                         fileName condition
# 1    RRCtlP1 RRCtlP1_QC.geneCounts.formatted.for.DESeq.txt.gz     RRCtl
# 2    RRCtlP2 RRCtlP2_QC.geneCounts.formatted.for.DESeq.txt.gz     RRCtl
# 3    RRCtlP3 RRCtlP3_QC.geneCounts.formatted.for.DESeq.txt.gz     RRCtl
# 4      RDrP1   RDrP1_QC.geneCounts.formatted.for.DESeq.txt.gz       RDr
# 5      RDrP2   RDrP2_QC.geneCounts.formatted.for.DESeq.txt.gz       RDr
# 6      RDrP3   RDrP3_QC.geneCounts.formatted.for.DESeq.txt.gz       RDr


deg.rnai <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable,
                                        directory = "./",
                                        design= ~ condition)

# colData(deg.rnai)
# rowData(deg.rnai)
rcounts <- assay(deg.rnai) 

nams <- rownames(rcounts)
idx <- grep("ENSRNA", nams)
deg.rnai <- deg.rnai[ -idx ]

# To remove genes with constant counts
rcounts <- assay(deg.rnai)   # update read counts
vars <- apply(rcounts, 1, var)
deg.rnai <- deg.rnai[which(vars > 0)]
rcounts <- assay(deg.rnai) # update read counts

factors <- calcNormFactors(object = rcounts)
rcounts <- t(t(rcounts) * factors)

head(rcounts)
## An experiment design is set.
colData <- colData(deg.rnai)
rcounts <-  round(rcounts)


Tomato_ITAG3.0_gene_models <- import("/data/users/xzy50/ITAG3.0_gene_models.gff")
head(Tomato_ITAG3.0_gene_models)
seqlevels(Tomato_ITAG3.0_gene_models)
gene = Tomato_ITAG3.0_gene_models[ Tomato_ITAG3.0_gene_models$type == "gene", c("ID","Alias","Note","Ontology_term")]
#seqlevels( gene ) <- c( "0","1","2","3","4","5","6","7","8","9","10","11","12" )
#seqlevels(gene, pruning.mode = "coarse") <- c("0","1","2","3","4","5","6","7","8","9","10","11","12" )
gene = sortBySeqnameAndStart(gene)

gene$ID <- substr(gene$ID, start =6 ,stop =22)
rcounts_gene <- data.frame(gene[match(row.names(rcounts),gene$ID)], name=rcounts)

gene_gf <- gene[match(row.names(rcounts),gene$ID)]
mcols(gene_gf) <- NULL

rcounts_gene_gr <- makeGRangesFromDataFrame(rcounts_gene, keep.extra.columns=TRUE,
                         start.field="start",
                         seqnames.field="seqnames",
                         end.field=c("end"),
                         strand.field="strand")

rcounts_gene_gr_2 <- rcounts_gene_gr[,5:10]
colnames(mcols(rcounts_gene_gr_2)) <- c("RRCtlP1","RRCtlP2","RRCtlP3","RDrP1","RDrP2","RDrP3")

grlist=list()
for(item in 1:26392){grlist[item]=GRanges(gene_gf)[item]}

names(grlist) <- gene_gf$ID
names(rowRanges(deg.rnai))
grlist_GRangesList <- GRangesList(grlist)
rowRanges(deg.rnai) <- grlist_GRangesList

# the FPKM values
fpkm(deg.rnai)



# DEGs = countTest2(DS = deg.rnai, num.cores = 4L, minCountPerIndv = 3, 
#                    countFilter = T, maxGrpCV = c(0.5, 0.5),
#                    FilterLog2FC = TRUE,Minlog2FC = 0.5, pvalCutOff = 0.05,
#                    MVrate = .95, verbose = FALSE, test = "Wald")

# parameters for 279 DEGs
DEGs = countTest2(DS = ds, num.cores = 24L, minCountPerIndv = 10, 
                  countFilter = T, maxGrpCV = c(1, 1),
                  FilterLog2FC = TRUE, Minlog2FC = 1, pvalCutOff = 0.05,
                  MVrate = .95, verbose = FALSE, test = "LRT")

DEGs
