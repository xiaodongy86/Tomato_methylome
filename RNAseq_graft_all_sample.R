## ##################################################################### ##
############################# ---RNA-Seq --- ##############################
## ##################################################################### ##
#
#### **** ---RNA_Seq data. RNAi NULL & WT ------**** ####
#
################################################ ##

library(DESeq2)
library(genefilter)
setwd("/data5/tomato_RNAseq/QoRTs_out")

sampleFiles <- list.files( pattern = ".geneCounts.formatted.for.DESeq.txt.gz" )
sampleName = sub("_QC.geneCounts.formatted.for.DESeq.txt.gz","", sampleFiles )

sampleCondition <- factor( c( "RDr", "RDr", "RDr",
                              "RRCtl", "RRCtl", "RRCtl",
                              "RutgersWt","RutgersWt","RutgersWt"),
                           levels = c("RDr", "RRCtl", "RutgersWt"))

sampleTable <- data.frame( sampleName = sampleName,
                           fileName = sampleFiles,
                           condition = sampleCondition )

#       sampleName                      fileName                     condition
# 1       RDrP1       RDrP1_QC.geneCounts.formatted.for.DESeq.txt.gz       RDr
# 2       RDrP2       RDrP2_QC.geneCounts.formatted.for.DESeq.txt.gz       RDr
# 3       RDrP3       RDrP3_QC.geneCounts.formatted.for.DESeq.txt.gz       RDr
# 4     RRCtlP1     RRCtlP1_QC.geneCounts.formatted.for.DESeq.txt.gz     RRCtl
# 5     RRCtlP2     RRCtlP2_QC.geneCounts.formatted.for.DESeq.txt.gz     RRCtl
# 6     RRCtlP3     RRCtlP3_QC.geneCounts.formatted.for.DESeq.txt.gz     RRCtl
# 7 RutgersWtP1 RutgersWtP1_QC.geneCounts.formatted.for.DESeq.txt.gz RutgersWt
# 8 RutgersWtP2 RutgersWtP2_QC.geneCounts.formatted.for.DESeq.txt.gz RutgersWt
# 9 RutgersWtP3 RutgersWtP3_QC.geneCounts.formatted.for.DESeq.txt.gz RutgersWt

deg.rnai <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable,
                                        directory = "./",
                                        design= ~ condition)
# Pre-filtering
# deg.rnai <- deg.rnai[ rowSums( counts( deg.rnai ) ) > 1, ]

######################################## # pheatmap
#### Differential expression analysis ####
######################################## # SummarizedExperiment
deg.rnai <- estimateSizeFactors( deg.rnai, type = "iterate", locfunc = shorth )

# deg.rnai <- estimateSizeFactors( deg.rnai )
deg.rnai <- estimateDispersions( deg.rnai, fitType = "local", maxit = 10^8 )
idx <- rowSums( counts( deg.rnai, normalized = TRUE ) >= 10 ) >= 3
deg.rnai = deg.rnai[ idx ]

deg.rnai <- nbinomWaldTest( deg.rnai, betaPrior = TRUE, useOptim = FALSE, maxit = 10^8 )

# plotDispEsts( deg.rnai )

# Extract results from a DESeq analysis with filter rv
rv <- rowVars( counts( deg.rnai, normalized = TRUE ) ) # row variance

DEG.MM <- results( deg.rnai, filter = rv , pAdjustMethod = "BH" ) #

DEG.MM = DEG.MM[ !is.na( DEG.MM$padj ), ]
DEG.MM = DEG.MM[ abs(DEG.MM$log2FoldChange) > 1, ] # 
DEG.MM$adj.pval =  p.adjust( DEG.MM$pvalue, method = "BH" )
DEG.MM = DEG.MM[ DEG.MM$adj.pval <= 0.05, ] # 1368
DEG.MM

# sum( abs(DEG.MM$log2FoldChange) > 1, na.rm = T ) # 157
# sum( abs(DEG.MM$log2FoldChange) > 0.5, na.rm = T ) # 864
# min( abs(DEG.MM$log2FoldChange) ) # 0.2818615
# summary(abs(DEG.MM$log2FoldChange) )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2819  0.4417  0.5772  0.6547  0.7685  3.0314 

save( deg.rnai, DEG.MM, file = "/JOB/Work/EpisTDNA/RNASeq/RNAiNull/DEGs/DEGs_WT_vs_MM_10-31-17.RData" )
write.csv( DEG.MM, file = "/JOB/Work/EpisTDNA/RNASeq/RNAiNull/DEGs/DEGs_WT_vs_MM_10-31-17.csv" )


