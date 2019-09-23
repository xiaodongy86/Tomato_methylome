
library(MethylIT)
library(DSS)
library(dplyr)
setwd("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor")

# 12-P1	Solanum lycopersicum	grafted Gen1 R/DR
# 12-P2	Solanum lycopersicum	grafted Gen1 R/DR
# 12-P3	Solanum lycopersicum	grafted Gen1 R/DR
# 13-P1	Solanum lycopersicum	grafted Gen1 R/R
# 13-P2	Solanum lycopersicum	grafted Gen1 R/R
# 13-P3	Solanum lycopersicum	grafted Gen1 R/R


tomato_bismark_wt_RDR_p1_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/12-P1/12-P1_CG.txt")
tomato_bismark_wt_RDR_p2_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/12-P2/12-P2_CG.txt")
tomato_bismark_wt_RDR_p3_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/12-P3/12-P3_CG.txt")


tomato_bismark_wt_RR_p1_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/13-P1/13-P1_CG.txt")
tomato_bismark_wt_RR_p2_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/13-P2/13-P2_CG.txt")
tomato_bismark_wt_RR_p3_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/13-P3/13-P3_CG.txt")


# bismark format
# chromosome number, genomic coordinate, strand, methylated count, unmethylated count, c-context , trinucleotide context
head(tomato_bismark_BE_P1_CG)

# DSS wants is:
# chromosome number, genomic coordinate, total number of reads, and number of reads showing methylation
# and colnames has to be 
# chr pos N X

colnames(tomato_bismark_BE_P1_CG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_BE_P2_CG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_BE_P3_CG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_GE_P1_CG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_GE_P2_CG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_GE_P3_CG) <- c("chr","pos","stand","X","uC","context","tri_context")









