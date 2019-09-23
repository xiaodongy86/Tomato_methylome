library(MethylIT)
library(DSS)
library(dplyr)

.libPaths <- c("/usr/share/R/library","/usr/lib64/R/library","/data/users/xzy50/R/x86_64-redhat-linux-gnu-library/3.6")

setwd("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor")


####################################################      WT gen1 vs gen2   ############################################################

#################                                              CG                                                 #####################

##### 1-P1 1-p2 WT gen1             2-P1  2-P2   gen2 

tomato_bismark_wt_gen1_p1_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/1-P1/1-P1_CG.txt")
tomato_bismark_wt_gen1_p2_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/1-P2/1-P2_CG.txt")

tomato_bismark_wt_gen2_p1_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/2-P1/2-P1_CG.txt")
tomato_bismark_wt_gen2_p2_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/2-P2/2-P2_CG.txt")

# bismark format
# chromosome number, genomic coordinate, strand, methylated count, unmethylated count, c-context , trinucleotide context
head(tomato_bismark_wt_gen1_p1_CG)

# DSS wants is:
# chromosome number, genomic coordinate, total number of reads, and number of reads showing methylation
# and colnames has to be 
# chr pos N X

colnames(tomato_bismark_wt_gen1_p1_CG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_wt_gen1_p2_CG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_wt_gen2_p1_CG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_wt_gen2_p2_CG) <- c("chr","pos","stand","X","uC","context","tri_context")

head(tomato_bismark_wt_gen1_p1_CG)

tomato_bismark_wt_gen1_p1_CG_A <- tomato_bismark_wt_gen1_p1_CG %>% mutate(N = X + uC)
tomato_bismark_wt_gen1_p1_CG_B <- tomato_bismark_wt_gen1_p1_CG_A[,c(1,2,8,4)]

tomato_bismark_wt_gen1_p2_CG_A <- tomato_bismark_wt_gen1_p2_CG %>% mutate(N = X + uC)
tomato_bismark_wt_gen1_p2_CG_B <- tomato_bismark_wt_gen1_p2_CG_A[,c(1,2,8,4)]

tomato_bismark_wt_gen2_p1_CG_A <- tomato_bismark_wt_gen2_p1_CG %>% mutate(N = X + uC)
tomato_bismark_wt_gen2_p1_CG_B <- tomato_bismark_wt_gen2_p1_CG_A[,c(1,2,8,4)]

tomato_bismark_wt_gen2_p2_CG_A <- tomato_bismark_wt_gen2_p2_CG %>% mutate(N = X + uC)
tomato_bismark_wt_gen2_p2_CG_B <- tomato_bismark_wt_gen2_p2_CG_A[,c(1,2,8,4)]

head(tomato_bismark_wt_gen2_p2_CG_B)

save(tomato_bismark_wt_gen1_p1_CG_B, tomato_bismark_wt_gen1_p2_CG_B, 
     file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/bismark_wt_gen1_CG_ready_for_DSS_9_23_2019.RData")

save(tomato_bismark_wt_gen2_p1_CG_B, tomato_bismark_wt_gen2_p2_CG_B, 
     file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/bismark_wt_gen2_CG_ready_for_DSS_9_23_2019.RData")



###############################################################################################################################################

####                                                             BS oject 

###############################################################################################################################################


##################        CG        ##########################
BSobj_CG_wt_gen1_vs_wt_gen2 <- makeBSseqData(list(tomato_bismark_wt_gen1_p1_CG_B,
                                                  tomato_bismark_wt_gen1_p2_CG_B,
                                                  tomato_bismark_wt_gen2_p1_CG_B,
                                                  tomato_bismark_wt_gen2_p2_CG_B), c("C1","C2","N1","N2"))
save(BSobj_CG_wt_gen1_vs_wt_gen2, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/BSobj_CG_wt_gen1_vs_wt_gen2_9_23_2019.RData")

dmlTest_CG_wt_gen1_vs_wt_gen2 <- DMLtest(BSobj_CG_wt_gen1_vs_wt_gen2, 
                                         group1=c("C1","C2"), group2=c("N1","N2"), smoothing=TRUE,smoothing.span=100)

save(dmlTest_CG_wt_gen1_vs_wt_gen2, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/dmlTest_CG_wt_gen1_vs_wt_gen2_9_23_2019.RData")




