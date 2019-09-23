library(MethylIT)
library(DSS)
library(dplyr)
setwd("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor")


####################################################      epiline   ############################################################

                          #################                 CG               #####################

##### BE-P1 BE-P2  BE-P3   GE-P1 GE-P2 GE-P3

tomato_bismark_BE_P1_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P1/BE-P1_CG.txt")
tomato_bismark_BE_P2_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P2/BE-P2_CG.txt")
tomato_bismark_BE_P3_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P3/BE-P3_CG.txt")

tomato_bismark_GE_P1_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P1/GE-P1_CG.txt")
tomato_bismark_GE_P2_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P2/GE-P2_CG.txt")
tomato_bismark_GE_P3_CG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P3/GE-P3_CG.txt")

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

head(tomato_bismark_BE_P1_CG)

tomato_bismark_BE_P1_CG_A <- tomato_bismark_BE_P1_CG %>% mutate(N = X + uC)
tomato_bismark_BE_P1_CG_B <- tomato_bismark_BE_P1_CG_A[,c(1,2,8,4)]

tomato_bismark_BE_P2_CG_A <- tomato_bismark_BE_P2_CG %>% mutate(N = X + uC)
tomato_bismark_BE_P2_CG_B <- tomato_bismark_BE_P2_CG_A[,c(1,2,8,4)]


tomato_bismark_BE_P3_CG_A <- tomato_bismark_BE_P3_CG %>% mutate(N = X + uC)
tomato_bismark_BE_P3_CG_B <- tomato_bismark_BE_P3_CG_A[,c(1,2,8,4)]

tomato_bismark_GE_P1_CG_A <- tomato_bismark_GE_P1_CG %>% mutate(N = X + uC)
tomato_bismark_GE_P1_CG_B <- tomato_bismark_GE_P1_CG_A[,c(1,2,8,4)]

tomato_bismark_GE_P2_CG_A <- tomato_bismark_GE_P2_CG %>% mutate(N = X + uC)
tomato_bismark_GE_P2_CG_B <- tomato_bismark_GE_P2_CG_A[,c(1,2,8,4)]

tomato_bismark_GE_P3_CG_A <- tomato_bismark_GE_P3_CG %>% mutate(N = X + uC)
tomato_bismark_GE_P3_CG_B <- tomato_bismark_GE_P3_CG_A[,c(1,2,8,4)]

save(tomato_bismark_BE_P1_CG_B, tomato_bismark_BE_P2_CG_B, tomato_bismark_BE_P3_CG_B, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/bismark_BE_CG_ready_for_DSS_5_6_2019.RData")
save(tomato_bismark_GE_P1_CG_B, tomato_bismark_GE_P2_CG_B, tomato_bismark_GE_P3_CG_B, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/bismark_GE_CG_ready_for_DSS_5_6_2019.RData")



                            #################                 CHG               #####################

##### BE-P1 BE-P2  BE-P3   GE-P1 GE-P2 GE-P3

tomato_bismark_BE_P1_CHG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P1/BE-P1_CHG.txt")
tomato_bismark_BE_P2_CHG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P2/BE-P2_CHG.txt")
tomato_bismark_BE_P3_CHG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P3/BE-P3_CHG.txt")

tomato_bismark_GE_P1_CHG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P1/GE-P1_CHG.txt")
tomato_bismark_GE_P2_CHG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P2/GE-P2_CHG.txt")
tomato_bismark_GE_P3_CHG <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P3/GE-P3_CHG.txt")

# bismark format
# chromosome number, genomic coordinate, strand, methylated count, unmethylated count, c-context , trinucleotide context
head(tomato_bismark_BE_P1_CHG)


# DSS wants is:
# chromosome number, genomic coordinate, total number of reads, and number of reads showing methylation
# and colnames has to be 
# chr pos N X


colnames(tomato_bismark_BE_P1_CHG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_BE_P2_CHG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_BE_P3_CHG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_GE_P1_CHG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_GE_P2_CHG) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_GE_P3_CHG) <- c("chr","pos","stand","X","uC","context","tri_context")

head(tomato_bismark_BE_P1_CHG)

tomato_bismark_BE_P1_CHG_A <- tomato_bismark_BE_P1_CHG %>% mutate(N = X + uC)
tomato_bismark_BE_P1_CHG_B <- tomato_bismark_BE_P1_CHG_A[,c(1,2,8,4)]

tomato_bismark_BE_P2_CHG_A <- tomato_bismark_BE_P2_CHG %>% mutate(N = X + uC)
tomato_bismark_BE_P2_CHG_B <- tomato_bismark_BE_P2_CHG_A[,c(1,2,8,4)]

tomato_bismark_BE_P3_CHG_A <- tomato_bismark_BE_P3_CHG %>% mutate(N = X + uC)
tomato_bismark_BE_P3_CHG_B <- tomato_bismark_BE_P3_CHG_A[,c(1,2,8,4)]

tomato_bismark_GE_P1_CHG_A <- tomato_bismark_GE_P1_CHG %>% mutate(N = X + uC)
tomato_bismark_GE_P1_CHG_B <- tomato_bismark_GE_P1_CHG_A[,c(1,2,8,4)]

tomato_bismark_GE_P2_CHG_A <- tomato_bismark_GE_P2_CHG %>% mutate(N = X + uC)
tomato_bismark_GE_P2_CHG_B <- tomato_bismark_GE_P2_CHG_A[,c(1,2,8,4)]

tomato_bismark_GE_P3_CHG_A <- tomato_bismark_GE_P3_CHG %>% mutate(N = X + uC)
tomato_bismark_GE_P3_CHG_B <- tomato_bismark_GE_P3_CHG_A[,c(1,2,8,4)]

save(tomato_bismark_BE_P1_CHG_B, tomato_bismark_BE_P2_CHG_B, tomato_bismark_BE_P3_CHG_B, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/bismark_BE_CHG_ready_for_DSS_5_6_2019.RData")
save(tomato_bismark_GE_P1_CHG_B, tomato_bismark_GE_P2_CHG_B, tomato_bismark_GE_P3_CHG_B, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/bismark_GE_CHG_ready_for_DSS_5_6_2019.RData")



                                     #################                 CHH               #####################

##### BE-P1 BE-P2  BE-P3   GE-P1 GE-P2 GE-P3

tomato_bismark_BE_P1_CHH <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P1/BE-P1_CHH.txt")
tomato_bismark_BE_P2_CHH <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P2/BE-P2_CHH.txt")
tomato_bismark_BE_P3_CHH <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P3/BE-P3_CHH.txt")

tomato_bismark_GE_P1_CHH <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P1/GE-P1_CHH.txt")
tomato_bismark_GE_P2_CHH <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P2/GE-P2_CHH.txt")
tomato_bismark_GE_P3_CHH <- read.delim("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P3/GE-P3_CHH.txt")

# bismark format
# chromosome number, genomic coordinate, strand, methylated count, unmethylated count, c-context , trinucleotide context
head(tomato_bismark_BE_P1_CHH)


# DSS wants is:
# chromosome number, genomic coordinate, total number of reads, and number of reads showing methylation
# and colnames has to be 
# chr pos N X


colnames(tomato_bismark_BE_P1_CHH) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_BE_P2_CHH) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_BE_P3_CHH) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_GE_P1_CHH) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_GE_P2_CHH) <- c("chr","pos","stand","X","uC","context","tri_context")
colnames(tomato_bismark_GE_P3_CHH) <- c("chr","pos","stand","X","uC","context","tri_context")

head(tomato_bismark_BE_P1_CHH)

tomato_bismark_BE_P1_CHH_A <- tomato_bismark_BE_P1_CHH %>% mutate(N = X + uC)
tomato_bismark_BE_P1_CHH_B <- tomato_bismark_BE_P1_CHH_A[,c(1,2,8,4)]

tomato_bismark_BE_P2_CHH_A <- tomato_bismark_BE_P2_CHH %>% mutate(N = X + uC)
tomato_bismark_BE_P2_CHH_B <- tomato_bismark_BE_P2_CHH_A[,c(1,2,8,4)]

tomato_bismark_BE_P3_CHH_A <- tomato_bismark_BE_P3_CHH %>% mutate(N = X + uC)
tomato_bismark_BE_P3_CHH_B <- tomato_bismark_BE_P3_CHH_A[,c(1,2,8,4)]

tomato_bismark_GE_P1_CHH_A <- tomato_bismark_GE_P1_CHH %>% mutate(N = X + uC)
tomato_bismark_GE_P1_CHH_B <- tomato_bismark_GE_P1_CHH_A[,c(1,2,8,4)]

tomato_bismark_GE_P2_CHH_A <- tomato_bismark_GE_P2_CHH %>% mutate(N = X + uC)
tomato_bismark_GE_P2_CHH_B <- tomato_bismark_GE_P2_CHH_A[,c(1,2,8,4)]

tomato_bismark_GE_P3_CHH_A <- tomato_bismark_GE_P3_CHH %>% mutate(N = X + uC)
tomato_bismark_GE_P3_CHH_B <- tomato_bismark_GE_P3_CHH_A[,c(1,2,8,4)]

save(tomato_bismark_BE_P1_CHH_B, tomato_bismark_BE_P2_CHH_B, tomato_bismark_BE_P3_CHH_B, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/bismark_BE_CHH_ready_for_DSS_5_6_2019.RData")
save(tomato_bismark_GE_P1_CHH_B, tomato_bismark_GE_P2_CHH_B, tomato_bismark_GE_P3_CHH_B, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/bismark_GE_CHH_ready_for_DSS_5_6_2019.RData")



###############################################################################################################################################
                            
####                                                             BS oject 

###############################################################################################################################################


                                          ##################        CG        ##########################
BSobj_CG_BE_vs_GE <- makeBSseqData(list(tomato_bismark_BE_P1_CG_B,
                                         tomato_bismark_BE_P2_CG_B,
                                         tomato_bismark_BE_P3_CG_B,
                                         tomato_bismark_GE_P1_CG_B,
                                         tomato_bismark_GE_P2_CG_B,
                                         tomato_bismark_GE_P3_CG_B),
                                    c("C1","C2","C3","N1","N2","N3"))
save(BSobj_CG_BE_vs_GE, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/BSobj_CG_BE_vs_GE_5_6_2019.RData")

dmlTest_CG_BE_vs_GE <- DMLtest(BSobj_CG_BE_vs_GE, 
                                group1=c("C1","C2","C3"), group2=c("N1","N2","N3"), smoothing=TRUE,smoothing.span=100)

save(dmlTest_CG_BE_vs_GE, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/dmlTest_CG_BE_vs_GE_5_6_2019.RData")

load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/dmlTest_CG_BE_vs_GE_5_6_2019.RData")

# test the DMR number for different cutoff 
cutOff_CG <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8)
for (i in cutOff_CG) {
dmrs_CG <- callDMR(dmlTest_CG_BE_vs_GE, delta= i, minlen=100, minCG=3, pct.sig=0.5, dis.merge=50, p.threshold=0.05)
print(paste("with cutoff",i,"The DMR number is",nrow(dmrs_CG)))
}

#  5_7_2019
# [1] "with cutoff 0.2 The DMR number is 676"
# [1] "with cutoff 0.3 The DMR number is 218"
# [1] "with cutoff 0.4 The DMR number is 59"
# [1] "with cutoff 0.5 The DMR number is 14"
# [1] "with cutoff 0.6 The DMR number is 4"
# no DMR for 0.7 and 0.8

### proceed with cutoff 0.6  for CG according to  The Plant Journal(2018)93,460-471
dmrs_CG_0.6_BE_GE <- callDMR(dmlTest_CHG_BE_vs_GE, delta= 0.6, minlen=100, minCG=3, pct.sig=0.5, dis.merge=50, p.threshold=0.05)
save(dmrs_CG_0.6_BE_GE, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/dmrs_CG_0.6_BE_GE_5_9_2019.RData")


                                         ##################        CHG        ##########################

BSobj_CHG_BE_vs_GE <- makeBSseqData(list(tomato_bismark_BE_P1_CHG_B,
                                        tomato_bismark_BE_P2_CHG_B,
                                        tomato_bismark_BE_P3_CHG_B,
                                        tomato_bismark_GE_P1_CHG_B,
                                        tomato_bismark_GE_P2_CHG_B,
                                        tomato_bismark_GE_P3_CHG_B),
                                   c("C1","C2","C3","N1","N2","N3"))
save(BSobj_CHG_BE_vs_GE, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/BSobj_CHG_BE_vs_GE_5_6_2019.RData")

dmlTest_CHG_BE_vs_GE <- DMLtest(BSobj_CHG_BE_vs_GE, 
                               group1=c("C1","C2","C3"), group2=c("N1","N2","N3"), smoothing=TRUE,smoothing.span=100)

save(dmlTest_CHG_BE_vs_GE, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/dmlTest_CHG_BE_vs_GE_5_6_2019.RData")

load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/dmlTest_CHG_BE_vs_GE_5_6_2019.RData")

# test the DMR number for different cutoff 
cutOff_CHG <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
for (i in cutOff_CHG) {
  dmrs_CHG <- callDMR(dmlTest_CHG_BE_vs_GE, delta= i, minlen=100, minCG=3, pct.sig=0.5, dis.merge=50, p.threshold=0.05)
  print(paste("with cutoff",i,"The DMR number is",nrow(dmrs_CHG)))
}

# [1] "with cutoff 0.1 The DMR number is 3498"
# [1] "with cutoff 0.2 The DMR number is 756"
# [1] "with cutoff 0.3 The DMR number is 216"
# [1] "with cutoff 0.4 The DMR number is 48"
# [1] "with cutoff 0.5 The DMR number is 11"
# [1] "with cutoff 0.6 The DMR number is 3"
# [1] "with cutoff 0.7 The DMR number is "
# [1] "with cutoff 0.8 The DMR number is "

### proceed with cutoff 0.4  for CHG according to  The Plant Journal(2018)93,460-471
dmrs_CHG_0.4_BE_GE <- callDMR(dmlTest_CHG_BE_vs_GE, delta= 0.4, minlen=100, minCG=3, pct.sig=0.5, dis.merge=50, p.threshold=0.05)
save(dmrs_CHG_0.4_BE_GE, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/dmrs_CHG_0.4_BE_GE_5_9_2019.RData")

                                      ##################        CHH        ##########################
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/bismark_BE_CHH_ready_for_DSS_5_6_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/bismark_GE_CHH_ready_for_DSS_5_6_2019.RData")

BSobj_CHH_BE_vs_GE <- makeBSseqData(list(tomato_bismark_BE_P1_CHH_B,
                                         tomato_bismark_BE_P2_CHH_B,
                                         tomato_bismark_BE_P3_CHH_B,
                                         tomato_bismark_GE_P1_CHH_B,
                                         tomato_bismark_GE_P2_CHH_B,
                                         tomato_bismark_GE_P3_CHH_B),
                                    c("C1","C2","C3","N1","N2","N3"))
save(BSobj_CHH_BE_vs_GE, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/BSobj_CHH_BE_vs_GE_5_6_2019.RData")

dmlTest_CHH_BE_vs_GE <- DMLtest(BSobj_CHH_BE_vs_GE, 
                                group1=c("C1","C2","C3"), group2=c("N1","N2","N3"), smoothing=TRUE,smoothing.span=100)

save(dmlTest_CHH_BE_vs_GE, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/dmlTest_CHH_BE_vs_GE_5_6_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/dmlTest_CHH_BE_vs_GE_5_6_2019.RData")

# test the DMR number for different cutoff 
cutOff_CHH <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
for (i in cutOff_CHH) {
  dmrs_CHH <- callDMR(dmlTest_CHH_BE_vs_GE, delta= i, minlen=100, minCG=3, pct.sig=0.5, dis.merge=50, p.threshold=0.05)
  print(paste("with cutoff",i,"The DMR number is",nrow(dmrs_CHH)))
}

#[1] "with cutoff 0.1 The DMR number is 2332"
#[1] "with cutoff 0.2 The DMR number is 77"
#[1] "with cutoff 0.3 The DMR number is 7"
#[1] "with cutoff 0.4 The DMR number is 2"

### proceed with cutoff 0.1  for CHH according to  The Plant Journal(2018)93,460-471

dmrs_CHH_0.1_BE_GE <- callDMR(dmlTest_CHH_BE_vs_GE, delta= 0.1, minlen=100, minCG=3, pct.sig=0.5, dis.merge=50, p.threshold=0.05)
save(dmrs_CHH_0.1_BE_GE, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/dmrs_CHH_0.1_BE_GE_5_9_2019.RData")


################################################################################
#
#                              tomato gene annotation 
#
################################################################################

Tomato_ITAG3.0_gene_models <- import("/data/users/xzy50/ITAG3.0_gene_models.gff")
head(Tomato_ITAG3.0_gene_models)
seqlevels(Tomato_ITAG3.0_gene_models)
gene = Tomato_ITAG3.0_gene_models[ Tomato_ITAG3.0_gene_models$type == "gene", c("ID","Alias","Note","Ontology_term")]
seqlevels( gene ) <- c( "0","1","2","3","4","5","6","7","8","9","10","11","12" )
seqlevels(gene, pruning.mode = "coarse") <- c("0","1","2","3","4","5","6","7","8","9","10","11","12" )
gene = sortBySeqnameAndStart(gene)


geneUpDown2kb <- GeneUpDownStream(gene, upstream = 2000, downstream = 2000)




dmrs_CG_0.6_BE_GE_gr <- makeGRangesFromDataFrame(dmrs_CG_0.6_BE_GE, keep.extra.columns = TRUE,
                                              seqnames.field=c("chr"),
                                              start.field="start",
                                              end.field=c("end"))

dmrs_CHG_0.4_BE_GE_gr <- makeGRangesFromDataFrame(dmrs_CHG_0.4_BE_GE, keep.extra.columns = TRUE,
                                                  seqnames.field=c("chr"),
                                                  start.field="start",
                                                  end.field=c("end"))

dmrs_CHH_0.1_BE_GE_gr <- makeGRangesFromDataFrame(dmrs_CHH_0.1_BE_GE, keep.extra.columns = TRUE,
                                                  seqnames.field=c("chr"),
                                                  start.field="start",
                                                  end.field=c("end"))

dmrs_BE_GE_CG_CHG_CHH <- rbind(dmrs_CG_0.6_BE_GE, dmrs_CHG_0.4_BE_GE, dmrs_CHH_0.1_BE_GE)

dmrs_BE_GE_CG_CHG_CHH_gr <- makeGRangesFromDataFrame(dmrs_BE_GE_CG_CHG_CHH, keep.extra.columns = TRUE,
                                                     seqnames.field=c("chr"),
                                                     start.field="start",
                                                     end.field=c("end"))
seqlevels(dmrs_BE_GE_CG_CHG_CHH_gr) <- c("0","1","2","3","4","5","6","7","8","9","10","11","12" )


head(dmrs_CHH_0.1_BE_GE)

seqlevels(dmrs_CHH_0.1_BE_GE_gr)
seqlevels(dmrs_CHH_0.1_BE_GE_gr) <- c("0","1","2","3","4","5","6","7","8","9","10","11","12" )

dmrs_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb  = subsetByOverlaps(geneUpDown2kb, dmrs_BE_GE_CG_CHG_CHH_gr, type = "any", ignore.strand = FALSE)
length(dmrs_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb)
unique(dmrs_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb)


dmrs_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376 <- substr(x = dmrs_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb$ID, start = 6, stop = 21)



#############################################################################################################
####################################            blastp with R           #####################################

#############################################################################################################


#' @rdname blastp
#' @title BLASTP programs search protein subjects using a protein query

#' @description This function is wrapping of BLASTP programs search protein

#'     subjects using a protein query. The function applies the command line

#'     version of the Basic Local Alignment Search Tool (BLAST), as provided by

#'     the \href{https://www.ncbi.nlm.nih.gov/books/NBK279671/}{NCBI}.

#'    

#'       The reason why the function is addressed to only BLASTP is because we

#'     are interested only in the codon sequences and, in general, the sequence

#'     alignment based on DNA sequence breaks the codons. The best accuracy is

#'     obtained by traslating the codon sequence into amino acid sequence.

#'    

#' @details BLAST command line sofware must be previously installed in the local

#'     computer where this function will be called, as explained at 

#'     \href{https://www.ncbi.nlm.nih.gov/books/NBK279671/}{NCBI}. 

#'    

#'     *Install BLAST in Ubuntu*

#'     This can be found in the repositories: "ncbi-blast+"

#'     To manage available BLAST databases, a subdirectory named db should be

#'     created. For the example installation, the following command creates such

#'     directory under "ncbi-blast+" directory: 

#'     \enumerate{

#'         \item _login as root_:

#'               sudo su root

#'         \item _Locate folder "ncbi-blast+" in your PC. In my case_

#'               _"/usr/share/doc/ncbi-blast+/". First update!_

#'              updatedb

#'              locate -b ncbi-blast

#'         \item _Change to directory "ncbi-blast+"_

#'               cd /usr/share/doc/ncbi-blast+/

#'         \item _Make directory "db"_

#'               mkdir ./db

#'         \item _Export BLASTDB and PATH for blast (if different from_ 

#'               _"/usr/bin/blastp")_

#'               export BLASTDB=”$/usr/share/doc/ncbi-blast+/db”

#'               BLASTDB=/usr/share/doc/ncbi-blast+/db'               

#'     }

#'    

#'     It is worthy to notice that a ncbi-blast database contains six files

#'     with the following structure:

#'     \enumerate{

#'         \item dataBaseName.phr

#'         \item dataBaseName.pin

#'         \item dataBaseName.pog

#'         \item dataBaseName.psd

#'         \item dataBaseName.psi

#'         \item dataBaseName.psq

#'     }

#' @param query.seq An AAStringSet object

#'     (\code{link[Biostrings]{XStringSet-class}}) from *Biostrings* package

#'     carrying the query amino acid sequence(s). This obejct can be created

#'     reading the fasta file of interest with function *readDNAStringSet* from

#'     *Biostrings* package.

#' @param dtb whole path to the database of sequences. For example, it would be

#'     dtb = "/Path_to_database/dataBaseName" (do not add any extension

#'     after the "dataBaseName", see details). The database must be

#'     created by function 'blastdbcmd', as specified for

#'     \href{https://www.ncbi.nlm.nih.gov/books/NBK279671/}{NCBI} blast. Default

#'     is NULL. If not provided, then function *blasp* will try to create it

#'     using the information provided in parameters *db.fa* and *dir.fa*.

#' @param db.fa,dir.fa Name and directory of the file containing the fasta

#'     sequence(s) to build a ncbi-blast database if dtb = NULL.

#' @param tmp A directory where to write a temporal file. 

#' @param maxTargetSeqs Integer value to pass to blastp. The maximum number of

#'     sequences to target in the database. Default is 2, i.e., for each query

#'     sequence, the two top sequences from the database with the minimum

#'     pairwise alignment *evalue* will be returned.

#' @param numcode The NCBI genetic code number for translation. By default the

#'     standard genetic code is used.

#' @param num.cores,tasks Parameters for parallel computation using

#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to

#'     use, i.e. at most how many child processes will be run simultaneously

#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job

#'     (only for Linux OS).

#' @importFrom seqinr syncodons

#' @importFrom Biostrings writeXStringSet

#' @importFrom BiocParallel MulticoreParam SnowParam bplapply

#' @export

#' @examples

#' # ## ============================= For Linux OS =============================

#' # ## ---------- Download the amino acid sequence from github ----------------
#'
#' 
#' 
library(seqinr)
library(Biostrings)
library(BiocParallel)


dir.create("/data/users/xzy50/TomatoBlast")
# download tomato CDS file ITAG3.0
download.file(url = "ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/annotation/ITAG3.0_release/ITAG3.0_proteins.fasta", 
              destfile = "/data/users/xzy50/TomatoBlast/ITAG3.0_proteins.fasta", method = "wget", quiet=FALSE)


download.file(url ="https://www.arabidopsis.org/download_files/Sequences/Araport11_blastsets/Araport11_genes.201606.pep.fasta.gz",
              destfile = "/data/users/xzy50/TomatoBlast/AT11_PEP.fasta", method = "wget", quiet=FALSE)

# ---------------- Set the information needed for blastp ------------------

# set up the reference 

# 
# Arabidopsis Information Portal
# 
# ```
# Arabidopsis thaliana Genome Annotation Official Release (Approved by NCBI GenBank)
# Version: Araport11
# Release date: June 2016
# ```
# 
# Highlights of the Araport11 Official Release
# 
# * 27,655 protein-coding genes
# * 5,178 non-coding genes
# * 3,901 transposable element genes
# * 952 pseudogenes
# * 508 novel transcribed regions
# * 111 upstream open reading frames
# * No changes made to the genome assembly
# 
# Changes compared to previous release (Apr 2016)
# 
# * Incorporates previously missed TAIR10 non-coding genes (of type "other RNA")
# 
# Files in this release:
#   
# Araport11_GFF3_genes_transposons.201606.gff3.gz   - Araport11 annotation in GFF3 format
# Araport11_GFF3_genes_transposons.201606.gtf.gz    - Araport11 annotation in GTF format
# Araport11_GFF3_genes_transposons.201606.stats.txt - Annotation summary statistics
# Araport11_genes.201606.cds.fasta.gz               - Coding sequences in FASTA format
# Araport11_genes.201606.pep.fasta.gz               - Protein translations in FASTA format
# Araport11_genes.201606.cdna.fasta.gz              - Transcript sequences in FASTA format
# README.201606.md                                  - This README file
# 
# 
# * * *
#   
#   If you have any questions regarding the data, please write to <mailto:araport@jcvi.org>
#   
#   [1]: http://bit.ly/aip-logo

db.fa <- "/Araport11_genes.201606.pep.fasta"
dir.fa <- "/data/users/xzy50/TomatoBlast"
dir.db <- "/data/users/xzy50/TomatoBlast"

# library(Biostrings)
# url1 <- paste0("https://github.com/genomaths/seqalignments/raw/master/","Pyrococcus/pep_seq_all/","Pyrococcus_abyssi_ge5.ASM19593v2.pep.all.fa")
# url2 <- paste0("https://github.com/genomaths/seqalignments/raw/master/","Pyrococcus/pep_seq_all/Pyrococcus_furiosus_dsm_3638.","ASM730v1.pep.all.fa")
# 
# setwd("/data/users/xzy50")
# dir.create("/data/users/xzy50/pyroc")
# outfile1 <- "/data/users/xzy50/pyroc/p_abiss.fa"
# outfile2 <- "/data/users/xzy50/pyroc/p_furiosus.fa"
# 
# download.file(url = url1, destfile = outfile1, method = "wget", quiet=TRUE)
# download.file(url = url2, destfile = outfile2, method = "wget", quiet=TRUE)
# 
# 
# 
# # ---------------- Set the information needed for blastp ------------------
# 
# # set up the reference 
# db.fa <- "/data/users/xzy50/pyroc/p_abiss.fa"
# dir.fa <- "/data/users/xzy50/pyroc"
# dir.db <- "/data/users/xzy50/pyroc"

tmp <- dir.db

### set up the query sequence 
file <- "/data/users/xzy50/TomatoBlast/ITAG3.0_proteins.fasta"
seqs <- readAAStringSet(filepath = file, format = "fasta")


#grep(DMG_13_vs_12_CG_208, names(seqs))
#?AAStringSet
keepers<-c()  
i=1
for(item in dmrs_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376) {
  tempy <- grep(item,names(seqs))
  if(is.integer(tempy)){
    keepers[i]<-tempy
    i=i+1
  }else{
    print("problem, is.integer is not an integer")
  }
}

query.seq <- seqs[keepers]

# query.seq <- seqs[1:20]
seq.name  <- substr(x = names(query.seq), start = 1, stop = 16)

seqsAT <- readAAStringSet(filepath = "/data/users/xzy50/TomatoBlast/AT11_PEP.fasta", format = "fasta")
query.seqsAT <- seqsAT[1:20]
names(seqsAT)

DMG_blastp <- blastp(query.seq=query.seq, dtb = NULL, seq.name = seq.name, db.fa = db.fa, 
                     dir.fa = dir.fa, tmp = tmp,  num.cores = 10L, tasks = 10L)

blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376 <- blastp(query.seq=query.seq, dtb = NULL, seq.name = seq.name, db.fa = db.fa, maxTargetSeqs = 1,
                     dir.fa = dir.fa, tmp = tmp,  num.cores = 30L, tasks = 30L)

save(blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376.RData")
write.csv(blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376.csv")


blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376$qseqid <- substr(x = blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376$qseqid , start = 1, stop = 16)
blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376$sseqid <- substr(x = blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376$sseqid , start = 1, stop = 9)



## ------------------------- delete 'pyroc' folder  -----------------------

# unlink(x = "/tmp/pyroc", recursive = TRUE)



blastp <- function(query.seq, dtb = NULL, seq.name = NULL, db.fa = NULL,
                   
                   dir.fa = NULL, tmp = getwd(), maxTargetSeqs = 2, numcode = 1,
                   
                   num.cores = 1L, tasks = 0L) {
  
  if (missing(query.seq)) stop("Provide query sequence(s) in FASTa format")
  
  if (is.null(db.fa) && is.null(dtb))
    
    stop("Provide sequence database or a fasta file to create it")
  
  
  
  if (!is.null(dtb)) {
    
    db <- suppressWarnings(try(system(paste0("blastdbcmd -db ", dtb,
                                             
                                             " -info"),
                                      
                                      intern = TRUE, ignore.stderr = TRUE),
                               
                               silent = TRUE))
    
    if (inherits(db, "try-error") || length(db) == 0)
      
      stop("The database provided is not valid")
    
  }
  
  
  
  if (is.null(dir.fa) && is.null(dtb)) stop("Provide the fasta file directory")
  
  else {
    
    if (is.null(dtb) && !is.null(db.fa)) {
      
      db <- suppressWarnings(try(system(paste0("blastdbcmd -db ", dir.db,
                                               
                                               db.fa, ".prot -info"),
                                        
                                        intern = TRUE,
                                        
                                        ignore.stderr = TRUE),
                                 
                                 silent = TRUE))
      
      if (!inherits(db, "try-error") && length(db) > 0)
        
        dtb <-paste0(dir.db, db.fa, ".prot")
      
    }
    
  }
  
  
  
  if (is.null(dtb) && !is.null(db.fa) && !is.null(dir.fa)) {
    
    newdb <- paste0("makeblastdb -in ", dir.fa, db.fa, " -dbtype prot ",
                    
                    "-parse_seqids -out ", dir.db, db.fa, ".prot ",
                    
                    "-title ", db.fa)
    
    system(newdb)
    
    dtb <- try(system(paste0("blastdbcmd -db ", dir.db, db.fa, ".prot -info")
                      
                      , intern = TRUE))
    
    if (inherits(dtb, "try-error") || length(dtb) == 0)
      
      stop("Database creation failed")
    
    dtb <-paste0(dir.db, db.fa, ".prot")
    
  }
  
  
  
  # === Auxiliar function to perform blastp through OS command  ===
  
  blast <- function(k, query, dtb, tmp, seq.name, maxTargetSeqs) {
    
    sname <- seq.name[k]
    
    writeXStringSet(query[k], filepath = paste0(tmp, "tmp", sname, ".fasta"))
    
    str1 = paste0("blastp -db ", dtb, " -query ", tmp, "tmp", sname, ".fasta")
    
    str2 = paste0("-out ", tmp, "tmp", sname,
                  
                  ".txt -outfmt '6 qseqid sseqid pident ",
                  
                  "qcovs bitscore score evalue' -max_target_seqs ",
                  
                  maxTargetSeqs)
    
    system(paste( str1, str2, sep = " " ))
    
    tmp1 <- try( read.delim(paste0(tmp, "tmp", sname,".txt"),
                            
                            header = FALSE ), silent = TRUE)
    
    if (!inherits(tmp1, "try-error")) {
      
      res <- DataFrame(tmp1)
      
      colnames(res) <- c("qseqid", "sseqid", "pident", "qcovs",
                         
                         "bitscore", "score", "evalue")
      
      file.remove(c(paste0(tmp, "tmp", sname,".txt"),
                    
                    paste0(tmp, "tmp", sname, ".fasta" ) ) )
      
      return(res)
      
    } else {
      
      file.remove(paste0(tmp, "tmp", sname, ".fasta" ))
      
      res = DataFrame(qseqid = NA, sseqid = NA, pident = NA,
                      
                      qcovs = NA, bitscore = NA, score = NA, evalue = NA)
      
    }
    
    return(res)
    
  }
  
  # ------------------------------------------------------------------------- #
  
  
  
  if (length(query.seq) > 1) {
    
    # Set parallel computation
    
    if (Sys.info()['sysname'] == "Linux") {
      
      bpparam <- MulticoreParam(workers=num.cores, tasks = tasks)
      
    } else bpparam <- SnowParam(workers = num.cores, type = "SOCK")
    
    num.seq <- 1:length(query.seq)
    
    res <- bplapply(num.seq, blast, query=query.seq, dtb=dtb, tmp=tmp,
                    
                    seq.name=seq.name, maxTargetSeqs=maxTargetSeqs,
                    
                    BPPARAM = bpparam)
    
    res <- do.call(rbind, res)
    
  } else res <- blast(1, query=query.seq, dtb=dtb, tmp=tmp,
                      
                      seq.name=seq.name, maxTargetSeqs)
  
  
  
  return(res)
  
}


######################################################################## #
#
##### -------------- NETWORK BASES ENRICHMENT ANALYSIS ---------- #######
#
######################################################################## #

library(EnrichmentBrowser)
library(GenomicRanges)
library(data.table )

# ============= KEGG datasets ============
kegg.gs <- getGenesets( "ath", db = "kegg" )
str(head(kegg.gs))

# GO terms of a selected ontology as gene sets
# ontology of interest (BP, MF or CC):
# BP: Biological Process
# MF: Molecular Function
# CC: Cellular Component
ath.bp.gs <- getGenesets(org = "ath", db = "go", go.onto = "BP", go.mode = "GO.db")
str(head(ath.bp.gs))

ath.mf.gs <- getGenesets(org = "ath", db = "go", go.onto = "MF", go.mode = "GO.db")
str(head(ath.mf.gs))

ath.cc.gs <- getGenesets(org = "ath", db = "go", go.onto = "CC", go.mode = "GO.db")
str(head(ath.cc.gs))

#### *** Compilation of a gene regulatory network from KEGG pathways (already
#### done) !!!
ath.grn <- compileGRN("ath")
head(ath.grn)
#      FROM        TO          TYPE
# [1,] "AT1G48970" "AT1G04170" "+" 
# [2,] "AT1G48970" "AT2G18720" "+" 
# [3,] "AT1G48970" "AT2G40290" "+" 
# [4,] "AT1G48970" "AT3G07920" "+" 
# [5,] "AT1G48970" "AT5G01940" "+" 
# [6,] "AT1G48970" "AT5G05470" "+"

#########################################################################################
#
# --------------- Network Enrichment Analysis Test . Based on GO -----------------------
#
#########################################################################################
# Compute NEAT (Signorelli et al., 2016) with small correction, a test for network enrichment
# analysis between/from a first list of sets ('A sets') and/to a second list of sets ('B sets').
# neat is the R package that implements NEAT, the Network Enrichment Analysis Test which is 
# presented in Signorelli, M., Vinciotti, V., Wit, E. C. (2016). NEAT: an efficient network 
# enrichment analysis test. BMC Bioinformatics, 17:352.

library( neat ) # https://cran.r-project.org/web/packages/neat/vignettes/neat.html
source("/data/R_functions/.neat.R" )
# ATH_GOs <- read.delim("/media/robersy/OS/Users/Robersy/Documents/Work/TAIR10_gff3/GOs/ATH_GO_GOSLIM.txt", header = FALSE)
# # ATH_GOs <- read.delim("C:/Users/Robersy/Documents/Work/TAIR10_gff3/GOs/ATH_GO_GOSLIM.txt", header = FALSE)
# ATH_GOs = data.table( ATH_GOs )
# colnames( ATH_GOs ) <- c( "locus", "TAIR.accession", "object.name", "relationship.type", "GO.term", "GO.ID",
#                           "TAIR.Keyword.ID", "Aspect", "GOslim term", "Evidence.code", "Evidence.description", 
#                           "Evidence with", "Reference", "Annotator", "Date.annotated" )
# 
# # *** To build custom annotations ****
# locus = ATH_GOs$locus
# GOs = paste0( ATH_GOs$GO.ID, "_", ATH_GOs$GO.term )
# GO2Genes = by( locus, GOs, function(s) unique( as.character(s) ), simplify = T )
# GO2Genes = GO2Genes[1:length(GO2Genes)]
#ATn = rownames( DMG )
ATn = blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376$sseqid
# ath.grn <- compile.grn.from.kegg( "/JOB/Work/KEGG/ath.zip" )
# head( ath.grn )
genesNet = unique( as.vector( ath.grn[ , 1:2 ] ) ) 

test = .neat( alist = list('AT.DEG' = ATn ), blist = ath.bp.gs, network = ath.grn[ , 1:2 ],
              nettype = 'directed', nodes = genesNet, alpha = 0.05 )
sum( test$pvalue <= 0.05 )
# ntest = test[ , ]
ntest = test[( test$nab > 0 & test$pvalue <= 0.05 ), ]
ntest = ntest[ (ntest$conclusion == "Overenrichment" ), ]

write.csv(ntest, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/NEAT_GO_blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376.csv")


GO.NetAB = ath.bp.gs[ match( as.character( ntest$B ), names( ath.bp.gs ) ) ]
GeneGO.Net = sapply( as.list( GO.NetAB ), function(s) s[ na.omit( match( ATn, s ) ) ] ) 

Netw = names( GeneGO.Net )
GenesGO.Net = c()
for( k in 1:length( GeneGO.Net ) ) {
  if( length( GeneGO.Net[[ k ]]) > 0 ) 
    GenesGO.Net = rbind( GenesGO.Net, data.frame( GeneID = GeneGO.Net[[ k ]], Network = Netw[ k ] ) )
}

GenesGO.Net = data.table( GenesGO.Net )
GenesGO.Nets = GenesGO.Net[ , list( Network = list( unique( as.character( Network ) ) ) ), by = GeneID ]
#DEGlist = rownames( DMG )
DEGlist = blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376$sseqid


### Anotation
AG_gff3 = import(con = paste0("ftp://ftp.ensemblgenomes.org/pub/release-38/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.38.gff3.gz"))
# or AG_gff3 = import("/data/TAIR10_gff3/Arabidopsis_thaliana.TAIR10.38.gff3.gz")
AG_anotation = AG_gff3[ AG_gff3$type == "gene",  c( "gene_id", "Name", "description" ) ]
seqlevels( AG_anotation ) <- c( "1","2","3","4","5","M","C" )
seqlevels(AG_anotation, pruning.mode = "coarse") <- c("1", "2", "3", "4", "5")

AG_anotation_df <- as.data.frame(mcols(AG_anotation))
GenesGO.Net_name <- data.frame(GenesGO.Net, name = AG_anotation_df[match(GenesGO.Net$GeneID, AG_anotation_df$gene_id),])
head(GenesGO.Net_name)
write.csv(GenesGO.Net_name_gen1, file = "2520_dmrs_genes_CG0.3_CHG_0.15_CHH0.15_gr_UpDown_1k_NEAT_GO_GENES_gen1.csv")

blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376

GenesGO.Net_name_tomato_ID <- data.frame(GenesGO.Net_name, name = blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376[match(GenesGO.Net_name$GeneID, blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376$sseqid),])
write.csv(GenesGO.Net_name_tomato_ID, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/DSS/Tomato_NEAT_GO_gebe_blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376.csv")

