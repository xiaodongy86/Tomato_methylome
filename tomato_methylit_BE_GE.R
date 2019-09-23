
#.libPaths(c("/data/users/xzy50/R/x86_64-redhat-linux-gnu-library/3.5/MethylIT_20181112","/usr/lib64/R/library","/usr/share/R/library","~/R/x86_64-redhat-linux-gnu-library/"))
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.5","/data/users/xzy50/R/x86_64-redhat-linux-gnu-library/3.5/MethylIT_20181112","/usr/lib64/R/library","/usr/share/R/library"))
library(MethylIT)
library(DSS)
library(dplyr)
setwd("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor")
.libPaths()
############################################################################################################
#                                          readCounts2GRangesList
#                                                    CG
#
############################################################################################################

samples_GE_CG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P1/GE-P1_CG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P2/GE-P2_CG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P3/GE-P3_CG.txt")
sample.id_GE_CG <- c("GE_P1","GE_P2","GE_P3")

GE_CG <- readCounts2GRangesList(filenames = samples_GE_CG, 
                                   sample.id = sample.id_GE_CG,
                                   columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                   verbose = TRUE)

save(GE_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE_CG_5_7_2019.RData")

samples_BE_CG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P1/BE-P1_CG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P2/BE-P2_CG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P3/BE-P3_CG.txt")
sample.id_BE_CG <- c("BE_P1","BE_P2","BE_P3")

BE_CG <- readCounts2GRangesList(filenames = samples_BE_CG, 
                                sample.id = sample.id_BE_CG,
                                columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                verbose = TRUE)

save(BE_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE_CG_5_7_2019.RData")


############################################################################################################
#                                          readCounts2GRangesList
#                                                    CHG
#
############################################################################################################


samples_GE_CHG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P1/GE-P1_CHG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P2/GE-P2_CHG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P3/GE-P3_CHG.txt")
sample.id_GE_CHG <- c("GE_P1","GE_P2","GE_P3")

GE_CHG <- readCounts2GRangesList(filenames = samples_GE_CHG, 
                                sample.id = sample.id_GE_CHG,
                                columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                verbose = TRUE)

save(GE_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE_CHG_5_8_2019.RData")

samples_BE_CHG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P1/BE-P1_CHG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P2/BE-P2_CHG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P3/BE-P3_CHG.txt")
sample.id_BE_CHG <- c("BE_P1","BE_P2","BE_P3")

BE_CHG <- readCounts2GRangesList(filenames = samples_BE_CHG, 
                                sample.id = sample.id_BE_CHG,
                                columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                verbose = TRUE)

save(BE_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE_CHG_5_8_2019.RData")


############################################################################################################
#                                          readCounts2GRangesList
#                                                    CHH
#
############################################################################################################


samples_GE_CHH <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P1/GE-P1_CHH.txt",
                    "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P2/GE-P2_CHH.txt",
                    "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE-P3/GE-P3_CHH.txt")
sample.id_GE_CHH <- c("GE_P1","GE_P2","GE_P3")

GE_CHH <- readCounts2GRangesList(filenames = samples_GE_CHH, 
                                 sample.id = sample.id_GE_CHH,
                                 columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                 verbose = TRUE)

save(GE_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE_CHH_5_8_2019.RData")


samples_BE_CHH <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P1/BE-P1_CHH.txt",
                    "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P2/BE-P2_CHH.txt",
                    "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE-P3/BE-P3_CHH.txt")
sample.id_BE_CHH <- c("BE_P1","BE_P2","BE_P3")

BE_CHH <- readCounts2GRangesList(filenames = samples_BE_CHH, 
                                 sample.id = sample.id_BE_CHH,
                                 columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                 verbose = TRUE)

save(BE_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE_CHH_5_8_2019.RData")


############################################################################################################
#                              BE three samples as reference (poolFromGRlist)
#                                                 CG CHG CHH
#
############################################################################################################
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE_CG_5_7_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE_CG_5_7_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE_CHG_5_8_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE_CHG_5_8_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/GE_CHH_5_8_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/BE_CHH_5_8_2019.RData")


GE_BE_CG <- c(GE_CG,BE_CG)
GE_BE_CHG <- c(GE_CHG,BE_CHG)
GE_BE_CHH <- c(GE_CHH,BE_CHH)

ref_BE_CG <- poolFromGRlist(list(BE_CG$BE_P1, BE_CG$BE_P2, BE_CG$BE_P3), stat = "sum", num.cores = 12L)
ref_BE_CHG <- poolFromGRlist(list(BE_CHG$BE_P1, BE_CHG$BE_P2, BE_CHG$BE_P3), stat = "sum", num.cores = 12L)
ref_BE_CHH <- poolFromGRlist(list(BE_CHH$BE_P1, BE_CHH$BE_P2, BE_CHH$BE_P3), stat = "sum", num.cores = 12L)

save(ref_BE_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_BE_CG.RData")
save(ref_BE_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_BE_CHG.RData")
save(ref_BE_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_BE_CHH.RData")


save(ref_BE_CG, ref_BE_CHG, ref_BE_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_BE_CG_CHG_CHH.RData")


############################################################################################################
#                                             estimateDivergence
#                                                    CG
#                                                 BE  GE 
############################################################################################################

HD_ref_BE_CG = estimateDivergence(ref = ref_BE_CG, 
                                       indiv = GE_BE_CG, 
                                       Bayesian = TRUE, 
                                       min.coverage = 4, 
                                       high.coverage = 300, 
                                       percentile = 0.999, 
                                       num.cores = 4L, tasks = 0L, verbose = FALSE )
save(HD_ref_BE_CG,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_BE_CG_5_8_2019.RData")

load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_BE_CG_5_8_2019.RData")

nlms_ref_BE_CG = nonlinearFitDist(HD_ref_BE_CG, column = 9, num.cores = 10L, verbose = TRUE)
PS_ref_BE_CG = getPotentialDIMP(LR = HD_ref_BE_CG, nlms = nlms_ref_BE_CG, div.col = 9, alpha = 0.05,
                            tv.col = 7, tv.cut = 0.2)
save(nlms_ref_BE_CG, PS_ref_BE_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_BE_CG_5_8_2019.RData")

############################################################################################################
#                                             estimateDivergence
#                                                    CHG
#                                                 BE  GE 
############################################################################################################

HD_ref_BE_CHG = estimateDivergence(ref = ref_BE_CHG, 
                                  indiv = GE_BE_CHG, 
                                  Bayesian = TRUE, 
                                  min.coverage = 4, 
                                  high.coverage = 300, 
                                  percentile = 0.999, 
                                  num.cores = 4L, tasks = 0L, verbose = FALSE )
save(HD_ref_BE_CHG,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_BE_CHG_5_8_2019.RData")

nlms_ref_BE_CHG = nonlinearFitDist(HD_ref_BE_CHG, column = 9, num.cores = 4L, verbose = TRUE)
PS_ref_BE_CHG = getPotentialDIMP(LR = HD_ref_BE_CHG, nlms = nlms_ref_BE_CHG, div.col = 9, alpha = 0.05,
                                tv.col = 7, tv.cut = 0.2)
save(nlms_ref_BE_CHG, PS_ref_BE_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_BE_CHG_5_8_2019.RData")


############################################################################################################
#                                             estimateDivergence
#                                                    CHH
#                                                 BE  GE 
############################################################################################################

#  BE_CHH
HD_ref_BE_CHH_BE = estimateDivergence(ref = ref_BE_CHH, 
                                   indiv = BE_CHH, 
                                   Bayesian = TRUE, 
                                   min.coverage = 4, 
                                   high.coverage = 300, 
                                   percentile = 0.999, 
                                   num.cores = 4L, tasks = 0L, verbose = FALSE )
save(HD_ref_BE_CHH_BE,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_BE_CHH_BE_5_8_2019.RData")

nlms_ref_BE_CHH_BE = nonlinearFitDist(HD_ref_BE_CHH_BE, column = 9, num.cores = 4L, verbose = TRUE)
PS_ref_BE_CHH_BE = getPotentialDIMP(LR = HD_ref_BE_CHH_BE, nlms = nlms_ref_BE_CHH_BE, div.col = 9, alpha = 0.05,
                                 tv.col = 7, tv.cut = 0.2)
save(nlms_ref_BE_CHH_BE, PS_ref_BE_CHH_BE, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_BE_CHH_BE_5_8_2019.RData")

#  GE_CHH
HD_ref_BE_CHH_GE = estimateDivergence(ref = ref_BE_CHH, 
                                      indiv = GE_CHH, 
                                      Bayesian = TRUE, 
                                      min.coverage = 4, 
                                      high.coverage = 300, 
                                      percentile = 0.999, 
                                      num.cores = 4L, tasks = 0L, verbose = FALSE )
save(HD_ref_BE_CHH_GE,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_BE_CHH_GE_5_8_2019.RData")

nlms_ref_BE_CHH_GE = nonlinearFitDist(HD_ref_BE_CHH_GE, column = 9, num.cores = 4L, verbose = TRUE)
PS_ref_BE_CHH_GE = getPotentialDIMP(LR = HD_ref_BE_CHH_GE, nlms = nlms_ref_BE_CHH_GE, div.col = 9, alpha = 0.05,
                                    tv.col = 7, tv.cut = 0.2)
save(nlms_ref_BE_CHH_GE, PS_ref_BE_CHH_GE, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_BE_CHH_GE_5_8_2019.RData")



load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_BE_CG_5_8_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_BE_CHG_5_8_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_BE_CHH_BE_5_8_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_BE_CHH_GE_5_8_2019.RData")

head(PS_ref_BE_CHH_GE)
PS_ref_BE_CHH <- c(PS_ref_BE_CHH_GE, PS_ref_BE_CHH_BE)


PS = getPotentialDIMP(LR = HD, nlms = nlms, div.col = 9, alpha = 0.05,
                      tv.col = 7, tv.cut = 0.2)


cutpoints = estimateCutPoint(LR = PS_ref_BE_CG, control.names = c("BE_P1", "BE_P2",  "BE_P3"),
                             treatment.names = c( "GE_P1", "GE_P2", "GE_P3"),
                             #treatment.names = c( "M_1_105", "M_1_168", "M_1_179", "M_1_27", "M_1_3"),
                             div.col = 9, verbose = T)
# > cutpoints
# $cutpoint
# BE_P1    BE_P2    BE_P3
# GE_P1 1.133861 1.133819 1.158631
# GE_P2 1.103369 1.100426 1.106281
# GE_P3 1.125184 1.106350 1.128637
# 
# $auc
# BE_P1     BE_P2     BE_P3
# GE_P1 0.9452266 0.9549155 0.8877164
# GE_P2 0.9303981 0.9452491 0.8566591
# GE_P3 0.9372825 0.9497046 0.8713461
# 
# $accuracy
# BE_P1     BE_P2     BE_P3
# GE_P1 0.9451488 0.9629636 0.8383582
# GE_P2 0.9458777 0.9652051 0.8160430
# GE_P3 0.9423486 0.9630204 0.8325256

cutpoints = estimateCutPoint(LR = PS_ref_BE_CG, control.names = c("BE_P1", "BE_P2",  "BE_P3"),
                             treatment.names = c( "GE_P1", "GE_P2", "GE_P3"),
                             #treatment.names = c( "M_1_105", "M_1_168", "M_1_179", "M_1_27", "M_1_3"),
                             div.col = 7L, verbose = T, absolute = T)
# $cutpoint
# BE_P1     BE_P2     BE_P3
# GE_P1 0.2212048 0.2098246 0.2065255
# GE_P2 0.2470760 0.2470760 0.2130542
# GE_P3 0.2470657 0.2470657 0.2100070
# 
# $auc
# BE_P1     BE_P2     BE_P3
# GE_P1 0.4197609 0.4067778 0.2482288
# GE_P2 0.4682549 0.4550395 0.2971489
# GE_P3 0.4355039 0.4223065 0.2617011
# 
# $accuracy
# BE_P1     BE_P2     BE_P3
# GE_P1 0.5994087 0.6177062 0.4078111
# GE_P2 0.6446291 0.6517606 0.4633759
# GE_P3 0.6112046 0.6180972 0.4356800

cutpoints = estimateCutPoint(LR = PS_ref_BE_CHG, control.names = c("BE_P1", "BE_P2",  "BE_P3"),
                             treatment.names = c( "GE_P1", "GE_P2", "GE_P3"),
                             #treatment.names = c( "M_1_105", "M_1_168", "M_1_179", "M_1_27", "M_1_3"),
                             div.col = 7L, verbose = T, absolute = T)
# > cutpoints
# $cutpoint
# BE_P1     BE_P2     BE_P3
# GE_P1 0.2828297 0.2785829 0.2100070
# GE_P2 0.3238157 0.2802260 0.2159443
# GE_P3 0.2828361 0.2802276 0.2137827
# 
# $auc
# BE_P1     BE_P2     BE_P3
# GE_P1 0.5141223 0.5087112 0.3726544
# GE_P2 0.5558524 0.5495562 0.4124599
# GE_P3 0.5364705 0.5306622 0.3945428
# 
# $accuracy
# BE_P1     BE_P2     BE_P3
# GE_P1 0.5792861 0.5882848 0.5901655
# GE_P2 0.5801345 0.6235429 0.5976832
# GE_P3 0.5992190 0.6067977 0.5987309

cutpoints = estimateCutPoint(LR = PS_ref_BE_CHH, control.names = c("BE_P1", "BE_P2",  "BE_P3"),
                             treatment.names = c( "GE_P1", "GE_P2", "GE_P3"),
                             #treatment.names = c( "M_1_105", "M_1_168", "M_1_179", "M_1_27", "M_1_3"),
                             div.col = 7L, verbose = T, absolute = T)
# $cutpoint
# BE_P1     BE_P2     BE_P3
# GE_P1 0.3000230 0.3000230 0.9782609
# GE_P2 0.3077176 0.4167852 0.3235460
# GE_P3 0.3077114 0.4167267 0.9782609
# 
# $auc
# BE_P1     BE_P2     BE_P3
# GE_P1 0.5538652 0.5455926 0.4686295
# GE_P2 0.5999885 0.5919347 0.5162147
# GE_P3 0.5709986 0.5628809 0.4867955
# 
# $accuracy
# BE_P1     BE_P2     BE_P3
# GE_P1 0.5685093 0.5676842 0.3541403
# GE_P2 0.6068254 0.5182049 0.5660600
# GE_P3 0.5781023 0.4819985 0.3382032


DIMPs_CG  = selectDIMP(PS_ref_BE_CG, div.col = 9, cutpoint = 1.133819 )
DIMPs_CHG = selectDIMP(PS_ref_BE_CHG, div.col = 9, cutpoint = 1.133819)
DIMPs_CHH = selectDIMP(PS_ref_BE_CHH, div.col = 9, cutpoint = 1.133819 )
save(DIMPs_CG, DIMPs_CHG, DIMPs_CHH, file = "DIMPs_BE_GE_5_23.RData" )

load("DIMPs_BE_GE_5_23.RData")

CG_groups = evaluateDIMPclass(LR = DIMPs_CG, control.names = c("BE_P1", "BE_P2", "BE_P3"),
                              treatment.names = c("GE_P1", "GE_P2", "GE_P3"),
                              column = c(hdiv = TRUE, TV = TRUE, wprob = TRUE, pos = TRUE),
                              classifier = "pca.qda",
                              n.pc = 4,
                              num.boot = 5,
                              mc.cores = 5,
                              center = TRUE, scale = TRUE,
                              output = "all", prop = 0.6)


# cutpoint = 1.133819
# Accuracy : 0.9689          
# 95% CI : (0.9686, 0.9692)
# No Information Rate : 0.8674          
# P-Value [Acc > NIR] : < 2.2e-16       
# 
# Kappa : 0.8501          
# 
# Mcnemar's Test P-Value : < 2.2e-16       
# 
# Sensitivity : 1.0000          
# Specificity : 0.7657          
# Pos Pred Value : 0.9654          
# Neg Pred Value : 1.0000          
# Prevalence : 0.8674          
# Detection Rate : 0.8674          
# Detection Prevalence : 0.8985          
# Balanced Accuracy : 0.8829          
# 
# 'Positive' Class : TT              
# 
# 
# $con.mat$FDR
# [1] 0.03457722



# Reference
# Prediction      CT      TT
# CT  119108       0
# TT   36442 1017489
# 
# Accuracy : 0.9689          
# 95% CI : (0.9686, 0.9692)
# No Information Rate : 0.8674          
# P-Value [Acc > NIR] : < 2.2e-16       
# 
# Kappa : 0.8501          
# 
# Mcnemar's Test P-Value : < 2.2e-16       
# $con.mat$FDR
# [1] 0.03457722


################################################################################
#
#                              tomato gene annotation 
#
################################################################################

Tomato_ITAG3.0_gene_models <- import("/data/users/xzy50/ITAG3.0_gene_models.gff")
head(Tomato_ITAG3.0_gene_models)
seqlevels(Tomato_ITAG3.0_gene_models)
gene = Tomato_ITAG3.0_gene_models[ Tomato_ITAG3.0_gene_models$type == "gene", c("ID","Alias","Note","Ontology_term")]
#seqlevels( gene ) <- c( "0","1","2","3","4","5","6","7","8","9","10","11","12" )
#seqlevels(gene, pruning.mode = "coarse") <- c("0","1","2","3","4","5","6","7","8","9","10","11","12" )
gene = sortBySeqnameAndStart(gene)

seqlevels(gene)
geneUpDown2kb <- GeneUpDownStream(gene, upstream = 2000, downstream = 2000)
seqlevels(DIMPs_CG)

DIMPs_GE_P1 = c(DIMPs_CG$GE_P1, DIMPs_CHG$GE_P1, DIMPs_CHH$GE_P1)
DIMPs_GE_P2 = c(DIMPs_CG$GE_P2, DIMPs_CHG$GE_P2, DIMPs_CHH$GE_P2)
DIMPs_GE_P3 = c(DIMPs_CG$GE_P3, DIMPs_CHG$GE_P3, DIMPs_CHH$GE_P3)
DIMPs_BE_P1 = c(DIMPs_CG$BE_P1, DIMPs_CHG$BE_P1, DIMPs_CHH$BE_P1)
DIMPs_BE_P2 = c(DIMPs_CG$BE_P2, DIMPs_CHG$BE_P2, DIMPs_CHH$BE_P2)
DIMPs_BE_P3 = c(DIMPs_CG$BE_P3, DIMPs_CHG$BE_P3, DIMPs_CHH$BE_P3)

DIMPs_GEP1_gene2kb = getDIMPatGenes(DIMPs_GE_P1, geneUpDown2kb)
DIMPs_GEP2_gene2kb = getDIMPatGenes(DIMPs_GE_P2, geneUpDown2kb)
DIMPs_GEP3_gene2kb = getDIMPatGenes(DIMPs_GE_P3, geneUpDown2kb)
DIMPs_BEP1_gene2kb = getDIMPatGenes(DIMPs_BE_P1, geneUpDown2kb)
DIMPs_BEP2_gene2kb = getDIMPatGenes(DIMPs_BE_P2, geneUpDown2kb)
DIMPs_BEP3_gene2kb = getDIMPatGenes(DIMPs_BE_P3, geneUpDown2kb)

Genes_DIMPs = uniqueGRanges(list( DIMPs_BEP1_gene2kb[, 2], DIMPs_BEP2_gene2kb[, 2], DIMPs_BEP3_gene2kb[, 2], 
                                  DIMPs_GEP1_gene2kb[, 2], DIMPs_GEP2_gene2kb[, 2], DIMPs_GEP3_gene2kb[, 2])
                            ,type = "equal", verbose = TRUE,ignore.strand = FALSE )

colnames( mcols(Genes_DIMPs)) <- c("BEP1","BEP2","BEP3","GEP1","GEP2","GEP3")

HITS = findOverlaps(geneUpDown2kb, Genes_DIMPs, type = "equal")
Genes_DIMPs_ID <- gene[queryHits(HITS), ]

dmps = data.frame( mcols( Genes_DIMPs ))
dmps = apply( dmps, 2, as.numeric )
rownames(dmps) <- as.vector(Genes_DIMPs_ID$Alias)

condition = data.frame(condition = factor(c("CT", "CT", "CT", 
                                            "TT", "TT", "TT"),
                                          levels = c("CT", "TT")))

rownames(condition) <- c("BEP1","BEP2","BEP3","GEP1","GEP2","GEP3")
DIMR <- DESeqDataSetFromMatrix(countData = dmps,
                                  colData = condition,
                                  design = formula( ~ condition ),
                                  rowRanges = Genes_DIMPs)
save(DIMR, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/DIMR_5_23_2019.RData")

####  5_16_2019   parameters for 4961 DMGs
DMGs = countTest2(DS = DIMR, num.cores = 4L, minCountPerIndv = 8, 
                  countFilter = TRUE, CountPerBp = 0.005, maxGrpCV = c(0.5, 0.5),
                  FilterLog2FC = TRUE,Minlog2FC = 1, pvalCutOff = 0.05,
                  MVrate = .95, verbose = FALSE)

write.csv(DMGs, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/DMGs_2kb_4961_5_16_2019.csv",row.names = ROWNAMES(DMGs), quote = FALSE)

#########################################################################################################################################################
#
#                                                                             blastp 
#########################################################################################################################################################

library(seqinr)
library(Biostrings)
library(BiocParallel)


dir.create("/data/users/xzy50/TomatoBlast")
# download tomato CDS file ITAG3.0
#download.file(url = "ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/annotation/ITAG3.0_release/ITAG3.0_proteins.fasta", 
#              destfile = "/data/users/xzy50/TomatoBlast/ITAG3.0_proteins.fasta", method = "wget", quiet=FALSE)


#download.file(url ="https://www.arabidopsis.org/download_files/Sequences/Araport11_blastsets/Araport11_genes.201606.pep.fasta.gz",
#              destfile = "/data/users/xzy50/TomatoBlast/AT11_PEP.fasta", method = "wget", quiet=FALSE)


db.fa <- "/Araport11_genes.201606.pep.fasta"
dir.fa <- "/data/users/xzy50/TomatoBlast"
dir.db <- "/data/users/xzy50/TomatoBlast"

tmp <- dir.db

### set up the query sequence 
file <- "/data/users/xzy50/TomatoBlast/ITAG3.0_proteins.fasta"
seqs <- readAAStringSet(filepath = file, format = "fasta")

#grep(DMG_13_vs_12_CG_208, names(seqs))
#?AAStringSet
DMGs_4961_BE_GE <- ROWNAMES(DMGs)

keepers<-c()  
i=1
for(item in DMGs_4961_BE_GE) {
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



blastp_DMGs_4961_BE_GE_geneUpDown2kb <- blastp(query.seq=query.seq, dtb = NULL, seq.name = seq.name, db.fa = db.fa, maxTargetSeqs = 1,
                                                        dir.fa = dir.fa, tmp = tmp,  num.cores = 30L, tasks = 30L)

blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit <- na.omit(blastp_DMGs_4961_BE_GE_geneUpDown2kb)
blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100 <- blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit[blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit$score > 100,]

blastpIDs <- unique(blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100$qseqid)
blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df <- as.data.frame(blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100)

unique_match_blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df <-as.data.frame(blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df[1,])  
i=1
for(item in blastpIDs){
  unique_match_blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df[i,] = blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df[blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df$qseqid==item,][1,]; i=i+1}
dim(unique_match_blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df)
#3464


unique_match_blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df$qseqid <- substr(x = unique_match_blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df$qseqid , start = 1, stop = 16)
unique_match_blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df$sseqid <- substr(x = unique_match_blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df$sseqid , start = 1, stop = 9)

write.csv(unique_match_blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/unique_match_blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df.csv")
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
ATn = unique_match_blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df$sseqid
# ath.grn <- compile.grn.from.kegg( "/JOB/Work/KEGG/ath.zip" )
# head( ath.grn )
genesNet = unique( as.vector( ath.grn[ , 1:2 ] ) ) 

test = .neat( alist = list('AT.DEG' = ATn ), blist = ath.bp.gs, network = ath.grn[ , 1:2 ],
              nettype = 'directed', nodes = genesNet, alpha = 0.05 )
sum( test$pvalue <= 0.05 )
# ntest = test[ , ]
ntest = test[( test$nab > 0 & test$pvalue <= 0.05 ), ]
ntest = ntest[ (ntest$conclusion == "Overenrichment" ), ]

write.csv(ntest, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/unique_match_blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_df_NEAT_GO.csv")


GO.NetAB = ath.bp.gs[ match( as.character( ntest$B ), names( ath.bp.gs ) ) ]
GeneGO.Net = sapply( as.list( GO.NetAB ), function(s) s[ na.omit( match( ATn, s ) ) ] ) 
 

# [1] "GO:0000160_phosphorelay_signal_transduction_system"                            "GO:0006351_transcription,_DNA-templated"                                      
# [3] "GO:0006355_regulation_of_transcription,_DNA-templated"                         "GO:0006457_protein_folding"                                                   
# [5] "GO:0006914_autophagy"                                                          "GO:0006952_defense_response"                                                  
# [7] "GO:0007623_circadian_rhythm"                                                   "GO:0009408_response_to_heat"                                                  
# [9] "GO:0009585_red,_far-red_light_phototransduction"                               "GO:0009626_plant-type_hypersensitive_response"                                
# [11] "GO:0009649_entrainment_of_circadian_clock"                                     "GO:0009664_plant-type_cell_wall_organization"                                 
# [13] "GO:0009735_response_to_cytokinin"                                              "GO:0009736_cytokinin-activated_signaling_pathway"                             
# [15] "GO:0009744_response_to_sucrose"                                                "GO:0009873_ethylene-activated_signaling_pathway"                              
# [17] "GO:0009886_post-embryonic_animal_morphogenesis"                                "GO:0009963_positive_regulation_of_flavonoid_biosynthetic_process"             
# [19] "GO:0010031_circumnutation"                                                     "GO:0010075_regulation_of_meristem_growth"                                     
# [21] "GO:0010103_stomatal_complex_morphogenesis"                                     "GO:0010114_response_to_red_light"                                             
# [23] "GO:0010440_stomatal_lineage_progression"                                       "GO:0010629_negative_regulation_of_gene_expression"                            
# [25] "GO:0042023_DNA_endoreduplication"                                              "GO:0042127_regulation_of_cell_proliferation"                                  
# [27] "GO:0045087_innate_immune_response"                                             "GO:0048316_seed_development"                                                  
# [29] "GO:0048443_stamen_development"                                                 "GO:0048507_meristem_development"                                              
# [31] "GO:0050821_protein_stabilization"                                              "GO:0051726_regulation_of_cell_cycle"                                          
# [33] "GO:0071281_cellular_response_to_iron_ion"                                      "GO:2000072_regulation_of_defense_response_to_fungus,_incompatible_interaction"

Netw = names( GeneGO.Net )
GenesGO.Net = c()
for( k in 1:length( GeneGO.Net ) ) {
  if( length( GeneGO.Net[[ k ]]) > 0 ) 
    GenesGO.Net = rbind( GenesGO.Net, data.frame( GeneID = GeneGO.Net[[ k ]], Network = Netw[ k ] ) )
}

GenesGO.Net = data.table( GenesGO.Net )
GenesGO.Nets = GenesGO.Net[ , list( Network = list( unique( as.character( Network ) ) ) ), by = GeneID ]
#DEGlist = rownames( DMG )
DEGlist = DMGs_3905_BE_GE_geneUpDown2kb_na.omit_score100


### Anotation
AG_gff3 = import(con = paste0("ftp://ftp.ensemblgenomes.org/pub/release-38/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.38.gff3.gz"))
# or AG_gff3 = import("/data/TAIR10_gff3/Arabidopsis_thaliana.TAIR10.38.gff3.gz")
AG_anotation = AG_gff3[ AG_gff3$type == "gene",  c( "gene_id", "Name", "description" ) ]
seqlevels( AG_anotation ) <- c( "1","2","3","4","5","M","C" )
seqlevels(AG_anotation, pruning.mode = "coarse") <- c("1", "2", "3", "4", "5")

AG_anotation_df <- as.data.frame(mcols(AG_anotation))
GenesGO.Net_name <- data.frame(GenesGO.Net, name = AG_anotation_df[match(GenesGO.Net$GeneID, AG_anotation_df$gene_id),])
head(GenesGO.Net_name)

#write.csv(GenesGO.Net_name, file = "/data5/tomato_RNAseq/GenesGO.Net_name_NEAT_GO_GENES_5_15_2019.csv")

blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100

GenesGO.Net_name_tomato_ID <- data.frame(GenesGO.Net_name, name = blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100[match(GenesGO.Net_name$name.gene_id, blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100$sseqid),])
dim(GenesGO.Net_name_tomato_ID)


GenesGO.Net_name_tomato_ID$name.qseqid <- substr(x = GenesGO.Net_name_tomato_ID$name.qseqid , start = 1, stop = 14)
GenesGO.Net_name_tomato_ID <- na.omit(GenesGO.Net_name_tomato_ID)

DMGs_GenesGO.Net_name_tomato_ID <- data.frame(GenesGO.Net_name_tomato_ID, name = DMGs[match(GenesGO.Net_name_tomato_ID$name.qseqid,  ROWNAMES(DMGs)),])

write.csv(DMGs_GenesGO.Net_name_tomato_ID, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/DMGs_GenesGO.Net_name_tomato_ID_BE_GE_5_16_2019.csv")


blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_name <- data.frame(blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100, name = AG_anotation_df[match(blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100$sseqid, AG_anotation_df$gene_id),])


blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_name$qseqid <- substr(x = blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_name$qseqid , start = 1, stop = 14)

blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_name_tomato_ID <- data.frame(blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_name, DMGS = DMGs[match(blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_name$qseqid, ROWNAMES(DMGs)),])



write.csv(blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_name_tomato_ID_unique, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_name_tomato_ID_unique_5_16_2019.csv")
dim(blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_name_tomato_ID)

blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_name_tomato_ID_unique <- unique(blastp_DMGs_4961_BE_GE_geneUpDown2kb_na.omit_score100_name_tomato_ID$qseqid)



########################################################################################################################################
#
#                                              make wig / bed for the IGB vasulization  
#
#
########################################################################################################################################
library(MethylIT)
library(rtracklayer)
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_BE_CG_5_8_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_BE_CHG_5_8_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_BE_CHH_GE_5_8_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_BE_CHH_BE_5_8_2019.RData")


HD_BE_P1 = c(HD_ref_BE_CG$BE_P1 , HD_ref_BE_CHG$BE_P1 , HD_ref_BE_CHH_BE$BE_P1 )
HD_BE_P2 = c(HD_ref_BE_CG$BE_P2 , HD_ref_BE_CHG$BE_P2 , HD_ref_BE_CHH_BE$BE_P2 )
HD_BE_P3 = c(HD_ref_BE_CG$BE_P3 , HD_ref_BE_CHG$BE_P3 , HD_ref_BE_CHH_BE$BE_P3 )

HD_GE_P1 = c(HD_ref_BE_CG$GE_P1 , HD_ref_BE_CHG$GE_P1 , HD_ref_BE_CHH_GE$GE_P1 )
HD_GE_P2 = c(HD_ref_BE_CG$GE_P2 , HD_ref_BE_CHG$GE_P2 , HD_ref_BE_CHH_GE$GE_P2 )
HD_GE_P3 = c(HD_ref_BE_CG$GE_P3 , HD_ref_BE_CHG$GE_P3 , HD_ref_BE_CHH_GE$GE_P3 )

HD_BE_P1_tv <- HD_BE_P1[,7]
HD_BE_P2_tv <- HD_BE_P2[,7]
HD_BE_P3_tv <- HD_BE_P3[,7]
HD_GE_P1_tv <- HD_GE_P1[,7]
HD_GE_P2_tv <- HD_GE_P2[,7]
HD_GE_P3_tv <- HD_GE_P3[,7]

colnames(mcols(HD_BE_P1_tv)) <- "score" 
colnames(mcols(HD_BE_P2_tv)) <- "score" 
colnames(mcols(HD_BE_P3_tv)) <- "score" 
colnames(mcols(HD_GE_P1_tv)) <- "score" 
colnames(mcols(HD_GE_P2_tv)) <- "score" 
colnames(mcols(HD_GE_P3_tv)) <- "score" 


export.wig(HD_BE_P1_tv, con = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/wig_file/HD_BE_P1_tv.wig", format = "wig")
export.wig(HD_BE_P2_tv, con = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/wig_file/HD_BE_P2_tv.wig", format = "wig")
export.wig(HD_BE_P3_tv, con = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/wig_file/HD_BE_P3_tv.wig", format = "wig")
export.wig(HD_GE_P1_tv, con = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/wig_file/HD_GE_P1_tv.wig", format = "wig")
export.wig(HD_GE_P2_tv, con = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/wig_file/HD_GE_P2_tv.wig", format = "wig")
export.wig(HD_GE_P3_tv, con = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/wig_file/HD_GE_P3_tv.wig", format = "wig")



HD_BE_P1_MethyLevel <- HD_BE_P1[,6]
HD_BE_P2_MethyLevel <- HD_BE_P2[,6]
HD_BE_P3_MethyLevel <- HD_BE_P3[,6]
HD_GE_P1_MethyLevel <- HD_GE_P1[,6]
HD_GE_P2_MethyLevel <- HD_GE_P2[,6]
HD_GE_P3_MethyLevel <- HD_GE_P3[,6]

colnames(mcols(HD_BE_P1_MethyLevel)) <- "score" 
colnames(mcols(HD_BE_P2_MethyLevel)) <- "score" 
colnames(mcols(HD_BE_P3_MethyLevel)) <- "score" 
colnames(mcols(HD_GE_P1_MethyLevel)) <- "score" 
colnames(mcols(HD_GE_P2_MethyLevel)) <- "score" 
colnames(mcols(HD_GE_P3_MethyLevel)) <- "score" 


export.wig(HD_BE_P1_MethyLevel, con = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/wig_file/HD_BE_P1_MethyLevel.wig", format = "wig")
export.wig(HD_BE_P2_MethyLevel, con = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/wig_file/HD_BE_P2_MethyLevel.wig", format = "wig")
export.wig(HD_BE_P3_MethyLevel, con = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/wig_file/HD_BE_P3_MethyLevel.wig", format = "wig")
export.wig(HD_GE_P1_MethyLevel, con = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/wig_file/HD_GE_P1_MethyLevel.wig", format = "wig")
export.wig(HD_GE_P2_MethyLevel, con = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/wig_file/HD_GE_P2_MethyLevel.wig", format = "wig")
export.wig(HD_GE_P3_MethyLevel, con = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/wig_file/HD_GE_P3_MethyLevel.wig", format = "wig")




