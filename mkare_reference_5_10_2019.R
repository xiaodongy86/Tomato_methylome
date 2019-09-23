#####################################################################################
#                    all four wt as reference
#####################################################################################
ref_all_four_wt_CG <- poolFromGRlist(list(WT_CG$WT1_P1, WT_CG$WT1_P2, WT_CG$WT2_P1,WT_CG$WT2_P2), stat = "sum", num.cores = 12L)

ref_all_four_wt_CHG <- poolFromGRlist(list(WT_CHG$WT1_P1, WT_CHG$WT1_P2, WT_CHG$WT2_P1,WT_CHG$WT2_P2), stat = "sum", num.cores = 12L)

ref_all_four_wt_CHH <- poolFromGRlist(list(WT_CHH$WT1_P1, WT_CHH$WT1_P2, WT_CHH$WT2_P1,WT_CHH$WT2_P2), stat = "sum", num.cores = 12L)

load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_minus_CHH_5_9_2019.RData")

save(ref_all_four_wt_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_all_four_wt_CG_5_10_2019.RData")
save(ref_all_four_wt_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_all_four_wt_CHG_5_10_2019.RData")
save(ref_all_four_wt_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_all_four_wt_CHH_5_10_2019.RData")


HD_WT_CG = estimateDivergence(ref = ref_all_four_wt_CG, 
                                  indiv = WT_CG, 
                                  Bayesian = TRUE, 
                                  min.coverage = 4, 
                                  high.coverage = 300, 
                                  percentile = 0.999, 
                                  num.cores = 4L, tasks = 0L, verbose = FALSE )
save(HD_WT_CG,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_all_four_wt_CG_5_10_2019.RData")

nlms_WT_CG = nonlinearFitDist(HD_WT_CG, column = 9, num.cores = 10L, verbose = TRUE)
PS_WT_CG = getPotentialDIMP(LR = HD_WT_CG, nlms = nlms_WT_CG, div.col = 9, alpha = 0.05,
                                tv.col = 7, tv.cut = 0.2)
save(nlms_WT_CG, PS_WT_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_WT_CG_5_10_2019.RData")
rm(HD_WT_CG)
rm(ref_all_four_wt_CG)
rm(WT_CG)

HD_WT_CHG = estimateDivergence(ref = ref_all_four_wt_CHG, 
                              indiv = WT_CHG, 
                              Bayesian = TRUE, 
                              min.coverage = 4, 
                              high.coverage = 300, 
                              percentile = 0.999, 
                              num.cores = 4L, tasks = 0L, verbose = FALSE )
save(HD_WT_CHG,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_all_four_wt_CHG_5_10_2019.RData")

nlms_WT_CHG = nonlinearFitDist(HD_WT_CHG, column = 9, num.cores = 10L, verbose = TRUE)
PS_WT_CHG = getPotentialDIMP(LR = HD_WT_CHG, nlms = nlms_WT_CHG, div.col = 9, alpha = 0.05,
                            tv.col = 7, tv.cut = 0.2)
save(nlms_WT_CHG, PS_WT_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_WT_CHG_5_10_2019.RData")
rm(HD_WT_CHG)
rm(ref_all_four_wt_CHG)
rm(WT_CHG)

