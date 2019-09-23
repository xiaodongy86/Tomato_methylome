library(MethylIT)

load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_all_four_wt_CG_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/WT_CG_5_9_2019.RData")

HD_WT_CG = estimateDivergence(ref = ref_all_four_wt_CG, 
                              indiv = WT_CG, 
                              Bayesian = TRUE, 
                              min.coverage = 4, 
                              high.coverage = 300, 
                              percentile = 0.999, 
                              num.cores = 2L, tasks = 0L, verbose = FALSE )
save(HD_WT_CG,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_all_four_wt_CG_5_10_2019.RData")

nlms_WT_CG = nonlinearFitDist(HD_WT_CG, column = 9, num.cores = 4L, verbose = TRUE)
PS_WT_CG = getPotentialDIMP(LR = HD_WT_CG, nlms = nlms_WT_CG, div.col = 9, alpha = 0.05,
                            tv.col = 7, tv.cut = 0.2)
save(nlms_WT_CG, PS_WT_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_WT_CG_5_10_2019.RData")
rm(HD_WT_CG)
rm(ref_all_four_wt_CG)
rm(WT_CG)


load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_all_four_wt_CHG_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/WT_CHG_5_9_2019.RData")

HD_WT_CHG = estimateDivergence(ref = ref_all_four_wt_CHG, 
                               indiv = WT_CHG, 
                               Bayesian = TRUE, 
                               min.coverage = 4, 
                               high.coverage = 300, 
                               percentile = 0.999, 
                               num.cores = 2L, tasks = 0L, verbose = FALSE )
save(HD_WT_CHG,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_all_four_wt_CHG_5_10_2019.RData")

nlms_WT_CHG = nonlinearFitDist(HD_WT_CHG, column = 9, num.cores = 4L, verbose = TRUE)
PS_WT_CHG = getPotentialDIMP(LR = HD_WT_CHG, nlms = nlms_WT_CHG, div.col = 9, alpha = 0.05,
                             tv.col = 7, tv.cut = 0.2)
save(nlms_WT_CHG, PS_WT_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_WT_CHG_5_10_2019.RData")
rm(HD_WT_CHG)
rm(ref_all_four_wt_CHG)
rm(WT_CHG)


load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_all_four_wt_CHH_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/WT_CHG_5_9_2019.RData")

HD_WT_CHH = estimateDivergence(ref = ref_all_four_wt_CHH, 
                               indiv = WT_CHH, 
                               Bayesian = TRUE, 
                               min.coverage = 4, 
                               high.coverage = 300, 
                               percentile = 0.999, 
                               num.cores = 2L, tasks = 0L, verbose = FALSE )
save(HD_WT_CHH,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_ref_all_four_wt_CHH_5_10_2019.RData")

nlms_WT_CHH = nonlinearFitDist(HD_WT_CHH, column = 9, num.cores = 4L, verbose = TRUE)
PS_WT_CHH = getPotentialDIMP(LR = HD_WT_CHH, nlms = nlms_WT_CHH, div.col = 9, alpha = 0.05,
                             tv.col = 7, tv.cut = 0.2)
save(nlms_WT_CHH, PS_WT_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_WT_CHH_5_10_2019.RData")
rm(HD_WT_CHH)
rm(ref_all_four_wt_CHH)
rm(WT_CHH)

