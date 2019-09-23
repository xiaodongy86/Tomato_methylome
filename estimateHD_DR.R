library(MethylIT)
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_all_four_wt_CG_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_minus_CG_5_9_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_plus_CG_5_9_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/dwarf_plus_CG_5_9_2019.RData")

DR_CG <- c(mild_plus_CG,mild_minus_CG,dwarf_plus_CG)

HD_dr_CG = estimateDivergence(ref = ref_all_four_wt_CG, 
                              indiv = DR_CG, 
                              Bayesian = TRUE, 
                              min.coverage = 4, 
                              high.coverage = 300, 
                              percentile = 0.999, 
                              num.cores = 2L, tasks = 0L, verbose = FALSE )
save(HD_dr_CG,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_DR_ref_all_four_wt_CG_5_10_2019.RData")

nlms_dr_CG = nonlinearFitDist(HD_dr_CG, column = 9, num.cores = 4L, verbose = TRUE)
PS_dr_CG = getPotentialDIMP(LR = HD_dr_CG, nlms = nlms_dr_CG, div.col = 9, alpha = 0.05,
                            tv.col = 7, tv.cut = 0.2)
save(nlms_dr_CG, PS_dr_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_dr_CG_5_10_2019.RData")

rm(HD_dr_CG)
rm(ref_all_four_wt_CG)
rm(DR_CG)
rm(mild_plus_CG)
rm(mild_minus_CG)
rm(dwarf_plus_CG)


load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_all_four_wt_CHG_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_minus_CHG_5_9_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_plus_CHG_5_9_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/dwarf_plus_CHG_5_9_2019.RData")

DR_CHG <- c(mild_plus_CHG,mild_minus_CHG,dwarf_plus_CHG)

HD_dr_CHG = estimateDivergence(ref = ref_all_four_wt_CHG, 
                              indiv = DR_CHG, 
                              Bayesian = TRUE, 
                              min.coverage = 4, 
                              high.coverage = 300, 
                              percentile = 0.999, 
                              num.cores = 2L, tasks = 0L, verbose = FALSE )
save(HD_dr_CHG,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_DR_ref_all_four_wt_CHG_5_10_2019.RData")

nlms_dr_CHG = nonlinearFitDist(HD_dr_CHG, column = 9, num.cores = 4L, verbose = TRUE)
PS_dr_CHG = getPotentialDIMP(LR = HD_dr_CHG, nlms = nlms_dr_CHG, div.col = 9, alpha = 0.05,
                            tv.col = 7, tv.cut = 0.2)
save(nlms_dr_CHG, PS_dr_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_dr_CHG_5_10_2019.RData")

rm(HD_dr_CHG)
rm(ref_all_four_wt_CHG)
rm(DR_CHG)
rm(mild_plus_CHG)
rm(mild_minus_CHG)
rm(dwarf_plus_CHG)


load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_all_four_wt_CHH_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_minus_CHH_5_9_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_plus_CHH_5_9_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/dwarf_plus_CHH_5_9_2019.RData")

DR_CHH <- c(mild_plus_CHH,mild_minus_CHH,dwarf_plus_CHH)

HD_dr_CHH = estimateDivergence(ref = ref_all_four_wt_CHH, 
                               indiv = DR_CHH, 
                               Bayesian = TRUE, 
                               min.coverage = 4, 
                               high.coverage = 300, 
                               percentile = 0.999, 
                               num.cores = 2L, tasks = 0L, verbose = FALSE )
save(HD_dr_CHH,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_DR_ref_all_four_wt_CHH_5_10_2019.RData")

nlms_dr_CHH = nonlinearFitDist(HD_dr_CHH, column = 9, num.cores = 4L, verbose = TRUE)
PS_dr_CHH = getPotentialDIMP(LR = HD_dr_CHH, nlms = nlms_dr_CHH, div.col = 9, alpha = 0.05,
                             tv.col = 7, tv.cut = 0.2)
save(nlms_dr_CHH, PS_dr_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_dr_CHH_5_10_2019.RData")

rm(HD_dr_CHH)
rm(ref_all_four_wt_CHH)
rm(DR_CHH)
rm(mild_plus_CHH)
rm(mild_minus_CHH)
rm(dwarf_plus_CHH)





