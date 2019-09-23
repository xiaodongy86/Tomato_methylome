

library(MethylIT)

############################################################################################################################
#                                                         CG                                                               #
############################################################################################################################

#########################################################################################################################
#
#                                     WT   gen1  gen2    1-P1  1-P2  2-P1   2-P2                                        #
#
#########################################################################################################################

samples_WT_CG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/1-P1/1-P1_CG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/1-P2/1-P2_CG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/2-P1/2-P1_CG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/2-P2/2-P2_CG.txt")

sample.id_WT_CG <- c("WT1_P1","WT1_P2","WT2_P1","WT2_P2")

WT_CG <- readCounts2GRangesList(filenames = samples_WT_CG, 
                                sample.id = sample.id_WT_CG,
                                columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                verbose = TRUE)


names(WT_CG) <- c("WT1_P1","WT1_P2","WT2_P1","WT2_P2")
save(WT_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/WT_CG_5_9_2019.RData")

#########################################################################################################################
#
#                                     DR   mild_minus   21-P1  21-P22  21-P32                                           #
#
#########################################################################################################################

samples_21_mild_minus_CG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P1/21-P1_CG.txt",
                              "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P22/21-P22_CG.txt",
                              "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P32/21-P32_CG.txt")

samples.id_21_mild_minus_CG <- c("Dr21_P1","Dr21_P22","Dr21_P32")

mild_minus_CG <- readCounts2GRangesList(filenames = samples_21_mild_minus_CG, 
                                sample.id = samples.id_21_mild_minus_CG,
                                columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                verbose = TRUE)

save(mild_minus_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_minus_CG_5_9_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_minus_CG_5_9_2019.RData")
names(mild_minus_CG) <- c("Dr21_P1","Dr21_P2","Dr21_P3")



#########################################################################################################################
#
#                                     DR   mild_plus  21-P3  21-P6  21-P34                                              #
#
#########################################################################################################################

samples_21_mild_plus_CG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P3/21-P3_CG.txt",
                              "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P6/21-P6_CG.txt",
                              "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P34/21-P34_CG.txt")

samples.id_21_mild_plus_CG <- c("Dr21_P3","Dr21_P6","Dr21_P34")

mild_plus_CG <- readCounts2GRangesList(filenames = samples_21_mild_plus_CG, 
                                        sample.id = samples.id_21_mild_plus_CG,
                                        columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                        verbose = TRUE)

save(mild_plus_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_plus_CG_5_9_2019.RData")

load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_plus_CG_5_9_2019.RData")


#########################################################################################################################
#
#                                     grafting  grafted gen1 R/DR  12-P1  12-P2  12-P3                                  #
#
#########################################################################################################################

samples_12_CG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/12-P1/12-P1_CG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/12-P2/12-P2_CG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/12-P3/12-P3_CG.txt")

samples.id_12_RDR_CG <- c("RDRP1","RDRP2","RDRP3")

RDR_12_CG <- readCounts2GRangesList(filenames = samples_12_CG, 
                                    sample.id = samples.id_12_RDR_CG,
                                    columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                    verbose = TRUE)

save(RDR_12_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/RDR_12_CG_5_10_2019.RData")



#########################################################################################################################
#
#                                     grafting  grafted gen1 R/R  13-P1  13-P2  13-P3                                    #
#
#########################################################################################################################

samples_13_CG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/13-P1/13-P1_CG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/13-P2/13-P2_CG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/13-P3/13-P3_CG.txt")

samples.id_13_RR_CG <- c("RRP1","RRP2","RRP3")

RR_13_CG <- readCounts2GRangesList(filenames = samples_13_CG, 
                                    sample.id = samples.id_13_RR_CG,
                                    columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                    verbose = TRUE)

save(RR_13_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/RR_13_CG_5_10_2019.RData")





#############################################################################################################################
#                                                         CHG
#############################################################################################################################

#########################################################################################################################
#
#                                     WT   gen1  gen2    1-P1  1-P2  2-P1   2-P2                                        #
#
#########################################################################################################################


samples_WT_CHG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/1-P1/1-P1_CHG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/1-P2/1-P2_CHG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/2-P1/2-P1_CHG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/2-P2/2-P2_CHG.txt")

sample.id_WT_CHG <- c("1_P1","1_P2","2_P1","2_P2")

WT_CHG <- readCounts2GRangesList(filenames = samples_WT_CHG, 
                                sample.id = sample.id_WT_CHG,
                                columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                verbose = TRUE)

save(WT_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/WT_CHG_5_9_2019.RData")


#########################################################################################################################
#
#                                     DR   mild_minus   21-P1  21-P22  21-P32                                           #
#
#########################################################################################################################

#  21-P1  21-P22  21-P32  
samples_21_mild_minus_CHG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P1/21-P1_CHG.txt",
                              "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P22/21-P22_CHG.txt",
                              "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P32/21-P32_CHG.txt")

samples.id_21_mild_minus_CHG <- c("Dr21_P1","Dr21_P22","Dr21_P32")

mild_minus_CHG <- readCounts2GRangesList(filenames = samples_21_mild_minus_CHG, 
                                        sample.id = samples.id_21_mild_minus_CHG,
                                        columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                        verbose = TRUE)

save(mild_minus_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_minus_CHG_5_9_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_minus_CHG_5_9_2019.RData")
names(mild_minus_CHG) <- c("Dr21_P1","Dr21_P22","Dr21_P32")


#########################################################################################################################
#
#                                     DR   mild_plus  21-P3  21-P6  21-P34                                              #
#
#########################################################################################################################

#  21-P3  21-P6  21-P34  
samples_21_mild_plus_CHG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P3/21-P3_CHG.txt",
                             "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P6/21-P6_CHG.txt",
                             "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P34/21-P34_CHG.txt")

samples.id_21_mild_plus_CHG <- c("Dr21_P3","Dr21_P6","Dr21_P34")

mild_plus_CHG <- readCounts2GRangesList(filenames = samples_21_mild_plus_CHG, 
                                       sample.id = samples.id_21_mild_plus_CHG,
                                       columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                       verbose = TRUE)

save(mild_plus_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_plus_CHG_5_9_2019.RData")

#########################################################################################################################
#
#                                     DR   dwarf_plus  21-P17  21-P19  21-P37                                          #
#
#########################################################################################################################
#  21-P17  21-P19  21-P37  
samples_21_dwarf_plus_CHG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P17/21-P17_CHG.txt",
                              "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P19/21-P19_CHG.txt",
                              "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P37/21-P37_CHG.txt")

samples.id_21_dwarf_plus_CHG <- c("Dr21_P17","Dr21_P19","Dr21_P37")

dwarf_plus_CHG <- readCounts2GRangesList(filenames = samples_21_dwarf_plus_CHG, 
                                        sample.id = samples.id_21_dwarf_plus_CHG,
                                        columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                        verbose = TRUE)

save(dwarf_plus_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/dwarf_plus_CHG_5_9_2019.RData")


#########################################################################################################################
#
#                                     grafting  grafted gen1 R/DR  12-P1  12-P2  12-P3                                  #
#
#########################################################################################################################

samples_12_CHG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/12-P1/12-P1_CHG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/12-P2/12-P2_CHG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/12-P3/12-P3_CHG.txt")

samples.id_12_RDR_CHG <- c("RDRP1","RDRP2","RDRP3")

RDR_12_CHG <- readCounts2GRangesList(filenames = samples_12_CHG, 
                                    sample.id = samples.id_12_RDR_CHG,
                                    columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                    verbose = TRUE)

save(RDR_12_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/RDR_12_CHG_5_10_2019.RData")



#########################################################################################################################
#
#                                     grafting  grafted gen1 R/R  13-P1  13-P2  13-P3                                    #
#
#########################################################################################################################

samples_13_CHG <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/13-P1/13-P1_CHG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/13-P2/13-P2_CHG.txt",
                   "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/13-P3/13-P3_CHG.txt")

samples.id_13_RR_CHG <- c("RRP1","RRP2","RRP3")

RR_13_CHG <- readCounts2GRangesList(filenames = samples_13_CHG, 
                                   sample.id = samples.id_13_RR_CHG,
                                   columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                   verbose = TRUE)

save(RR_13_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/RR_13_CHG_5_10_2019.RData")


#############################################################################################################################
#                                                         CHH
#############################################################################################################################

#########################################################################################################################
#
#                                     WT   gen1  gen2    1-P1  1-P2  2-P1   2-P2                                        #
#
#########################################################################################################################


samples_WT_CHH <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/1-P1/1-P1_CHH.txt",
                    "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/1-P2/1-P2_CHH.txt",
                    "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/2-P1/2-P1_CHH.txt",
                    "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/2-P2/2-P2_CHH.txt")

sample.id_WT_CHH <- c("WT1_P1","WT1_P2","WT2_P1","WT2_P2")

WT_CHH <- readCounts2GRangesList(filenames = samples_WT_CHH, 
                                 sample.id = sample.id_WT_CHH,
                                 columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                 verbose = TRUE)

save(WT_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/WT_CHH_5_9_2019.RData")


#########################################################################################################################
#
#                                     DR   mild_minus   21-P1  21-P22  21-P32                                           #
#
#########################################################################################################################

#  21-P1  21-P22  21-P32  
samples_21_mild_minus_CHH <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P1/21-P1_CHH.txt",
                               "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P22/21-P22_CHH.txt",
                               "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P32/21-P32_CHH.txt")

samples.id_21_mild_minus_CHH <- c("Dr21_P1","Dr21_P22","Dr21_P32")

mild_minus_CHH <- readCounts2GRangesList(filenames = samples_21_mild_minus_CHH, 
                                         sample.id = samples.id_21_mild_minus_CHH,
                                         columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                         verbose = TRUE)
save(mild_minus_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_minus_CHH_5_9_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_minus_CHH_5_9_2019.RData")
names(mild_minus_CHH) <- c("Dr21_P1","Dr21_P22","Dr21_P32")

########################################################################################################################
#
#                                     DR   mild_plus  21-P3  21-P6  21-P34                                              #
#
#########################################################################################################################

#  21-P3  21-P6  21-P34  
samples_21_mild_plus_CHH <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P3/21-P3_CHH.txt",
                              "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P6/21-P6_CHH.txt",
                              "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P34/21-P34_CHH.txt")

samples.id_21_mild_plus_CHH <- c("Dr21_P3","Dr21_P6","Dr21_P34")

mild_plus_CHH <- readCounts2GRangesList(filenames = samples_21_mild_plus_CHH, 
                                        sample.id = samples.id_21_mild_plus_CHH,
                                        columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                        verbose = TRUE)

save(mild_plus_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/mild_plus_CHH_5_9_2019.RData")

#########################################################################################################################
#
#                                     DR   dwarf_plus  21-P17  21-P19  21-P37                                          #
#
#########################################################################################################################

#  21-P17  21-P19  21-P37  
samples_21_dwarf_plus_CHH <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P17/21-P17_CHH.txt",
                               "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P19/21-P19_CHH.txt",
                               "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/21-P37/21-P37_CHH.txt")

samples.id_21_dwarf_plus_CHH <- c("Dr21_P17","Dr21_P19","Dr21_P37")

dwarf_plus_CHH <- readCounts2GRangesList(filenames = samples_21_dwarf_plus_CHH, 
                                         sample.id = samples.id_21_dwarf_plus_CHH,
                                         columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                         verbose = TRUE)

save(dwarf_plus_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/dwarf_plus_CHH_5_9_2019.RData")

#########################################################################################################################
#
#                                     grafting  grafted gen1 R/DR  12-P1  12-P2  12-P3                                  #
#
#########################################################################################################################

samples_12_CHH <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/12-P1/12-P1_CHH.txt",
                    "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/12-P2/12-P2_CHH.txt",
                    "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/12-P3/12-P3_CHH.txt")

samples.id_12_RDR_CHH <- c("RDRP1","RDRP2","RDRP3")

RDR_12_CHH <- readCounts2GRangesList(filenames = samples_12_CHH, 
                                     sample.id = samples.id_12_RDR_CHH,
                                     columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                     verbose = TRUE)

save(RDR_12_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/RDR_12_CHH_5_10_2019.RData")



#########################################################################################################################
#
#                                     grafting  grafted gen1 R/R  13-P1  13-P2  13-P3                                    #
#
#########################################################################################################################

samples_13_CHH <- c("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/13-P1/13-P1_CHH.txt",
                    "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/13-P2/13-P2_CHH.txt",
                    "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/13-P3/13-P3_CHH.txt")

samples.id_13_RR_CHH <- c("RRP1","RRP2","RRP3")

RR_13_CHH <- readCounts2GRangesList(filenames = samples_13_CHH, 
                                    sample.id = samples.id_13_RR_CHH,
                                    columns = c(seqnames = 1, start = 2, strand = 3, mC = 4, uC = 5),
                                    verbose = TRUE)

save(RR_13_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/RR_13_CHH_5_10_2019.RData")





