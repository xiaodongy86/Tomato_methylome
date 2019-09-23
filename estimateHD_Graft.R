library(MethylIT)
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_all_four_wt_CG_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/RDR_12_CG_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/RR_13_CG_5_10_2019.RData")

graft_CG <- c(RDR_12_CG,RR_13_CG)
ref_RR_CG <- poolFromGRlist(list(RR_13_CG$RRP1, RR_13_CG$RRP2, RR_13_CG$RRP3), stat = "sum", num.cores = 12L)

HD_graft_CG = estimateDivergence(ref = ref_RR_CG, 
                              indiv = graft_CG, 
                              Bayesian = TRUE, 
                              min.coverage = 4, 
                              high.coverage = 300, 
                              percentile = 0.999, 
                              num.cores = 1L, tasks = 0L, verbose = FALSE )
save(HD_graft_CG,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_graft_CG_5_10_2019.RData")

nlms_graft_CG = nonlinearFitDist(HD_graft_CG, column = 9, num.cores = 1L, verbose = TRUE)
PS_graft_CG = getPotentialDIMP(LR = HD_graft_CG, nlms = nlms_graft_CG, div.col = 9, alpha = 0.05,
                            tv.col = 7, tv.cut = 0.2)
save(nlms_graft_CG, PS_graft_CG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_graft_CG_5_10_2019.RData")
rm(HD_graft_CG)
rm(ref_RR_CG)
rm(graft_CG)

load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_all_four_wt_CHG_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/RDR_12_CHG_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/RR_13_CHG_5_10_2019.RData")

graft_CHG <- c(RDR_12_CHG,RR_13_CHG)
ref_RR_CHG <- poolFromGRlist(list(RR_13_CHG$RRP1, RR_13_CHG$RRP2, RR_13_CHG$RRP3), stat = "sum", num.cores = 12L)

HD_graft_CHG = estimateDivergence(ref = ref_RR_CHG, 
                                 indiv = graft_CHG, 
                                 Bayesian = TRUE, 
                                 min.coverage = 4, 
                                 high.coverage = 300, 
                                 percentile = 0.999, 
                                 num.cores = 1L, tasks = 0L, verbose = FALSE )
save(HD_graft_CHG,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_graft_CHG_5_10_2019.RData")

nlms_graft_CHG = nonlinearFitDist(HD_graft_CHG, column = 9, num.cores = 1L, verbose = TRUE)
PS_graft_CHG = getPotentialDIMP(LR = HD_graft_CHG, nlms = nlms_graft_CHG, div.col = 9, alpha = 0.05,
                               tv.col = 7, tv.cut = 0.2)
save(nlms_graft_CHG, PS_graft_CHG, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_graft_CHG_5_10_2019.RData")
rm(HD_graft_CHG)
rm(ref_RR_CHG)
rm(graft_CHG)

load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/ref_all_four_wt_CHH_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/RDR_12_CHH_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/RR_13_CHH_5_10_2019.RData")

graft_CHH <- c(RDR_12_CHH,RR_13_CHH)
ref_RR_CHH <- poolFromGRlist(list(RR_13_CHH$RRP1, RR_13_CHH$RRP2, RR_13_CHH$RRP3), stat = "sum", num.cores = 12L)

HD_graft_CHH = estimateDivergence(ref = ref_RR_CHH, 
                                  indiv = graft_CHH, 
                                  Bayesian = TRUE, 
                                  min.coverage = 4, 
                                  high.coverage = 300, 
                                  percentile = 0.999, 
                                  num.cores = 8L, tasks = 4L, verbose = FALSE )
save(HD_graft_CHH,  file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/HD_graft_CHH_5_15_2019.RData")

nlms_graft_CHH = nonlinearFitDist(HD_graft_CHH, column = 9, num.cores = 4L, verbose = TRUE)
PS_graft_CHH = getPotentialDIMP(LR = HD_graft_CHH, nlms = nlms_graft_CHH, div.col = 9, alpha = 0.05,
                                tv.col = 7, tv.cut = 0.2)
save(nlms_graft_CHH, PS_graft_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_graft_CHH_5_15_2019.RData")
rm(HD_graft_CHH)
rm(ref_RR_CHH)
rm(graft_CHH)


load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_graft_CG_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_graft_CHG_5_10_2019.RData")
load("/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/PS_nlms_graft_CHH_5_15_2019.RData")


cutpoints = estimateCutPoint(LR = PS_graft_CG, control.names = c("RRP1", "RRP2",  "RRP3"),
                             treatment.names = c( "RDRP1", "RDRP2", "RDRP3"),
                             #treatment.names = c( "M_1_105", "M_1_168", "M_1_179", "M_1_27", "M_1_3"),
                             div.col = 9, verbose = T)

# > cutpoints
# RRP1     RRP2     RRP3
# RDRP1 1.384532 1.384814 1.384465
# RDRP2 1.133270 1.140125 1.287662
# RDRP3 1.383042 1.383042 1.383032
# 
# $auc
# RRP1      RRP2      RRP3
# RDRP1 0.9796265 0.9728403 0.9758743
# RDRP2 0.9674036 0.9587680 0.9297028
# RDRP3 0.9745073 0.9665647 0.9680937
# 
# $accuracy
# RRP1      RRP2      RRP3
# RDRP1 0.9791494 0.9752350 0.9598436
# RDRP2 0.9595165 0.9541404 0.8576724
# RDRP3 0.9771465 0.9728274 0.9564647

DIMPs_CG  = selectDIMP(PS_graft_CG, div.col = 9, cutpoint = 1.5 )
DIMPs_CHG = selectDIMP(PS_graft_CHG, div.col = 9, cutpoint = 1.5)
DIMPs_CHH = selectDIMP(PS_graft_CHH, div.col = 9, cutpoint = 1.5 )
save(DIMPs_CG, DIMPs_CHG, DIMPs_CHH, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/DIMPs_graft_5_24_2019.RData")


CG_groups = evaluateDIMPclass(LR = DIMPs_CG, control.names = c("RRP1", "RRP2", "RRP3"),
                              treatment.names = c("RDRP1", "RDRP2", "RDRP3"),
                              column = c(hdiv = TRUE, TV = TRUE, wprob = TRUE, pos = TRUE),
                              classifier = "pca.qda",
                              n.pc = 4,
                              num.boot = 5,
                              mc.cores = 5,
                              center = TRUE, scale = TRUE,
                              output = "all", prop = 0.6)


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

# set.seed(11133)
# list1 <- sample(gene,1630)
# list2 <- sample(gene,1710)
# IDX <- findOverlaps(list1,list2, type = c("equal"))


seqlevels(gene)
geneUpDown2kb <- GeneUpDownStream(gene, upstream = 2000, downstream = 2000)
seqlevels(DIMPs_CG)

DIMPs_RDRP1 = c(DIMPs_CG$RDRP1, DIMPs_CHG$RDRP1, DIMPs_CHH$RDRP1)
DIMPs_RDRP2 = c(DIMPs_CG$RDRP2, DIMPs_CHG$RDRP2, DIMPs_CHH$RDRP2)
DIMPs_RDRP3 = c(DIMPs_CG$RDRP3, DIMPs_CHG$RDRP3, DIMPs_CHH$RDRP3)
DIMPs_RRP1 = c(DIMPs_CG$RRP1, DIMPs_CHG$RRP1, DIMPs_CHH$RRP1)
DIMPs_RRP2 = c(DIMPs_CG$RRP2, DIMPs_CHG$RRP2, DIMPs_CHH$RRP2)
DIMPs_RRP3 = c(DIMPs_CG$RRP3, DIMPs_CHG$RRP3, DIMPs_CHH$RRP3)

DIMPs_RDRP1_gene2kb = getDIMPatGenes(DIMPs_RDRP1, geneUpDown2kb)
DIMPs_RDRP2_gene2kb = getDIMPatGenes(DIMPs_RDRP2, geneUpDown2kb)
DIMPs_RDRP3_gene2kb = getDIMPatGenes(DIMPs_RDRP3, geneUpDown2kb)
DIMPs_RRP1_gene2kb = getDIMPatGenes(DIMPs_RRP1, geneUpDown2kb)
DIMPs_RRP2_gene2kb = getDIMPatGenes(DIMPs_RRP2, geneUpDown2kb)
DIMPs_RRP3_gene2kb = getDIMPatGenes(DIMPs_RRP3, geneUpDown2kb)

Genes_DIMPs = uniqueGRanges(list(DIMPs_RRP1_gene2kb[, 2], DIMPs_RRP2_gene2kb[, 2], DIMPs_RRP3_gene2kb[, 2],
                                 DIMPs_RDRP1_gene2kb[, 2], DIMPs_RDRP2_gene2kb[, 2], DIMPs_RDRP3_gene2kb[, 2])
                            ,type = "equal", verbose = TRUE,ignore.strand = FALSE )

colnames( mcols(Genes_DIMPs)) <- c("RRP1","RRP2","RRP3","RDRP1","RDRP2","RDRP3")

HITS = findOverlaps(geneUpDown2kb, Genes_DIMPs, type = "equal")
Genes_DIMPs_ID <- gene[queryHits(HITS), ]

dmps = data.frame( mcols( Genes_DIMPs ))
dmps = apply( dmps, 2, as.numeric )
rownames(dmps) <- as.vector(Genes_DIMPs_ID$Alias)

condition = data.frame(condition = factor(c("CT", "CT", "CT", 
                                            "TT", "TT", "TT"),
                                          levels = c("CT", "TT")))

rownames(condition) <- c("RRP1","RRP2","RRP3","RDRP1","RDRP2","RDRP3")
DIMR <- DESeqDataSetFromMatrix(countData = dmps,
                               colData = condition,
                               design = formula( ~ condition ),
                               rowRanges = Genes_DIMPs)
save(DIMR, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/DIMR_graft_5_24_2019.RData")

####  5_24_2019   parameters for 1710 DMGs

DMGs = countTest2(DS = DIMR, num.cores = 4L, minCountPerIndv = 8, 
                  countFilter = TRUE, CountPerBp = 0.005, maxGrpCV = c(0.5, 0.5),
                  FilterLog2FC = TRUE,Minlog2FC = 1, pvalCutOff = 0.05,
                  MVrate = .95, verbose = FALSE)

write.csv(DMGs, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/DMGs_2kb_1710_5_24_2019.csv",row.names = ROWNAMES(DMGs), quote = FALSE)


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
DMGs_1710 <- ROWNAMES(DMGs)

keepers<-c()  
i=1
for(item in DMGs_1710) {
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



blastp_DMGs_1710_geneUpDown2kb_graft <- blastp(query.seq=query.seq, dtb = NULL, seq.name = seq.name, db.fa = db.fa, maxTargetSeqs = 1,
                                               dir.fa = dir.fa, tmp = tmp,  num.cores = 30L, tasks = 30L)

blastp_DMGs_1710_geneUpDown2kb_graft_na.omit <- na.omit(blastp_DMGs_1710_geneUpDown2kb_graft)
blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100 <- blastp_DMGs_1710_geneUpDown2kb_graft_na.omit[blastp_DMGs_1710_geneUpDown2kb_graft_na.omit$score > 100,]

blastpIDs <- unique(blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100$qseqid)
blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df <- as.data.frame(blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100)

unique_match_blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df <-as.data.frame(blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df[1,])  
i=1
for(item in blastpIDs){
  unique_match_blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df[i,] = blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df[blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df$qseqid==item,][1,]; i=i+1}
dim(unique_match_blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df)
#1157


unique_match_blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df$qseqid <- substr(x = unique_match_blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df$qseqid , start = 1, stop = 16)
unique_match_blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df$sseqid <- substr(x = unique_match_blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df$sseqid , start = 1, stop = 9)

write.csv(unique_match_blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/unique_match_blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df.csv")



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
ATn = unique_match_blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df$sseqid
# ath.grn <- compile.grn.from.kegg( "/JOB/Work/KEGG/ath.zip" )
# head( ath.grn )
genesNet = unique( as.vector( ath.grn[ , 1:2 ] ) ) 

test = .neat( alist = list('AT.DMG' = ATn ), blist = ath.bp.gs, network = ath.grn[ , 1:2 ],
              nettype = 'directed', nodes = genesNet, alpha = 0.05 )
sum( test$pvalue <= 0.05 )
# ntest = test[ , ]
ntest = test[( test$nab > 0 & test$pvalue <= 0.05 ), ]
ntest = ntest[ (ntest$conclusion == "Overenrichment" ), ]

write.csv(ntest, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/unique_match_blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_df_NEAT_GO.csv")


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

blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100

GenesGO.Net_name_tomato_ID <- data.frame(GenesGO.Net_name, name = blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100[match(GenesGO.Net_name$name.gene_id, blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100$sseqid),])
dim(GenesGO.Net_name_tomato_ID)


GenesGO.Net_name_tomato_ID$name.qseqid <- substr(x = GenesGO.Net_name_tomato_ID$name.qseqid , start = 1, stop = 14)
GenesGO.Net_name_tomato_ID <- na.omit(GenesGO.Net_name_tomato_ID)

DMGs_GenesGO.Net_name_tomato_ID <- data.frame(GenesGO.Net_name_tomato_ID, name = DMGs[match(GenesGO.Net_name_tomato_ID$name.qseqid,  ROWNAMES(DMGs)),])

write.csv(DMGs_GenesGO.Net_name_tomato_ID, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/DMGs_GenesGO.Net_name_tomato_ID_BE_GE_5_16_2019.csv")


blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_name <- data.frame(blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100, name = AG_anotation_df[match(blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100$sseqid, AG_anotation_df$gene_id),])


blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_name$qseqid <- substr(x = blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_name$qseqid , start = 1, stop = 14)

blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_name_tomato_ID <- data.frame(blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_name, DMGS = DMGs[match(blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_name$qseqid, ROWNAMES(DMGs)),])



write.csv(blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_name_tomato_ID_unique, file = "/data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/methyl_extractor/blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_name_tomato_ID_unique_5_16_2019.csv")
dim(blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_name_tomato_ID)

blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_name_tomato_ID_unique <- unique(blastp_DMGs_1710_geneUpDown2kb_graft_na.omit_score100_name_tomato_ID$qseqid)



