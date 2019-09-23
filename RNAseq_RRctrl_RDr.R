## ##################################################################### ##
############################# ---RNA-Seq --- ##############################
## ##################################################################### ##
#
#### **** ---RNA_Seq data. RNAi NULL & WT ------**** ####
#
################################################ ##
install.packages("hR") 
library(hR)
library(DESeq2)
library(genefilter)
setwd("/data5/tomato_RNAseq/QoRTs_out")

sampleFiles <- c("RRCtlP1_QC.geneCounts.formatted.for.DESeq.txt.gz",
                 "RRCtlP2_QC.geneCounts.formatted.for.DESeq.txt.gz",
                 "RRCtlP3_QC.geneCounts.formatted.for.DESeq.txt.gz",
                 "RDrP1_QC.geneCounts.formatted.for.DESeq.txt.gz",
                 "RDrP2_QC.geneCounts.formatted.for.DESeq.txt.gz",   
                 "RDrP3_QC.geneCounts.formatted.for.DESeq.txt.gz")

sampleName <- c("RRCtlP1","RRCtlP2","RRCtlP3","RDrP1","RDrP2","RDrP3")


sampleCondition <- factor( c( "RRCtl", "RRCtl", "RRCtl",
                              "RDr", "RDr", "RDr"),
                           levels = c("RRCtl", "RDr"))

sampleTable <- data.frame( sampleName = sampleName,
                           fileName = sampleFiles,
                           condition = sampleCondition )
# sampleName                                         fileName condition
# 1    RRCtlP1 RRCtlP1_QC.geneCounts.formatted.for.DESeq.txt.gz     RRCtl
# 2    RRCtlP2 RRCtlP2_QC.geneCounts.formatted.for.DESeq.txt.gz     RRCtl
# 3    RRCtlP3 RRCtlP3_QC.geneCounts.formatted.for.DESeq.txt.gz     RRCtl
# 4      RDrP1   RDrP1_QC.geneCounts.formatted.for.DESeq.txt.gz       RDr
# 5      RDrP2   RDrP2_QC.geneCounts.formatted.for.DESeq.txt.gz       RDr
# 6      RDrP3   RDrP3_QC.geneCounts.formatted.for.DESeq.txt.gz       RDr


deg.rnai <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable,
                                        directory = "./",
                                        design= ~ condition)

assay(deg.rnai)
colData(deg.rnai)
rowData(deg.rnai)







# DEGs = countTest2(DS = deg.rnai, num.cores = 4L, minCountPerIndv = 3, 
#                    countFilter = T, maxGrpCV = c(0.5, 0.5),
#                    FilterLog2FC = TRUE,Minlog2FC = 0.5, pvalCutOff = 0.05,
#                    MVrate = .95, verbose = FALSE, test = "Wald")

# parameters for 1026 DEGs
DEGs = countTest2(DS = deg.rnai, num.cores = 4L, minCountPerIndv = 20, 
                  countFilter = T, maxGrpCV = c(0.5, 0.5),
                  FilterLog2FC = TRUE,Minlog2FC = 0.5, pvalCutOff = 0.05,
                  MVrate = .95, verbose = FALSE, test = "Wald")

##############################################################################################

#                                            DESeq
##############################################################################################

# Pre-filtering
# deg.rnai <- deg.rnai[ rowSums( counts( deg.rnai ) ) > 1, ]

######################################## # pheatmap
#### Differential expression analysis ####
######################################## # SummarizedExperiment
deg.rnai <- estimateSizeFactors( deg.rnai, type = "iterate", locfunc = shorth )

#deg.rnai_fpkm <- fpkm(deg.rnai,robust = TRUE)

# deg.rnai <- estimateSizeFactors( deg.rnai )
deg.rnai <- estimateDispersions( deg.rnai, fitType = "local", maxit = 10^8 )
idx <- rowSums( counts( deg.rnai, normalized = TRUE ) >= 10 ) >= 3
deg.rnai = deg.rnai[ idx ]

deg.rnai_test <- nbinomWaldTest( deg.rnai, betaPrior = TRUE, useOptim = FALSE, maxit = 10^8 )

deg.rnai_test <- DESeq(deg.rnai)
deg.rnai_test_result <- results(deg.rnai_test)

deg.rnai_test_result
deg.rnai_test_result_1 = deg.rnai_test_result[ abs(deg.rnai_test_result$log2FoldChange) > 0.5, ]
deg.rnai_test_result_1 <- na.omit(deg.rnai_test_result_1)
deg.rnai_test_result_2 = deg.rnai_test_result_1[ abs(deg.rnai_test_result_1$padj <= 0.1), ]
unique(deg.rnai_test_result_2)

plotDispEsts( deg.rnai_test )
# Extract results from a DESeq analysis with filter rv
rv <- rowVars( counts( deg.rnai_test, normalized = TRUE ) ) # row variance

deg.rnai_test_1 <- results( deg.rnai_test, filter = rv , pAdjustMethod = "BH" ) #
deg.rnai_test_2 = deg.rnai_test_1[ !is.na( deg.rnai_test_1$padj ), ]
deg.rnai_test_3 = deg.rnai_test_2[ deg.rnai_test_2$padj <= 0.1, ]
deg.rnai_test_4 = deg.rnai_test_3[ abs(deg.rnai_test_3$log2FoldChange) > 0.5, ]

deg.rnai_test_4 = deg.rnai_test_3[ deg.rnai_test_3$padj <= 0.05, ]


DEG.MM = DEG.MM[ abs(DEG.MM$log2FoldChange) > 0.5, ] # 
DEG.MM$adj.pval =  p.adjust( DEG.MM$pvalue, method = "BH" )
#DEG.MM = DEG.MM[ DEG.MM$adj.pval <= 0.05, ] # 1368
DEG.MM = DEG.MM[ DEG.MM$adj.pval <= 0.1, ]
DEG.MM

# sum( abs(DEG.MM$log2FoldChange) > 1, na.rm = T ) # 157
# sum( abs(DEG.MM$log2FoldChange) > 0.5, na.rm = T ) # 864
# min( abs(DEG.MM$log2FoldChange) ) # 0.2818615
# summary(abs(DEG.MM$log2FoldChange) )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2819  0.4417  0.5772  0.6547  0.7685  3.0314 

save( deg.rnai, DEG.MM, file = "/JOB/Work/EpisTDNA/RNASeq/RNAiNull/DEGs/DEGs_WT_vs_MM_10-31-17.RData" )
write.csv( DEG.MM, file = "/JOB/Work/EpisTDNA/RNASeq/RNAiNull/DEGs/DEGs_WT_vs_MM_10-31-17.csv" )



Tomato_ITAG3.0_gene_models <- import("/data/users/xzy50/ITAG3.0_gene_models.gff")
head(Tomato_ITAG3.0_gene_models)
seqlevels(Tomato_ITAG3.0_gene_models)
gene = Tomato_ITAG3.0_gene_models[ Tomato_ITAG3.0_gene_models$type == "gene", c("ID","Alias","Note","Ontology_term")]
seqlevels( gene ) <- c( "0","1","2","3","4","5","6","7","8","9","10","11","12" )
seqlevels(gene, pruning.mode = "coarse") <- c("0","1","2","3","4","5","6","7","8","9","10","11","12" )
gene = sortBySeqnameAndStart(gene)

library(seqinr)
library(Biostrings)
library(BiocParallel)
db.fa <- "/Araport11_genes.201606.pep.fasta"
dir.fa <- "/data/users/xzy50/TomatoBlast"
dir.db <- "/data/users/xzy50/TomatoBlast"

### set up the query sequence 
file <- "/data/users/xzy50/TomatoBlast/ITAG3.0_proteins.fasta"
seqs <- readAAStringSet(filepath = file, format = "fasta")

DEGs_RR_RDR <- rownames(DEGs)

tmp <- dir.db
#grep(DMG_13_vs_12_CG_208, names(seqs))
#?AAStringSet
keepers<-c()  
i=1
for(item in DEGs_RR_RDR) {
  tempy <- grep(item,names(seqs))
  if(is.integer(tempy)){
    keepers[i]<-tempy
    i=i+1
  }else{
    print("problem, is.integer is not an integer")
  }
}

query.seq <- seqs[na.omit(keepers)]

# query.seq <- seqs[1:20]
seq.name  <- substr(x = names(query.seq), start = 1, stop = 16)


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


blastp_DEGs_RR_RDR <- blastp(query.seq=query.seq, dtb = NULL, seq.name = seq.name, db.fa = db.fa, maxTargetSeqs = 1,
                                                        dir.fa = dir.fa, tmp = tmp,  num.cores = 30L, tasks = 30L)
sum(is.na(blastp_DEGs_RR_RDR))
blastp_DEGs_RR_RDR_na.omit <-na.omit(blastp_DEGs_RR_RDR)

blastp_DEGs_RR_RDR_na.omit_score_100 = blastp_DEGs_RR_RDR_na.omit[blastp_DEGs_RR_RDR_na.omit$score >100, ] 



blastp_DEGs_RR_RDR_na.omit_score_100$qseqid <- substr(x = blastp_DEGs_RR_RDR_na.omit_score_100$qseqid , start = 1, stop = 16)
DEGs_RR_RDR_na.omit_score_100_total_1096<- substr(x = blastp_DEGs_RR_RDR_na.omit_score_100$sseqid , start = 1, stop = 9)
blastp_DEGs_RR_RDR_na.omit_score_100$sseqid <- substr(x = blastp_DEGs_RR_RDR_na.omit_score_100$sseqid , start = 1, stop = 9)

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
ATn = DEGs_RR_RDR_na.omit_score_100_total_1096
# ath.grn <- compile.grn.from.kegg( "/JOB/Work/KEGG/ath.zip" )
# head( ath.grn )
genesNet = unique( as.vector( ath.grn[ , 1:2 ] ) ) 

test = .neat( alist = list('AT.DEG' = ATn ), blist = ath.bp.gs, network = ath.grn[ , 1:2 ],
              nettype = 'directed', nodes = genesNet, alpha = 0.05 )
sum( test$pvalue <= 0.05 )
# ntest = test[ , ]
ntest = test[( test$nab > 0 & test$pvalue <= 0.05 ), ]
ntest = ntest[ (ntest$conclusion == "Overenrichment" ), ]

write.csv(ntest, file = "/data5/tomato_RNAseq/NEAT_GO_blastp_DEGs_RR_RDR_na.omit_score_100_total_1096.csv")


GO.NetAB = ath.bp.gs[ match( as.character( ntest$B ), names( ath.bp.gs ) ) ]
GeneGO.Net = sapply( as.list( GO.NetAB ), function(s) s[ na.omit( match( ATn, s ) ) ] ) 

# [1] "GO:0000303_response_to_superoxide"                                                                       
# [2] "GO:0006355_regulation_of_transcription,_DNA-templated"                                                   
# [3] "GO:0006470_protein_dephosphorylation"                                                                    
# [4] "GO:0008219_cell_death"                                                                                   
# [5] "GO:0009410_response_to_xenobiotic_stimulus"                                                              
# [6] "GO:0009414_response_to_water_deprivation"                                                                
# [7] "GO:0009610_response_to_symbiotic_fungus"                                                                 
# [8] "GO:0009625_response_to_insect"                                                                           
# [9] "GO:0009723_response_to_ethylene"                                                                         
# [10] "GO:0009737_response_to_abscisic_acid"                                                                    
# [11] "GO:0009738_abscisic_acid-activated_signaling_pathway"                                                    
# [12] "GO:0009755_hormone-mediated_signaling_pathway"                                                           
# [13] "GO:0009788_negative_regulation_of_abscisic_acid-activated_signaling_pathway"                             
# [14] "GO:0009871_jasmonic_acid_and_ethylene-dependent_systemic_resistance,_ethylene_mediated_signaling_pathway"
# [15] "GO:0009873_ethylene-activated_signaling_pathway"                                                         
# [16] "GO:0009909_regulation_of_flower_development"                                                             
# [17] "GO:0009933_meristem_structural_organization"                                                             
# [18] "GO:0010082_regulation_of_root_meristem_growth"                                                           
# [19] "GO:0010105_negative_regulation_of_ethylene-activated_signaling_pathway"                                  
# [20] "GO:0010182_sugar_mediated_signaling_pathway"                                                             
# [21] "GO:0016567_protein_ubiquitination"                                                                       
# [22] "GO:0019915_lipid_storage"                                                                                
# [23] "GO:0030968_endoplasmic_reticulum_unfolded_protein_response"                                              
# [24] "GO:0042538_hyperosmotic_salinity_response"                                                               
# [25] "GO:0045893_positive_regulation_of_transcription,_DNA-templated"                                          
# [26] "GO:0050826_response_to_freezing"                                                                         
# [27] "GO:0071281_cellular_response_to_iron_ion"                                                                
# [28] "GO:0071368_cellular_response_to_cytokinin_stimulus"                                                      
# [29] "GO:0080113_regulation_of_seed_growth"     

Netw = names( GeneGO.Net )
GenesGO.Net = c()
for( k in 1:length( GeneGO.Net ) ) {
  if( length( GeneGO.Net[[ k ]]) > 0 ) 
    GenesGO.Net = rbind( GenesGO.Net, data.frame( GeneID = GeneGO.Net[[ k ]], Network = Netw[ k ] ) )
}

GenesGO.Net = data.table( GenesGO.Net )
GenesGO.Nets = GenesGO.Net[ , list( Network = list( unique( as.character( Network ) ) ) ), by = GeneID ]
#DEGlist = rownames( DMG )
DEGlist = DEGs_RR_RDR_na.omit_score_100_total_1096


### Anotation
AG_gff3 = import(con = paste0("ftp://ftp.ensemblgenomes.org/pub/release-38/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.38.gff3.gz"))
# or AG_gff3 = import("/data/TAIR10_gff3/Arabidopsis_thaliana.TAIR10.38.gff3.gz")
AG_anotation = AG_gff3[ AG_gff3$type == "gene",  c( "gene_id", "Name", "description" ) ]
seqlevels( AG_anotation ) <- c( "1","2","3","4","5","M","C" )
seqlevels(AG_anotation, pruning.mode = "coarse") <- c("1", "2", "3", "4", "5")

AG_anotation_df <- as.data.frame(mcols(AG_anotation))
GenesGO.Net_name <- data.frame(GenesGO.Net, name = AG_anotation_df[match(GenesGO.Net$GeneID, AG_anotation_df$gene_id),])
head(GenesGO.Net_name)
write.csv(GenesGO.Net_name, file = "/data5/tomato_RNAseq/GenesGO.Net_name_NEAT_GO_GENES_5_15_2019.csv")

blastp_BE_GE_CG_CHG_CHH_gr_geneUpDown2kb_1376

GenesGO.Net_name_tomato_ID <- data.frame(GenesGO.Net_name, name = blastp_DEGs_RR_RDR_na.omit_score_100[match(GenesGO.Net_name$name.gene_id, blastp_DEGs_RR_RDR_na.omit_score_100$sseqid),])
dim(GenesGO.Net_name_tomato_ID)

DEGs_GenesGO.Net_name_tomato_ID <- data.frame(GenesGO.Net_name_tomato_ID, name = DEGs[match(GenesGO.Net_name_tomato_ID$name.qseqid, rownames(DEGs)),])
write.csv(DEGs_GenesGO.Net_name_tomato_ID, file = "/data5/tomato_RNAseq/DEGs_GenesGO.Net_name_tomato_ID_countest2_FC0.5_5_10_2019.csv")




blastp_DEGs_RR_RDR_na.omit_score_100_name <- data.frame(blastp_DEGs_RR_RDR_na.omit_score_100, name = AG_anotation_df[match(blastp_DEGs_RR_RDR_na.omit_score_100$sseqid, AG_anotation_df$gene_id),])
blastp_DEGs_RR_RDR_na.omit_score_100_name_tomato_ID <- data.frame(blastp_DEGs_RR_RDR_na.omit_score_100_name, name = DEGs[match(blastp_DEGs_RR_RDR_na.omit_score_100_name$qseqid, rownames(DEGs)),])
write.csv(blastp_DEGs_RR_RDR_na.omit_score_100_name_tomato_ID, file = "/data5/tomato_RNAseq/DEGs_total_1096_name_tomato_ID_countest2_FC0.5_5_10_2019.csv")
dim(blastp_DEGs_RR_RDR_na.omit_score_100_name_tomato_ID)

head(GenesGO.Net_name_tomato_ID)

