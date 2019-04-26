```
cd /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2
```
ls |sed s'#\(.*\)#ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/\1/\1_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/\1/\1_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/\1#'> tomato_bismark_bowtie2.sh
```
set up ts two job a time
``````
 ts -S 2
``````
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/12-P1/12-P1_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/12-P1/12-P1_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/12-P1
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/12-P2/12-P2_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/12-P2/12-P2_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/12-P2
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/12-P3/12-P3_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/12-P3/12-P3_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/12-P3
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/13-P1/13-P1_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/13-P1/13-P1_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/13-P1
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/13-P2/13-P2_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/13-P2/13-P2_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/13-P2
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/13-P3/13-P3_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/13-P3/13-P3_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/13-P3
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/1-P1/1-P1_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/1-P1/1-P1_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/1-P1
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/1-P2/1-P2_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/1-P2/1-P2_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/1-P2
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P1/21-P1_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P1/21-P1_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/21-P1
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P17/21-P17_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P17/21-P17_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/21-P17
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P19/21-P19_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P19/21-P19_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/21-P19
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P22/21-P22_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P22/21-P22_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/21-P22
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P3/21-P3_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P3/21-P3_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/21-P3
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P32/21-P32_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P32/21-P32_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/21-P32
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P34/21-P34_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P34/21-P34_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/21-P34
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P37/21-P37_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P37/21-P37_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/21-P37
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P6/21-P6_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/21-P6/21-P6_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/21-P6
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/2-P1/2-P1_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/2-P1/2-P1_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/2-P1
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/2-P2/2-P2_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/2-P2/2-P2_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/2-P2
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/BE-P1/BE-P1_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/BE-P1/BE-P1_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/BE-P1
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/BE-P2/BE-P2_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/BE-P2/BE-P2_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/BE-P2
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/BE-P3/BE-P3_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/BE-P3/BE-P3_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/BE-P3
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/GE-P1/GE-P1_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/GE-P1/GE-P1_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/GE-P1
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/GE-P2/GE-P2_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/GE-P2/GE-P2_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/GE-P2
ts /usr/local/Bismark/bismark --bam --multicore 4 --bowtie2 -p 2 --phred33-quals --genome_folder /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/ftp.solgenomics.net/genomes/Solanum_lycopersicum/assembly/build_3.00 -1 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/GE-P3/GE-P3_1_val_1.fq.gz -2 /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/trim_glore/GE-P3/GE-P3_2_val_2.fq.gz -o /data5/F15FTSUSAT0747_TOMrcwM/YXD_Bismark_pipeline/bismark-bt2/GE-P3
```