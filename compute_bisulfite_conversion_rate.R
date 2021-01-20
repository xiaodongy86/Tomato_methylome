```r
library(dplyr)

#NO: This code should work but read.delim() causes chloroplast/mitochondria rows to NOT be there.
# cov22.cov <- read.delim("/data/experiments/F15FTSUSAT0821_ARAhlaM/bismarkCytosineReport/22.cov")
# cov23.cov <- read.delim("/data/experiments/F15FTSUSAT0821_ARAhlaM/bismarkCytosineReport/23.cov")
# cov24.cov <- read.delim("/data/experiments/F15FTSUSAT0821_ARAhlaM/bismarkCytosineReport/24.cov")
# cov31.cov <- read.delim("/data/experiments/F15FTSUSAT0821_ARAhlaM/bismarkCytosineReport/31.cov")
# cov32.cov <- read.delim("/data/experiments/F15FTSUSAT0821_ARAhlaM/bismarkCytosineReport/32.cov")
# cov33.cov <- read.delim("/data/experiments/F15FTSUSAT0821_ARAhlaM/bismarkCytosineReport/33.cov")

cov22 <- read.csv("/data/experiments/F15FTSUSAT0821_ARAhlaM/bismarkCytosineReport/22.cov", sep="\t")
cov23 <- read.csv("/data/experiments/F15FTSUSAT0821_ARAhlaM/bismarkCytosineReport/23.cov", sep="\t")
cov24 <- read.csv("/data/experiments/F15FTSUSAT0821_ARAhlaM/bismarkCytosineReport/24.cov", sep="\t")
cov31 <- read.csv("/data/experiments/F15FTSUSAT0821_ARAhlaM/bismarkCytosineReport/31.cov", sep="\t")
cov32 <- read.csv("/data/experiments/F15FTSUSAT0821_ARAhlaM/bismarkCytosineReport/32.cov", sep="\t")
cov33 <- read.csv("/data/experiments/F15FTSUSAT0821_ARAhlaM/bismarkCytosineReport/33.cov", sep="\t")

#colnames(cov22) <- c("chr", "start", "stop", "str", "ctx", "seq", "mC", "uC")
colnames(cov22) <- c("chr", "pos", "str", "mC", "uC", "ctx", "seq")
colnames(cov23) <- c("chr", "pos", "str", "mC", "uC", "ctx", "seq")
colnames(cov24) <- c("chr", "pos", "str", "mC", "uC", "ctx", "seq")
colnames(cov31) <- c("chr", "pos", "str", "mC", "uC", "ctx", "seq")
colnames(cov32) <- c("chr", "pos", "str", "mC", "uC", "ctx", "seq")
colnames(cov33) <- c("chr", "pos", "str", "mC", "uC", "ctx", "seq")

cov22_methylationLevel <- cov22 %>% select(c(chr,mC,uC)) %>% filter(chr=="chloroplast") %>% select(c("mC","uC")) %>% summarise_all(funs(sum)) %>% mutate(pctMethylated=uC/(mC+uC))
cov23_methylationLevel <- cov23 %>% select(c(chr,mC,uC)) %>% filter(chr=="chloroplast") %>% select(c("mC","uC")) %>% summarise_all(funs(sum)) %>% mutate(pctMethylated=uC/(mC+uC))
cov24_methylationLevel <- cov24 %>% select(c(chr,mC,uC)) %>% filter(chr=="chloroplast") %>% select(c("mC","uC")) %>% summarise_all(funs(sum)) %>% mutate(pctMethylated=uC/(mC+uC))
cov31_methylationLevel <- cov31 %>% select(c(chr,mC,uC)) %>% filter(chr=="chloroplast") %>% select(c("mC","uC")) %>% summarise_all(funs(sum)) %>% mutate(pctMethylated=uC/(mC+uC))
cov32_methylationLevel <- cov32 %>% select(c(chr,mC,uC)) %>% filter(chr=="chloroplast") %>% select(c("mC","uC")) %>% summarise_all(funs(sum)) %>% mutate(pctMethylated=uC/(mC+uC))
cov33_methylationLevel <- cov33 %>% select(c(chr,mC,uC)) %>% filter(chr=="chloroplast") %>% select(c("mC","uC")) %>% summarise_all(funs(sum)) %>% mutate(pctMethylated=uC/(mC+uC))

print(paste("cov22_methylationLevel = ",cov22_methylationLevel$pctMethylated))
print(paste("cov23_methylationLevel = ",cov23_methylationLevel$pctMethylated))
print(paste("cov24_methylationLevel = ",cov24_methylationLevel$pctMethylated))
print(paste("cov31_methylationLevel = ",cov31_methylationLevel$pctMethylated))
print(paste("cov32_methylationLevel = ",cov32_methylationLevel$pctMethylated))
print(paste("cov33_methylationLevel = ",cov33_methylationLevel$pctMethylated))

# cov22_methylationLevel =  0.996404953684659
# cov23_methylationLevel =  0.995368106880399
# cov24_methylationLevel =  0.996418570946231
# cov31_methylationLevel =  0.995241354025608
# cov32_methylationLevel =  0.996284067456537
# cov33_methylationLevel =  0.99674703801451
```r
