library(ggplot2)
#setwd("/Users/mackenzie/Documents")
SS <- read.csv("/Users/mackenzie/Rproject/tomato_alignment_statisctics.csv")
head(SS)

class(SS$C_methylated_in_CpG)
#   CG 
CG_all<- ggplot(SS,  aes(x  =  TYPE,  y  =  CG,  fill  =  TYPE))  + 
  scale_fill_manual(values=c("#b3cde0", "#6497b1", "#005b96","#d35b5a","#fd7675","#9cb807","#41d608","#51ccf7")) +
  stat_boxplot(geom ='errorbar', width=0.5) + 
  geom_boxplot(outlier.size = 0.1, width=0.5)+ 
  #facet_wrap(~gen, nrow = 1)+ 
  #geom_jitter(width = 0.05)+
  scale_x_discrete(limits = c("WT","Dr_minus_mild","Dr_plus_dwarf","Dr_plus_mild","grafted_Gen1_RR","grafted_Gen1_RDR","F2_epi_nor","F2_epi_high"))+
  xlab("epi genometype")  +      
  ylab("Cytopsine methylated in CG (%) ")  +  
  ggtitle("CG")+  
  theme(plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, colour = "black", face = "bold"),
        axis.text.x =  element_text(angle  =  45,  hjust  =  1, face = "bold", size = 10),
        axis.text.y =  element_text(face = "bold", size = 10))  
ggsave(CG_all, width=12, height=5, file="/Users/mackenzie/Rproject/tomato_CG_all_5_20_2019.tiff", dpi = 1200)


#   CG 
CHG_all<- ggplot(SS,  aes(x  =  TYPE,  y  =  CHG,  fill  =  TYPE))  + 
  scale_fill_manual(values=c("#b3cde0", "#6497b1", "#005b96","#d35b5a","#fd7675","#9cb807","#41d608","#51ccf7")) +
  stat_boxplot(geom ='errorbar', width=0.5) + 
  geom_boxplot(outlier.size = 0.1, width=0.5)+ 
  #facet_wrap(~gen, nrow = 1)+ 
  #geom_jitter(width = 0.05)+
  scale_x_discrete(limits = c("WT","Dr_minus_mild","Dr_plus_dwarf","Dr_plus_mild","grafted_Gen1_RR","grafted_Gen1_RDR","F2_epi_nor","F2_epi_high"))+
  xlab("epi genometype")  +      
  ylab("Cytopsine methylated in CHG (%) ")  +  
  ggtitle("CHG")+  
  theme(plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, colour = "black", face = "bold"),
        axis.text.x =  element_text(angle  =  45,  hjust  =  1, face = "bold", size = 10),
        axis.text.y =  element_text(face = "bold", size = 10))  
ggsave(CHG_all, width=12, height=5, file="/Users/mackenzie/Rproject/tomato_CHG_all_5_20_2019.tiff", dpi = 1200)


#   CHH 
CHH_all<- ggplot(SS,  aes(x  =  TYPE,  y  =  CHH,  fill  =  TYPE))  + 
  scale_fill_manual(values=c("#b3cde0", "#6497b1", "#005b96","#d35b5a","#fd7675","#9cb807","#41d608","#51ccf7")) +
  stat_boxplot(geom ='errorbar', width=0.5) + 
  geom_boxplot(outlier.size = 0.1, width=0.5)+ 
  #facet_wrap(~gen, nrow = 1)+ 
  #geom_jitter(width = 0.05)+
  scale_x_discrete(limits = c("WT","Dr_minus_mild","Dr_plus_dwarf","Dr_plus_mild","grafted_Gen1_RR","grafted_Gen1_RDR","F2_epi_nor","F2_epi_high"))+
  xlab("epi genometype")  +      
  ylab("Cytopsine methylated in CHH (%) ")  +  
  ggtitle("CHH")+  
  theme(plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, colour = "black", face = "bold"),
        axis.text.x =  element_text(angle  =  45,  hjust  =  1, face = "bold", size = 10),
        axis.text.y =  element_text(face = "bold", size = 10))  
ggsave(CHH_all, width=12, height=5, file="/Users/mackenzie/Rproject/tomato_CHH_all_5_20_2019.tiff", dpi = 1200)

