#Figure 1
library(dplyr)
library(ggplot2)
library(MASS)
library(grid)
library(gridExtra)
library(ppcor)
library(ggpubr)
library(noise)
library(sm)
library(class)
library(gam)
library(FNN)
library('GEOquery')
#load data
C57_genes<-data.frame(readRDS("Data/mouse.c57.counts.rds"))
Cas_genes<-data.frame(readRDS("Data/mouse.cast.counts.rds"))
Mouse_expression<-data.frame(readRDS("Data/mouse.rpkm.3digits.rds"))
Mouse_total_counts<-data.frame(readRDS("Data/mouse.counts.rds"))
Mouse_genes_all<-read.table(file="Data/Mouse_genes_all.txt",sep="\t",header=TRUE)
C57_genes<-C57_genes[,61:120]
Cas_genes<-Cas_genes[,61:120]
C57_high<-rowMeans(C57_genes)
Cas_high<-rowMeans(Cas_genes)
Mouse_expression<-Mouse_expression[,61:120]
Mouse_total_counts<-Mouse_total_counts[,61:120]
spike_in_total<-as.numeric(as.character(colSums(Mouse_total_counts[35209:35280,])))
conversion_factor<-spike_in_total/max(spike_in_total)
for(i in 1:60){
  C57_genes[,i]<-C57_genes[,i]/conversion_factor[i]
  Cas_genes[,i]<-Cas_genes[,i]/conversion_factor[i]
}
focus_gene_index<-((C57_high>5)&(Cas_high>5))
C57_focus_genes<-C57_genes[focus_gene_index,]
Cas_focus_genes<-Cas_genes[focus_gene_index,]
Focus_expression<-rowMeans(Mouse_expression[focus_gene_index,])
Focus_genes<-names(Focus_expression)
Extrinsic_noise<-numeric(5566)
Intrinsic_noise<-numeric(5566)
read_counts<-numeric(5566)
ks_p_value<-numeric(5566)
for(i in 1:5566){
  scaling_factor<-mean(as.numeric(as.character((C57_focus_genes[i,]))))*mean(as.numeric(as.character((Cas_focus_genes[i,]))))
  
  Extrinsic_noise[i]<-(computeExtrinsicNoise(as.numeric(as.character(C57_focus_genes[i,])),
                                             as.numeric(as.character(Cas_focus_genes[i,])))$unbiased)/scaling_factor
  Intrinsic_noise[i]<-(computeIntrinsicNoise(as.numeric(as.character(C57_focus_genes[i,])),
                                             as.numeric(as.character(Cas_focus_genes[i,])))$unbiasedGeneral)/scaling_factor
  read_counts[i]<-(mean(as.numeric(as.character(C57_focus_genes[i,])))+mean(as.numeric(as.character(Cas_focus_genes[i,]))))/2
  ks_p_value[i]<-ks.test(as.numeric(as.character(C57_focus_genes[i,])),as.numeric(as.character(Cas_focus_genes[i,])))$p.value
}
filter_by_dataAvailable<-((!is.na(Intrinsic_noise))&(!is.na(Extrinsic_noise)))
Focus_expression<-Focus_expression[filter_by_dataAvailable]
Focus_genes<-Focus_genes[filter_by_dataAvailable]
Intrinsic_noise<-Intrinsic_noise[filter_by_dataAvailable]
Extrinsic_noise<-Extrinsic_noise[filter_by_dataAvailable]
read_counts<-read_counts[filter_by_dataAvailable]
ks_p_value<-ks_p_value[filter_by_dataAvailable]
Gene_noise_table<-data.frame(Focus_genes,Focus_expression,Intrinsic_noise,Extrinsic_noise,read_counts,ks_p_value)
Chr_index<-match(Gene_noise_table$Focus_genes,Mouse_genes_all$Gene.name)
Gene_noise_table$chr<-Mouse_genes_all$Chromosome.scaffold.name[Chr_index]
Gene_noise_table<-Gene_noise_table%>%
  filter(chr!="X"&chr!="3"&chr!="4")
Gene_noise_table$noise<-Gene_noise_table$Intrinsic_noise+Gene_noise_table$Extrinsic_noise
Gene_noise_table$log_intrinsic_noise<-log(Gene_noise_table$Intrinsic_noise)
Gene_noise_table$log_readcounts<-log(Gene_noise_table$read_counts)
Gene_noise_table$log_exp<-log(Gene_noise_table$Focus_expression)
Gene_noise_table$adjust_ksp<-p.adjust(Gene_noise_table$ks_p_value,method="fdr")
Gene_noise_table<-Gene_noise_table%>%
  filter(adjust_ksp>0.05)
Gene_noise_table_for_cl7<-read.table(file="Data/Gene_noise_table_all_cells",sep="\t",header=TRUE)
Gene_noise_table_for_cl7_1<-read.table(file="Data/Gene_noise_table_first_30_cells",sep="\t",header=TRUE)
Gene_noise_table_for_cl7_2<-read.table(file="Data/Gene_noise_table_last_30_cells",sep="\t",header=TRUE)
Gene_match_index<-match(Gene_noise_table_for_cl7_1$Genes,Gene_noise_table_for_cl7_2$Genes)
Gene_noise_table_for_cl7_1$cl7_2_intrinsic_dm<-Gene_noise_table_for_cl7_2$intrinsic_residual[Gene_match_index]
Gene_noise_table_for_cl7_1$cl7_2_extrinsic_dm<-Gene_noise_table_for_cl7_2$extrinsic_residual[Gene_match_index]
Gene_noise_table_for_common<-Gene_noise_table_for_cl7_1%>%
  filter(!is.na(cl7_2_intrinsic_dm))
Gene_noise_raw1<-read.table(file="Data/cl7_raw_noise_first_30",sep="\t",header=TRUE)
Gene_noise_raw2<-read.table(file="Data/cl7_raw_noise_last_30",sep="\t",header=TRUE)
Gene_match_indexRaw<-match(Gene_noise_raw1$Focus_genes,Gene_noise_raw2$Focus_genes)
Gene_noise_raw1$cl7_2_intrinsic<-Gene_noise_raw2$Intrinsic_noise[Gene_match_indexRaw]
Gene_noise_raw1$cl7_2_extrinsic<-Gene_noise_raw2$Extrinsic_noise[Gene_match_indexRaw]
Gene_noise_table_raw_for_common<-Gene_noise_raw1%>%
  filter(!is.na(cl7_2_extrinsic))
#Fig.1b
ggplot(Gene_noise_table_raw_for_common,aes(x=log(Intrinsic_noise/2),y=log(cl7_2_intrinsic/2)))+
  geom_point(alpha=0.1,size=0.5)+
  geom_abline(slope=1,color="#D55E00",size=0.5,intercept=0)+
  xlab(expression('ln '*eta[int]^{2}~"(subsample 1)"))+
  ylab(expression('ln '*eta[int]^{2}~"(subsample 2)"))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black",margin=margin(t=2,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black",angle=90,hjust = 0.5,margin=margin(t=0,r=2,b=0,l=0)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(limits=c(-6.5,3),breaks=c(-6,-4,-2,0,2))+
  scale_y_continuous(limits=c(-6.5,3),breaks=c(-6,-4,-2,0,2))

#Fig.1c
Gene_noise_table_raw_for_common$log_extrinsic<-log(Gene_noise_table_raw_for_common$Extrinsic_noise-min(Gene_noise_table_raw_for_common$Extrinsic_noise)+0.1)
Gene_noise_table_raw_for_common$log_extrinsic_2<-log(Gene_noise_table_raw_for_common$cl7_2_extrinsic-min(Gene_noise_table_raw_for_common$cl7_2_extrinsic)+0.1)
ggplot(Gene_noise_table_raw_for_common,aes(x=log_extrinsic,y=log_extrinsic_2))+
  geom_point(alpha=0.1,size=0.5)+
  geom_abline(slope=1,color="#D55E00",size=0.5,intercept=0)+
  xlab(expression('ln '*eta[ext]^{2}~"(subsample 1)"))+
  ylab(expression('ln '*eta[ext]^{2}~"(subsample 2)"))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black",margin=margin(t=2,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black",angle=90,hjust = 0.5,margin=margin(t=0,r=2,b=0,l=0)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(limits=c(-3.3,3.3),breaks=c(-3,-2,-1,0,1,2,3))+
  scale_y_continuous(limits=c(-3.3,3.3),breaks=c(-3,-2,-1,0,1,2,3))

#Fig.1d
ggplot(Gene_noise_table,aes(x=log_exp,y=log_intrinsic_noise-log(2)))+
  geom_point(size=0.5,alpha=0.1)+
  xlab(label="ln(mean expression in RPKM)")+
  ylab(expression('ln '*eta[int]^{2}))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black",hjust=1))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black",angle=90,hjust = 0.6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Fig.1e
Gene_noise_table$extrinsic_non_negative<-log(Gene_noise_table$Extrinsic_noise-min(Gene_noise_table$Extrinsic_noise)+0.1)
ggplot(Gene_noise_table,aes(x=log_exp,y=extrinsic_non_negative))+
  geom_point(size=0.5,alpha=0.1)+
  xlab(label="ln(mean expression in RPKM)")+
  ylab(expression('ln '*eta[ext]^{2}))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black",hjust=1))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black",angle=90,hjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Fig.1f
ggplot(Gene_noise_table_for_common,aes(x=intrinsic_residual,y=cl7_2_intrinsic_dm))+
  geom_point(size=0.5,alpha=0.1)+
  geom_abline(slope=1,color="#D55E00",size=0.5,intercept=0)+
  xlab(expression(italic(D)[int]~"(subsample 1)"))+
  ylab(expression(italic(D)[int]~"(subsample 2)"))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black",margin=margin(t=2,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black",angle=90,hjust = 0.5,margin=margin(t=0,r=2,b=0,l=0)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(limits=c(-4000,4000),breaks=c(-3000,0,3000))+
  scale_y_continuous(limits=c(-4000,4000),breaks=c(-3000,0,3000))+
  coord_fixed()

#Fig.1g
ggplot(Gene_noise_table_for_common,aes(x=extrinsic_residual,y=cl7_2_extrinsic_dm))+
  geom_point(alpha=0.1,size=0.5)+
  geom_abline(slope=1,color="#D55E00",size=0.5,intercept=0)+
  xlab(expression(italic(D)[ext]~"(subsample 1)"))+
  ylab(expression(italic(D)[ext]~"(subsample 2)"))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black",margin=margin(t=2,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black",angle=90,hjust = 0.5,margin=margin(t=0,r=2,b=0,l=0)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_continuous(limits=c(-4000,4000),breaks=c(-3000,0,3000))+
  scale_y_continuous(limits=c(-4000,4000),breaks=c(-3000,0,3000))+
  coord_fixed()

#Fig.1h
ggplot(Gene_noise_table_for_cl7,aes(x=extrinsic_residual,y=intrinsic_residual))+
  geom_point(alpha=0.1,size=0.5)+
  xlab(expression(italic(D)[ext]))+
  ylab(expression(italic(D)[int]))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black",margin=margin(t=2,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black",angle=90,hjust = 0.5,margin=margin(t=0,r=2,b=0,l=0)))+
  scale_x_continuous(breaks=c(-2000,0,2000))+
  scale_y_continuous(breaks=c(-2000,0,2000))+
  theme(aspect.ratio=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
