library(dplyr)
library(ggplot2)
library(ppcor)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggExtra)
Gene_noise_raw<-read.table(file="Data/cl7_noise_raw",sep="\t",header=TRUE)
Gene_noise_table<-read.table(file="Data/Gene_noise_table_all_cells",sep="\t",header=TRUE)
Gene_noise_table$intrinsic_noise<-Gene_noise_raw$Intrinsic_noise
Gene_noise_table$extrinsic_noise<-Gene_noise_raw$Extrinsic_noise
Gene_expression_table<-read.table(file="Data/cl7_expression",sep="\t",header=TRUE)
noise_expression_index<-match(Gene_noise_table$Genes,Gene_expression_table$Genes)
Gene_noise_table$expression<-Gene_expression_table$rpkm[noise_expression_index]
Mouse_genes_all<-read.table(file="Data/Mouse_genes_all.txt",sep="\t",header=TRUE)
Mouse_genes_all$TSS<-ifelse(Mouse_genes_all$Strand==1,Mouse_genes_all$Gene.start..bp.,Mouse_genes_all$Gene.end..bp.)
Gene_infor_index<-match(Gene_noise_table$Genes,Mouse_genes_all$Gene.name)
Gene_noise_table$chr<-Mouse_genes_all$Chromosome.scaffold.name[Gene_infor_index]
Gene_noise_table$TSS<-Mouse_genes_all$TSS[Gene_infor_index]
Gene_noise_table$gene_end<-Mouse_genes_all$Gene.end..bp.[Gene_infor_index]
Gene_noise_table$gene_start<-Mouse_genes_all$Gene.start..bp.[Gene_infor_index]
Mouse_expression_all<-read.table(file="Data/mouse_rpkm.txt",sep="\t",header=TRUE,quote="",fill=FALSE)
exp_name_index<-match(rownames(Mouse_expression_all),Mouse_genes_all$Gene.stable.ID)
Mouse_expression_all$geneName<-Mouse_genes_all$Gene.name[exp_name_index]
Mouse_expression_all$mean_exp<-rowMeans(Mouse_expression_all[,1:39])
noise_allexp_index<-match(Gene_noise_table$Genes,Mouse_expression_all$geneName)
Gene_noise_table$exp_all<-Mouse_expression_all$mean_exp[noise_allexp_index]
TATA_box<-read.table(file="Data/mouse_TATA_all.bed",sep="\t")

TATA_box$genes = unlist(lapply(TATA_box$V4, function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][1]))
TATA_index<-match(Gene_noise_table$Genes,TATA_box$genes)
TATA_box$motif<-1
Gene_noise_table$TATA_motif<-TATA_box$motif[TATA_index]
Gene_noise_table$TATA_motif[is.na(Gene_noise_table$TATA_motif)]<-0
#Figure 2a
TATA_intrinsic_tab<-data.frame(c(Gene_noise_table$intrinsic_residual,Gene_noise_table$intrinsic_residual_controlEx),
                               c(Gene_noise_table$TATA_motif,Gene_noise_table$TATA_motif),
                               c(rep("Dint",length(Gene_noise_table$Genes)),rep("Dint_c",length(Gene_noise_table$Genes))))
names(TATA_intrinsic_tab)<-c("int","is_TATA_motif","i_or_c")
TATA_intrinsic_tab$is_TATA_motif[TATA_intrinsic_tab$is_TATA_motif==1]<-"With TATA box"
TATA_intrinsic_tab$is_TATA_motif[TATA_intrinsic_tab$is_TATA_motif==0]<-"Without TATA box"
TATA_intrinsic_tab$f12<-interaction(TATA_intrinsic_tab$is_TATA_motif,TATA_intrinsic_tab$i_or_c)
ggplot(TATA_intrinsic_tab,aes(y=int,x=f12,fill=is_TATA_motif))+
  geom_boxplot(outlier.size=0.1,fatten=0.5)+
  scale_fill_manual(values=c("#CC79A7", "#56B4E9"))+
  xlab(label="")+
  ylab(label="Intrinsic noise")+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme_linedraw()+
  theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12,family="Times New Roman",color="black"))+
  theme(legend.direction = "vertical")+
  theme(legend.margin=margin(-0.5))+
  removeGridX()+
  theme(axis.text.x=element_blank())+
  theme(axis.ticks.x=element_blank())+
  theme(axis.title.x=element_blank())+
  scale_y_continuous(limits=c(-3000,6000),breaks=c(-2000,0,2000,4000))+
  theme(legend.position = "none")

#Fig.2b
TATA_extrinsic_tab<-data.frame(c(Gene_noise_table$extrinsic_residual,Gene_noise_table$extrinsic_residual_controlIn),
                               c(Gene_noise_table$TATA_motif,Gene_noise_table$TATA_motif),
                               c(rep("Dext",length(Gene_noise_table$Genes)),rep("Dext_c",length(Gene_noise_table$Genes))))
names(TATA_extrinsic_tab)<-c("ext","is_TATA_motif","e_or_c")
TATA_extrinsic_tab$is_TATA_motif[TATA_extrinsic_tab$is_TATA_motif==1]<-"With TATA box"
TATA_extrinsic_tab$is_TATA_motif[TATA_extrinsic_tab$is_TATA_motif==0]<-"Without TATA box"
TATA_extrinsic_tab$f12<-interaction(TATA_extrinsic_tab$is_TATA_motif,TATA_extrinsic_tab$e_or_c)
ggplot(TATA_extrinsic_tab,aes(y=ext,x=f12,fill=is_TATA_motif))+
  geom_boxplot(outlier.size=0.1,fatten=0.5)+
  scale_fill_manual(values=c("#CC79A7", "#56B4E9"))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black"))+
  xlab(label="")+
  ylab(label="Extrinsic noise")+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme_linedraw()+
  theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12,family="Times New Roman",color="black"))+
  theme(legend.direction = "vertical")+
  theme(legend.margin=margin(-0.5))+
  removeGridX()+
  theme(axis.text.x=element_blank())+
  theme(axis.ticks.x=element_blank())+
  theme(axis.title.x=element_blank())+
  scale_y_continuous(limits=c(-3000,6000),breaks=c(-2000,0,2000,4000))+
  theme(legend.position = "none")

#Fig. 2c-f
TF_target_datmiRNA<-read.table(file="Data/Mouse_regulatory_interaction.csv",sep=",",header=TRUE)
miRNA_target_dat<-TF_target_datmiRNA[(grepl("miR",as.character(TF_target_datmiRNA$regulator_symbol))),]
target_number_dat<-miRNA_target_dat%>%
  group_by(target_symbol)%>%
  dplyr::summarise(target_number=length(target_symbol))
noise_target_index<-match(Gene_noise_table$Genes,target_number_dat$target_symbol)
Gene_noise_table$target_number<-target_number_dat$target_number[noise_target_index]
Gene_noise_table$target_number[is.na(Gene_noise_table$target_number)]<-0
#Fig.2c
miRNA_intrinsic_tab<-data.frame(c(Gene_noise_table$intrinsic_residual,Gene_noise_table$intrinsic_residual_controlEx),
                                c(Gene_noise_table$target_number,Gene_noise_table$target_number),
                                c(rep("Dint",length(Gene_noise_table$Genes)),rep("Dint_c",length(Gene_noise_table$Genes))))
names(miRNA_intrinsic_tab)<-c("int","target_by_miRNA","i_or_c")
miRNA_intrinsic_tab$target_by_miRNA[miRNA_intrinsic_tab$target_by_miRNA>0]<-"Targeted by miRNA"
miRNA_intrinsic_tab$target_by_miRNA[miRNA_intrinsic_tab$target_by_miRNA==0]<-"Not targeted by miRNA"
miRNA_intrinsic_tab$target_by_miRNA<-factor(miRNA_intrinsic_tab$target_by_miRNA,
                                            levels = c('Targeted by miRNA','Not targeted by miRNA'),ordered = TRUE)
miRNA_intrinsic_tab$f12<-interaction(miRNA_intrinsic_tab$target_by_miRNA,miRNA_intrinsic_tab$i_or_c)
ggplot(miRNA_intrinsic_tab,aes(y=int,x=f12,fill=target_by_miRNA))+
  geom_boxplot(outlier.size=0.1,fatten=0.5)+
  scale_fill_manual(values=c("#009E73","#F0E442"))+
  xlab(label="")+
  ylab(label="Intrinsic noise")+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme_linedraw()+
  theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12,family="Times New Roman",color="black"))+
  theme(legend.direction = "vertical")+
  theme(legend.margin=margin(-0.5))+
  removeGridX()+
  theme(axis.text.x=element_blank())+
  theme(axis.ticks.x=element_blank())+
  theme(axis.title.x=element_blank())+
  scale_y_continuous(limits=c(-3000,6000),breaks=c(-2000,0,2000,4000))

#Fig.2d
ggplot(Gene_noise_table,aes(x=target_number,y=intrinsic_residual))+
  geom_point(size=0.5,alpha=0.1)+
  geom_smooth(method="lm",se=FALSE,size=0.5,color="#0072B2")+
  xlab(label="Number of miRNA species")+
  ylab(expression(italic(D)[int]))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.y=element_text(angle=90,size=12,family="Times New Roman",color="black",hjust=0.5, vjust=0.5))+
  theme_linedraw()

#Fig.2e
ggplot(Gene_noise_table,aes(x=target_number,y=intrinsic_residual_controlEx))+
  geom_point(size=0.5,alpha=0.1)+
  geom_smooth(method="lm",se=FALSE,size=0.5,color="#0072B2")+
  xlab(label="Number of miRNA species")+
  ylab(expression(italic(D)*{"'"}[int]))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black",angle=90,hjust = 0.5,vjust=0.5))+
  theme_linedraw()

#Fig.2f
miRNA_extrinsic_tab<-data.frame(c(Gene_noise_table$extrinsic_residual,Gene_noise_table$extrinsic_residual_controlIn),
                                c(Gene_noise_table$target_number,Gene_noise_table$target_number),
                                c(rep("Dext",length(Gene_noise_table$Genes)),rep("Dext_c",length(Gene_noise_table$Genes))))
names(miRNA_extrinsic_tab)<-c("ext","target_by_miRNA","e_or_c")
miRNA_extrinsic_tab$target_by_miRNA[miRNA_extrinsic_tab$target_by_miRNA>0]<-"Targeted by miRNA"
miRNA_extrinsic_tab$target_by_miRNA[miRNA_extrinsic_tab$target_by_miRNA==0]<-"Not targeted by miRNA"
miRNA_extrinsic_tab$target_by_miRNA<-factor(miRNA_extrinsic_tab$target_by_miRNA,
                                            levels = c('Targeted by miRNA','Not targeted by miRNA'),ordered = TRUE)
miRNA_extrinsic_tab$f12<-interaction(miRNA_extrinsic_tab$target_by_miRNA,miRNA_extrinsic_tab$e_or_c)
ggplot(miRNA_extrinsic_tab,aes(y=ext,x=f12,fill=target_by_miRNA))+
  geom_boxplot(outlier.size=0.1,fatten=0.5)+
  scale_fill_manual(values=c("#009E73","#F0E442"))+
  xlab(label="")+
  ylab(label="Extrinsic noise")+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme_linedraw()+
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  #       panel.background = element_blank(),axis.line = element_line(colour = "black"),panel.border=element_blank())+
  theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=12,family="Times New Roman",color="black"))+
  theme(legend.direction = "vertical")+
  theme(legend.margin=margin(-0.5))+
  removeGridX()+
  theme(axis.text.x=element_blank())+
  theme(axis.ticks.x=element_blank())+
  theme(axis.title.x=element_blank())+
  scale_y_continuous(limits=c(-3000,6000),breaks=c(-2000,0,2000,4000))

#Fig.2g~j
TF_target_dat<-read.table(file="Data/Mouse_regulatory_interaction.csv",sep=",",header=TRUE)
TF_high_quality<-TF_target_dat%>%
  filter(confidence=="High"|confidence=="Medium"|confidence=="Low")
target_match_index<-match(TF_high_quality$target_symbol,Gene_noise_table$Genes)
TF_high_quality$target_extrinsic_noise<-Gene_noise_table$extrinsic_residual[target_match_index]
#for D'ext, switch extrinsic_residual to extrinsic_residual_controlIn
TF_high_quality$target_intrinsic_noise<-Gene_noise_table$intrinsic_residual[target_match_index]
#for D'int, switch extrinsic_residual to intrinsic_residual_controlEx
TF_match_index<-match(TF_high_quality$regulator_symbol,Gene_noise_table$Genes)
TF_high_quality$TF_intrinsic_noise<-Gene_noise_table$intrinsic_noise[TF_match_index]
TF_high_quality$TF_extrinsic_noise<-Gene_noise_table$extrinsic_noise[TF_match_index]
TF_high_quality_with_noise<-TF_high_quality%>%
  filter((!is.na(TF_extrinsic_noise))&(!is.na(target_extrinsic_noise)))
TF_high_quality_with_noise<-TF_high_quality_with_noise[!(TF_high_quality_with_noise$target_symbol%in%TF_high_quality_with_noise$regulator_symbol),]
TF_high_quality_with_noise<-TF_high_quality_with_noise[as.character(TF_high_quality_with_noise$regulator_symbol)!=as.character(TF_high_quality_with_noise$target_symbol),]
TF_high_quality_with_noise$TF_total<-TF_high_quality_with_noise$TF_intrinsic_noise+TF_high_quality_with_noise$TF_extrinsic_noise
Target_list<-names(table(as.character(TF_high_quality_with_noise$target_symbol)))
Target_keep_index<-match(Target_list,TF_high_quality_with_noise$target_symbol)
TF_high_quality_with_noise<-TF_high_quality_with_noise[Target_keep_index,]
TF_noise_and_target_noise<-TF_high_quality_with_noise%>%
  group_by(regulator_symbol)%>%
  dplyr::summarise(TF_i_noise=mean(TF_intrinsic_noise),TF_e_noise=mean(TF_extrinsic_noise),
                   Target_i_noise=mean(target_intrinsic_noise),Target_e_noise=mean(target_extrinsic_noise))
TF_exp_index<-match(TF_noise_and_target_noise$regulator_symbol,Gene_expression_table$Genes)
TF_noise_and_target_noise$TF_exp<-Gene_expression_table$rpkm[TF_exp_index]

#Fig.2g~h
ggplot(TF_noise_and_target_noise,aes(x=TF_i_noise/2+TF_e_noise,y=Target_i_noise))+
  geom_point(size=0.5)+
  geom_smooth(method="lm",se=FALSE,size=0.5,color="#0072B2")+
  xlab(expression(atop(italic(Trans)*"-regulator noise",{"("}*eta[int]^{2}*{"+"}*eta[ext]^{2}*{")"})))+
  ylab(expression(atop("Mean target noise",{"("}*italic(D)[int]*{")"})))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black",angle=90,hjust = 0.5,vjust=0.5))+
  theme_linedraw()+
  scale_y_continuous(limits=c(-3000,6000),breaks=c(-2000,0,2000,4000))

#Fig.2i~j
ggplot(TF_noise_and_target_noise,aes(x=TF_i_noise/2+TF_e_noise,y=Target_e_noise))+
  geom_point(size=0.5)+
  geom_smooth(method="lm",se=FALSE,size=0.5,color="#0072B2")+
  xlab(expression(atop(italic(Trans)*"-regulator noise",{"("}*eta[int]^{2}*{"+"}*eta[ext]^{2}*{")"})))+
  ylab(expression(atop("Mean target noise",{"("}*italic(D)[ext]*{")"})))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black",angle=90,hjust = 0.5,vjust=0.5))+
  theme_linedraw()+
  scale_y_continuous(limits=c(-3000,6000),breaks=c(-2000,0,2000,4000))

#Fig.2k~l
TF_high_quality_for_shuffle<-TF_high_quality%>%
  filter(!is.na(target_extrinsic_noise))
TF_high_quality_for_shuffle<-TF_high_quality_for_shuffle[!TF_high_quality_for_shuffle$target_symbol%in%TF_high_quality_for_shuffle$regulator_symbol,]
regulator_list_unique<-names(table(as.character(TF_high_quality_for_shuffle$regulator_symbol)))
regulator_list_unique<-regulator_list_unique[table(as.character(TF_high_quality_for_shuffle$regulator_symbol))>1]
TF_high_quality_for_shuffle<-TF_high_quality_for_shuffle[TF_high_quality_for_shuffle$regulator_symbol%in%regulator_list_unique,]
target_list_unique<-unique(as.character(TF_high_quality_for_shuffle$target_symbol))
target_unique_index<-match(target_list_unique,TF_high_quality_for_shuffle$target_symbol)
target_intrinsic_noise<-as.numeric(as.character(TF_high_quality_for_shuffle$target_intrinsic_noise[target_unique_index]))
target_extrinsic_noise<-as.numeric(as.character(TF_high_quality_for_shuffle$target_extrinsic_noise[target_unique_index]))
target_index_back<-match(TF_high_quality_for_shuffle$target_symbol,target_list_unique)
sd_extrinsic_shuffled<-numeric(10000)
sd_intrinsic_shuffled<-numeric(10000)
set.seed(7)
for(i in 1:10000){
  shuffled_target_index<-sample(c(1:length(target_list_unique)))
  shuffled_index<-shuffled_target_index[target_index_back]
  TF_high_quality_for_shuffle$target_shuffle_intrinsic<-target_intrinsic_noise[shuffled_index]
  TF_high_quality_for_shuffle$target_shuffle_extrinsic<-target_extrinsic_noise[shuffled_index]
  shuffled_results<-TF_high_quality_for_shuffle%>%
    group_by(regulator_symbol)%>%
    dplyr::summarise(sd_extrinsic=sd(target_shuffle_extrinsic),sd_intrinsic=sd(target_shuffle_intrinsic))
  sd_extrinsic_shuffled[i]<-median(shuffled_results$sd_extrinsic)
  sd_intrinsic_shuffled[i]<-median(shuffled_results$sd_intrinsic)
}
Target_noise_observed<-TF_high_quality_for_shuffle%>%
  group_by(regulator_symbol)%>%
  dplyr::summarise(sd_extrinsic=sd(target_extrinsic_noise),sd_intrinsic=sd(target_intrinsic_noise))
dat_shuffled<-data.frame(sd_extrinsic_shuffled,sd_intrinsic_shuffled)
write.table(dat_shuffled,file="cl7_shuffled_tab",row.names=FALSE,col.names=TRUE,sep="\t")
#Fig.2k
dat_shuffled<-read.table(file="cl7_shuffled_tab",sep="\t",header=TRUE)
tiff(file="Figure2i.tiff",width=2.57,height=1.94,units="in",res=300)
ggplot(dat_shuffled, aes(x=sd_intrinsic_shuffled)) +
  theme_bw() +
  geom_histogram(colour="grey", fill="grey70",aes(y=..count../sum(..count..)),bins=100)+
  geom_segment(aes(x = 586.8035, y = 0.005, xend =586.8035, yend = 0.001), size=0.5,arrow = arrow(length = unit(0.1, "cm")),color="red")+
  ylab(label="Frequency")+
  xlab(expression(atop("Median within-group", {"standard deviation"}~"in"~italic(D)[int])))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.y=element_text(size=12,angle=90,vjust=0,family="Times New Roman",color="black",margin=margin(t=0,r=8,b=0,l=0)))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black"))+
  scale_y_continuous(expand=c(0,0))
dev.off()


#Fig.2l
ggplot(dat_shuffled, aes(x=sd_extrinsic_shuffled)) +
  theme_bw() +
  geom_histogram(colour="grey", fill="grey70",aes(y=..count../sum(..count..)),bins=100)+
  geom_segment(aes(x = 1073.151, y = 0.005, xend =1073.151, yend = 0.001), size=0.5,arrow = arrow(length = unit(0.1, "cm")),color="red")+
  ylab(label="Frequency")+
  xlab(expression(atop("Median within-group", {"standard deviation"}~"in"~italic(D)[ext])))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.y=element_text(size=12,angle=90,vjust=0,family="Times New Roman",color="black",margin=margin(t=0,r=8,b=0,l=0)))+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.x=element_text(size=12,family="Times New Roman",color="black"))+
  scale_x_discrete(limits=c(1070,1100,1130))+expand_limits(x=c(1060,1160))+
  scale_y_continuous(expand=c(0,0))
