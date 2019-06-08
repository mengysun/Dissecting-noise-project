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
#Figure 3
Mouse_GO_term<-read.table(file="Data/GO_term_name.txt",sep="\t",header=TRUE,quote="", fill=FALSE)
Mouse_Mito<-Mouse_GO_term%>%
  filter(GO.term.name=="mitochondrion")
Mito_index<-match(Gene_noise_table$Genes,Mouse_Mito$Gene.name)
Gene_noise_table$Mito<-(!is.na(Mito_index))
Gene_noise_Mito<-Gene_noise_table%>%
  filter(Mito)
Gene_noise_nonMito<-Gene_noise_table%>%
  filter(!Mito)

#Fig.3b
Mito_intrinsic_tab<-data.frame(c(Gene_noise_table$intrinsic_residual,Gene_noise_table$intrinsic_residual_controlEx),
                               c(Gene_noise_table$Mito,Gene_noise_table$Mito),
                               c(rep("Dint",length(Gene_noise_table$Genes)),rep("Dint_c",length(Gene_noise_table$Genes))))
names(Mito_intrinsic_tab)<-c("int","is_Mito","i_or_c")
Mito_intrinsic_tab$is_Mito[Mito_intrinsic_tab$is_Mito==TRUE]<-"Mitochondrial genes"
Mito_intrinsic_tab$is_Mito[Mito_intrinsic_tab$is_Mito==FALSE]<-"Non-mitochondrial genes"
Mito_intrinsic_tab$is_Mito<-factor(Mito_intrinsic_tab$is_Mito,
                                   levels = c("Mitochondrial genes","Non-mitochondrial genes"),ordered = TRUE)
Mito_intrinsic_tab$f12<-interaction(Mito_intrinsic_tab$is_Mito,Mito_intrinsic_tab$i_or_c)
ggplot(Mito_intrinsic_tab,aes(y=int,x=f12,fill=is_Mito))+
  scale_fill_manual(values=c("#CC79A7", "#56B4E9"))+
  geom_boxplot(outlier.size=0.1,fatten=0.5)+
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

#Fig.3a
Mito_extrinsic_tab<-data.frame(c(Gene_noise_table$extrinsic_residual,Gene_noise_table$extrinsic_residual_controlIn),
                               c(Gene_noise_table$Mito,Gene_noise_table$Mito),
                               c(rep("Dext",length(Gene_noise_table$Genes)),rep("Dext_c",length(Gene_noise_table$Genes))))
names(Mito_extrinsic_tab)<-c("ext","is_Mito","e_or_c")
Mito_extrinsic_tab$is_Mito[Mito_extrinsic_tab$is_Mito==TRUE]<-"Mitochondrial genes"
Mito_extrinsic_tab$is_Mito[Mito_extrinsic_tab$is_Mito==FALSE]<-"Non-mitochondrial genes"
Mito_extrinsic_tab$is_Mito<-factor(Mito_extrinsic_tab$is_Mito,
                                   levels = c("Mitochondrial genes","Non-mitochondrial genes"),ordered = TRUE)
Mito_extrinsic_tab$f12<-interaction(Mito_extrinsic_tab$is_Mito,Mito_extrinsic_tab$e_or_c)
ggplot(Mito_extrinsic_tab,aes(y=ext,x=f12,fill=is_Mito))+
  scale_fill_manual(values=c("#CC79A7", "#56B4E9"))+
  geom_boxplot(outlier.size=0.1,fatten=0.5)+
  xlab(label="")+
  ylab(label="Extrinsic noise")+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
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

#Figure 3c~d
#========================================================================================================================
Gene_noise_raw<-read.table(file="Data/cl7_noise_raw",sep="\t",header=TRUE)
Mouse_genes_all<-read.table(file="Data/Mouse_protein_coding.txt",sep="\t",header=TRUE)
Mouse_genes_all<-Mouse_genes_all%>%
  filter(Chromosome.scaffold.name!="MT")
unique_genes<-unique(as.character(Mouse_genes_all$Gene.name))
gene_unique_index<-match(unique_genes,Mouse_genes_all$Gene.name)
Mouse_genes_all<-Mouse_genes_all[gene_unique_index,]
Mouse_GO_term<-read.table(file="Data/GO_term_name.txt",sep="\t",header=TRUE,quote="", fill=FALSE)
Mouse_expression_all<-read.table(file="Data/mouse_rpkm.txt",sep="\t",header=TRUE,quote="",fill=FALSE)
Mouse_expression_all$mean_exp<-apply(Mouse_expression_all[,1:39],1,mean)
gene_exp_index<-match(Mouse_genes_all$Gene.stable.ID,rownames(Mouse_expression_all))
Mouse_genes_all$mean_exp<-Mouse_expression_all$mean_exp[gene_exp_index]
Mouse_genes_all<-Mouse_genes_all%>%
  filter(!is.na(mean_exp))
#Genomic features
#TATA box
TATA_box<-read.table(file="Data/mouse_TATA_all.bed",sep="\t")
TATA_box$genes = unlist(lapply(TATA_box$V4, function (x) strsplit(as.character(x), "_", fixed=TRUE)[[1]][1]))
TATA_index<-match(Mouse_genes_all$Gene.name,TATA_box$genes)
Mouse_genes_all$TATA<-(!is.na(TATA_index))


#RegNetwork data
TF_target_dat<-read.table(file="Data/Mouse_regulatory_interaction.csv",sep=",",header=TRUE)
miRNA_target_dat<-TF_target_dat[(grepl("miR",as.character(TF_target_dat$regulator_symbol))),]
target_number_dat<-miRNA_target_dat%>%
  group_by(target_symbol)%>%
  dplyr::summarise(target_number=length(target_symbol))
noise_target_index<-match(Mouse_genes_all$Gene.name,target_number_dat$target_symbol)
Mouse_genes_all$target_number<-target_number_dat$target_number[noise_target_index]
Mouse_genes_all$target_number[is.na(Mouse_genes_all$target_number)]<-0
Mouse_GO_term<-read.table(file="Data/GO_term_name.txt",sep="\t",header=TRUE,quote="", fill=FALSE)
Mouse_Mito<-Mouse_GO_term%>%
  filter(GO.term.name=="mitochondrion")
Mito_index<-match(Mouse_genes_all$Gene.name,Mouse_Mito$Gene.name)
Mouse_genes_all$Mito<-(!is.na(Mito_index))
Mito_genes<-Mouse_genes_all%>%
  filter(Mito)
Mito_control<-Mouse_genes_all%>%
  filter(!Mito)
bins<-quantile(Mito_genes$mean_exp,prob=c(1:51)/51)
Mito_stratified<-split(Mito_control,cut(Mito_control$mean_exp,breaks=bins))
genes_each_bin<-numeric(50)
for(i in 1:50){
  genes_each_bin[i]<-length(Mito_stratified[[i]]$Gene.name)
}
sam_per_bin<-min(genes_each_bin)
stratified_Mito_control<-Mito_control[FALSE,]
names(stratified_Mito_control)<-names(Mito_control)
set.seed(8)
for(i in 1:50){
  sub_tab_index<-sample(c(1:genes_each_bin[i]),sam_per_bin)
  sub_tab<-Mito_stratified[[i]][sub_tab_index,]
  stratified_Mito_control<-rbind(stratified_Mito_control,sub_tab)
}


#Figure 3c
sum(stratified_Mito_control$TATA)
344/2850
sum(Mito_genes$TATA)
125/1603
df <- data.frame(genes=c("Mitochondria_genes", "Control"),
                 ratio=c(125/1603,344/2850))
df$bar_order <- factor(df$genes, as.character(df$genes))
ggplot(data=df, aes(x=bar_order, y=ratio,fill=bar_order)) +
  geom_bar(stat="identity",width=0.3)+
  scale_fill_manual(values=c("#CC79A7", "#F0E442"))+
  xlab(label="")+
  ylab(label="Fraction of genes with TATA-box")+
  theme(axis.text.x = element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.y = element_text(size=12,angle=90,vjust = 0.5,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme_linedraw()+
  scale_y_continuous(limits=c(0,0.15),expand=c(0,0))+
  scale_x_discrete(labels=c("Mitochondrial\n genes","Non-mitochondrial \n genes (stratified)"))+
  scale_x_discrete(breaks = NULL)+
  theme(legend.position = "none")

#Figure 3d
df <- data.frame(genes=c("Mitochondria_genes", "Control"),
                 ratio=c(1058/1603,2085/2850))
df$bar_order <- factor(df$genes, as.character(df$genes))
ggplot(data=df, aes(x=bar_order, y=ratio,fill=bar_order)) +
  geom_bar(stat="identity",width=0.3)+
  scale_fill_manual(values=c("#CC79A7", "#F0E442"))+
  xlab(label="")+
  ylab(label="Fraction of genes targeted by miRNA")+
  theme(axis.text.x = element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.y = element_text(size=12,angle=90,vjust = 0.5,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme_linedraw()+
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  # theme(axis.ticks=element_line(size=1.5))+
  # theme(axis.ticks.length=unit(5,"mm"))+
  scale_y_continuous(limits=c(0,1),expand=c(0,0))+
  scale_x_discrete(labels=c("Mitochondrial\n genes","Non-mitochondrial \n genes (stratified)"))+
  scale_x_discrete(breaks = NULL)+
  theme(legend.position = "none")

