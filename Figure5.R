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
#Figure 5
#figure 5 a~b
Mouse_GO_term<-read.table(file="Data/GO_term_name.txt",sep="\t",header=TRUE,quote="", fill=FALSE)
Mouse_cell_cycle<-Mouse_GO_term%>%
  filter(GO.term.name=="cell cycle")
cell_cycle_index<-match(Gene_noise_table$Genes,Mouse_cell_cycle$Gene.name)
Gene_noise_table$cell_cycle<-(!is.na(cell_cycle_index))
Gene_noise_cell_cycle<-Gene_noise_table%>%
  filter(cell_cycle==TRUE)
Gene_noise_non_cell_cycle<-Gene_noise_table%>%
  filter(cell_cycle!=TRUE)
#Figure 5a
CellCycle_intrinsic_tab<-data.frame(c(Gene_noise_table$intrinsic_residual,Gene_noise_table$intrinsic_residual_controlEx),
                                    c(Gene_noise_table$cell_cycle,Gene_noise_table$cell_cycle),
                                    c(rep("Dint",length(Gene_noise_table$Genes)),rep("Dint_c",length(Gene_noise_table$Genes))))
names(CellCycle_intrinsic_tab)<-c("int","is_cell_cycle","i_or_c")
CellCycle_intrinsic_tab$is_cell_cycle[CellCycle_intrinsic_tab$is_cell_cycle==TRUE]<-"Cell-cycle genes"
CellCycle_intrinsic_tab$is_cell_cycle[CellCycle_intrinsic_tab$is_cell_cycle==FALSE]<-"Non-cell-cycle genes"
CellCycle_intrinsic_tab$is_cell_cycle<-factor(CellCycle_intrinsic_tab$is_cell_cycle,
                                              levels = c("Cell-cycle genes","Non-cell-cycle genes"),ordered = TRUE)
CellCycle_intrinsic_tab$f12<-interaction(CellCycle_intrinsic_tab$is_cell_cycle,CellCycle_intrinsic_tab$i_or_c)
ggplot(CellCycle_intrinsic_tab,aes(y=int,x=f12,fill=is_cell_cycle))+
  scale_fill_manual(values=c("#CC79A7", "#56B4E9"))+
  geom_boxplot(outlier.size=0.1,fatten=0.5)+
  xlab(label="")+
  ylab(label="Intrinsic noise")+
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

#figure 5b
CellCycle_extrinsic_tab<-data.frame(c(Gene_noise_table$extrinsic_residual,Gene_noise_table$extrinsic_residual_controlIn),
                                    c(Gene_noise_table$cell_cycle,Gene_noise_table$cell_cycle),
                                    c(rep("Dext",length(Gene_noise_table$Genes)),rep("ext_c",length(Gene_noise_table$Genes))))
names(CellCycle_extrinsic_tab)<-c("ext","is_cell_cycle","e_or_c")
CellCycle_extrinsic_tab$is_cell_cycle[CellCycle_extrinsic_tab$is_cell_cycle==TRUE]<-"Cell-cycle genes"
CellCycle_extrinsic_tab$is_cell_cycle[CellCycle_extrinsic_tab$is_cell_cycle==FALSE]<-"Non-cell-cycle genes"
CellCycle_extrinsic_tab$is_cell_cycle<-factor(CellCycle_extrinsic_tab$is_cell_cycle,
                                              levels = c("Cell-cycle genes","Non-cell-cycle genes"),ordered = TRUE)
CellCycle_extrinsic_tab$f12<-interaction(CellCycle_extrinsic_tab$is_cell_cycle,CellCycle_extrinsic_tab$e_or_c)
ggplot(CellCycle_extrinsic_tab,aes(y=ext,x=f12,fill=is_cell_cycle))+
  scale_fill_manual(values=c("#CC79A7", "#56B4E9"))+
  geom_boxplot(outlier.size=0.1,fatten=0.5)+
  xlab(label="")+
  ylab(label="Extrinsic noise")+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black"))+
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
  scale_y_continuous(limits=c(-3000,6000),breaks=c(-2000,0,2000,4000))

#figure5 c~d
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

Mouse_cell_cycle<-Mouse_GO_term%>%
  filter(GO.term.name=="cell cycle")
cell_cycle_index<-match(Mouse_genes_all$Gene.name,Mouse_cell_cycle$Gene.name)
Mouse_genes_all$cell_cycle<-(!is.na(cell_cycle_index))

cycle_genes<-Mouse_genes_all%>%
  filter(cell_cycle)
cycle_control<-Mouse_genes_all%>%
  filter(!cell_cycle)
bins<-quantile(cycle_genes$mean_exp,prob=c(1:51)/51)
cycle_stratified<-split(cycle_control,cut(cycle_control$mean_exp,breaks=bins))
genes_each_bin<-numeric(50)
for(i in 1:50){
  genes_each_bin[i]<-length(cycle_stratified[[i]]$Gene.name)
}
sam_per_bin<-min(genes_each_bin)
stratified_cycle_control<-cycle_control[FALSE,]
names(stratified_cycle_control)<-names(cycle_control)
set.seed(8)
for(i in 1:50){
  sub_tab_index<-sample(c(1:genes_each_bin[i]),sam_per_bin)
  sub_tab<-cycle_stratified[[i]][sub_tab_index,]
  stratified_cycle_control<-rbind(stratified_cycle_control,sub_tab)
}


df <- data.frame(genes=c("cycle_genes", "Control"),
                 ratio=c(438/570,4006/5600))
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
  scale_x_discrete(labels=c("Cell-cycle\n genes","Non-cell-cycle \n genes (stratified)"))+
  scale_x_discrete(breaks=NULL)+
  theme(legend.position = "none")

#Figure 5d
dat_cycle<-rbind(stratified_cycle_control,cycle_genes)
dat_cycle$is_cycle<-ifelse(dat_cycle$cell_cycle,"cycle genes","Background")
colorder<-c("cycle genes","Background")
ggplot(dat_cycle,aes(factor(is_cycle),target_number,fill=is_cycle))+
  geom_boxplot(outlier.shape = NA,fatten=0.5)+
  scale_fill_manual(values=c("#F0E442","#CC79A7"))+
  xlab(label="")+
  ylab(label="Number of miRNA species")+
  theme(axis.text.x=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.text.y=element_text(size=12,family="Times New Roman",color="black"))+
  theme(axis.title.y=element_text(size=12,family="Times New Roman",color="black",angle=90,vjust = 0.5))+
  theme_linedraw()+
  coord_cartesian(ylim=quantile(dat_cycle$target_number,c(0.1,0.99)))+
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  # theme(axis.ticks=element_line(size=1.5))+
  # theme(axis.ticks.length=unit(5,"mm"))+
  scale_x_discrete(limits=colorder,labels=c("Cell-cycle\n genes","Non-cell-cycle\n genes(stratified)"))+
  scale_x_discrete(limits=colorder,breaks=NULL)+
  theme(legend.position="none")
