#
remove(list = ls(all = TRUE))
setwd("H:/AAAAA20180917")

# Peripheral blood metabolome of all samples
Abundance_PBmet <- read.csv("metabolome/Metabolome_PB1.csv",header=T,row.names=1)

# In total, 141 metabolites were taken into consideration; 
# To identify the significantly differential metabolites using T-test
data = Abundance_PBmet
tt=NULL
for (i in 1:141){
  tt$mean_HC[i] = mean(as.numeric(data[1:18,i]))
  tt$sd_HC[i] = sd(as.numeric(data[1:18,i]))
  tt$mean_MM[i] = mean(as.numeric(data[19:37,i]))
  tt$sd_MM[i] = sd(as.numeric(data[19:37,i]))
  # log2FC=log2(mean(as.numeric(data[19:37,i]))/mean(as.numeric(data[1:18,i]))) # log2FC
  # tt$log2FC[i]=log2FC
  
  my_test<- t.test(as.numeric(as.numeric(data[19:37,i]),data[1:18,i]),alternative="two.sided",
                   paired=FALSE,var.equal=F)
  tt$p_value[i] = my_test$p.value
}
T_result=as.data.frame(tt)
rownames(T_result)=colnames(data)
T_result=T_result[order(rownames(T_result)),]
T_result1 = T_result[which(T_result$p_value <0.05),] #  & abs(T_result$log2FC) >1



# Read the data of VIP derived from SIMCA
# To identify metabolites with VIP>1 and p_value<0.05 (ttest)
VIP_PBmet <- read.table("metabolome/VIP__OPLS-DA_PBmet_Par.txt",sep='\t',header=T,row.names=1)
VIP_PBmet1 <- VIP_PBmet[which(VIP_PBmet$M1.VIP.1. >1),]
# VIP_PBmet2 <- VIP_PBmet1[which(VIP_PBmet1$M1.VIP.1. > VIP_PBmet1$X1.89456..M1.VIP.1.cvSE),] # 20
much = intersect(rownames(VIP_PBmet1),rownames(T_result1))
Result = cbind.data.frame(T_result1[much,],VIP_PBmet1[much,])


Abund_PBmet <- read.csv("metabolome/Metabolome_PB2_name.csv",header=T,row.names=4)
Diff_mets_name = Abund_PBmet[rownames(Result),]
RESULT = cbind.data.frame(Result,Diff_mets_name)
write.csv(RESULT,"metabolome/Diff_Met/Diff_PBmet.csv") #



# Heatmap, 26 differential metabolites
library(pheatmap)
Abundance_PBmet <- read.csv("metabolome/Metabolome_PB1.csv",header=T,row.names=1)
Abundance_PBmet <- t(Abundance_PBmet)
Diff_mets <- read.csv("metabolome/Diff_Met/Diff_PBmet.csv",header=T,row.names=1)
Abundance_Diff_mets <- Abundance_PBmet[rownames(Diff_mets),]
data = t(scale(t(Abundance_Diff_mets),center=T,scale=T))
rownames(data) <- Diff_mets$Name

annotation_col = data.frame(Group=as.factor(c(rep("HC",18),rep("MM",19))))
rownames(annotation_col) = colnames(data)

annotation_row = data.frame(Class=as.factor(c(rep("Amino_acids",13),rep("Carbohydrates",6),
                                              rep("Fatty_acids",2),rep("Organic_acids",5))))
                                            
rownames(annotation_row) = rownames(data)

ann_color = list(Group=c(MM="red",HC="blue"),
                 Class=c(Amino_acids="Yellow1",Carbohydrates="Purple1",
                         Fatty_acids="VioletRed1",Organic_acids="Green1"))


c1=colorRampPalette(c("blue","white"),bias=1)(-1*min(data)*100)
c2=colorRampPalette(c("white","red"),bias=1)(max(data)*100)

pheatmap(data,cellwidth=5,cellheight=8,border_color=NA,cluster_col=F,cluster_rows=TRUE,
         fontname="Times New Roman",fontsize=8,treeheight_row=20,cutree_row=2,cutree_cols=2,
         show_colnames=F,fontsize_row=8,fontsize_col=8,
         annotation_col = annotation_col,annotation_row = annotation_row,annotation_colors=ann_color,
         color=c(c1,c2),fontsize_number=8,scale="none",tl.srt=45,
         filename = "metabolome/Diff_Met/Diff_PBmet.pdf",width=8,height =5)



















