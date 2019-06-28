
# Spearman's correlation between differential species and metabolites
# Heatmap
remove(list = ls(all = TRUE))
setwd("H:/AAAAA20180917")

Abundance_S=read.csv("microbiota/SpeciesRarefied.csv",header=T,row.names=1)# 37个样本，有37个门
# Noise removal function - note this can be increased to focus on the organisms
# most relevant to distance and PCoA calculation (e.g. set to 1% to get a few species)
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)}
Abundance_S <- noise.removal(Abundance_S, percent=0.01) # 去除noise（cutoff = 0.001%）
A=c()
for (i in 1:nrow(Abundance_S)){
  if (length(which(Abundance_S[i,]==0))/length(Abundance_S[i,]) <0.5)
    A = c(A,i)}
Abundance_S[A,] -> Abundance_S1 #################


diffSpecies=read.csv("Diff_species_ko2/Diff_species_DESeq2_STRAINs.csv",header=T,row.names=1) # 36
rownames(diffSpecies) <- diffSpecies$Species

diffSpecies_MMenriched=diffSpecies[which(diffSpecies$log2FoldChange >1),] # 20
Abundance_diffSpecies_MMenriched=Abundance_S1[rownames(diffSpecies_MMenriched),]
diffSpecies_HCenriched=diffSpecies[which(diffSpecies$log2FoldChange < -1),] # 16
Abundance_diffSpecies_HCenriched=Abundance_S1[rownames(diffSpecies_HCenriched),]
diffSpecies=rbind(diffSpecies_MMenriched,diffSpecies_HCenriched) # 28MM+16HC
Abundance_diffSpecies = Abundance_S1[rownames(diffSpecies),]

data_Sp = t(Abundance_diffSpecies)
data_Sp=data_Sp[order(rownames(data_Sp)),]
data_Sp1=data_Sp[,c(1,2,3,4,5,6,7,8,11,12,13,28,29,33)]


# Read the abundance of the differential metabolites in peripheral blood, PBmet
Abundance_diffPBmet=read.csv("Diff_Met/Abundance_diffPBmet.csv",header=T,row.names=1)
Abundance_diffPBmet=Abundance_diffPBmet[order(Abundance_diffPBmet$Class),]
Abundance_diffPBmet=Abundance_diffPBmet[,-1]

# Calculate the spearman's correlation between MM-enriched(20) species, HC-riched(16) species and PBmet(26)
# Heatmap
data_PBmet = t(Abundance_diffPBmet)
# data_Sp = t(Abundance_diffSpecies)

R=data.frame();P=data.frame()
library("psych")
for (i in 1:ncol(data_PBmet)){
  for (j in 1:ncol(data_Sp1)){
    res <- corr.test(data_PBmet[,i],data_Sp1[,j],method="spearman",adjust="BH")
    R[i,j]=res$r
    P[i,j]=res$p
   
  }
}
dimnames(R) = list(colnames(data_PBmet),colnames(data_Sp1))
dimnames(P) = list(colnames(data_PBmet),colnames(data_Sp1))

# annotation_col1 = data.frame(Family=as.factor(rep(c("Enterobacteriaceae","Streptococcaceae","unclassified_Clostridiales",
#                                                     "Prevotellaceae","Coriobacteriaceae","Bifidobacteriaceae","Lachnospiraceae",
#                                                     "Clostridiaceae","Peptostreptococcaceae","Fusobacteriaceae"),c(6,8,1,3,1,4,5,5,2,1))))
# 
# annotation_col2 = data.frame(Phylum=as.factor(rep(c("Proteobacteria","Firmicutes","Bacteroidetes","Actinobacteria",
#                                                     "Firmicutes","Fusobacteria"),c(6,9,3,5,12,1))))

annotation_col1 = data.frame(Family=as.factor(rep(c("Enterobacteriaceae","Streptococcaceae","Lachnospiraceae","Clostridiaceae"),c(6,5,1,2))))
annotation_col2 = data.frame(Phylum=as.factor(rep(c("Proteobacteria","Firmicutes"),c(6,8))))
annotation_col = cbind(annotation_col1,annotation_col2)
rownames(annotation_col) = colnames(R)

annotation_row = data.frame(Class=as.factor(c(rep("Amino_acids",13),rep("Carbohydrates",6),
                                              rep("Fatty_acids",2),rep("Organic_acids",5))))
rownames(annotation_row) = rownames(R)

ann_colors = list(Class=c(Amino_acids="Yellow1",Carbohydrates="Purple1",Fatty_acids="VioletRed1",Organic_acids="Green1"),
                  Phylum=c(Proteobacteria="#D55E00",Firmicutes="#E69F00",Bacteroidetes="#56B4E9",
                           Actinobacteria="#009E73",Fusobacteria="#0072B2"),
                  Family=c(Enterobacteriaceae="#D02090",Streptococcaceae="#D2691E",unclassified_Clostridiales="#B0B0B0",
                           Prevotellaceae="#CDCD00",Coriobacteriaceae="#CDB5CD",Bifidobacteriaceae="#CD950C",
                           Lachnospiraceae="#C1FFC1",Clostridiaceae="#C0FF3E",
                           Peptostreptococcaceae="#BA55D3",Fusobacteriaceae="#CDC673"))


c1=colorRampPalette(c("blue","white"))(-1*min(R)*1000)
c2=colorRampPalette(c("white","red"))(max(R)*1000)

library(pheatmap)
A=pheatmap(R,cluster_col=T,cluster_rows=T,treeheight_row=35,treeheight_col=35,
           cellwidth=12,cellheight=8,color=c(c1,c2),legend=T,angle_col=315,
           annotation_row=annotation_row,annotation_col=annotation_col,annotation_colors=ann_colors,
           cutree_rows=1,cutree_cols=1,number_color="black",
           fontsize_row=8,fontsize_col=8,fontsize=8,fontsize_number=7,cex=1,
           clustering_distance_rows="correlation",clustering_distance_cols="correlation",
           display_numbers=matrix(ifelse(P<0.001,"#",ifelse(P<0.01,"+",ifelse(P<0.05,"*",""))),nrow(P)),
           filename = "Interaction/Sp_PBmet_2_1.pdf",width=10,height=6)





