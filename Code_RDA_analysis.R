

## Redundancy analysis，RDA (Different species and different metabolites)
remove(list = ls(all = TRUE))
setwd("H:/AAAAA20180917")

# Read species abundance
Abundance_S=read.csv("microbiota/SpeciesRarefied.csv",header=T,row.names=1)# 37个样本，有37个门
# Noise removal function - note this can be increased to focus on the organisms
# most relevant to distance and PCoA calculation (e.g. set to 1% to get a few species)
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)}
Abundance_S <- noise.removal(Abundance_S, percent=0.01) # remove noise（cutoff = 0.001%）
A=c()
for (i in 1:nrow(Abundance_S)){
  if (length(which(Abundance_S[i,]==0))/length(Abundance_S[i,]) <0.5) #
    A = c(A,i)}
Abundance_S[A,] -> Abundance_S1 #################

diffSpecies=read.csv("Diff_species_ko2/Diff_species_DESeq2_STRAINs.csv",header=T,row.names=1) # 36
rownames(diffSpecies) <- diffSpecies$Species

diffSpecies_MMenriched=diffSpecies[which(diffSpecies$log2FoldChange >1),] # 20
Abundance_diffSpecies_MMenriched=Abundance_S1[rownames(diffSpecies_MMenriched),]
diffSpecies_HCenriched=diffSpecies[which(diffSpecies$log2FoldChange < -1),] # 16
Abundance_diffSpecies_HCenriched=Abundance_S1[rownames(diffSpecies_HCenriched),]
diffSpecies=rbind(diffSpecies_MMenriched,diffSpecies_HCenriched) # 28MM+17HC
Abundance_diffSpecies = Abundance_S1[rownames(diffSpecies),]
data_Sp = t(Abundance_diffSpecies)
data_Sp=data_Sp[order(rownames(data_Sp)),]

data_Sp1=data_Sp[,c(1,2,3,4,5,6,7,8,11,12,13,28,29,33)]

# Read PBmet abundance
Abundance_diffPBmet=read.csv("metabolome/Diff_Met/Abundance_diffPBmet_123.csv",header=T,row.names=1)
Abundance_diffPBmet=Abundance_diffPBmet[,-1]
data_PBmet = t(Abundance_diffPBmet)
data_PBmet=data_PBmet[order(rownames(data_PBmet)),]

data_PBmet1=data_PBmet[,c(2,3,4,6,8,9,13,14,15,16,17,19,25)]

# CCA  
library(vegan)
vare.cca <- cca(data_Sp1,data_PBmet1)
vare.cca
plot(vare.cca)






