

### 筛选差异ko (level 3/4)
remove(list = ls(all = TRUE))
setwd("H:/AAAAA20180917")

Abundance_ko0<-read.csv("metacv/Abundance_ko_level4_Rarefied.csv",header=T,row.names=1)
# Abundance_ko<-read.csv("metacv/Abundance_eggNOG_Rarefied.csv",header=T,row.names=1)

# Noise removal function - note this can be increased to focus on the organisms
# most relevant to distance and PCoA calculation (e.g. set to 1% to get a few species)
noise.removal <- function(dataframe, percent=0.0001, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)}
Abundance_ko <- noise.removal(Abundance_ko0, percent=0.001) # remove noise（cutoff = 0.001%）

A=c()
for (i in 1:nrow(Abundance_ko)){
  if (length(which(Abundance_ko[i,]==0))/length(Abundance_ko[i,]) <0.5) #
    A = c(A,i)}
Abundance_ko[A,] -> Abundance_ko1 # 

# To identify differential KOs between the two groups using DESeq2

library(DESeq2)
data=Abundance_ko1 ##### read count
group=as.factor(c(rep("HC",18),rep("MM",19)))
conditions=data.frame(colnames(data),group)
dds <- DESeqDataSetFromMatrix(data, colData=conditions, design= ~group)
dds <- dds[ rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res <- results(dds)
deg2=as.data.frame(res)
result_DESeq2=deg2[which(deg2$padj<0.05 & abs(deg2$log2FoldChange)>1),] # 
# result_DESeq2_abundance = Abundance_ko1[rownames(result_DESeq2),]
# write.csv(result_DESeq2,file="20181209_ko/Diff_ko_DESeq2_Rarefied.csv") # 





### Calculate the beta-diversity between the two groups of all samples (Function, KO)
## PCOA

remove(list = ls(all = TRUE))
setwd("H:/AAAAA20180917")
library(ade4)
library(vegan)
Abundance_ko<-read.csv("metacv/Abundance_ko_level4_Rarefied.csv",header=T,row.names=1)
data = t(Abundance_ko)
gp=as.factor(c(rep("HC",18),rep("MM",19)))
dis<-vegdist(wisconsin(sqrt(data)),method="bray")

# anosim(dis,gp,permutations=999,distance="bray") # 0.002, Analysis of similarities (ANOSIM)
adonis2(dis~gp,gp,permutations=999,method="bray") # p=0.003 **
# plot PCoA
# library(GUniFrac) # Unifrac 
library(ape) # pcoa
library(ggplot2)
# dis<-vegdist(wisconsin(sqrt(data)),method="bray") # 计算 Bray–Curtis dissimilarity
PCOA <- pcoa(dis, correction="none",rn=NULL) #
result <-PCOA$values[,"Relative_eig"]  
pro1 = as.numeric(sprintf("%.3f",result[1]))*100  
pro2 = as.numeric(sprintf("%.3f",result[2]))*100  
x = PCOA$vectors  
sample_names = rownames(x)  
pc = as.data.frame(PCOA$vectors)  
pc$names = sample_names  
legend_title = ""
Group=as.factor(c(rep("HC",18),rep("MM",19)))
shape <- c("HC"=16, "MM"=17) #
color <- c("HC"='blue', "MM"='red') #
pc$group = Group
xlab=paste("PCoA1(",pro1,"%)",sep="")   
ylab=paste("PCoA2(",pro2,"%)",sep="")  
pca=ggplot(pc,aes(Axis.1,Axis.2)) +  # ggplot 
  geom_point(size=2,aes(color=group,shape=group)) +   
  #  geom_text(aes(label=names),size=4,vjust=-1) +  
  labs(x=xlab,y=ylab,title=" ",color=legend_title,shape=legend_title) +   
  geom_hline(yintercept=0,linetype=4,color="grey") +
  geom_vline(xintercept=0,linetype=4,color="grey") +
  scale_shape_manual(values=shape) +  
  scale_color_manual(values=color) +  
  theme_bw()

png(filename="DIFF_species_ko2/pcoa_KOs_Rarefied.tiff",width=2000,height=2000,res=400)
pdf(file="DIFF_species_ko2/pcoa_KOs_Rarefied.pdf",width=5,height=4) #
plot(pca)  
dev.off()



