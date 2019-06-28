

remove(list = ls(all = TRUE))
setwd("H:/AAAAA20180917")

# Noise removal function - note this can be increased to focus on the organisms
# most relevant to distance and PCoA calculation (e.g. set to 1% to get a few species)

Abundance_S=read.csv("microbiota/SpeciesRarefied.csv",header=T,row.names=1)# 37 samples
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)}
Abundance_S <- noise.removal(Abundance_S, percent=0.01) # 去除noise（cutoff = 0.001%）

A=c()
for (i in 1:nrow(Abundance_S)){
  if (length(which(Abundance_S[i,]==0))/length(Abundance_S[i,]) <0.5) # 控制样本中0的个数
    A = c(A,i)}
Abundance_S[A,] -> Abundance_S1 ################# 160


# To identify differential species between the two groups using DESeq2
library(DESeq2)
data=Abundance_S1 ##### only accept read count
group=as.factor(c(rep("HC",18),rep("MM",19)))
conditions=data.frame(colnames(Abundance_S),group)
dds <- DESeqDataSetFromMatrix(Abundance_S, colData=conditions, design= ~group)
dds <- dds[ rowSums(counts(dds)) > 10, ]
dds2 <- DESeq(dds)
res<-results(dds2,contrast=c("group","MM","HC")) # 
deg2=as.data.frame(res)
result_DESeq2=deg2[which(deg2$padj<0.01 & abs(deg2$log2FoldChange)>1),] #
result_DESeq2=result_DESeq2[order(result_DESeq2$log2FoldChange,decreasing=T),]
write.csv(result_DESeq2,file="Diff_species_ko2/Diff_species_DESeq2.csv") # 




## Heatmap, 36 differential species were identified and focused, including 20 MM-enriched species and 16 HC-enriched species
result_DESeq2=read.csv("DIFF_species_ko2/Diff_species_DESeq2_STRAINs.csv",header=T,row.names=4)

x=c(1:36)
y=result_DESeq2$log2FoldChange
se=result_DESeq2$lfcSE
yplus=y-se
yminus=y+se
FC=cbind.data.frame(x,y,yplus,yminus)


pdf(file="DIFF_species_ko2/diff_error_2_1111.pdf",width=5,height=3) #
library(Hmisc)
plot(x,c(rep(0,36)),type="l",lty=4,ylim=c(-6,6),xlim=c(1,36),ylab=expression("log"[2]*"(Abundance)"))
errbar(x,y,yplus,yminus,add=T)
dev.off()


diff_sp = rownames(result_DESeq2)
Abund_diff_sp0 = Abundance_S1[diff_sp,]
Abund_diff_sp = log10(Abund_diff_sp0) #
Abund_diff_sp[Abund_diff_sp== -Inf] <- 0

Abund_diff_sp1 = t(scale(t(Abund_diff_sp),center=T,scale=T))

# display_numbers=matrix(ifelse(P<0.001,"***",ifelse(P<0.01,"**",ifelse(P<0.05,"*",""))),nrow(P)),

annotation_col = data.frame(Group=as.factor(c(rep("HC",18),rep("MM",19))))
rownames(annotation_col) = colnames(Abund_diff_sp1)

annotation_row1 = data.frame(Family=as.factor(rep(c("Enterobacteriaceae","Streptococcaceae","unclassified_Clostridiales",
                                            "Prevotellaceae","Coriobacteriaceae","Bifidobacteriaceae","Lachnospiraceae",
                                  "Clostridiaceae","Peptostreptococcaceae","Fusobacteriaceae"),c(6,8,1,3,1,4,5,5,2,1))))

annotation_row2 = data.frame(Phylum=as.factor(rep(c("Proteobacteria","Firmicutes","Bacteroidetes","Actinobacteria",
                                                   "Firmicutes","Fusobacteria"),c(6,9,3,5,12,1))))

annotation_row = cbind(annotation_row1,annotation_row2)
rownames(annotation_row) = rownames(Abund_diff_sp1)

ann_colors = list(Group=c(HC="blue2",MM="red2"),
                  Phylum=c(Proteobacteria="#D55E00",Firmicutes="#E69F00",Bacteroidetes="#56B4E9",
                           Actinobacteria="#009E73",Fusobacteria="#0072B2"),
                  Family=c(Enterobacteriaceae="#D02090",Streptococcaceae="#D2691E",unclassified_Clostridiales="#B0B0B0",
                           Prevotellaceae="#CDCD00",Coriobacteriaceae="#CDB5CD",Bifidobacteriaceae="#CD950C",
                           Lachnospiraceae="#C1FFC1",Clostridiaceae="#C0FF3E",
                           Peptostreptococcaceae="#BA55D3",Fusobacteriaceae="#CDC673"))


c1=colorRampPalette(c("blue","white"))(-min(Abund_diff_sp1)*1000)
c2=colorRampPalette(c("white","red"))(max(Abund_diff_sp1)*1000)

library(pheatmap)
pheatmap(Abund_diff_sp1,cellwidth=4,cellheight=6.5,border_color=NA,cluster_rows=F,cluster_col=F,
        treeheight_row=50,cutree_row=2,cutree_cols=1,scale="none",show_rownames=T,show_colnames=F,
        fontsize=6,fontsize_row=6,fontsize_col=6,fontsize_number=6,clustering_distance_rows="correlation",
        color=c(c1,c2),annotation_col = annotation_col,annotation_row = annotation_row,annotation_colors=ann_colors,
        filename = "Diff_species_ko2/diff_1111_111.pdf",width=10,height=6)



          



