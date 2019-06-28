
# Step1: Analysis of data at phylum level，Top4 plylum
remove(list = ls(all = TRUE))
setwd("H:/AAAAA20180917")

Abundance_P=read.csv("microbiota/PhylumRarefied.csv",header=T,row.names=1) # 37 samples

# The top plylum in HC and MM group
data0 <- Abundance_P
gp=as.factor(c(rep("HC",18),rep("MM",19))) # 
data1=data0[,gp=="HC"]; data2=data0[,gp=="MM"]
S1=rowMeans(data1); S2=rowMeans(data2)
SS=cbind(S1,S2);colnames(SS)=c("HC","MM")
SS=SS[order(SS[,1],decreasing =T),] # 
SSS=rbind(SS[1:4,1:2],apply(SS[5:nrow(SS),],2,sum))
rownames(SSS)=c("Bacteroidetes","Firmicutes","Proteobacteria","Actinobacteria","Others") # The top 4 plylum with the most abundance

SSS1 <- apply(SSS,2,function(x) x/sum(x)) # 
rn=dim(SSS1)[1]
# png(filename = "phylum_genus/Phylum.tiff",width=500,height=500,res=80)
pdf(file="phylum_genus/Phylum.pdf",width=4,height=4) #
par(mar=c(8,5,3,2))
barplot(as.matrix(SSS1*100),width=1,col=rainbow(rn),beside=FALSE,ylim =c(0,100), ylab="Percentage (%)")
legend(x=2.3,y=-38,rownames(SSS1),fill=rainbow(rn),adj=0,pt.lwd=0.5,text.width=0.5,
       xjust=1,yjust=0,cex=0.6,xpd=T,ncol=3)
dev.off()


# Top4 plylum, boxplot
top4p <- c("Bacteroidetes","Firmicutes","Proteobacteria","Actinobacteria")
Ab_top4p <- Abundance_P[top4p,]
Ab_top4p <- log10(Ab_top4p)


# Calculate the significance using two-tailed Wilcox rank sum test
T_result=NULL
data=t(Abundance_P)
for(i in 1:ncol(data)){
  wil=cbind(value=data[,i],group=as.factor(c(rep("HC",18),rep("MM",19))))
  my_test<- wilcox.test(value~group,wil,alternative ="two.sided",paired=FALSE,exact=FALSE)
  logFC=log2(mean(as.numeric(wil[19:37]))/mean(as.numeric(wil[1:18]))) # log2(MM/HC)
  T_result$p.value[i] = my_test$p.value
  T_result$logfc[i]=logFC
}
tt=as.data.frame(T_result)
rownames(tt)=colnames(data)
tt=tt[order(tt$p.value),]
fdr.w <- p.adjust(as.numeric(tt[,1]),method="BH",length(tt[,1]))
tt$fdr = fdr.w


P_top4p = tt[rownames(Ab_top4p),]$fdr
mark1<- symnum(P_top4p[1],cutpoints=c(0, 0.001,0.01,0.05, 1),symbols=c("***","**","*","."))
mark2<- symnum(P_top4p[2],cutpoints=c(0, 0.001,0.01,0.05, 1),symbols=c("***","**","*","."))
mark3<- symnum(P_top4p[3],cutpoints=c(0, 0.001,0.01,0.05, 1),symbols=c("***","**","*","."))
mark4<- symnum(P_top4p[4],cutpoints=c(0, 0.001,0.01,0.05, 1),symbols=c("***","**","*","."))


# png(filename = "phylum_genus/Top4_Phylum.tiff",width=2000,height=2000,res=400)
pdf(file="phylum_genus/Top4_Phylum.pdf",width=4,height=4) #
par(mar=c(4,1,1,7))
data = t(Ab_top4p)
dataset <- data.frame(value=data[,1], group=as.factor(c(rep("HC",18),rep("MM",19))))
boxplot(value~group,data=dataset,at = c(1,1.4), boxwex = 0.3,col=c("blue","red"),horizontal =T,yaxt='n',
        notch=F,xlim=c(0.8,4.6),ylim=c(3,7.2), xlab="lg(Abundance)",names = NA)

dataset <- data.frame(value = data[,2], group = factor(c(rep("HC",18),rep("MM",19))))
boxplot(value~group,data = dataset,at = c(2,2.4), boxwex = 0.3,col=c("blue","red"),horizontal =T,yaxt='n',
        add=T,names = NA)

dataset <- data.frame(value = data[,3], group = factor(c(rep("HC",18),rep("MM",19))))
boxplot(value~group,data = dataset,at = c(3,3.4), boxwex = 0.3,col=c("blue","red"),horizontal =T,yaxt='n',
        add=T,names = NA)

dataset <- data.frame(value = data[,4], group = factor(c(rep("HC",18),rep("MM",19))))
boxplot(value~group,data = dataset,at = c(4,4.4), boxwex = 0.3,col=c("blue","red"),horizontal =T,yaxt='n',
        add=T,names = NA)

axis(4,c(1.2,2.2,3.2,4.2),c("Bacteroidetes","Firmicutes","Proteobacteria","Actinobacteria"),
                                   tick = T,las =2,col.axis = "black",cex.axis=0.85)
text(7.1,1.2,mark1,col = "black",cex = 2.5)
text(7.1,2.2,mark2,col = "black",cex = 2.5)
text(7.1,3.2,mark3,col = "black",cex = 2.5)
text(7.1,4.2,mark4,col = "black",cex = 2.5)

dev.off()





# Step2: Analysis of data at phylum level，Top4 plylum
remove(list = ls(all = TRUE))
setwd("H:/AAAAA20180917")

Abundance_G=read.csv("microbiota/GenusRarefied.csv",header = T,row.names = 1) # 

# Noise removal function - note this can be increased to focus on the organisms
# most relevant to distance and PCoA calculation (e.g. set to 1% to get a few species)
noise.removal <- function(dataframe, percent=0.001, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)}
Abundance_G1 <- noise.removal(Abundance_G, percent=0.001) # remove noise（cutoff = 0.001%）


# To get the top 30 genus (HC + MM)，draw boxplot
data=Abundance_G1  ######
data1=apply(data,1,mean)
data1 = data.frame(data1)
data1 = cbind.data.frame(data1,data1)
data2 = data1[order(data1$data1,decreasing=T),]
data2 <- apply(data2,2,function(x) x/sum(x)) #
data2 = data.frame(data2)
data2$order = c(1:nrow(data2))

# sum(data2$data1[1:30]) # 0.9587929
Data30g = data[rownames(data2)[1:30],]
write.csv(Data30g,"phylum_genus/Abundance_top30_genus.csv")


# Calculate the significance using two-tailed Wilcox rank sum test
data=t(Abundance_G1)
T_result=NULL
for(i in 1:ncol(data)){
  wil=cbind(value=data[,i],gp=as.factor(c(rep("HC",18),rep("MM",19))))
  my_test<- wilcox.test(value~gp,wil,alternative ="two.sided",paired=FALSE,exact=FALSE)
  logFC=log2(mean(as.numeric(wil[19:38]))/mean(as.numeric(wil[1:18]))) # log2(MM/HC)
  T_result$p.value[i] = my_test$p.value
  T_result$logfc[i]=logFC}

tt=as.data.frame(T_result)
rownames(tt)=colnames(data)
tt=tt[order(tt$p.value),]
fdr.w <- p.adjust(as.numeric(tt[,1]),method="BH",length(tt[,1]))
tt$fdr = fdr.w



Ab_top30g = t(log10(Data30g))
P_top30g = tt[colnames(Ab_top30g),]$p.value

# png(filename = "phylum_genus/Top30_Genus_pvalue.tiff",width = 3000,height = 4000,res = 400)
pdf(file="phylum_genus/Top30_Genus_pvalue.pdf",width=3.9,height=5.5) #
par(mar=c(4,1,1,8))

dataset <- data.frame(value=Ab_top30g[,1],group=factor(c(rep("HC",18),rep("MM",19))))
boxplot(value~group,data=dataset,at=c(1,1.4),boxwex=0.3,col=c("blue","red"),horizontal=T,notch=F,yaxt='n',
        xlim=c(1.6,29.8),ylim=c(1.5,7.1),xlab="lg (Abundance)",names=NA) # expression("log"[10]*"(Abundance)")

mark<- symnum(P_top30g[1],cutpoints=c(0, 0.001,0.01,0.05, 1),symbols=c("***","**","*","."))
text(7,1.2,mark,col = "black",cex = 1.2)

for (i in 2:30){
  dataset <- data.frame(value = Ab_top30g[,i], group=factor(c(rep("HC",18),rep("MM",19))))
  boxplot(value~group,data = dataset,at = c(i,i+ 0.4), boxwex = 0.3,col=c("blue","red"),horizontal =T,yaxt='n',
          add=T,names = NA)
  mark<- symnum(P_top30g[i],cutpoints=c(0, 0.001,0.01,0.05, 1),symbols=c("***","**","*","."))
  text(7,i+0.2,mark,col = "black",cex = 1.2)
}

axis(4,c(1.2:30.2),colnames(Ab_top30g),tick = T,las =2,col.axis = "black",cex.axis=0.9)
dev.off()



# install.packages("Hmisc")
# library(Hmisc)
# Calculate the spearman's correlation between the top 30 genus(30)
data = Data30g # Abundance_30g, for spearman, cytoscape
library("psych")
my_data = t(data)
result <- corr.test(as.matrix(my_data),method="spearman",adjust="BH")

# A simple function to format the correlation matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
Result <- flattenCorrMatrix(result$r, result$p)

Result1 = Result[which(abs(Result$cor) >0.5 & Result$p <0.01),]
Result1$log10p <- log10(Result1$p)* -1
write.csv(Result1,"phylum_genus/top30g_spearman_P0.01_cytoscape.csv") # for cytoscape 





