
# 
remove(list = ls(all = TRUE))
setwd("H:/AAAAA20180917")

# Calculate the alpha-diversity of all samples at species level
# Shannon, Simpson, Fisher diversity indices, Species richness(S)and Pielou's evenness(J)

# Read the abundance of all samples
Abundance_S=read.csv("microbiota/SpeciesRarefied.csv",header=T,row.names=1)# 37 samples

library(vegan)
Species = t(Abundance_S) ######

H <- data.frame(diversity(Species))
simpson <- data.frame(diversity(Species, "simpson")) 
shannon<-data.frame(diversity(Species, "shannon"))
invsimp <- data.frame(diversity(Species, "inv"))
alpha <- data.frame(fisher.alpha(Species))
S <- data.frame(specnumber(Species))
J <- data.frame(H/log(S))
Diversity<-cbind(simpson, shannon, invsimp, alpha, S, J)
Diversity$Sample<-row.names(Diversity)
colnames(Diversity)<-c("Simpson","Shannon","InvSimpson","Alpha","SpeciesNo","Evenness","Sample")
Diversity<-Diversity[,c(7,1,2,3,4,5,6)]
Diversity$Class<-c(rep("HC",18),rep("MM",19))
write.csv(Diversity,"Species analysis/SampleDiversity.csv")

div = Diversity[,c("Simpson","Shannon","Alpha","SpeciesNo","Evenness")]
div$class = as.factor(Diversity$Class)

pdf(file="Species analysis/Sample_diversity.pdf",width=12,height=12) #
#tiff(file="Species analysis/Sample_diversity.tiff",width=2000,height=2000,res=200) #
pairs(div[1:5],pch=21,bg=c("red","blue")[unclass(div$class)])
dev.off()
#dev.off()


# Draw violin plot on the basis of shannon index
library("vioplot")
wil_test_shannon<-wilcox.test(Diversity$Shannon[1:18],Diversity$Shannon[19:37],alternative="less",paired=FALSE)
mark_shannon<- symnum(wil_test_shannon$p.value,cutpoints=c(0, 0.001,0.01,0.05, 1),symbols=c("***","**","*",".")) # 0.04595955
# png(filename = "Species analysis/Shannon-Wiener_species_violin.tiff",width = 1000,height = 2000,res = 400)
pdf(file="Species analysis/Shannon-Wiener_species_violin.pdf",width=2.5,height=5) #
x=1:1;y=1:1
plot(x,y,xlim=c(-0.5,1.5),xlab="",cex.main=0.9,ylab="Shannon-Wiener index",main="",
     ylim=c(1.0,4.5),pch=21,col='white',xaxt="n")

for(i in 1:2){
  HC = Diversity$Shannon[1:18]
  MM = Diversity$Shannon[19:37]
  vioplot(add = T,col='blue',HC,at=3*(i-1),lty=1,border="blue",rectCol="white",colMed="black")
  vioplot(add = T,col='red',MM,at=3*(i-1)+1,lty=1,border="red",rectCol="white",colMed="black")
  text(0.5,4.4,mark_shannon,col="black",cex=2)
}
axis(side=1, at=c(0,1), labels = F, tick = TRUE)
text(c(0,1),rep(par("usr")[3]-0.3,8),xpd=NA,cex=0.8,labels=c('HC','MM'))
text(0.5,4.8,xpd=NA,cex=1.1,labels=c('Species α-diversity'))
dev.off()





# Calculate the alpha-diversity of all samples at genus level
# shannon index
Abundance_G=read.csv("microbiota/GenusRarefied.csv",header = T,row.names = 1)#
Shannon_g<-data.frame(diversity(t(Abundance_G), "shannon")) # shannon index
wil_test_shannon_g<-wilcox.test(Shannon_g[1:18,1],Shannon_g[19:37,1],alternative="less",paired=FALSE)
mark_shannon<- symnum(wil_test_shannon_g$p.value,cutpoints=c(0, 0.001,0.01,0.05, 1),symbols=c("***","**","*",".")) # 0.04303628

# Draw violin plot on the basis of shannon index 
# png(filename = "Species analysis/Shannon-Wiener_genus_violin.tiff",width=1000,height=2000,res=400)
pdf(file="Species analysis/Shannon-Wiener_genus_violin.pdf",width=2.5,height=5) #
x=1:1;y=1:1
plot(x,y,xlim=c(-0.5,1.5),xlab="",cex.main=0.9,ylab="Shannon-Wiener index",main="",
                   ylim=c(0.2,3.6),pch=21,col='white',xaxt="n")
     
for(i in 1:2){
  HC = Shannon_g[1:18,1]
  MM = Shannon_g[19:37,1]
  vioplot(add = T,col='blue',HC,at=3*(i-1),lty=1,border="blue",rectCol="white",colMed="black")
  vioplot(add = T,col='red',MM,at=3*(i-1)+1,lty=1,border="red",rectCol="white",colMed="black")
  text(0.5,3.5,mark_shannon,col="black",cex=2)
}
axis(side=1, at=c(0,1), labels = F, tick = TRUE)
text(c(0,1),rep(par("usr")[3]-0.3,8),xpd=NA,cex=0.8,labels=c('HC','MM'))
text(0.5,3.9,xpd=NA,cex=1.1,labels=c('Genera α-diversity'))
dev.off()





# Rarefaction curves of all samples
# png(filename="Species analysis/Rarefaction_curves_species.tiff",width=2000,height=2000,res=400)
pdf(file="Species analysis/Rarefaction_curves_species.pdf",width=3,height=3,pointsize=6) #矢量图，pdf格式
raremax <- min(rowSums(Species))
A=rarecurve(Species[1:37,],sample=raremax,step=10000,col=c(rep("blue",18),rep("red",19)),
                xlab="Reads count",ylab="Species",label=F,cex=0.5,lwd=0.8)
dev.off()






# ade4
# library(ade4)
# library(vegan)
# Calculate the beta-diversity between the two groups of all samples

Abundance_S=read.csv("microbiota/SpeciesRarefied.csv",header=T,row.names=1)# 37 samples
data = t(Abundance_S)
gp=as.factor(c(rep("HC",18),rep("MM",19)))
dis<-vegdist(wisconsin(sqrt(data)),method="bray")

# anosim(dis,gp,permutations=999,distance="bray") # 0.001, Analysis of similarities (ANOSIM)
adonis2(dis~gp,gp,permutations=999,method="bray") # R2=2.2149, 0.001 ***

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
pdf(file="Species analysis/pcoa_species.pdf",width=5,height=5) #
# png(filename="Species analysis/pcoa_species.tiff",width=2000,height=2000,res=400)
plot(pca)  
dev.off()





