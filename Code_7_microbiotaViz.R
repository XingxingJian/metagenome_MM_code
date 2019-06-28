
remove(list = ls(all = TRUE))
setwd("H:/AAAAA20180917")

# Install R package microbiomeViz
# options(download.file.method="libcurl")
# devtools::install_github("lch14forever/microbiomeViz") #
library(microbiomeViz)
# 
allData0 <- read.csv("microbiota/microbiomeViz.csv",header=T,row.names=1)

# Noise removal function - note this can be increased to focus on the organisms
# most relevant to distance and PCoA calculation (e.g. set to 1% to get a few species)
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)}
allData <- noise.removal(allData0, percent=0.001) # remove noise（cutoff = 0.001%）

A=c()
for (i in 1:nrow(allData)){
  if (length(which(allData[i,]==0))/length(allData[i,]) <0.5) # 控制样本中0的个数
    A = c(A,i)}
allData[A,] -> allData1 ########
ALLDATA <- allData1

for (i in 1:37){ALLDATA[,i]=ALLDATA[,i]/ALLDATA[1,i]}
write.csv(ALLDATA,"microbiota/microbiomeViz_100%.csv")


################# run microbiomeViz ############################

df <- read.csv("microbiota/microbiomeViz_100%.csv",header=T)
# Calculate the mean to render the node size
dat <- data.frame(V1=df[,1], V2=rowMeans(df[,-1]), stringsAsFactors=FALSE)
dat[,2] <- dat[,2]*100

# To generate tree skeletons on the basis of OTU and abundance 
tree <- parseMetaphlanTSV(dat,node.size.offset=3,node.size.scale=1.5)
p <- tree.backbone(tree,size=0.71,layout="circular",shape=21,
                   fill="white",color="black")

p

# annotation
lefse_lists = data.frame(node=c("s__Enterobacter_cloacae","s__Citrobacter_freundii",
              "s__Klebsiella_aerogenes","s__Klebsiella_pneumoniae","s__Klebsiella_variicola",
              "s__Streptococcus_oralis","s__Streptococcus_gordonii","s__Streptococcus_salivarius",
              "s__Streptococcus_mitis","s__Streptococcus_pneumoniae","s__Raoultella_ornithinolytica",
              "g__Klebsiella","g__Streptococcus","f__Enterobacteriaceae","f__Streptococcaceae","g__Clostridium",
              "s__Anaerostipes_hadrus","s__Clostridium_saccharobutylicum","s__Clostridium_butyricum"),
                         color=c(rep('red',15), rep('blue',4)),
                         stringsAsFactors = FALSE)

p <- clade.anno(p,lefse_lists,alpha=0.35,anno.x=0,anno.y=40)
p















# # example:
# df<-read.table("http://bailab.genetics.ac.cn/markdown/R/microbiomeViz/merged_abundance_table.txt",
#                head=TRUE, stringsAsFactors = FALSE)
# 
# # 计算均值用于呈现结点大小
# dat <- data.frame(V1=df[,1], V2=rowMeans(df[,-1]), stringsAsFactors = FALSE)
# 
# # 用物种和丰度生成树骨架
# tr <- parseMetaphlanTSV(dat, node.size.offset=2, node.size.scale=0.8)
# p <- tree.backbone(tr, size=0.5)
# p
# 
# # 差异物种注释
# # 读取需要颜色标注的差异物种列表，本质上是两列和颜色对应表
# lefse_lists = data.frame(node=c('s__Haemophilus_parainfluenzae','p__Proteobacteria',
#                                 'f__Veillonellaceae','o__Selenomonadales',
#                                 'c__Negativicutes', 's__Streptococcus_parasanguinis',
#                                 'p__Firmicutes','f__Streptococcaceae',
#                                 'g__Streptococcus','o__Lactobacillales',
#                                 'c__Bacilli','s__Streptococcus_mitis'),
#                          color=c(rep('darkgreen',6), rep('red',6)),
#                          stringsAsFactors = FALSE)
# 
# # 注释树
# p <- clade.anno(p, lefse_lists, alpha=0.3)
# p




