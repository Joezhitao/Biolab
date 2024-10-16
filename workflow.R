## Bio_worflow

## TCGA数据整理_TPM文件
#加载R包
library(rjson)
library(tidyverse)
#设置路径
setwd("E:/CRC_model")
#CART文件名字
cart_file <- "gdc_download_20240814_054206.318162"
#json文件名字
json_file <- "metadata.cart.2024-08-14.json"
#读入meta.data文件
json <- jsonlite::fromJSON(json_file)
#View(json)
#获取样本名称及文件名称
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
#sample_id[1:10]
file_sample <- data.frame(sample_id,file_name=json$file_name)  
#View(file_sample)
#获取counts所在位置
count_file <- list.files(paste0(cart_file,'/',sep=''),
                         pattern = '*.tsv',recursive = TRUE)
#count_file[1:10]
#获取每个文件名称
count_file_name <- strsplit(count_file,split='/')
#count_file_name[1:10]
count_file_name <- sapply(count_file_name,function(x){x[2]})
#count_file_name[1:10]
#构建一个空的数据框
matrix = data.frame(matrix(nrow=60660,ncol=0))
#逐个读取及合并
for (i in 1:length(count_file)){
  path = paste0(cart_file,'//',count_file[i])   #Counts文件夹名
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  #3：unstranded,counts；4：stranded_first；5：stranded_second；6：tpm_unstranded；7：fpkm_unstranded；8：fpkm_uq_unstranded
  #data <- data[3]
  data <- data[6]
  #data <- data[7]
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}
#转化为gene_symbol
path = paste0(cart_file,'//',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
#gene_name[1:10]
matrix0 <- cbind(gene_name,matrix)
#获取基因类型
gene_type <- data[-c(1:6),2]
#gene_type[1:10]
matrix0 <- cbind(gene_type,matrix0)
#将gene_name列去除重复的基因，保留基因最大表达量结果,min,mean
matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)
table(gene_name)
#保留mRNA
matrix0 <- subset(x = matrix0, gene_type == "protein_coding")
#table(gene_type)
#将gene_name列设为行名，并转化为导出格式
rownames(matrix0) <- matrix0[,1]
matrix0 <- matrix0[,-c(1,2)]
matrix1 = data.frame(ID=rownames(matrix0),matrix0)
colnames(matrix1) = gsub('[.]', '-', colnames(matrix1))
#导出
write.table(matrix1,'TCGA_CRC_TPM.txt', sep="\t", quote=F, row.names = F)
#write.table(matrix1,'TCGA_LUSC_count.txt', sep="\t", quote=F, row.names = F)

## TCGA数据整理_COUNT文件
#加载R包
library(rjson)
library(tidyverse)
#设置路径
setwd("E:/CRC_model")
#CART文件名字
cart_file <- "gdc_download_20240814_054206.318162"
#json文件名字
json_file <- "metadata.cart.2024-08-14.json"
#读入meta.data文件
#出现报错，注意上面的包是否加载，工作目录下是否存在这个文件，文件是否损坏（先跑示例数据）
json <- jsonlite::fromJSON(json_file)
#View(json)
#获取样本名称及文件名称
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
#sample_id[1:10]
file_sample <- data.frame(sample_id,file_name=json$file_name)  
#View(file_sample)
#获取counts所在位置
count_file <- list.files(paste0(cart_file,'/',sep=''),
                         pattern = '*.tsv',recursive = TRUE)
#count_file[1:10]
#获取每个文件名称
count_file_name <- strsplit(count_file,split='/')
#count_file_name[1:10]
count_file_name <- sapply(count_file_name,function(x){x[2]})
#count_file_name[1:10]
#构建一个空的数据框
matrix = data.frame(matrix(nrow=60660,ncol=0))
#逐个读取及合并
#出现报错，注意解压方式（仔细看视频），文件名修改，文件后缀（解压软件）
for (i in 1:length(count_file)){
  path = paste0(cart_file,'//',count_file[i])   #Counts文件夹名
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  #3：unstranded,counts；4：stranded_first；5：stranded_second；6：tpm_unstranded；7：fpkm_unstranded；8：fpkm_uq_unstranded
  #data <- data[3]
  data <- data[3]
  #data <- data[7]
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}
#转化为gene_symbol
path = paste0(cart_file,'//',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
#gene_name[1:10]
matrix0 <- cbind(gene_name,matrix)
#获取基因类型
gene_type <- data[-c(1:6),2]
#gene_type[1:10]
matrix0 <- cbind(gene_type,matrix0)
#将gene_name列去除重复的基因，保留基因最大表达量结果,min,mean
matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)
table(gene_name)
#保留mRNA
matrix0 <- subset(x = matrix0, gene_type == "protein_coding")
#table(gene_type)
#将gene_name列设为行名，并转化为导出格式
rownames(matrix0) <- matrix0[,1]
matrix0 <- matrix0[,-c(1,2)]
matrix1 = data.frame(ID=rownames(matrix0),matrix0)
colnames(matrix1) = gsub('[.]', '-', colnames(matrix1))
#导出
#write.table(matrix1,'TCGA_CORD_TPM.txt', sep="\t", quote=F, row.names = F)
write.table(matrix1,'TCGA_CRC_count.txt', sep="\t", quote=F, row.names = F)

## 获取数据基因
## 通过KEGGREST包获取部分相关基因
setwd("E:/基因建模/gene/KEGG")

library(KEGGREST)

# 定义路径
pathways <- c('hsa04370','hsa04151','hsa04010','hsa04330','hsa04350','hsa04310')

# 初始化一个空的数据框用于存储所有基因信息
all_genes <- data.frame()

# 遍历每个通路
for (pathway_id in pathways) {
  # 获取通路信息
  pathway_info <- keggGet(pathway_id)
  
  # 提取基因信息
  genes <- unlist(lapply(gs[[1]]$GENE, function(x) strsplit(x, ";")[[1]]))
  
  genelist <- genes[which((seq_along(genes) %% 3) == 2)]
  
  # 将基因信息转换为数据框
  genelist <-  data.frame(genelist)
  
  # 合并到总的数据框中
  all_genes <- rbind(all_genes, genelist)
}

# 去除重复的基因
all_genes <- unique(all_genes)

# 输出到CSV文件
write.table(all_genes, "genelist.csv", row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

## TCGA临床信息整理
#加载包
library("XML")
library("methods")
library("tinyarray")
#设置路径
setwd("E:/CRC_model")
#设置一个目录
dir="E:/CRC_model/gdc_download_20240814_142022.604604"
#获取每个样本所在位置
all_fiels=list.files(path = dir ,pattern='*.xml$',recursive=T)##导入文件
cl = lapply(all_fiels, function(x){
  #读取 XML 文件  
  result <- xmlParse(file = file.path(dir,x)) 
  #获取根元素
  rootnode <- xmlRoot(result)  
  #转化为数据框
  xmldataframe <- xmlToDataFrame( rootnode[2] ) 
  #转置
  return(t(xmldataframe)) })

#由于读入的时候自动读入了一些NA值，导致存在不同的列
#因此获取共有列来保持列数相同
n = intersect_all(lapply(cl, rownames))
#提取共有列
cl3 = lapply(cl, function(x){x[n,]})
#对每个列表进行合并
clinical <- t(do.call(cbind,cl3))
#导出
# 将第12列(TCGA样本编号列)设置为行名
rownames(clinical) <- clinical[, 12]
# 从数据框中移除第十二列
#clinical <- clinical[, -12]
#导出
write.table(clinical,file="clinical.txt",sep="\t",quote=F,row.names = F)

## TCGA差异基因分析_limma
#加载包
library(limma)
library(pheatmap)
library(stringr)
library(edgeR)
library(ggrepel)
#设置工作目录
setwd("E:/CRC_model")
#读取输入文件
data=read.table("TCGA_CRC_count.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
#去除低表达的基因
data=data[rowMeans(data)>1,]
#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#将2变为1
group=gsub("2", "1", group)
#获取正常及肿瘤组样本数目
conNum=length(group[group==1])
treatNum=length(group[group==0])
Type=c(rep(1,conNum), rep(2,treatNum))
#根据正常和肿瘤排序
data1 = data[,group == 1]
data2 = data[,group == 0]
data = cbind(data1,data2)
#设置分组信息
Type = factor(Type)
design <- model.matrix(~0+Type)
rownames(design) = colnames(data)
colnames(design) <- levels(Type)
#创建对象
DGElist <- DGEList( counts = data, group = Type )
#过滤掉cpm小于等于1的基因
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2 # 自定义
table(keep_gene)
#生成DGEList对象
DGElist <- DGElist[keep_gene, , keep.lib.sizes = FALSE ]
#计算比例因子以将原始库大小转换为有效库大小。
DGElist <- calcNormFactors( DGElist )
#转换logCPM
v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")
#非线性最小二乘法
fit <- lmFit(v, design)
colnames(design)=c("normal","tumor")
#构建比较矩阵
cont.matrix <- makeContrasts(contrasts = c('tumor-normal'), levels = design)
#线性拟合模型构建
fit2 <- contrasts.fit(fit, cont.matrix)
#用经验贝叶斯调整t-test中方差的部分
fit2 <- eBayes(fit2)
nrDEG_limma_voom = topTable(fit2, coef = 'tumor-normal', n = Inf)
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
#保存全部计算的基因
nrDEG_limma_voom$gene <- rownames(nrDEG_limma_voom)

write.table(nrDEG_limma_voom,file="TCGA.all.limma.txt",sep="\t",row.names=F,quote=F)
#筛选标准
padj = 0.05
logFC = 1
#筛选
outDiff = nrDEG_limma_voom[(nrDEG_limma_voom$adj.P.Val < padj &
                              (nrDEG_limma_voom$logFC>logFC | nrDEG_limma_voom$logFC<(-logFC))),]
outDiff = outDiff[order(outDiff$logFC),]
write.table(data.frame(ID=rownames(outDiff),outDiff),file="TCGA.diff.limma.txt",sep="\t",row.names=F,quote=F)
nrDEG_limma_voom$logFC=as.numeric(as.vector(nrDEG_limma_voom$logFC))
nrDEG_limma_voom$change = ifelse(nrDEG_limma_voom$adj.P.Val < 0.05 & abs(nrDEG_limma_voom$logFC) >= 1, 
                                 ifelse(nrDEG_limma_voom$logFC> 1 ,'Up','Down'),
                                 'Stable')
#热图
#设置展示基因的数目
geneNum=50
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=as.vector(outDiff[,7])
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=log2(data[hmGene,]+0.01)
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
#pdf(file="heatmap_CRC_Top50.pdf", width=10, height=8)
png("heatmap_CRC_limma_Top50.png", width=10, height=8, units = "in", res = 800)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5,
         fontsize_col=8)
dev.off()
#火山图2,可以显示希望标记的基因,例：TUBA1C
geneList0 <- c('NELL2',"NKX2-2","HPCAL4")
for (i in geneList0) {
  geneList <- nrDEG_limma_voom[i,]
  rownames(nrDEG_limma_voom) <- nrDEG_limma_voom$gene
  sum(is.na(nrDEG_limma_voom$logFC))  # 检查 logFC 列中有多少个 NA 值
  sum(is.na(nrDEG_limma_voom$adj.P.Val))    # 检查 adj.P.Val 列中有多少个 NA 值
  sum(is.na(nrDEG_limma_voom$change)) # 检查 change 列中有多少个 NA 值
  nrDEG_limma_voom$logFC[is.na(nrDEG_limma_voom$logFC)] <- 0  # 将 logFC 列中的 NA 值替换为 0
  nrDEG_limma_voom$adj.P.Val[is.na(nrDEG_limma_voom$adj.P.Val)] <- 1      # 将 adj.P.Val 列中的 NA 值替换为 1
  nrDEG_limma_voom <- na.omit(nrDEG_limma_voom)
  p <- ggplot(# 数据、映射、颜色
    nrDEG_limma_voom, aes(x = logFC, y = -log10(adj.P.Val), colour=change)) +
    geom_point(alpha=0.5, size=3.5) +
    scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+
    #突出表示差异基因
    geom_point(data=geneList,aes(x = logFC, y = -log10(adj.P.Val)),colour="yellow",size=3.5)+
    #辅助线
    geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.8) +
    labs(x="log2(fold change)",y="-log10 (p-value)")+   # 坐标轴# 坐标轴和图标题title="Volcano plot",
    theme_bw()+    #去除背景色
    theme(panel.grid = element_blank())+  #去除网格线
    #xlim(-2, 2)+   #设置坐标轴范围
    #图例
    theme(plot.title = element_text(hjust = 0.5,size=24), 
          legend.position="bottom", 
          legend.title = element_blank(),
          legend.text=element_text(size=18),
          legend.key.size = unit(1, 'cm'),
          legend.background = element_rect(fill="gray90", linetype="solid",colour ="gray"),
          axis.title.x =element_text(size=18), 
          axis.title.y=element_text(size=18),
          axis.text=element_text(size=14,face = "bold"))
  #标记出选定基因的label
  geneList1 <- nrDEG_limma_voom[rownames(nrDEG_limma_voom) %in% i,]
  geneList1 <- subset(geneList1, select = -change)
  geneList1$label <- rownames(geneList1)
  #pdf("hotplot.pdf", width = 8, height = 10)
  volplot <- p + geom_label_repel(data = geneList1, 
                                  aes(x = logFC, y = -log10(adj.P.Val), label = label),
                                  size = 4,color="black",
                                  box.padding = unit(0.4, "lines"), 
                                  segment.color = "black",   #连线的颜色
                                  segment.size = 0.4,  #连线的粗细
  )
  png(paste("volcano_limma_",i,".png"), width = 8, height = 10, units = "in", res = 800)
  print(volplot)
  dev.off()
}

## VENN图
#加载
library(VennDiagram)
#设置工作目录
setwd("E:/CRC_model")
#读取
DEG <- c("Wilcoxon","limma","edgeR","DESeq2")
for (i in DEG) {
  data=read.table(paste("TCGA.diff.",i,".txt",sep=""), header=F, sep="\t", check.names=F)
  DEG = data[,1]
  data=read.table("ARG_gene.txt", header=F, sep="\t", check.names=F)
  ARG = data[,1]
  
  #画图
  venn.diagram(x = list('ARG' = ARG,i = DEG),
               filename = paste("VN_",i,".png",sep = ""),fill = c("darkorange1","green"))
  #交集基因
  intersectGenes = intersect(DEG,ARG)
  
  write.table(file= paste("intersectGenes_",i,".txt",sep = "") , intersectGenes, sep="\t", quote=F, col.names=F, row.names=F)
}

## Lasso回归 变量初步筛选，处理多重共线性问题
#加载包
library("glmnet")
library("survival")
# 设置工作目录
base_dir <- "E:/CRC_model"
setwd(base_dir)
#读入生存数据
cli=read.table("time_CRC.txt", header=T, sep="\t", check.names=F, row.names=1)
#生存时间以年为单位
cli$time=cli$time/365
#读入
data=read.table("TCGA_CRC_count.txt", header=T, sep="\t", check.names=F,row.names = 1)
# 读取并合并交集基因
#DEG_methods <- c("Wilcoxon", "limma", "edgeR", "DESeq2")
# 读取只有一列且有表头的文件
DEG_data <- read.table("intersectGenes_DESeq2.txt", header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#仅保留肿瘤样本
data = data[,group == 0]
# 从表达数据中提取交集基因
data = data[DEG_data$ID,]
#转置
data=t(data)
#样本名仅保留前12字符
rownames(data)=substr(rownames(data),1,12)
#将.改为-
rownames(data) = gsub('[.]', '-', rownames(data))
#获取共同样本
sameSample=intersect(row.names(data), row.names(cli))
#提取共同样本
data=data[sameSample,]
cli=cli[sameSample,]
#合并
rt=cbind(cli, data)
#设置随机种子，K折交叉验证，每次的数据都是随机的，随机数种子一致，就结果给固定住。
set.seed(123456)
#构建lasso回归模型
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$time, rt$status))
fit=glmnet(x, y, family="cox", nfolds = 10)
#c-index
cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C", nfolds = 10)
#pdf("lasso.c-index.pdf")
png(file="lasso.c-index_limma.png", width=8, height=8, units="in", res=800)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()
#deviance图形
cvfit=cv.glmnet(x, y, family="cox", type.measure = "deviance", nfolds = 10)
#pdf("lasso.cvfit.pdf")
png(file="lasso.cvfit_limma.png", width=8, height=8, units="in", res=800)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()
#pdf("lasso.lambda.pdf")
png(file="lasso.lambda_limma.png", width=8, height=8, units="in", res=800)
plot(fit, xvar="lambda", label=T)
abline(v=log(cvfit$lambda.min), lty="dashed")
dev.off()
#输出lasso显著基因表达量
coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
row.names(coef)[index]
actCoef
index
lassoSigExp=rt[,c("time", "status", lassoGene)]
lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
write.table(lassoSigExpOut,file="lasso.SigExp_limma.txt",sep="\t",row.names=F,quote=F)

## Boruta 算法进一步筛选变量
# Joe_2024_10_12
rm(list=ls())

# 加载包
library(tidyverse)
library(survival)
library(Boruta)

# 设置工作目录
base_dir <- "E:/CRC_model"
setwd(base_dir)

# 读取数据
sadata <- read.table("lasso.SigExp_limma.txt", header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
#sadata <- readr::read_csv("F:/Mach_learn_data/sadata.csv")
colnames(sadata)

# 修正变量类型
# 分类变量转factor
#for(i in c(3)){
#  sadata[[i]] <- factor(sadata[[i]])
#}

# 剔除变量-无关变量
#sadata$pid <- NULL
# 剔除样本-含有缺失值的样本，时间非正的样本
# sadata <- na.omit(sadata) # 剔除任意变量下有缺失值的样本
sadata <- sadata %>%
  # na.omit() %>%
  # drop_na(age) %>%  # 剔除指定变量下有缺失值的样本
  filter(time > 0) # 剔除时间非正的样本

#sadata删除id那一列
sadata <- sadata[,-1]

# 数据概况
skimr::skim(sadata)

# 感兴趣的时间节点
range(unique(sadata$time))
itps <- c(1 * c(1, 3, 5))
itps
table(cut(sadata$time, c(0, itps, Inf)))


###################################################

# 数据拆分构建任务对象
set.seed(42)
datasplit <- rsample::initial_split(
  sadata, prop = 0.6, strata = time, breaks = 10
)
traindata <- rsample::training(datasplit) %>%
  sample_n(nrow(.))
testdata <- rsample::testing(datasplit) %>%
  sample_n(nrow(.))

# 拆分得到的数据的生存曲线比较
sadata2 <- sadata
sadata2$set <- "Test"
sadata2$set[datasplit$in_id] <- "Train"
sadata2$set <- factor(sadata2$set)
sfit_set <- survfit(Surv(time, status) ~ set, data=sadata2)
survminer::ggsurvplot(
  sfit_set,
  pval=TRUE, 
  pval.coord = c(0.1, 0.5),
  risk.table=TRUE,
  ggtheme = survminer::theme_survminer() +
    theme(text = element_text(family = "serif")),
  font.family = "serif"
)

###################################################

# 数据预处理
datarecipe_coxph <- 
  recipes::recipe(time + status ~ ., traindata) %>%
  recipes::prep()
datarecipe_coxph

# 按方处理训练集和测试集
traindata2 <- 
  recipes::bake(datarecipe_coxph, new_data = NULL) %>%
  dplyr::select(time, status, everything())
testdata2 <- 
  recipes::bake(datarecipe_coxph, new_data = testdata) %>%
  dplyr::select(time, status, everything())

# 特征筛选
set.seed(42)
result_boruta <- Boruta(
  x = traindata2[,-c(1,2)],
  y = Surv(traindata2$time, traindata2$status), 
  doTrace = 1,
  maxRuns = 100
)
# 筛选结果概况
result_boruta

# 对于存疑变量进一步确认
result_boruta <- TentativeRoughFix(result_boruta)
# 具体结果
attStats(result_boruta)
# 图示
par(las = 3, mar = c(8, 4, 4, 2))
plot(result_boruta, xlab = "", main = "Boruta算法筛选结果")
legend("topleft", 
       legend = c("confirmed", "rejected"),
       pch = 22,
       pt.bg = c("green", "red"))

# 筛选得到的变量
getSelectedAttributes(result_boruta)

# 筛选之后的数据
traindata <- traindata2 %>%
  dplyr::select(time, status, 
                getSelectedAttributes(result_boruta))
testdata <- testdata2 %>%
  dplyr::select(time, status, 
                getSelectedAttributes(result_boruta))

##模型1 coxph
library(tidyverse)
library(survival)
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
source("H:/Biolab/Biolab/tidyfuncs4sa.R")
# file.choose()

# 性能指标
measure_sa <- msrs("surv.cindex")
measure_sa

# 训练任务对象
task_train <- as_task_surv(
  traindata, 
  time = "time",
  event = "status", 
  type = "right"
)
task_train
# 测试任务对象
task_test <- as_task_surv(
  testdata, 
  time = "time",
  event = "status", 
  type = "right"
)
task_test

###################################################

# coxph模型
# https://mlr3proba.mlr-org.com/reference/mlr_learners_surv.coxph.html

# 模型设定
learner_coxph <- lrn("surv.coxph")
learner_coxph$id <- "coxph"
learner_coxph

# 模型训练
set.seed(42)
learner_coxph$train(task_train)
learner_coxph

# 模型概况
learner_coxph$model
summary(learner_coxph$model)
# 森林图
survminer::ggforest(
  learner_coxph$model, 
  data = task_train$data()
)


###################################################

# 预测训练集
predtrain_coxph <- learner_coxph$predict(task_train)
predprobtrain_coxph <- predprob(
  pred = predtrain_coxph, 
  preddata = traindata, 
  etime = "time",
  estatus = "status",
  model = "coxph", 
  dataset = "train", 
  timepoints =itps
)

# 性能指标
predtrain_coxph$score(measure_sa)
cindex_bootci(learner_coxph, traindata)

evaltrain_coxph <- eval4sa(
  predprob = predprobtrain_coxph,
  preddata = traindata,
  etime = "time",
  estatus = "status",
  model = "coxph",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltrain_coxph$auc
evaltrain_coxph$rocplot
evaltrain_coxph$brierscore
evaltrain_coxph$brierscoretest
evaltrain_coxph$calibrationplot
evaltrain_coxph$riskplot

sadca(
  predprob = predprobtrain_coxph,
  preddata = traindata,
  etime = "time",
  estatus = "status",
  model = "coxph",
  dataset = "train",
  timepoints = itps,
  timepoint = 1, 
  xrange = 0:100 / 100
)

# 预测测试集
predtest_coxph <- learner_coxph$predict(task_test)
predprobtest_coxph <- predprob(
  pred = predtest_coxph, 
  preddata = testdata, 
  etime = "time",
  estatus = "status",
  model = "coxph", 
  dataset = "test", 
  timepoints =itps
)
# 性能指标
predtest_coxph$score(measure_sa)
cindex_bootci(learner_coxph, testdata)

evaltest_coxph <- eval4sa(
  predprob = predprobtest_coxph,
  preddata = testdata,
  etime = "time",
  estatus = "status",
  model = "coxph",
  dataset = "test",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltest_coxph$auc
evaltest_coxph$rocplot
evaltest_coxph$brierscore
evaltest_coxph$brierscoretest
evaltest_coxph$calibrationplot
evaltest_coxph$riskplot

sadca(
  predprob = predprobtest_coxph,
  preddata = testdata,
  etime = "time",
  estatus = "status",
  model = "coxph",
  dataset = "test",
  timepoints = itps,
  timepoint = 1, 
  xrange = 0:100 / 100
)

# 保存结果用于比较
save(predprobtrain_coxph,
     evaltrain_coxph,
     predprobtest_coxph,
     evaltest_coxph,
     file = "F:/Mach_learn_data/mlr_model/coxph.RData")

# 保存结果用于shiny网页计算器
traindata4sadata <- traindata
learner_coxph4sadata <- learner_coxph
save(traindata4sadata,
     learner_coxph4sadata,
     file = "F:/Mach_learn_data/shiny/coxph.RData")

#############################################

#############################################

# 模型解释
# 原始训练样本中自变量和因变量
traindatax <- traindata[, task_train$feature_names]
catvars <- 
  colnames(traindatax)[sapply(traindatax, is.factor)]
convars <- setdiff(colnames(traindatax), catvars)
traindatay <- survival::Surv(
  time = traindata$time, 
  event = traindata$status
)

# 解释器——基于训练集，可以不指定时间点
exper_coxph <- survex::explain(
  learner_coxph$model, 
  data = traindatax,
  y = traindatay,
  times = itps
)

# 变量重要性
viplot(exper_coxph, output_type = "risk")
viplot(exper_coxph, output_type = "survival")

# 偏依赖图
pdplot(exper_coxph, convars, output_type = "survival")
pdplot(exper_coxph, convars, output_type = "chf")
pdplot(exper_coxph, convars, output_type = "risk")


# 单样本预测分解
shap4one(exper_coxph, traindatax[1,], output_type = "survival")
shap4one(exper_coxph, traindatax[1,], output_type = "risk")

# 全局shap
shap4all_coxph <- survex::model_survshap(
  exper_coxph,
  new_observation = traindatax,
  y_true = traindatay,
  N = 100,
  calculation_method = "kernelshap",
  aggregation_method = "integral",
  output_type = "risk"
)
plot(shap4all_coxph)
plot(shap4all_coxph, geom = "beeswarm")
plot(shap4all_coxph, geom = "profile", variable = "GAL")
plot(shap4all_coxph, geom = "profile", variable = "grade")

# (抽样)汇总计算shap
sumshap_coxph <- sumshap(
  exper_coxph, 
  traindatax, 
  catvars, 
  convars, 
  sampleN=200
)
sumshap_coxph$shapvipplot
sumshap_coxph$shapplotd
sumshap_coxph$shapplotc

##模型2 ctree
# 加载包
# library(devtools)
# devtools::install_github("mlr-org/mlr3proba")
# devtools::install_github("mlr-org/mlr3extralearners@*release")


##模型3 xgboost
library(tidyverse)
library(survival)
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
source("H:/Biolab/Biolab/tidyfuncs4sa.R")
# file.choose()
# 性能指标
measure_sa <- msrs("surv.cindex")
measure_sa

# 训练任务对象
task_train <- as_task_surv(
  traindata, 
  time = "time",
  event = "status", 
  type = "right"
)
task_train
# 测试任务对象
task_test <- as_task_surv(
  testdata, 
  time = "time",
  event = "status", 
  type = "right"
)
task_test

###################################################

# xgboost模型
# https://mlr3extralearners.mlr-org.com/reference/mlr_learners_surv.xgboost.cox.html

# 设定模型
learner_xgboost <- as_learner(
  po("encode", method = "treatment") %>>%
    lrn(
      "surv.xgboost.cox",
      nrounds = to_tune(100, 500), 
      max_depth = to_tune(1, 5),
      eta = to_tune(1e-4, 1)
    )
)
learner_xgboost$id <- "xgboost"
learner_xgboost

# 多核并行
future::plan("multisession")

# 超参数调优设定
set.seed(42)
tune_xgboost <- tune(
  tuner = tnr("random_search"),
  task = task_train, 
  learner = learner_xgboost,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("evals", n_evals = 20)
)
tune_xgboost
tune_xgboost$archive$data %>%
  as.data.frame() %>%
  select(1:4) %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~surv.cindex, 
                colorscale = 'Jet', 
                showscale = T),
    dimensions = list(
      list(label = 'nrounds', values = ~surv.xgboost.cox.nrounds),
      list(label = 'max_depth', values = ~surv.xgboost.cox.max_depth),
      list(label = 'eta', values = ~surv.xgboost.cox.eta)
    )
  ) %>%
  plotly::layout(title = "xgboost HPO Guided by C-Index",
                 font = list(family = "serif"))

# 训练
learner_xgboost$param_set$values <-
  tune_xgboost$result_learner_param_vals
set.seed(42)
learner_xgboost$train(task_train)

# 模型概况
learner_xgboost

###################################################

# 预测训练集
predtrain_xgboost <- learner_xgboost$predict(task_train)
predprobtrain_xgboost <- predprob(
  pred = predtrain_xgboost, 
  preddata = traindata, 
  etime = "time",
  estatus = "status",
  model = "xgboost", 
  dataset = "train", 
  timepoints =itps
)

# 性能指标
predtrain_xgboost$score(measure_sa)
cindex_bootci(learner_xgboost, traindata)

evaltrain_xgboost <- eval4sa(
  predprob = predprobtrain_xgboost,
  preddata = traindata,
  etime = "time",
  estatus = "status",
  model = "xgboost",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltrain_xgboost$auc
evaltrain_xgboost$rocplot
evaltrain_xgboost$brierscore
evaltrain_xgboost$brierscoretest
evaltrain_xgboost$calibrationplot
evaltrain_xgboost$riskplot

sadca(
  predprob = predprobtrain_xgboost,
  preddata = traindata,
  etime = "time",
  estatus = "status",
  model = "xgboost",
  dataset = "train",
  timepoints = itps,
  timepoint = 1, 
  xrange = 0:100 / 100
)

# 预测测试集
predtest_xgboost <- learner_xgboost$predict(task_test)
predprobtest_xgboost <- predprob(
  pred = predtest_xgboost, 
  preddata = testdata, 
  etime = "time",
  estatus = "status",
  model = "xgboost", 
  dataset = "test", 
  timepoints =itps
)
# 性能指标
predtest_xgboost$score(measure_sa)
cindex_bootci(learner_xgboost, testdata)

evaltest_xgboost <- eval4sa(
  predprob = predprobtest_xgboost,
  preddata = testdata,
  etime = "time",
  estatus = "status",
  model = "xgboost",
  dataset = "test",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltest_xgboost$auc
evaltest_xgboost$rocplot
evaltest_xgboost$brierscore
evaltest_xgboost$brierscoretest
evaltest_xgboost$calibrationplot
evaltest_xgboost$riskplot

sadca(
  predprob = predprobtest_xgboost,
  preddata = testdata,
  etime = "time",
  estatus = "status",
  model = "xgboost",
  dataset = "test",
  timepoints = itps,
  timepoint = 1, 
  xrange = 0:100 / 100
)

# 保存结果用于比较
save(predprobtrain_xgboost,
     evaltrain_xgboost,
     predprobtest_xgboost,
     evaltest_xgboost,
     file = "F:/Mach_learn_data/mlr_model/xgboost.RData")

# 保存结果用于shiny网页计算器
traindata4sadata <- traindata
learner_xgboost4sadata <- learner_xgboost
save(traindata4sadata,
     learner_xgboost4sadata,
     file = "F:/Mach_learn_data/shiny/xgboost.RData")

#############################################


#############################################

# 模型解释
# 原始训练样本中自变量和因变量
traindatax <- traindata[, task_train$feature_names]
catvars <- 
  colnames(traindatax)[sapply(traindatax, is.factor)]
convars <- setdiff(colnames(traindatax), catvars)
traindatay <- survival::Surv(
  time = traindata$time, 
  event = traindata$status
)

# 解释器——基于训练集，可以不指定时间点
exper_xgboost <- survex::explain_survival(
  learner_xgboost, 
  data = traindatax,
  y = traindatay,
  predict_function = risk_pred,
  predict_survival_function = surv_pred,
  predict_cumulative_hazard_function = chf_pred,
  label = "xgboost",
  times = itps
)

# 变量重要性
viplot(exper_xgboost, output_type = "risk")
viplot(exper_xgboost, output_type = "survival")

# 偏依赖图
pdplot(exper_xgboost, convars, output_type = "survival")
pdplot(exper_xgboost, convars, output_type = "chf")
pdplot(exper_xgboost, convars, output_type = "risk")


# 单样本预测分解
shap4one(exper_xgboost, traindatax[1,], output_type = "survival")
shap4one(exper_xgboost, traindatax[1,], output_type = "risk")

# 全局shap
shap4all_xgboost <- survex::model_survshap(
  exper_xgboost,
  new_observation = traindatax,
  y_true = traindatay,
  N = 100,
  calculation_method = "kernelshap",
  aggregation_method = "integral",
  output_type = "risk"
)
plot(shap4all_xgboost)
plot(shap4all_xgboost, geom = "beeswarm")
plot(shap4all_xgboost, geom = "profile", variable = "nodes")
plot(shap4all_xgboost, geom = "profile", variable = "grade")

# (抽样)汇总计算shap
sumshap_xgboost <- sumshap(
  exper_xgboost, 
  traindatax, 
  catvars, 
  convars, 
  sampleN=200
)
sumshap_xgboost$shapvipplot
sumshap_xgboost$shapplotd
sumshap_xgboost$shapplotc
