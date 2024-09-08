
#加载R包
library(maftools)
library(dplyr)
library(mclust)
library(NMF)
library(pheatmap)
library(barplot3d)
library(openxlsx)
library(TRONCO)
library(Matrix)

#目录
setwd("E:/CRC/genetic mutation/")
#合并所有数据
unnecessary_files <- list.files(pattern = '*.gz.parcel',recursive = TRUE)
files <- list.files(pattern = '*.gz',recursive = TRUE)
#去掉files中包含unnecessary_files的文件
files <- files[!files %in% unnecessary_files]

all_mut <- data.frame()
for (file in files) {
  mut <- read.delim(file,skip = 7, header = T, fill = TRUE,sep = "\t")
  all_mut <- rbind(all_mut,mut)
}
#仅保留前12个字符
all_mut$Tumor_Sample_Barcode = substr(all_mut$Tumor_Sample_Barcode,1,12)

#数据读入
all_mut <- read.maf(all_mut)

#数据整理
a <- all_mut@data %>%
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>%
  as.data.frame() %>%
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))

#提取基因
gene <- as.character(unique(a$Hugo_Symbol))
#提取样本
sample <- as.character(unique(a$Tumor_Sample_Barcode))

#创建data.frame
mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))
#将信息填入
for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}

#创建data.frame
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))
#将信息填入
for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

#所有样本突变情况汇总/排序
gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))

colnames(gene_count)[1] = "Gene"
colnames(gene_count)[2] = "Num"
write.table(gene_count,'geneMut.txt', sep="\t", quote=F, row.names = F)
write.mafSummary(maf = all_mut,basename = "input")

#绘制瀑布图oncoplot
#missense_mutation:错义突变 frame_shift_del：移码缺失突变 
#nonsense_mutation：无义突变 frame_shift_ins：移码插入突变 
#splice_site：剪接位点 in_frame_ins：inframe插入 
#in_frame_del：inframe缺失 translation_start_site:转录起始位点 
#nonstop_mutation：终止密码子突变 
png(file="plotmafSummary.png", width=15, height=15, units = "in", res = 800)
plotmafSummary(maf = all_mut, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

png(file="oncoplot.png", width=15, height=25, units = "in", res = 800)
oncoplot(maf = all_mut,
         top = 100, #显示前30个的突变基因信息
         fontSize = 0.6, #设置字体大小
         showTumorSampleBarcodes = F) #不显示病人信息
dev.off()
#通路上的突变
png(file="oncoplot_pathway.png", width=8, height=8, units = "in", res = 800)
OncogenicPathways(maf = all_mut)
dev.off()
#绘制特定通路上的突变基因
png(file="PlotOncogenicPathways.png", width=8, height=8, units = "in", res = 800)
PlotOncogenicPathways(maf = all_mut,pathways="RTK-RAS") 
dev.off()

#计算tmb值
tmb_table = tmb(maf = all_mut,logScale = F)
#保留需要tmb值信息
tmb_table = tmb_table[,c(1,3)]
tmb_table = as.data.frame(tmb_table)
tmb_table[,1]=substr(tmb_table[,1],1,12)
colnames(tmb_table)
tmb_table <- aggregate( . ~ Tumor_Sample_Barcode,data=tmb_table, max)
colnames(tmb_table)[1] = "id"
colnames(tmb_table)[2] = "TMB"
write.table(tmb_table,'TMB.txt', sep="\t", quote=F, row.names = F)

#过渡和颠倒
laml.titv = titv(maf = all_mut, plot = FALSE, useSyn = TRUE)
#plot titv summary
png(file="plotTiTv.png", width=8, height=8, units = "in", res = 800)
plotTiTv(res = laml.titv)
dev.off()

#绘制lollipop图
png(file="lollipopPlot.png", width=25, height=6, units = "in", res = 800)
lollipopPlot(
  maf = all_mut,
  gene = 'TTN',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  refSeqID = 'NM_003319'  # 替换为你需要的转录本
)
dev.off()

#蛋白质结构域
png(file="plotProtein.png", width=8, height=8, units = "in", res = 800)
plotProtein(gene = "TP53", refSeqID = "NM_000546")
dev.off()

#降雨图
png(file="rainfallPlot.png", width=8, height=8, units = "in", res = 800)
rainfallPlot(maf = all_mut, detectChangePoints = TRUE, pointSize = 0.4)
dev.off()

#将突变负荷与 TCGA 队列进行比较
png(file="tcgaCompare.png", width=18, height=18, units = "in", res = 800)
laml.mutload = tcgaCompare(maf = all_mut, cohortName = 'Example-LAML', logscale = TRUE, capture_size = 50)
dev.off()

#躯体相互作用
png(file="plotOncogenicSomaticInteractions.png", width=8, height=8, units = "in", res = 800)
somaticInteractions(maf = all_mut, top = 25, pvalue = c(0.05, 0.1))
dev.off()

#基于位置聚类检测癌症驱动基因
laml.sig = oncodrive(maf = all_mut, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
head(laml.sig)
#绘制结果
png(file="plotOncodrive.png", width=8, height=8, units = "in", res = 800)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
dev.off()

#药物-基因相互作用
png(file="drugInteractions.png", width=8, height=8, units = "in", res = 800)
dgi = drugInteractions(maf = all_mut, fontSize = 0.75)
dev.off()
#变异签名
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
laml.tnm = trinucleotideMatrix(maf = all_mut, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
#APOBEC富集和非富集样品之间的差异
plotApobecDiff(tnm = laml.tnm, maf = laml, pVal = 0.2)
