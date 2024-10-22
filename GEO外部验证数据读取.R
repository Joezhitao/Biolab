
#加载
library(GEOquery)
rm(list = ls())
#设置工作目录
setwd("E:/CRC_model/GEO")

#下载矩阵，目录下已存在，自动读取
gset <- getGEO("GSE71187",destdir = "E:/CRC_model/GEO",AnnotGPL = F,getGPL = F) 

#获取表达矩阵
dat=exprs(gset[[1]])

#如果数据已经经log后，显示log2 transform not needed；
#如果数据没有尽行log，需要log程序按照代码尽行log转换，完成后显示log2 transform finished。
ex <- dat
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))

LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NA
dat <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

#获取临床信息
pd=pData(gset[[1]])

# 将 'survival (year):ch1' 列重命名为 'time'
names(pd)[names(pd) == 'survival (year):ch1'] <- 'time'
# 将 'survival (year):ch1' 列重命名为 'time'
names(pd)[names(pd) == 'events (death = 1, alive = 0):ch1'] <- 'status'

# 将 "NA" 字符串转换为真正的 NA 值
pd$time[pd$time == "NA"] <- NA

# 删除 'time' 列中包含 NA 值的行
pd <- pd[!is.na(pd$time), ]

# 查看处理后的结果
head(pd$time)
nrow(pd)
colnames(pd)
# 只保留 'time' 和 'status' 列
pd <- pd[, c('time', 'status')]

#输出临床信息
write.csv(pd,'clinical_GSE74777.csv',row.names = TRUE)

#读入txt注释文件
gpl=read.table("GPL6480-9577.txt",
               header = TRUE,fill = T,sep = "\t",
               comment.char = "#",
               stringsAsFactors = FALSE,
               quote = "")

#查看
#View(gpl)
#colnames(gpl)

#提取探针ID及基因symbol
ids=gpl[,c("ID","GENE_SYMBOL")]

#修改列名
colnames(ids)=c('probe_id','symbol')

#获取基因symbol
library(stringr)  
#针对Gene Symbol列，使用str_split函数，以//为分隔符，将Gene Symbol列分割为两列，取第二列
ids$symbol=trimws(str_split(ids$symbol,'//',simplify = T)[,2])

#去掉没有注释symbol的探针
ids=ids[ids$symbol != '',]
ids=ids[ids$symbol != '---',]
ids=ids[ids$symbol != '--- ',]

#由于部分探针ID是由全数字组成的，会导致R自动将其识别为数值型
#因此，这里需要将探针的id改为字符型
rownames(dat) = as.character(rownames(dat))
ids$probe_id = as.character(ids$probe_id)

#平台文件的ID和矩阵中的ID匹配。%in%用于判断是否匹配。
ids=ids[ids$probe_id %in% rownames(dat),]

#获取匹配的表达数据
dat=dat[ids$probe_id,]

#查看探针名与注释文件名是否完全一致
table(rownames(dat) == ids$probe_id)

#合并
dat <- cbind(ids,dat)

#删除重复基因，保留最大值
dat <- aggregate( . ~ symbol,data=dat, max)

#一定注意看，是否需要删除几行
View(dat)
#dat = dat[-1,]

#转化行名
rownames(dat) <- dat[,1]

#删除第一二列
dat <- dat[,-c(1,2)]

#导出
write.table(data.frame(ID=rownames(dat),dat),file="GSE71187.txt", sep="\t", quote=F, row.names = F)

# 转置表达矩阵
dat_t <- t(dat)

# 确保行名（样本名）是字符型
rownames(dat_t) <- as.character(rownames(dat_t))

# 确保pd的行名（样本名）是字符型
rownames(pd) <- as.character(rownames(pd))

# 找出在临床数据中存在的样本
common_samples <- intersect(rownames(pd), rownames(dat_t))

# 只保留这些共同的样本
pd_filtered <- pd[common_samples, ]
dat_t_filtered <- dat_t[common_samples, ]

# 合并数据
merged_data <- cbind(pd_filtered, dat_t_filtered)

# 查看结果
dim(merged_data)
head(merged_data[, 1:10])  # 显示前10列，包括临床数据和一些基因表达数据

# 保存合并后的数据
write.table(merged_data, file="merged_expression_clinical.txt", sep="\t", quote=FALSE, row.names=TRUE)


#保留差异基因列
cgene <- getSelectedAttributes(result_boruta)
# 定义要保留的列，使用 cgene 和临床数据列
keep_columns <- c("time", "status", cgene)

# 检查所有需要的列是否存在
missing_columns <- setdiff(keep_columns, colnames(merged_data))
if (length(missing_columns) > 0) {
  warning(paste("以下列不存在于数据中:", paste(missing_columns, collapse = ", ")))
  keep_columns <- intersect(keep_columns, colnames(merged_data))
}

# 只保留指定的列
final_data <- merged_data[, keep_columns]

# 查看结果
dim(final_data)
head(final_data)
