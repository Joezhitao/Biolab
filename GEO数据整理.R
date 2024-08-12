#加载
library(GEOquery)

#设置工作目录
setwd("E:/CORD/GEO/GSE13067")

#下载矩阵，目录下已存在，自动读取
gset <- getGEO("GSE13067",destdir = "E:/CORD/GEO/GSE13067",AnnotGPL = F,getGPL = F) 

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

#输出临床信息
write.csv(pd,'clinical_GSE13067.csv',row.names = TRUE)

#读入txt注释文件
gpl=read.table("GPL570-55999.txt",
               header = TRUE,fill = T,sep = "\t",
               comment.char = "#",
               stringsAsFactors = FALSE,
               quote = "")
gpl$Gene.Symbol
#查看
#View(gpl)
#colnames(gpl)


#提取探针ID及基因symbol
ids=gpl[,c("ID","Gene.Symbol")]

#修改列名
colnames(ids)=c('probe_id','symbol')

#获取基因symbol
library(stringr)  
#[]中的数字，以及///的数量，根据表格中gene symbol的位置决定
ids$symbol=trimws(str_split(ids$symbol,'///',simplify = T)[,1])

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
write.table(data.frame(ID=rownames(dat),dat),file="GSE13067.txt", sep="\t", quote=F, row.names = F)

