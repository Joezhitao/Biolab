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

