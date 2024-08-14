#加载R包
library(rjson)
library(tidyverse)

#设置路径
setwd("E:/CORD/TCGAdata")

#读入meta.data文件
#出现报错，注意上面的包是否加载，工作目录下是否存在这个文件，文件是否损坏（先跑示例数据）
json <- jsonlite::fromJSON("metadata.cart.2024-08-11.json")
#View(json)

#获取样本名称及文件名称
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
#sample_id[1:10]
file_sample <- data.frame(sample_id,file_name=json$file_name)  
#View(file_sample)

#获取counts所在位置
count_file <- list.files('gdc_download_20240811_161933.816528/',
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
  path = paste0('gdc_download_20240811_161933.816528//',count_file[i])   #Counts文件夹名
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
path = paste0('gdc_download_20240811_161933.816528//',count_file[1])
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
write.table(matrix1,'TCGA_CORD_TPM.txt', sep="\t", quote=F, row.names = F)
#write.table(matrix1,'TCGA_LUSC_count.txt', sep="\t", quote=F, row.names = F)

################################################################################
#提取count数据
#设置路径
setwd("E:/CORD/TCGAdata")

#读入meta.data文件
#出现报错，注意上面的包是否加载，工作目录下是否存在这个文件，文件是否损坏（先跑示例数据）
json <- jsonlite::fromJSON("metadata.cart.2024-08-11.json")
#View(json)

#获取样本名称及文件名称
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
#sample_id[1:10]
file_sample <- data.frame(sample_id,file_name=json$file_name)  
#View(file_sample)

#获取counts所在位置
count_file <- list.files('gdc_download_20240811_161933.816528/',
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
  path = paste0('gdc_download_20240811_161933.816528//',count_file[i])   #Counts文件夹名
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
path = paste0('gdc_download_20240811_161933.816528//',count_file[1])
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
write.table(matrix1,'TCGA_CORD_count.txt', sep="\t", quote=F, row.names = F)