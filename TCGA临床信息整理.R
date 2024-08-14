#加载包
library("XML")
library("methods")
library("tinyarray")

#设置工作路径
setwd("E:/CORD/clinic")

#设置一个目录
dir="E:/CORD/clinic/gdc_download_20240814_142022.604604"      

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
write.table(clinical,file="clinical.txt",sep="\t",quote=F,row.names = F)  
