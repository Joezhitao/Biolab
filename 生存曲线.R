#加载
library(survival)
library(survminer)

#工作路径
setwd("E:/CORD/SURV")

#读入生存数据
cli=read.table("time_CORD.txt", header=T, sep="\t", check.names=F, row.names=1)
#生存时间以年为单位
cli$time=cli$time/365

#读取输入文件
data=read.table("TCGA_CORD_count.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#仅保留肿瘤样本
data = data[,group == 0]

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

#根据TMEM196表达量的中位值，把样品分为两组
group = ifelse(rt[,"TMEM196"]>quantile(rt[,"TMEM196"], seq(0,1,1/2))[2],"High","Low")
#group=ifelse(rt[,"ABCA3"]>quantile(rt[,"ABCA3"], seq(0,1,1/3))[3],"High","Low")
rt[,"group"] = group
length = length(levels(factor(group)))

#生存分析
diff = survdiff(Surv(time, state) ~group,data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
pValue=paste0("p=",sprintf("%.04f",pValue))
fit <- survfit(Surv(time, state) ~group,data = rt)

#绘制生存曲线
bioCol=c("Firebrick3","MediumSeaGreen","#6E568C","#223D6C")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="TMEM196 expression",
                   legend.labs=levels(factor(rt[,"group"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 2,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)

#输出图形
pdf(file="survival.pdf", width=6.5, height=6.25, onefile=FALSE)
print(surPlot)
dev.off()

#for循环批量输出基因数据
#读取cox差异分析基因文件
genetable=read.table("uniCox.txt", header=T, sep="\t", check.names=F,row.names = 1)
gene <- rownames(genetable)
#去除gene中重复基因
gene=unique(gene)

for (i in gene) {
  #设置保存路径
  filepath <- paste0("E:/CORD/SURV/genesurv/")
  
  #读入生存数据
  cli=read.table("time_CORD.txt", header=T, sep="\t", check.names=F, row.names=1)
  #生存时间以年为单位
  cli$time=cli$time/365
  
  #读取输入文件
  data=read.table("TCGA_CORD_count.txt", header=T, sep="\t", check.names=F,row.names = 1)
  #转化为matrix
  dimnames=list(rownames(data), colnames(data))
  data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
  
  #正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
  group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
  group=sapply(strsplit(group,""), "[", 1)
  #仅保留肿瘤样本
  data = data[,group == 0]
  
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
  
  #根据i基因表达量的中位值，把样品分为两组
  group = ifelse(rt[,i]>quantile(rt[,i], seq(0,1,1/2))[2],"High","Low")
  #group=ifelse(rt[,"ABCA3"]>quantile(rt[,"ABCA3"], seq(0,1,1/3))[3],"High","Low")
  rt[,"group"] = group
  length = length(levels(factor(group)))
  
  #生存分析
  diff = survdiff(Surv(time, state) ~group,data = rt)
  pValue=1-pchisq(diff$chisq, df=length-1)
  pValue=paste0("p=",sprintf("%.04f",pValue))
  fit <- survfit(Surv(time, state) ~group,data = rt)
  
  #绘制生存曲线
  bioCol=c("Firebrick3","MediumSeaGreen","#6E568C","#223D6C")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title=paste0(i," expression"),
                     legend.labs=levels(factor(rt[,"group"])),
                     legend = c(0.8, 0.8),
                     font.legend=10,
                     xlab="Time(years)",
                     break.time.by = 2,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table=T,
                     cumevents=F,
                     risk.table.height=.25)
  
  #输出图形
  pdf(file=paste(filepath,i,"_survial.pdf",sep = ""), width=6.5, height=6.25, onefile=FALSE)
  print(surPlot)
  dev.off()
  
}


#for循环只输出p值小于0.05的生存曲线
for (i in gene) {
  #设置保存路径
  filepath <- paste0("E:/CORD/SURV/genesurv_pvalue/")
  
  #读入生存数据
  cli=read.table("time_CORD.txt", header=T, sep="\t", check.names=F, row.names=1)
  #生存时间以年为单位
  cli$time=cli$time/365
  
  #读取输入文件
  data=read.table("TCGA_CORD_count.txt", header=T, sep="\t", check.names=F,row.names = 1)
  #转化为matrix
  dimnames=list(rownames(data), colnames(data))
  data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
  
  #正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
  group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
  group=sapply(strsplit(group,""), "[", 1)
  #仅保留肿瘤样本
  data = data[,group == 0]
  
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
  
  #根据i基因表达量的中位值，把样品分为两组
  group = ifelse(rt[,i]>quantile(rt[,i], seq(0,1,1/2))[2],"High","Low")
  rt[,"group"] = group
  length = length(levels(factor(group)))
  
  #生存分析
  diff = survdiff(Surv(time, state) ~group,data = rt)
  pValue=1-pchisq(diff$chisq, df=length-1)
  pValueFormatted=paste0("p=",sprintf("%.04f",pValue))
  fit <- survfit(Surv(time, state) ~group,data = rt)
  
  # 仅在pValue小于0.05时绘制和保存生存曲线
  if (pValue < 0.05) {
    #绘制生存曲线
    bioCol=c("Firebrick3","MediumSeaGreen","#6E568C","#223D6C")
    bioCol=bioCol[1:length]
    surPlot=ggsurvplot(fit, 
                       data=rt,
                       conf.int=F,
                       pval=pValueFormatted,
                       pval.size=6,
                       legend.title=paste0(i," expression"),
                       legend.labs=levels(factor(rt[,"group"])),
                       legend = c(0.8, 0.8),
                       font.legend=10,
                       xlab="Time(years)",
                       break.time.by = 2,
                       palette = bioCol,
                       surv.median.line = "hv",
                       risk.table=T,
                       cumevents=F,
                       risk.table.height=.25)
    
    #输出图形
    pdf(file=paste(filepath,i,"_survial.pdf",sep = ""), width=6.5, height=6.25, onefile=FALSE)
    print(surPlot)
    dev.off()
  }
}
