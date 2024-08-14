#加载
library(survival)
library(survminer)
library(xlsx)

#工作路径
setwd("E:/CORD/COX")

#读取表达文件
data=read.table("TCGA_CORD_count.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

#读取差异基因文件
data2=read.table("TCGA.diff.limma.txt", header=T, sep="\t", check.names=F,row.names = 1)
#提取差异基因的表达矩阵
data = data[rownames(data2),]

#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#仅保留肿瘤样本
data = data[,group == 0]

#转置
data=t(data)

#样本名仅保留前12字符
rownames(data)=substr(rownames(data),1,12)

#读入生存数据
cli <- read.xlsx("time_CORD.xlsx",sheetIndex = 1,header = T,encoding = "UTF-8")
#write.table(cli,'time_CORD.txt', sep="\t", quote=F, row.names = F)
# 将首列设置为行名
row.names(cli) <- cli[[1]]  # 将首列的值赋给行名
cli <- cli[-1]               # 删除首列

#cli=read.table("time_CORD.txt", header=T, sep="\t", check.names=F, row.names=1)
#生存时间以年为单位
cli$time=cli$time/365

#获取共同样本
sameSample=intersect(row.names(data), row.names(cli))
#提取共同样本
data=data[sameSample,]
cli=cli[sameSample,]
#合并
rt=cbind(cli, data)

#0.05,0.01,0.001
p.value = 0.01
#查找预后相关的基因
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox分析
  cox <- coxph(Surv(time, state) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP < p.value){
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

#输出单因素的结果
write.table(outTab, file="uniCox.txt", sep="\t", row.names=F, quote=F)


#画图
#读取输入文件
rt <- read.table("uniCox.txt",header=T,sep="\t",check.names=F,row.names=1)
#仅仅做前20个,实际运行不运行这一行
#rt = rt[1:20,]
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#绘制森林图
pdf(file="uniCoxforest.pdf", width = 12,height = height = nrow(rt)/13+15)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "red", "blue")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
axis(1)
dev.off()
