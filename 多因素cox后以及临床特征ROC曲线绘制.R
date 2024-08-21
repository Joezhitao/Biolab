#加载包
library(pROC)
library(ggplot2)
library(survival)
library(survminer)
library(timeROC)

#设置工作目录
setwd("E:/CORD/ROC")

#读取风险输入文件
risk=read.table("risk.txt", header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("time", "state", "riskScore")]

#读取临床数据文件
cli=read.table("clinical.txt", header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#定义颜色
bioCol=c("DarkOrchid","Orange2","MediumSeaGreen","NavyBlue","#8B668B","#FF4500","#135612","#561214")
#ROC分析
roc1 <- roc(rt$state ~ rt$riskScore)

# ROC曲线
#$pdf(file=paste0("ROC.riskscore.pdf"), width=5, height=5)
png(file=paste0("ROC.riskscore.png"), width=5, height=5, units="in", res=600)
plot(roc1, print.auc=TRUE, col=bioCol, legacy.axes=T)
dev.off()

#?timeROC
#T：事件时间
#delta：事件状态. 删失数据编码为0.
#marker ：计算ROC的biomaker，默认是marker值越大，事件越可能发生；反之的话，前面加-号。
#other_markers：矩阵形式输入，可多个marker，类似协变量. 默认值other_markers=NULL.
#cause：所关心的事件结局。没有竞争风险（Without competing risks）中，必须是非删失数据的编码方式，一般为1。
#存在竞争风险（With competing risks）中，和所关心的事件结局一致，通常为1 or 2.
#weighting：计算方法，默认是weighting="marginal"，KM模型；
#weighting="cox" 和weighting="aalen"分别为COX模型和additive Aalen 模型。
#times：想计算的ROC曲线的时间节点。
#ROC：默认值ROC = TRUE，保存sensitivities 和 specificties值。
#iid: 默认值iid = FALSE。iid = TRUE 才会保存置信区间，但是样本量大了后，耗时耗资源。

#绘制1 3 5年的ROC曲线
ROC_rt=timeROC(T=risk$time,delta=risk$state,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
#pdf(file="ROC.all.pdf", width=5, height=5)
png(file=paste0("ROC.all.png"), width=5, height=5, units="in", res=600)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=4)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=4)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=4)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       #设置分组信息，也可删除
       col=bioCol[1:3], lwd=4, bty = 'n',title = "All set")
dev.off()

#绘制临床的ROC曲线
predictTime=3     #定义预测年限
aucText=c()
#pdf(file="cliROC.all.pdf", width=5.5, height=5.5)
png(file=paste0("cliROC.all.png"), width=5.5, height=5.5, units="in", res=600)
#绘制风险得分的ROC曲线
i=3
ROC_rt=timeROC(T=risk$time,
               delta=risk$state,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=4)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)
#对临床数据进行循环，绘制临床数据的ROC曲线
for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$time,
                 delta=rt$state,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=4, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
#绘制图例，得到ROC曲线下的面积
#设置分组信息，也可删除
legend("bottomright", aucText,lwd=4,bty="n",col=bioCol[1:(ncol(rt)-1)],title = "All set")
dev.off()