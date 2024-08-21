#加载包
library(limma)
library(ggplot2)
library(scatterplot3d)

#工作路径
setwd("E:/CORD/PCA")

#读入
rt=read.table("risk.txt", header=T, sep="\t", check.names=F,row.names = 1)

#获取分组信息
risk=as.vector(rt$risk)

#删除生存时间及状态,,风险得分及分组信息
data=rt[,3:(ncol(rt)-2)]

#PCA分析
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)

#绘制图形
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], risk=risk)
PCA.mean=aggregate(PCA[,1:2], list(risk=PCA$risk), mean)

#画图2d
#pdf(file="PCA.2d.pdf", width=5.5, height=4.75)
png("PCA.2d.png", width=5.5, height=4.75, units = 'in', res = 800)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color=risk, shape=risk)) +
  scale_colour_manual(name="risk", values =c("DarkOrchid","Orange2"))+
  theme_bw()+
  labs(title ="PCA")+
  theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))+
  annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$risk, cex=7)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#画图3d
color=ifelse(risk=="high","DarkOrchid","Orange2")
#pdf(file="PCA.3d.pdf", width=7, height=7)
png("PCA.3d.png", width=8, height=8, units = 'in', res = 800)
par(oma=c(1,1,2.5,1))
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color, angle=60)
legend("top", legend = c("High risk","Low risk"),pch = 16, inset = -0.2, box.col="white", xpd = TRUE, horiz = TRUE,col=c("DarkOrchid","Orange2"))
dev.off()
