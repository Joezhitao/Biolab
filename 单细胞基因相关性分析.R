library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

filepath <- paste("E:/B_group/data_backup/B_group_sample",".RDS", sep = "")
pbmc <- readRDS(filepath)

cluster.averages <- AverageExpression(pbmc, group.by="group")

expr<-as.matrix(cluster.averages[["RNA"]])

expr <- log2(expr+1)
dim(expr)
expr[1:4,1:4]
dat<-as.data.frame(t(expr))
dat[1:4,1:4]
library(ggstatsplot)
png("E:/B_group/横向比较/Ppara_Ccnd1.png",width = 12,height = 12,units = "in",res = 800)
ggscatterstats(dat,
               y =Ppara,
               x =Ccnd1,
               type = "pearson",
               centrality.para = "mean",
               margins = "both",
               xfill = "#009E73",
               yfill = "#D55E00",
               marginal.type = "histogram",
               title = "Relationship between Ppara and Ccnd1")
dev.off()

y <-as.numeric(dat[,"CDK6"])
colnames <-colnames(dat)
cor_data_df <- data.frame(colnames)
for (i in 1:length(colnames)){
  test <-cor.test(as.numeric(dat[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

write.csv(cor_data_df,file = "D:\\分析文件存储\\COR_GENE.csv")
