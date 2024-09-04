
library(Seurat)
library(dplyr)
library(patchwork)

filepath = "E:/B组数据备份(4.29)/单细胞结果/横向比较/组间/"
genepath = 'E:/B组数据备份(4.29)/横向比较基因/'

filepath <- paste("E:/B组数据备份(4.29)/鉴定后总样本备份/",'all_sample_decontX035_pc30',".RDS", sep = "")
#saveRDS(pbmc, filepath)
pbmc <- readRDS(filepath)

s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(pbmc, features = c("Pcna", "Top2a", "Mcm5", "Mki67"), ncol = 2)

pbmc <- ScaleData(pbmc, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(pbmc))

FeaturePlot(pbmc, features = c("Hp"))
