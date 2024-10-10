setwd("E:/基因建模/gene/KEGG")

library(KEGGREST)

# 定义路径
pathways <- c('hsa04370','hsa04151','hsa04010','hsa04330','hsa04350','hsa04310')

# 初始化一个空的数据框用于存储所有基因信息
all_genes <- data.frame()

# 遍历每个通路
for (pathway_id in pathways) {
  # 获取通路信息
  pathway_info <- keggGet(pathway_id)
  
  # 提取基因信息
  genes <- unlist(lapply(gs[[1]]$GENE, function(x) strsplit(x, ";")[[1]]))
  
  genelist <- genes[which((seq_along(genes) %% 3) == 2)]
  
  # 将基因信息转换为数据框
  genelist <-  data.frame(genelist)
  
  # 合并到总的数据框中
  all_genes <- rbind(all_genes, genelist)
}

# 去除重复的基因
all_genes <- unique(all_genes)

# 输出到CSV文件
write.table(all_genes, "genelist.csv", row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
