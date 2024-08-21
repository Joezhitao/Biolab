# 加载必需的库
library(ggplot2)
library(pheatmap)
library(reshape2)
library(Hmisc)

# 设置工作目录
setwd("E:/CORD/DEG")

# 读取数据
df <- read.table("TCGA_CORD_TPM.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)

# 转置数据框
df_transposed <- t(df)

# 将转置后的数据框转换为数据框格式，并将行名设置为第一列
df_transposed <- as.data.frame(df_transposed)
df_transposed$Sample <- rownames(df_transposed)
rownames(df_transposed) <- NULL

# 将所有基因表达值进行 log2(TPM + 1) 转换
df_transposed[,-ncol(df_transposed)] <- log2(df_transposed[,-ncol(df_transposed)] + 1)

# 计算 NELL2 的中位数
nell2_median <- median(df_transposed$NELL2, na.rm = TRUE)

# 根据 NELL2 的表达量将样本分为高表达组和低表达组
df_transposed$ExpressionGroup <- ifelse(df_transposed$NELL2 > nell2_median, "high", "low")

# 查看结果
print(df_transposed)

# 计算 NELL2 与其他基因的 Pearson 相关性
correlations <- sapply(df_transposed[, -which(names(df_transposed) %in% c("Sample", "ExpressionGroup", "NELL2"))], 
                       function(x) cor(x, df_transposed$NELL2, use = "complete.obs"))

# 找到相关性绝对值最强的 10 个基因
top_genes <- names(sort(abs(correlations), decreasing = TRUE))[1:10]

# 创建一个新的数据框，仅保留 NELL2 和相关性最强的 10 个基因
df_top_genes <- df_transposed[, c("NELL2", top_genes, "Sample", "ExpressionGroup")]

# 查看结果
print(df_top_genes)

# 加载必要的包
library(pheatmap)

# 按照 NELL2 的表达量对样本进行排序
sorted_df <- df_top_genes[order(df_top_genes$NELL2), ]

# 提取用于热图的数据（移除 Sample 和 ExpressionGroup 列）
heatmap_data <- sorted_df[, -which(names(sorted_df) %in% c("Sample", "ExpressionGroup"))]

# 转置数据框以便样本在列上，基因在行上
heatmap_data_t <- t(heatmap_data)

# 计算 p 值和 R 值
cor_results <- apply(heatmap_data_t[-1, ], 1, function(x) {
  test <- cor.test(x, heatmap_data_t["NELL2", ])
  list(p_value = test$p.value, r_value = test$estimate)
})

# 创建基因名称向量，附加 p 值和 R 值
gene_names <- sapply(names(cor_results), function(name) {
  p_val <- cor_results[[name]]$p_value
  r_val <- cor_results[[name]]$r_value
  sprintf("%s\np=%.3f R=%.3f", name, p_val, r_val)
})

# 绘制热图
pheatmap(
  heatmap_data_t,
  cluster_cols = FALSE,  # 不对列（样本）进行聚类
  cluster_rows = TRUE,   # 对行（基因）进行聚类
  labels_row = c("NELL2", gene_names),  # 行标签
  show_rownames = TRUE,
  show_colnames = FALSE,
  scale = "row"  # 按行标准化
)
#######################################################################################

# 加载必要的包
library(pheatmap)
library(ggplot2)
library(cowplot)
library(grid)  # 用于单位设置

# 按照 NELL2 的表达量对样本进行排序
sorted_df <- df_top_genes[order(df_top_genes$NELL2), ]

# 提取用于热图的数据（移除 Sample、ExpressionGroup 和 NELL2 列）
heatmap_data <- sorted_df[, -which(names(sorted_df) %in% c("Sample", "ExpressionGroup", "NELL2"))]

# 转置数据框以便样本在列上，基因在行上
heatmap_data_t <- t(heatmap_data)

# 计算 p 值和 R 值
cor_results <- apply(heatmap_data_t, 1, function(x) {
  test <- cor.test(x, sorted_df$NELL2)
  list(p_value = test$p.value, r_value = test$estimate)
})

# 根据 p 值生成星号
get_stars <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

# 创建基因名称向量，附加星号和 R 值
gene_names <- sapply(names(cor_results), function(name) {
  p_val <- cor_results[[name]]$p_value
  r_val <- cor_results[[name]]$r_value
  stars <- get_stars(p_val)
  sprintf("%s %s R=%.3f", name, stars, r_val)
})

# 绘制热图
heatmap <- pheatmap(
  heatmap_data_t,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  labels_row = gene_names,
  show_rownames = TRUE,
  show_colnames = FALSE,
  scale = "row",
  silent = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  fontsize_row = 8  # 调整行名字体大小
)

# 绘制折线图
line_plot <- ggplot(line_data, aes(x = Sample, y = NELL2, fill = Group)) +
  geom_area(alpha = 0.5) +
  geom_line() +
  scale_fill_manual(values = c("high" = "red", "low" = "blue")) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),  # 移除网格线
    axis.text.x = element_blank(),  # 隐藏 x 轴刻度标签
    axis.ticks.x = element_blank(),  # 隐藏 x 轴刻度线
    axis.ticks.y = element_line(color = "black"),  # 显示 y 轴刻度线
    axis.ticks.length.y = unit(0.2, "cm"),  # 设置刻度线长度
    axis.line = element_line(color = "black")  # 显示轴线
  ) +
  labs(y = "NELL2 Log2(TPM+1)", x = NULL) +
  scale_x_continuous(expand = c(0, 0)) +  # 去掉x轴的空白
  scale_y_continuous(breaks = seq(0, ceiling(max(line_data$NELL2)), by = 2))  # 设置 y 轴间隔为 2

# 绘制热图，添加渐变色图的名字
heatmap <- pheatmap(
  heatmap_data_t,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  labels_row = gene_names,
  show_rownames = TRUE,
  show_colnames = FALSE,
  scale = "row",
  silent = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  fontsize_row = 8,  # 调整行名字体大小
  legend = TRUE,  # 显示图例
  legend_labels = "z-score"  # 添加渐变色图名字
)

# 保存为 PNG 并合并图形
png("E:/CORD/DEG/heatmap2.png", width = 10, height = 10, units = "in", res = 800)
plot_grid(
  line_plot, 
  heatmap$gtable, 
  ncol = 1, 
  align = "v", 
  axis = "tb",  # 指定在顶部和底部对齐
  rel_heights = c(1, 4)
)
dev.off()

# 合并热图和折线图
png("E:/CORD/DEG/heatmap.png", width = 8, height = 8, units = "in", res = 800)
dev.off()