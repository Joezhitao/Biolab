library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)

# 设置工作目录
base_dir <- "E:/CRC"
setwd(base_dir)

# 读入生存数据
cli = read.table("time_CRC.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cli$time = cli$time / 365

# 读取输入文件
data = read.table("TCGA_CRC_count.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
dimnames = list(rownames(data), colnames(data))
data = matrix(as.numeric(as.matrix(data)), nrow = nrow(data), dimnames = dimnames)

# 正常和肿瘤数目
group = sapply(strsplit(colnames(data), "\\-"), "[", 4)
group = sapply(strsplit(group, ""), "[", 1)
data = data[, group == 0]

# 转置
data = t(data)
rownames(data) = substr(rownames(data), 1, 12)
rownames(data) = gsub('[.]', '-', rownames(data))

# 获取共同样本
sameSample = intersect(row.names(data), row.names(cli))
data = data[sameSample, ]
cli = cli[sameSample, ]
rt = cbind(cli, data)

# 根据KRAS和TBC1D16表达量的中位值，把样品分为四组
kras_median = quantile(rt[,"KRAS"], 0.5)
tbc1d16_median = quantile(rt[,"TBC1D16"], 0.5)

group = ifelse(rt[,"KRAS"] > kras_median & rt[,"TBC1D16"] <= tbc1d16_median, "KRAS High, TBC1D16 Low",
               ifelse(rt[,"KRAS"] > kras_median & rt[,"TBC1D16"] > tbc1d16_median, "KRAS High, TBC1D16 High",
                      ifelse(rt[,"KRAS"] <= kras_median & rt[,"TBC1D16"] <= tbc1d16_median, "KRAS Low, TBC1D16 Low",
                             "KRAS Low, TBC1D16 High")))
rt[,"group"] = group

# 生存分析
diff = survdiff(Surv(time, status) ~ group, data = rt)
pValue = 1 - pchisq(diff$chisq, df = length(unique(group)) - 1)
pValue = paste0("p=", sprintf("%.04f", pValue))
fit <- survfit(Surv(time, status) ~ group, data = rt)

# 绘制生存曲线
bioCol = c("Firebrick3", "MediumSeaGreen", "#6E568C", "#223D6C")
bioCol = bioCol[1:length(unique(group))]
surPlot = ggsurvplot(fit, 
                     data = rt,
                     conf.int = FALSE,
                     pval = pValue,
                     pval.size = 6,
                     legend.title = "Expression Groups",
                     legend.labs = levels(factor(rt[,"group"])),
                     legend = c(0.8, 0.8),
                     font.legend = 10,
                     xlab = "Time(years)",
                     break.time.by = 2,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table = TRUE,
                     cumevents = FALSE,
                     risk.table.height = .25)

pdf(file = "survival_KRAS_TBC1D16_four_groups.pdf", width = 6.5, height = 6.25, onefile = FALSE)
print(surPlot)
dev.off()

# 去除没有TNM分期和BMI的样本
rt_tnm = rt[!is.na(rt$TNM), ]
rt_bmi = rt[!is.na(rt$BMI), ]

# 根据BMI分组
rt_bmi$BMI_group = ifelse(rt_bmi$BMI < 23, "Low", "High")

# 定义基因列表，包括KRAS、TBC1D16和GEF相关基因
genes_of_interest <- c("KRAS", "TBC1D16", "SOS1", "SOS2", "VAV1", "VAV2", "VAV3", "RAPGEF1", "RANBP9", "RANBP10")

# 绘制箱线图
pdf(file = "expression_boxplots_TNM_BMI.pdf", width = 8, height = 6)

for (gene in genes_of_interest) {
  if (gene %in% colnames(rt_tnm)) {
    p_tnm <- ggplot(rt_tnm, aes(x = as.factor(TNM), y = .data[[gene]], fill = as.factor(TNM))) +
      geom_boxplot() +
      labs(title = paste(gene, "Expression by TNM Stage"), x = "TNM Stage", y = paste(gene, "Expression")) +
      theme_minimal() +
      scale_fill_manual(values = c("skyblue", "orange", "green", "purple")) +
      stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"), c("1", "4"),
                                            c("2", "3"), c("2", "4"), c("3", "4")),
                         method = "wilcox.test", label = "p.signif")
    print(p_tnm)
  }
  
  if (gene %in% colnames(rt_bmi)) {
    p_bmi <- ggplot(rt_bmi, aes(x = BMI_group, y = .data[[gene]], fill = BMI_group)) +
      geom_boxplot() +
      labs(title = paste(gene, "Expression by BMI Group"), x = "BMI Group", y = paste(gene, "Expression")) +
      theme_minimal() +
      scale_fill_manual(values = c("lightblue", "pink")) +
      stat_compare_means(method = "wilcox.test", label = "p.signif")
    print(p_bmi)
  }
}

# 绘制BMI分组的生存曲线
fit_bmi <- survfit(Surv(time, status) ~ BMI_group, data = rt_bmi)
surPlot_bmi = ggsurvplot(fit_bmi, 
                         data = rt_bmi,
                         conf.int = FALSE,
                         pval = TRUE,
                         pval.size = 6,
                         legend.title = "BMI Group",
                         legend.labs = levels(factor(rt_bmi$BMI_group)),
                         legend = c(0.8, 0.8),
                         font.legend = 10,
                         xlab = "Time(years)",
                         break.time.by = 2,
                         palette = c("lightblue", "pink"),
                         surv.median.line = "hv",
                         risk.table = TRUE,
                         cumevents = FALSE,
                         risk.table.height = .25)

print(surPlot_bmi)
dev.off()
