# 加载必要的库
library(survival)
library(survminer)
library(xlsx)
library(VennDiagram)
library(tidyverse)
library(Boruta)
library(glmnet)

# 设置工作目录
setwd("E:/CRC_model")

DEG <- c("Wilcoxon", "limma", "edgeR", "DESeq2")

for (method in DEG) {
  
  # 读取差异基因
  data <- read.table(paste("TCGA.diff.", method, ".txt", sep=""), header=F, sep="\t", check.names=F)
  DEG_genes <- data[,1]
  
  # 读取ARG基因
  data <- read.table("ARG_gene.txt", header=F, sep="\t", check.names=F)
  ARG <- data[,1]
  
  # 画Venn图
  venn.diagram(x = list('ARG' = ARG, method = DEG_genes),
               filename = paste("VN_", method, ".png", sep = ""),
               fill = c("darkorange1", "green"))
  
  # 交集基因
  intersectGenes <- intersect(DEG_genes, ARG)
  write.table(file = paste("intersectGenes_", method, ".txt", sep = ""), 
              intersectGenes, sep="\t", quote=F, col.names=F, row.names=F)
  
  # 读取TCGA数据
  TCGA <- read.table("E:/CRC_model/PC/TCGA_cli.txt", header=T, sep="\t", check.names=F, row.names=1)
  
  # 选择需要的列
  keep_columns <- c("time", "status", intersectGenes)
  keep_columns <- intersect(keep_columns, colnames(TCGA))
  final_data <- TCGA[, keep_columns]
  
  # 单因素Cox分析
  p.value = 0.05
  outTab <- data.frame()
  for (i in colnames(final_data[,3:ncol(final_data)])) {
    cox <- coxph(Surv(time, status) ~ final_data[,i], data = final_data)
    coxSummary = summary(cox)
    coxP = coxSummary$coefficients[,"Pr(>|z|)"]
    if (coxP < p.value) {
      outTab <- rbind(outTab,
                      cbind(id=i,
                            HR=coxSummary$conf.int[,"exp(coef)"],
                            HR.95L=coxSummary$conf.int[,"lower .95"],
                            HR.95H=coxSummary$conf.int[,"upper .95"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
      )
    }
  }
  
  # 输出单因素Cox结果
  write.table(outTab, file=paste("uniCox_", method, ".txt", sep=""), sep="\t", row.names=F, quote=F)
  
  if (nrow(outTab) > 0) {
    # 多因素Cox分析
    multiCox <- coxph(Surv(time, status) ~ ., data = final_data[,c("time","status",outTab$id)])
    multiCox <- step(multiCox, direction = "both")
    multiCoxSum <- summary(multiCox)
    
    # 输出多因素Cox结果
    outTab <- data.frame(
      coef=multiCoxSum$coefficients[,"coef"],
      HR=multiCoxSum$conf.int[,"exp(coef)"],
      HR.95L=multiCoxSum$conf.int[,"lower .95"],
      HR.95H=multiCoxSum$conf.int[,"upper .95"],
      pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"]
    )
    outTab <- cbind(id=row.names(outTab), outTab)
    write.table(outTab, file=paste("multiCox_", method, ".txt", sep=""), sep="\t", row.names=F, quote=F)
    
    # 特征选择
    traindata <- final_data
    set.seed(42)
    result_boruta <- Boruta(
      x = traindata[,-c(1,2)],
      y = Surv(traindata$time, traindata$status), 
      doTrace = 1,
      maxRuns = 100
    )
    result_boruta <- TentativeRoughFix(result_boruta)
    selected_features <- getSelectedAttributes(result_boruta)
    traindata <- traindata[, c("time", "status", selected_features)]
    if (length(selected_features) > 1) {
      # LASSO
      x <- as.matrix(traindata[, selected_features])
      y <- Surv(traindata$time, traindata$status)
      
      cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C", nfolds = 10, intercept = FALSE)
      
      coefficient <- coef(cvfit, s = cvfit$lambda.min)
      Active.Index <- which(as.numeric(coefficient) != 0)
      active.coefficients <- as.numeric(coefficient)[Active.Index]
      sig_multi_cox <- rownames(coefficient)[Active.Index]
      
      coef_dataframe <- data.frame("features" = sig_multi_cox, 'coef' = active.coefficients)
      
      # 使用之前的方法计算风险评分
      traindata$risk_score <- rowSums(as.matrix(traindata[, sig_multi_cox, drop = FALSE]) * active.coefficients)
      
      # 读取测试集数据
      testdata <- read.table("E:/CRC_model/PC/GEO_cli.txt", header=T, sep="\t", check.names=F, row.names=1)
      testdata <- testdata[, c("time", "status", sig_multi_cox)]
      testdata$risk_score <- rowSums(as.matrix(testdata[, sig_multi_cox, drop = FALSE]) * active.coefficients)
      
      # 寻找最佳截断点
      tryCatch({
        cutpoint <- surv_cutpoint(traindata, time = "time", event = "status", variables = "risk_score")
        cutpoint <- cutpoint$cutpoint$cutpoint
      }, error = function(e) {
        cat("Error in finding cutpoint:", conditionMessage(e), "\n")
        cat("Using median as cutpoint.\n")
        cutpoint <- median(traindata$risk_score)
      })
      
      # 分组
      traindata$risk_group <- ifelse(traindata$risk_score > cutpoint, "High", "Low")
      testdata$risk_group <- ifelse(testdata$risk_score > cutpoint, "High", "Low")
      
      # 生存分析
      surv_train <- survfit(Surv(time, status) ~ risk_group, data = traindata)
      surv_test <- survfit(Surv(time, status) ~ risk_group, data = testdata)
      
      # 绘制生存曲线
      p1 <- ggsurvplot(
        surv_train, 
        data = traindata, 
        pval = TRUE, 
        title = paste("Training Set Survival Curve -", method),
        risk.table = TRUE,       
        risk.table.height = 0.25 
      )
      
      p2 <- ggsurvplot(
        surv_test, 
        data = testdata, 
        pval = TRUE, 
        title = paste("Test Set Survival Curve -", method),
        risk.table = TRUE,       
        risk.table.height = 0.25 
      )
      
      # 检查p值并保存图片
      if (p1$plot$labels$pval < 0.05 & p2$plot$labels$pval < 0.05) {
        pdf(paste("significant_survival_curves_", method, ".pdf", sep=""), width = 12, height = 10)
        print(p1)
        print(p2)
        dev.off()
      }
      
      # 输出结果
      cat("Method:", method, "\n")
      cat("Number of selected features:", length(sig_multi_cox), "\n")
      cat("Selected features:", paste(sig_multi_cox, collapse = ", "), "\n\n")
    } else {
      cat("Less than 2 features were selected by Boruta for method:", method, ". Skipping this method.\n\n")
    }
  } else {
    cat("No genes passed univariate Cox analysis for method:", method, ". Skipping this method.\n\n")
  }
}

# 显示警告信息
warnings()
