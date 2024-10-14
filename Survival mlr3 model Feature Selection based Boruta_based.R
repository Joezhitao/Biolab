# Joe_2024_10_12
# 加载包
library(tidyverse)
library(survival)
library(Boruta)

# 读取数据
sadata <- readr::read_csv("F:/Mach_learn_data/sadata.csv")
colnames(sadata)

# 修正变量类型
# 分类变量转factor
for(i in c(3, 5, 9)){
  sadata[[i]] <- factor(sadata[[i]])
}

# 剔除变量-无关变量
sadata$pid <- NULL
# 剔除样本-含有缺失值的样本，时间非正的样本
# sadata <- na.omit(sadata) # 剔除任意变量下有缺失值的样本
sadata <- sadata %>%
  # na.omit() %>%
  # drop_na(age) %>%  # 剔除指定变量下有缺失值的样本
  filter(rfstime > 0) # 剔除时间非正的样本

# 数据概况
skimr::skim(sadata)

# 感兴趣的时间节点
range(unique(sadata$rfstime))
itps <- c(365 * c(1, 3, 5))
itps
table(cut(sadata$rfstime, c(0, itps, Inf)))


###################################################

# 数据拆分构建任务对象
set.seed(42)
datasplit <- rsample::initial_split(
  sadata, prop = 0.8, strata = rfstime, breaks = 10
)
traindata <- rsample::training(datasplit) %>%
  sample_n(nrow(.))
testdata <- rsample::testing(datasplit) %>%
  sample_n(nrow(.))

# 拆分得到的数据的生存曲线比较
sadata2 <- sadata
sadata2$set <- "Test"
sadata2$set[datasplit$in_id] <- "Train"
sadata2$set <- factor(sadata2$set)
sfit_set <- survfit(Surv(rfstime, status) ~ set, data=sadata2)
survminer::ggsurvplot(
  sfit_set,
  pval=TRUE, 
  pval.coord = c(0.1, 0.5),
  risk.table=TRUE,
  ggtheme = survminer::theme_survminer() +
    theme(text = element_text(family = "serif")),
  font.family = "serif"
)

###################################################

# 数据预处理
datarecipe_coxph <- 
  recipes::recipe(rfstime + status ~ ., traindata) %>%
  recipes::prep()
datarecipe_coxph

# 按方处理训练集和测试集
traindata2 <- 
  recipes::bake(datarecipe_coxph, new_data = NULL) %>%
  dplyr::select(rfstime, status, everything())
testdata2 <- 
  recipes::bake(datarecipe_coxph, new_data = testdata) %>%
  dplyr::select(rfstime, status, everything())

# 特征筛选
set.seed(42)
result_boruta <- Boruta(
  x = traindata2[,-c(1,2)],
  y = Surv(traindata2$rfstime, traindata2$status), 
  doTrace = 1,
  maxRuns = 100
)
# 筛选结果概况
result_boruta

# 对于存疑变量进一步确认
result_boruta <- TentativeRoughFix(result_boruta)
# 具体结果
attStats(result_boruta)
# 图示
par(las = 3, mar = c(8, 4, 4, 2))
plot(result_boruta, xlab = "", main = "Boruta算法筛选结果")
legend("topleft", 
       legend = c("confirmed", "rejected"),
       pch = 22,
       pt.bg = c("green", "red"))

# 筛选得到的变量
getSelectedAttributes(result_boruta)

# 筛选之后的数据
traindata3 <- traindata2 %>%
  dplyr::select(rfstime, status, 
                getSelectedAttributes(result_boruta))
testdata3 <- testdata2 %>%
  dplyr::select(rfstime, status, 
                getSelectedAttributes(result_boruta))












