# R语言机器学习-特征筛选-无监督筛选+基于随机森林的筛选

# 读取数据
simdata <- readr::read_csv("F:/Mach_learn_data/simdata.csv")
colnames(simdata) 
# 修正变量类型
# 将分类变量转换为factor
for(i in c(7, 8, 20, 27, 28)){ 
  simdata[[i]] <- factor(simdata[[i]])
}
skimr::skim(simdata)
DataExplorer::plot_correlation(
  simdata,
  cor_args = list(use = "pairwise.complete.obs"),
  title = "所有变量的相关性",
  ggtheme = ggplot2::theme_minimal()
)

# 数据拆分
set.seed(42)
datasplit <- 
  rsample::initial_split(simdata, prop = 0.75, strata = Cl)
traindata <- rsample::training(datasplit)
testdata <- rsample::testing(datasplit)

#################################################################

# 无监督筛选
library(recipes)
# 配方
datarecipe <- recipe(Y+Cl+Cl3 ~ ., traindata) %>%
  # 主观剔除，样本id之类的变量或者主观分析要剔除的变量
  # step_rm(id) %>%
  # 按缺失比例剔除变量，threshold设定缺失比例界限
  step_filter_missing(all_predictors(), threshold = 0) %>%
  # 剔除零方差变量，所有样本取值一样的变量
  step_zv(all_predictors()) %>%
  # 分类独热编码，后续部分算法要求所有自变量为数值
  step_dummy(all_nominal_predictors()) %>%
  # 剔除近零方差变量，取值水平少，个别取值水平占绝大多数的变量
  step_nzv(all_predictors(), freq_cut = 95/5, unique_cut = 10) %>%
  # 剔除高度相关的变量
  step_corr(all_predictors()) %>%
  prep()
# 以上所述每一步处理了哪些变量
datarecipe


# 处理
trainx <- bake(datarecipe, new_data = traindata) %>%
  select(-Y, -Cl, -Cl3)
trainy <- traindata$Y
traincl <- traindata$Cl
traincl3 <- traindata$Cl3

testx <- bake(datarecipe, new_data = testdata) %>%
  select(-Y, -Cl, -Cl3)
testy <- testdata$Y
testcl <- testdata$Cl
testcl3 <- testdata$Cl3

#################################################################

# 有监督筛选
# 选其一或者多种结合均可

#################################

# 基于随机森林模型的筛选

###################################   1  回归或者分类问题均可

library(VSURF)
set.seed(42)
result_vsurf <- VSURF(
  x = trainx,
  y = trainy, # 回归
  # y = traincl, # 二分类
  # y = traincl3, # 多分类
  ##############################################
  # 以下参数一般保持默认
  # eliminate irrelevant variables
  ntree.thres = 500,
  nfor.thres = 20,
  nmin = 1,
  
  # select variables related to y for interpretation purpose
  ntree.interp = 100,
  nfor.interp = 10,
  nsd = 1,
  
  # eliminating redundancy
  ntree.pred = 100,
  nfor.pred = 10,
  nmj = 1,
  # 以上参数一般保持默认
  ##############################################
  
  # 随机森林模型引擎
  RFimplem = "randomForest"  # 如果数据量大速度慢换成ranger
)
# 筛选结果概况
summary(result_vsurf)
# plot(result_vsurf) # 以变量个数为横轴
par(las = 3) # 让下图中的变量名称竖直显示
plot(result_vsurf, var.names = TRUE) # 以变量名称为横轴，注意省略问题

# 分解结果
# 筛选得到的相关变量
result_vsurf$varselect.thres
colnames(trainx)[result_vsurf$varselect.thres]
# 图示
plot(result_vsurf, step = "thres")
plot(result_vsurf, step = "thres", imp.sd = F, var.names = TRUE)
plot(result_vsurf, step = "thres", imp.mean = F, var.names = TRUE)

# 筛选得到的解释变量
result_vsurf$varselect.interp
colnames(trainx)[result_vsurf$varselect.interp]
# 图示
plot(result_vsurf, step = "interp")
plot(result_vsurf, step = "interp", var.names = TRUE)

# 筛选得到的预测变量
result_vsurf$varselect.pred
colnames(trainx)[result_vsurf$varselect.pred]
# 图示
plot(result_vsurf, step = "pred")
plot(result_vsurf, step = "pred", var.names = TRUE)

par(las=1)

###################################   2   回归或者分类问题均可

library(Boruta)
set.seed(42)
result_boruta <- Boruta(
  x = trainx,
  y = trainy, # 回归
  # y = traincl, # 二分类
  # y = traincl3, # 多分类
  doTrace = 1,
  maxRuns = 100 # 当有部分变量未得到判定时，增加该参数值
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

###################################   3 只适用于分类问题

library(varSelRF)
set.seed(42)
result_varselrf <- varSelRF(
  xdata = trainx,
  # Class = traincl, # 二分类
  Class = traincl3, # 多分类
  c.sd = 1, # 1倍标准差规则，设定为0则采用最优规则
  vars.drop.frac = 0.2, # 每步变量筛选删减的变量比例
  whole.range = T # 变量删减至只剩下两个变量
)
# 筛选结果概况
result_varselrf

# 最终选定的变量
result_varselrf$selected.vars

# 图示过程
library(tidyverse)
history_varselrf <- result_varselrf$selec.history %>%
  mutate(`OOB+sd` = OOB + sd.OOB,
         `OOB-sd` = OOB - sd.OOB)
history_varselrf %>%
  pivot_longer(cols = c(3,5,6)) %>%
  ggplot(aes(x = Number.Variables,
             y = value,
             color = name)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = history_varselrf$Number.Variables) +
  labs(color = "", y = "OOB error", x = "Number of Variables",
       title = "varSelRF算法筛选结果") +
  theme_bw() +
  theme(legend.position = "bottom")

history_varselrf %>%
  ggplot(aes(x = Number.Variables,
             y = OOB)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = `OOB-sd`, ymax = `OOB+sd`)) +
  scale_x_continuous(breaks = history_varselrf$Number.Variables) +
  labs(color = "", y = "OOB error", x = "Number of Variables",
       title = "varSelRF算法筛选结果") +
  theme_bw()

# 所有自变量的最初重要性，scale=T
importance(result_varselrf$firstForest, scale = T) %>%
  as.data.frame() %>%
  rownames_to_column("x") %>%
  arrange(MeanDecreaseAccuracy) %>%
  mutate(x = as_factor(x)) %>%
  ggplot(aes(x = x, y = MeanDecreaseAccuracy)) +
  geom_col(fill = "skyblue") +
  labs(x = "") +
  coord_flip() +
  theme_bw()


