# R语言机器学习-特征筛选-无监督筛选+RFE

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

# Recursive Feature Elimination 递归特征消除 RFE

# 基于caret包的rfe函数

################################# 回归

# 基于随机森林
library(caret)
set.seed(42)
result_rfe_rf <- rfe(
  x = as.data.frame(trainx),
  y = trainy,
  size = 2:ncol(trainx), # 指定的保留的自变量的个数
  metric = "RMSE",
  maximize = F,
  rfeControl = rfeControl(
    functions = rfFuncs,
    method = "cv",
    number = 5,
    verbose = T
  )
)
result_rfe_rf$fit
# 筛选得到的变量
result_rfe_rf$optVariables
# 图示
result_rfe_rf$results %>%
  ggplot(aes(x = Variables, y = RMSE)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = length(result_rfe_rf$optVariables),
             color = "red",
             linetype = "dashed") +
  scale_x_continuous(breaks = result_rfe_rf$results$Variables) +
  labs(title = "RF-RFE结果") +
  theme_classic()
# 变量重要性
result_rfe_rf$variables %>%
  group_by(var) %>%
  summarise(Importance = mean(Overall)) %>%
  arrange(Importance) %>%
  mutate(var = forcats::as_factor(var),
         sel = ifelse(var %in% result_rfe_rf$optVariables, "S", "D")) %>%
  ggplot(aes(x = Importance, y = var, fill = sel)) +
  geom_col(show.legend = F) +
  labs(y = "", title = "所有变量的重要性排序") +
  theme_minimal()

# # 预测
# predtest_rfe_rf <- predict(result_rfe_rf, testx)

# 类似于rfFuncs的函数还有lmFuncs, treebagFuncs, nbFuncs
# # 也可以采用caretFuncs，然后自定义train method

################################# 二分类


# 基于支持向量机
library(MLmetrics)
myrfefuncs <- caretFuncs
myrfefuncs$summary <- 
  function(data, lev = levels(data$obs), model = NULL) {
    c(
      twoClassSummary(data = data, lev = levels(data$obs), model),
      prSummary(data = data, lev = levels(data$obs), model),
      defaultSummary(data = data, lev = levels(data$obs), model)
    )
  }

library(caret)
set.seed(42)
result_rfe_svm <- rfe(
  x = as.data.frame(trainx),
  y = traincl,
  size = 2:ncol(trainx), # 指定的保留的自变量的个数
  metric = "ROC",
  maximize = T,
  rfeControl = rfeControl(
    functions = myrfefuncs,
    method = "cv",
    number = 5,
    verbose = T
  ),
  method = 'svmRadial',
  trControl = trainControl(
    method = "cv",
    number = 3,
    summaryFunction = twoClassSummary,
    classProbs = TRUE
  ),
  prob.model = T
)
result_rfe_svm
# 筛选得到的变量
result_rfe_svm$optVariables
# 图示
result_rfe_svm$results %>%
  ggplot(aes(x = Variables, y = ROC)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = length(result_rfe_svm$optVariables),
             color = "red",
             linetype = "dashed") +
  scale_x_continuous(breaks = result_rfe_svm$results$Variables) +
  labs(title = "SVM-RFE结果") +
  theme_classic()
# 变量重要性
result_rfe_svm$variables %>%
  group_by(var) %>%
  summarise(Importance = mean(Overall)) %>%
  arrange(Importance) %>%
  mutate(var = forcats::as_factor(var),
         sel = ifelse(var %in% result_rfe_svm$optVariables, "S", "D")) %>%
  ggplot(aes(x = Importance, y = var, fill = sel)) +
  geom_col(show.legend = F) +
  labs(y = "", title = "所有变量的重要性排序") +
  theme_minimal()

# # 预测
# predtest_rfe_svm <-  predict(result_rfe_svm, testx)

################################# 多分类

# 基于treebag
library(caret)
set.seed(42)
result_rfe_treebag <- rfe(
  x = as.data.frame(trainx),
  y = traincl3,
  size = 2:ncol(trainx), # 指定的保留的自变量的个数
  rfeControl = rfeControl(
    functions = treebagFuncs,
    method = "cv",
    number = 5,
    verbose = T
  )
)
result_rfe_treebag
# 筛选得到的变量
result_rfe_treebag$optVariables
# 图示
result_rfe_treebag$results %>%
  ggplot(aes(x = Variables, y = Accuracy)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = length(result_rfe_treebag$optVariables),
             color = "red",
             linetype = "dashed") +
  scale_x_continuous(breaks = result_rfe_treebag$results$Variables) +
  labs(title = "Treebag-RFE结果") +
  theme_classic()
# 变量重要性
result_rfe_treebag$variables %>%
  group_by(var) %>%
  summarise(Importance = mean(Overall)) %>%
  arrange(Importance) %>%
  mutate(var = forcats::as_factor(var),
         sel = ifelse(var %in% result_rfe_treebag$optVariables, "S", "D")) %>%
  ggplot(aes(x = Importance, y = var, fill = sel)) +
  geom_col(show.legend = F) +
  labs(y = "", title = "所有变量的重要性排序") +
  theme_minimal()

# # 预测
# predtest_rfe_treebag <-  predict(result_rfe_treebag, testx)