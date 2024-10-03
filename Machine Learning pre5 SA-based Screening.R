# R语言机器学习-特征筛选-无监督筛选+基于SA的筛选

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


# 注册多核并行
library(doParallel)
cl <- makePSOCKcluster(3)
registerDoParallel(cl)

############################################################

# Simulated annealing (SA)--模拟退火

################################ 回归

# 基于随机森林

library(caret)
set.seed(42)
result_sa_rf <- safs(
  x = as.data.frame(trainx),
  y = trainy,
  iters = 100, # 迭代次数
  safsControl = safsControl(
    method = "cv",
    number = 5,
    improve = 5,  # 重启的间隔迭代次数
    functions = rfSA,
    verbose = T
  )
)
result_sa_rf
# 筛选得到的变量
result_sa_rf$optVariables
# plot(result_sa_rf) +
#   theme_bw()

# # 预测
# predtest_sa_rf <-  predict(result_sa_rf, testx)


# 类似于rfSA的函数还有treebagSA
# 也可以采用caretSA，然后自定义train method

################################ 二分类

# 基于支持向量机
library(MLmetrics)
mysafuncs <- caretSA
mysafuncs$fitness_extern <- 
  function(data, lev = levels(data$obs), model = NULL) {
    c(
      twoClassSummary(data = data, lev = levels(data$obs), model),
      prSummary(data = data, lev = levels(data$obs), model),
      defaultSummary(data = data, lev = levels(data$obs), model)
    )
  }

library(caret)
set.seed(42)
result_sa_svm <- safs(
  x = as.data.frame(trainx),
  y = traincl,
  iters = 100, # 迭代次数
  safsControl = safsControl(
    method = "cv",
    metric = c(internal = "ROC", external = "ROC"),
    maximize = c(internal = TRUE, external = TRUE),
    number = 5,
    improve = 5,  # 重启的间隔迭代次数
    functions = mysafuncs,
    verbose = T
  ),
  
  method = "svmRadial",
  tuneLength = 3,
  trControl = trainControl(
    method = "cv",
    number = 5,
    summaryFunction = twoClassSummary,
    classProbs = TRUE
  ),
  metric = "ROC"
)
result_sa_svm
result_sa_svm$optVariables
# plot(result_sa) +
#   theme_bw()

# # 预测
# predtest_sa_svm <-  predict(result_sa_svm, testx)


################################ 多分类

# 基于treebag

library(caret)
set.seed(42)
result_sa_treebag <- safs(
  x = as.data.frame(trainx),
  y = traincl3,
  iters = 100, # 迭代次数
  safsControl = safsControl(
    method = "cv",
    number = 5,
    improve = 5,  # 重启的间隔迭代次数
    functions = treebagSA,
    verbose = T
  )
)
result_sa_treebag
# 筛选得到的变量
result_sa_treebag$optVariables
# plot(result_sa_treebag) +
#   theme_bw()

# # 预测
# predtest_sa_treebag <-  predict(result_sa_treebag, testx)