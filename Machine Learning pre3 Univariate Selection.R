# R语言机器学习-特征筛选-无监督筛选+单变量筛选

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
set.seed(1220)
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

# 单变量筛选

# 基于caret包的sbf函数

################################# 回归

# 以随机森林为后建模型
library(caret)
set.seed(1220)
result_sbf_rf <- sbf(
  x = as.data.frame(trainx),
  y = trainy,
  sbfControl = sbfControl(
    functions = rfSBF, # 这个模型是单变量筛选之后的建模
    method = "cv", # 交叉验证
    number = 5, # 交叉验证次数
    saveDetails = T 
  )
)
result_sbf_rf
# 以后建模型性能比较得到的筛选变量
result_sbf_rf$optVariables

# # 预测
# predtest_sbf_rf <- predict(result_sbf_rf, testx)

# 将过滤结果求交集得到的筛选变量
# 此时可以不用关注sbfControl中的function
result_sbf_rf2 <- result_sbf_rf$variables
names(result_sbf_rf2) <- 
  paste0("Set", 1:length(result_sbf_rf2))
RVenn::overlap(RVenn::Venn(result_sbf_rf2))
VennDiagram::venn.diagram(result_sbf_rf2, 
                          "result_sbf_rf.png",
                          imagetype = "png",
                          main = "单变量筛选结果汇总",
                          main.pos = c(0.5, 0.95),
                          margin = 0.2,
                          fill = rainbow(length(result_sbf_rf2)))
RVenn::setmap(RVenn::Venn(result_sbf_rf2), 
              element_clustering = FALSE, 
              set_clustering = FALSE,
              title = "单变量筛选结果汇总")

# 类似于rfSBF的函数还有lmSBF, treebagSBF, ldaSBF, nbSBF
# # 也可以采用caretSBF，然后自定义train method


################################# 二分类

# 以SVM为后建模型
library(caret)
set.seed(1220)
result_sbf_svm <- sbf(
  x = as.data.frame(trainx),
  y = traincl,
  sbfControl = sbfControl(
    functions = caretSBF, # 这个模型是单变量筛选之后的建模
    method = "cv", # 交叉验证
    number = 5, # 交叉验证次数
    saveDetails = T 
  ),
  method = "svmRadial",  # 高斯核支持向量机
  trControl = trainControl(
    method = "cv",
    number = 3,
    summaryFunction = twoClassSummary,
    classProbs = T
  ),
  tuneLength = 5,
  metric = "ROC",
  prob.model = T
)
result_sbf_svm$fit$finalModel
# 以后建模型性能比较得到的筛选变量
result_sbf_svm$optVariables

# # 预测
# predtest_sbf_svm <- predict(result_sbf_svm, testx)

# 将过滤结果求交集得到的筛选变量
result_sbf_svm2 <- result_sbf_svm$variables
names(result_sbf_svm2) <- 
  paste0("Set", 1:length(result_sbf_svm2))
RVenn::overlap(RVenn::Venn(result_sbf_svm2))
VennDiagram::venn.diagram(result_sbf_svm2, 
                          "result_sbf_svm.png",
                          imagetype = "png",
                          main = "单变量筛选结果汇总",
                          main.pos = c(0.5, 0.95),
                          margin = 0.2,
                          fill = rainbow(length(result_sbf_svm2)))
RVenn::setmap(RVenn::Venn(result_sbf_svm2), 
              element_clustering = FALSE, 
              set_clustering = FALSE,
              title = "单变量筛选结果汇总")


################################# 多分类

# 以决策树为后建模型
library(caret)
set.seed(1220)
result_sbf_dt <- sbf(
  x = as.data.frame(trainx),
  y = traincl3,
  sbfControl = sbfControl(
    functions = caretSBF, # 这个模型是单变量筛选之后的建模
    method = "cv", # 交叉验证
    number = 5, # 交叉验证次数
    saveDetails = T
  ),
  method = "rpart",  # 决策树
  trControl = trainControl(
    method = "cv",
    number = 3
  ),
  tuneLength = 5
)
result_sbf_dt
# 以后建模型性能比较得到的筛选变量
result_sbf_dt$optVariables

# # 预测
# predtest_sbf_dt <- predict(result_sbf_dt, testx)

# 将过滤结果求交集得到的筛选变量
result_sbf_dt2 <- result_sbf_dt$variables
names(result_sbf_dt2) <- 
  paste0("Set", 1:length(result_sbf_dt2))
RVenn::overlap(RVenn::Venn(result_sbf_dt2))
VennDiagram::venn.diagram(result_sbf_dt2, 
                          "result_sbf_dt.png",
                          imagetype = "png",
                          main = "单变量筛选结果汇总",
                          main.pos = c(0.5, 0.95),
                          margin = 0.2,
                          fill = rainbow(length(result_sbf_dt2)))
RVenn::setmap(RVenn::Venn(result_sbf_dt2), 
              element_clustering = FALSE, 
              set_clustering = FALSE,
              title = "单变量筛选结果汇总")


# 基于tidymodels框架的colino包提供了一些特征筛选的step函数可以集成到recipes中