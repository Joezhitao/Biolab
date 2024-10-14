#Joe_2024_10_12
# 加载包
# library(devtools)
# devtools::install_github("mlr-org/mlr3proba")
# devtools::install_github("mlr-org/mlr3extralearners@*release")
library(tidyverse)
library(survival)
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
source("H:/Biolab/Biolab/tidyfuncs4sa.R")
# file.choose()

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

# 性能指标
measure_sa <- msrs("surv.cindex")
measure_sa

###################################################

# 数据拆分构建任务对象
set.seed(42)
datasplit <- rsample::initial_split(
  sadata, prop = 0.8, strata = status
)
traindata <- rsample::training(datasplit) %>%
  select(rfstime, status, everything()) %>%
  sample_n(nrow(.))
testdata <- rsample::testing(datasplit) %>%
  select(rfstime, status, everything()) %>%
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

# 训练任务对象
task_train <- as_task_surv(
  traindata, 
  time = "rfstime",
  event = "status", 
  type = "right"
)
task_train
# 测试任务对象
task_test <- as_task_surv(
  testdata, 
  time = "rfstime",
  event = "status", 
  type = "right"
)
task_test

###################################################

# coxph模型
# https://mlr3proba.mlr-org.com/reference/mlr_learners_surv.coxph.html

# 模型设定
learner_coxph <- lrn("surv.coxph")
learner_coxph$id <- "coxph"
learner_coxph

# 模型训练
set.seed(42)
learner_coxph$train(task_train)
learner_coxph

# 模型概况
learner_coxph$model
summary(learner_coxph$model)
# 森林图
survminer::ggforest(
  learner_coxph$model, 
  data = task_train$data()
)


###################################################

# 预测训练集
predtrain_coxph <- learner_coxph$predict(task_train)
predprobtrain_coxph <- predprob(
  pred = predtrain_coxph, 
  preddata = traindata, 
  etime = "rfstime",
  estatus = "status",
  model = "coxph", 
  dataset = "train", 
  timepoints =itps
)

# 性能指标
predtrain_coxph$score(measure_sa)
cindex_bootci(learner_coxph, traindata)

evaltrain_coxph <- eval4sa(
  predprob = predprobtrain_coxph,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "coxph",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltrain_coxph$auc
evaltrain_coxph$rocplot
evaltrain_coxph$brierscore
evaltrain_coxph$brierscoretest
evaltrain_coxph$calibrationplot
evaltrain_coxph$riskplot

sadca(
  predprob = predprobtrain_coxph,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "coxph",
  dataset = "train",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 预测测试集
predtest_coxph <- learner_coxph$predict(task_test)
predprobtest_coxph <- predprob(
  pred = predtest_coxph, 
  preddata = testdata, 
  etime = "rfstime",
  estatus = "status",
  model = "coxph", 
  dataset = "test", 
  timepoints =itps
)
# 性能指标
predtest_coxph$score(measure_sa)
cindex_bootci(learner_coxph, testdata)

evaltest_coxph <- eval4sa(
  predprob = predprobtest_coxph,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "coxph",
  dataset = "test",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltest_coxph$auc
evaltest_coxph$rocplot
evaltest_coxph$brierscore
evaltest_coxph$brierscoretest
evaltest_coxph$calibrationplot
evaltest_coxph$riskplot

sadca(
  predprob = predprobtest_coxph,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "coxph",
  dataset = "test",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 保存结果用于比较
save(predprobtrain_coxph,
     evaltrain_coxph,
     predprobtest_coxph,
     evaltest_coxph,
     file = "F:/Mach_learn_data/mlr_model/coxph.RData")

# 保存结果用于shiny网页计算器
traindata4sadata <- traindata
learner_coxph4sadata <- learner_coxph
save(traindata4sadata,
     learner_coxph4sadata,
     file = "F:/Mach_learn_data/shiny/coxph.RData")

#############################################

# 外部验证
# 读取数据
newdata <- readr::read_csv("F:/Mach_learn_data/new_sadata.csv")
# 分类变量转factor
for(i in c(3, 5, 9)){
  newdata[[i]] <- factor(newdata[[i]])
}
# 剔除变量-无关变量
newdata$pid <- NULL
# 剔除样本-含有缺失值的样本，时间非正的样本
# newdata <- na.omit(newdata) # 剔除任意变量下有缺失值的样本
newdata <- newdata %>%
  # na.omit() %>%
  # drop_na(age) %>%  # 剔除指定变量下有缺失值的样本
  filter(rfstime > 0) # 剔除时间非正的样本
# 数据概况
skimr::skim(newdata)

# 预测
prednew_coxph <- learner_coxph$predict_newdata(newdata)
predprobnew_coxph <- predprob(
  pred = prednew_coxph, 
  preddata = newdata, 
  etime = "rfstime",
  estatus = "status",
  model = "coxph", 
  dataset = "new", 
  timepoints =itps
)

# 性能指标
prednew_coxph$score(measure_sa)
cindex_bootci(learner_coxph, newdata)
evalnew_coxph <- eval4sa(
  predprob = predprobnew_coxph,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "coxph",
  dataset = "new",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evalnew_coxph$auc
evalnew_coxph$rocplot
evalnew_coxph$brierscore
evalnew_coxph$brierscoretest
evalnew_coxph$calibrationplot
evalnew_coxph$riskplot
sadca(
  predprob = predprobnew_coxph,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "coxph",
  dataset = "new",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

#############################################

# 模型解释
# 原始训练样本中自变量和因变量
traindatax <- traindata[, task_train$feature_names]
catvars <- 
  colnames(traindatax)[sapply(traindatax, is.factor)]
convars <- setdiff(colnames(traindatax), catvars)
traindatay <- survival::Surv(
  time = traindata$rfstime, 
  event = traindata$status
)

# 解释器——基于训练集，可以不指定时间点
exper_coxph <- survex::explain(
  learner_coxph$model, 
  data = traindatax,
  y = traindatay,
  times = itps
)

# 变量重要性
viplot(exper_coxph, output_type = "risk")
viplot(exper_coxph, output_type = "survival")

# 偏依赖图
pdplot(exper_coxph, convars, output_type = "survival")
pdplot(exper_coxph, convars, output_type = "chf")
pdplot(exper_coxph, convars, output_type = "risk")
pdplot(exper_coxph, catvars, output_type = "survival")
pdplot(exper_coxph, catvars, output_type = "chf")
pdplot(exper_coxph, catvars, output_type = "risk")


# 单样本预测分解
shap4one(exper_coxph, traindatax[1,], output_type = "survival")
shap4one(exper_coxph, traindatax[1,], output_type = "risk")

# 全局shap
shap4all_coxph <- survex::model_survshap(
  exper_coxph,
  new_observation = traindatax,
  y_true = traindatay,
  N = 100,
  calculation_method = "kernelshap",
  aggregation_method = "integral",
  output_type = "risk"
)
plot(shap4all_coxph)
plot(shap4all_coxph, geom = "beeswarm")
plot(shap4all_coxph, geom = "profile", variable = "nodes")
plot(shap4all_coxph, geom = "profile", variable = "grade")

# (抽样)汇总计算shap
sumshap_coxph <- sumshap(
  exper_coxph, 
  traindatax, 
  catvars, 
  convars, 
  sampleN=200
)
sumshap_coxph$shapvipplot
sumshap_coxph$shapplotd
sumshap_coxph$shapplotc

