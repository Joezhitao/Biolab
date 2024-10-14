##Joe_2024_10_12
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

# enet模型
# https://mlr3extralearners.mlr-org.com/reference/mlr_learners_surv.cv_glmnet.html

# 模型设定
learner_enet <- as_learner(
  po("encode", method = "treatment") %>>%
    lrn(
      "surv.cv_glmnet",
      type.measure = "C",
      alpha = to_tune(0.05, 0.95)
    )
)
learner_enet$id <- "enet"
learner_enet

# 多核并行
future::plan("multisession")

# 超参数调优设定
set.seed(42)
tune_enet <- tune(
  tuner = tnr(
    "grid_search", 
    resolution = 10
  ),
  task = task_train,
  learner = learner_enet,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)
tune_enet

# 训练最终模型
learner_enet$param_set$values <-
  tune_enet$result_learner_param_vals
set.seed(42)
learner_enet$train(task_train)
learner_enet

# 模型概况
learner_enet$base_learner()$model$model
coef(learner_enet$base_learner()$model$model)
plot(learner_enet$base_learner()$model$model)

###################################################

# 预测训练集
predtrain_enet <- learner_enet$predict(task_train)
predprobtrain_enet <- predprob(
  pred = predtrain_enet, 
  preddata = traindata, 
  etime = "rfstime",
  estatus = "status",
  model = "enet", 
  dataset = "train", 
  timepoints =itps
)

# 性能指标
predtrain_enet$score(measure_sa)
cindex_bootci(learner_enet, traindata)

evaltrain_enet <- eval4sa(
  predprob = predprobtrain_enet,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "enet",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltrain_enet$auc
evaltrain_enet$rocplot
evaltrain_enet$brierscore
evaltrain_enet$brierscoretest
evaltrain_enet$calibrationplot
evaltrain_enet$riskplot

sadca(
  predprob = predprobtrain_enet,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "enet",
  dataset = "train",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 预测测试集
predtest_enet <- learner_enet$predict(task_test)
predprobtest_enet <- predprob(
  pred = predtest_enet, 
  preddata = testdata, 
  etime = "rfstime",
  estatus = "status",
  model = "enet", 
  dataset = "test", 
  timepoints =itps
)
# 性能指标
predtest_enet$score(measure_sa)
cindex_bootci(learner_enet, testdata)

evaltest_enet <- eval4sa(
  predprob = predprobtest_enet,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "enet",
  dataset = "test",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltest_enet$auc
evaltest_enet$rocplot
evaltest_enet$brierscore
evaltest_enet$brierscoretest
evaltest_enet$calibrationplot
evaltest_enet$riskplot

sadca(
  predprob = predprobtest_enet,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "enet",
  dataset = "test",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 保存结果用于比较
save(predprobtrain_enet,
     evaltrain_enet,
     predprobtest_enet,
     evaltest_enet,
     file = "F:/Mach_learn_data/mlr_model/enet.RData")

# 保存结果用于shiny网页计算器
traindata4sadata <- traindata
learner_enet4sadata <- learner_enet
save(traindata4sadata,
     learner_enet4sadata,
     file = "F:/Mach_learn_data/shiny/enet.RData")

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
prednew_enet <- learner_enet$predict_newdata(newdata)
predprobnew_enet <- predprob(
  pred = prednew_enet, 
  preddata = newdata, 
  etime = "rfstime",
  estatus = "status",
  model = "enet", 
  dataset = "new", 
  timepoints =itps
)

# 性能指标
prednew_enet$score(measure_sa)
cindex_bootci(learner_enet, newdata)
evalnew_enet <- eval4sa(
  predprob = predprobnew_enet,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "enet",
  dataset = "new",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evalnew_enet$auc
evalnew_enet$rocplot
evalnew_enet$brierscore
evalnew_enet$brierscoretest
evalnew_enet$calibrationplot
evalnew_enet$riskplot
sadca(
  predprob = predprobnew_enet,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "enet",
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
exper_enet <- survex::explain_survival(
  learner_enet, 
  data = traindatax,
  y = traindatay,
  predict_function = risk_pred,
  predict_survival_function = surv_pred,
  predict_cumulative_hazard_function = chf_pred,
  label = "enet",
  times = itps
)

# 变量重要性
viplot(exper_enet, output_type = "risk")
viplot(exper_enet, output_type = "survival")

# 偏依赖图
pdplot(exper_enet, convars, output_type = "survival")
pdplot(exper_enet, convars, output_type = "chf")
pdplot(exper_enet, convars, output_type = "risk")
pdplot(exper_enet, catvars, output_type = "survival")
pdplot(exper_enet, catvars, output_type = "chf")
pdplot(exper_enet, catvars, output_type = "risk")


# 单样本预测分解
shap4one(exper_enet, traindatax[1,], output_type = "survival")
shap4one(exper_enet, traindatax[1,], output_type = "risk")

# 全局shap
shap4all_enet <- survex::model_survshap(
  exper_enet,
  new_observation = traindatax,
  y_true = traindatay,
  N = 100,
  calculation_method = "kernelshap",
  aggregation_method = "integral",
  output_type = "risk"
)
plot(shap4all_enet)
plot(shap4all_enet, geom = "beeswarm")
plot(shap4all_enet, geom = "profile", variable = "nodes")
plot(shap4all_enet, geom = "profile", variable = "grade")

# (抽样)汇总计算shap
sumshap_enet <- sumshap(
  exper_enet, 
  traindatax, 
  catvars, 
  convars, 
  sampleN=200
)
sumshap_enet$shapvipplot
sumshap_enet$shapplotd
sumshap_enet$shapplotc
