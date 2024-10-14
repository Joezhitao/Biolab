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

# xgboost模型
# https://mlr3extralearners.mlr-org.com/reference/mlr_learners_surv.xgboost.cox.html

# 设定模型
learner_xgboost <- as_learner(
  po("encode", method = "treatment") %>>%
    lrn(
      "surv.xgboost.cox",
      nrounds = to_tune(100, 500), 
      max_depth = to_tune(1, 5),
      eta = to_tune(1e-4, 1)
    )
)
learner_xgboost$id <- "xgboost"
learner_xgboost

# 多核并行
future::plan("multisession")

# 超参数调优设定
set.seed(42)
tune_xgboost <- tune(
  tuner = tnr("random_search"),
  task = task_train, 
  learner = learner_xgboost,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("evals", n_evals = 20)
)
tune_xgboost
tune_xgboost$archive$data %>%
  as.data.frame() %>%
  select(1:4) %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~surv.cindex, 
                colorscale = 'Jet', 
                showscale = T),
    dimensions = list(
      list(label = 'nrounds', values = ~surv.xgboost.cox.nrounds),
      list(label = 'max_depth', values = ~surv.xgboost.cox.max_depth),
      list(label = 'eta', values = ~surv.xgboost.cox.eta)
    )
  ) %>%
  plotly::layout(title = "xgboost HPO Guided by C-Index",
                 font = list(family = "serif"))

# 训练
learner_xgboost$param_set$values <-
  tune_xgboost$result_learner_param_vals
set.seed(42)
learner_xgboost$train(task_train)

# 模型概况
learner_xgboost

###################################################

# 预测训练集
predtrain_xgboost <- learner_xgboost$predict(task_train)
predprobtrain_xgboost <- predprob(
  pred = predtrain_xgboost, 
  preddata = traindata, 
  etime = "rfstime",
  estatus = "status",
  model = "xgboost", 
  dataset = "train", 
  timepoints =itps
)

# 性能指标
predtrain_xgboost$score(measure_sa)
cindex_bootci(learner_xgboost, traindata)

evaltrain_xgboost <- eval4sa(
  predprob = predprobtrain_xgboost,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "xgboost",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltrain_xgboost$auc
evaltrain_xgboost$rocplot
evaltrain_xgboost$brierscore
evaltrain_xgboost$brierscoretest
evaltrain_xgboost$calibrationplot
evaltrain_xgboost$riskplot

sadca(
  predprob = predprobtrain_xgboost,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "xgboost",
  dataset = "train",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 预测测试集
predtest_xgboost <- learner_xgboost$predict(task_test)
predprobtest_xgboost <- predprob(
  pred = predtest_xgboost, 
  preddata = testdata, 
  etime = "rfstime",
  estatus = "status",
  model = "xgboost", 
  dataset = "test", 
  timepoints =itps
)
# 性能指标
predtest_xgboost$score(measure_sa)
cindex_bootci(learner_xgboost, testdata)

evaltest_xgboost <- eval4sa(
  predprob = predprobtest_xgboost,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "xgboost",
  dataset = "test",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltest_xgboost$auc
evaltest_xgboost$rocplot
evaltest_xgboost$brierscore
evaltest_xgboost$brierscoretest
evaltest_xgboost$calibrationplot
evaltest_xgboost$riskplot

sadca(
  predprob = predprobtest_xgboost,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "xgboost",
  dataset = "test",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 保存结果用于比较
save(predprobtrain_xgboost,
     evaltrain_xgboost,
     predprobtest_xgboost,
     evaltest_xgboost,
     file = "F:/Mach_learn_data/mlr_model/xgboost.RData")

# 保存结果用于shiny网页计算器
traindata4sadata <- traindata
learner_xgboost4sadata <- learner_xgboost
save(traindata4sadata,
     learner_xgboost4sadata,
     file = "F:/Mach_learn_data/shiny/xgboost.RData")

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
prednew_xgboost <- learner_xgboost$predict_newdata(newdata)
predprobnew_xgboost <- predprob(
  pred = prednew_xgboost, 
  preddata = newdata, 
  etime = "rfstime",
  estatus = "status",
  model = "xgboost", 
  dataset = "new", 
  timepoints =itps
)

# 性能指标
prednew_xgboost$score(measure_sa)
cindex_bootci(learner_xgboost, newdata)
evalnew_xgboost <- eval4sa(
  predprob = predprobnew_xgboost,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "xgboost",
  dataset = "new",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evalnew_xgboost$auc
evalnew_xgboost$rocplot
evalnew_xgboost$brierscore
evalnew_xgboost$brierscoretest
evalnew_xgboost$calibrationplot
evalnew_xgboost$riskplot
sadca(
  predprob = predprobnew_xgboost,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "xgboost",
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
exper_xgboost <- survex::explain_survival(
  learner_xgboost, 
  data = traindatax,
  y = traindatay,
  predict_function = risk_pred,
  predict_survival_function = surv_pred,
  predict_cumulative_hazard_function = chf_pred,
  label = "xgboost",
  times = itps
)

# 变量重要性
viplot(exper_xgboost, output_type = "risk")
viplot(exper_xgboost, output_type = "survival")

# 偏依赖图
pdplot(exper_xgboost, convars, output_type = "survival")
pdplot(exper_xgboost, convars, output_type = "chf")
pdplot(exper_xgboost, convars, output_type = "risk")
pdplot(exper_xgboost, catvars, output_type = "survival")
pdplot(exper_xgboost, catvars, output_type = "chf")
pdplot(exper_xgboost, catvars, output_type = "risk")


# 单样本预测分解
shap4one(exper_xgboost, traindatax[1,], output_type = "survival")
shap4one(exper_xgboost, traindatax[1,], output_type = "risk")

# 全局shap
shap4all_xgboost <- survex::model_survshap(
  exper_xgboost,
  new_observation = traindatax,
  y_true = traindatay,
  N = 100,
  calculation_method = "kernelshap",
  aggregation_method = "integral",
  output_type = "risk"
)
plot(shap4all_xgboost)
plot(shap4all_xgboost, geom = "beeswarm")
plot(shap4all_xgboost, geom = "profile", variable = "nodes")
plot(shap4all_xgboost, geom = "profile", variable = "grade")

# (抽样)汇总计算shap
sumshap_xgboost <- sumshap(
  exper_xgboost, 
  traindatax, 
  catvars, 
  convars, 
  sampleN=200
)
sumshap_xgboost$shapvipplot
sumshap_xgboost$shapplotd
sumshap_xgboost$shapplotc
