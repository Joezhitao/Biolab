#Joe 2024-10-14
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

# 神经网络
# https://mlr3extralearners.mlr-org.com/reference/mlr_learners_surv.deepsurv.html

# 设定神经网络模型
learner_nn <- as_learner(
  po("encode", method = "treatment") %>>%
    lrn(
      "surv.deepsurv",
      optimizer = "adam",
      frac = 0.2,
      early_stopping = TRUE,
      patience = 20,
      batch_size = 32,
      epochs = 500
    )
)
learner_nn$id <- "nn"
learner_nn

# 超参数寻优空间
search_space_nn <- ps(
  num = p_int(lower = 1, upper = 3),
  nodes = p_int(lower = 5, upper = 13),
  surv.deepsurv.learning_rate = p_dbl(lower = 0, upper = 0.1),
  surv.deepsurv.dropout = p_dbl(lower = 0, upper = 0.5),
  surv.deepsurv.weight_decay = p_dbl(lower = 0, upper = 0.5)
)
search_space_nn$extra_trafo <- function(x, param_set) {
  x$surv.deepsurv.num_nodes = rep(x$nodes, x$num)
  x$nodes = NULL
  x$num = NULL
  return(x)
}
search_space_nn

# 多核并行
future::plan("multisession")

# 超参数调优设定
set.seed(42)
tune_nn <- tune(
  tuner = tnr("random_search"),
  task = task_train,
  search_space = search_space_nn,
  learner = learner_nn,
  resampling = rsmp("cv", folds = 3),
  measure = msr("surv.cindex"),
  terminator = trm("evals", n_evals = 20)
)
# tune_nn
tune_nn$archive$data %>%
  as.data.frame() %>%
  select(1:6) %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~surv.cindex, 
                colorscale = 'Jet', 
                showscale = T),
    dimensions = list(
      list(label = 'num', values = ~num),
      list(label = 'nodes', values = ~nodes),
      list(label = 'learning_rate', values = ~surv.deepsurv.learning_rate),
      list(label = 'dropout', values = ~surv.deepsurv.dropout),
      list(label = 'weight_decay', values = ~surv.deepsurv.weight_decay)
    )
  ) %>%
  plotly::layout(title = "NN HPO Guided by C-Index",
                 font = list(family = "serif"))

# 训练
learner_nn$param_set$values <-
  tune_nn$result_learner_param_vals
set.seed(42)
learner_nn$train(task_train)
learner_nn

# 模型概况
learner_nn
learner_nn$model

###################################################

# 预测训练集
predtrain_nn <- learner_nn$predict(task_train)
predprobtrain_nn <- predprob(
  pred = predtrain_nn, 
  preddata = traindata, 
  etime = "rfstime",
  estatus = "status",
  model = "nn", 
  dataset = "train", 
  timepoints =itps
)

# 性能指标
predtrain_nn$score(measure_sa)
cindex_bootci(learner_nn, traindata)

evaltrain_nn <- eval4sa(
  predprob = predprobtrain_nn,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "nn",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltrain_nn$auc
evaltrain_nn$rocplot
evaltrain_nn$brierscore
evaltrain_nn$brierscoretest
evaltrain_nn$calibrationplot
evaltrain_nn$riskplot

sadca(
  predprob = predprobtrain_nn,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "nn",
  dataset = "train",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 预测测试集
predtest_nn <- learner_nn$predict(task_test)
predprobtest_nn <- predprob(
  pred = predtest_nn, 
  preddata = testdata, 
  etime = "rfstime",
  estatus = "status",
  model = "nn", 
  dataset = "test", 
  timepoints =itps
)
# 性能指标
predtest_nn$score(measure_sa)
cindex_bootci(learner_nn, testdata)

evaltest_nn <- eval4sa(
  predprob = predprobtest_nn,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "nn",
  dataset = "test",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltest_nn$auc
evaltest_nn$rocplot
evaltest_nn$brierscore
evaltest_nn$brierscoretest
evaltest_nn$calibrationplot
evaltest_nn$riskplot

sadca(
  predprob = predprobtest_nn,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "nn",
  dataset = "test",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 保存结果用于比较
save(predprobtrain_nn,
     evaltrain_nn,
     predprobtest_nn,
     evaltest_nn,
     file = "F:/Mach_learn_data/nn.RData")

# 保存结果用于shiny网页计算器
traindata4sadata <- traindata
learner_nn4sadata <- learner_nn
save(traindata4sadata,
     learner_nn4sadata,
     file = "F:/Mach_learn_data/mlr_model/shiny_nn.RData")

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
prednew_nn <- learner_nn$predict_newdata(newdata)
predprobnew_nn <- predprob(
  pred = prednew_nn, 
  preddata = newdata, 
  etime = "rfstime",
  estatus = "status",
  model = "nn", 
  dataset = "new", 
  timepoints =itps
)

# 性能指标
prednew_nn$score(measure_sa)
cindex_bootci(learner_nn, newdata)
evalnew_nn <- eval4sa(
  predprob = predprobnew_nn,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "nn",
  dataset = "new",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evalnew_nn$auc
evalnew_nn$rocplot
evalnew_nn$brierscore
evalnew_nn$brierscoretest
evalnew_nn$calibrationplot
evalnew_nn$riskplot
sadca(
  predprob = predprobnew_nn,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "nn",
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
exper_nn <- survex::explain_survival(
  learner_nn, 
  data = traindatax,
  y = traindatay,
  predict_function = risk_pred,
  predict_survival_function = surv_pred,
  predict_cumulative_hazard_function = chf_pred,
  label = "nn",
  times = itps
)

# 变量重要性
viplot(exper_nn, output_type = "risk")
viplot(exper_nn, output_type = "survival")

# 偏依赖图
pdplot(exper_nn, convars, output_type = "survival")
pdplot(exper_nn, convars, output_type = "chf")
pdplot(exper_nn, convars, output_type = "risk")
pdplot(exper_nn, catvars, output_type = "survival")
pdplot(exper_nn, catvars, output_type = "chf")
pdplot(exper_nn, catvars, output_type = "risk")


# 单样本预测分解
shap4one(exper_nn, traindatax[1,], output_type = "survival")
shap4one(exper_nn, traindatax[1,], output_type = "risk")

# 全局shap
shap4all_nn <- survex::model_survshap(
  exper_nn,
  new_observation = traindatax,
  y_true = traindatay,
  N = 100,
  calculation_method = "kernelshap",
  aggregation_method = "integral",
  output_type = "risk"
)
plot(shap4all_nn)
plot(shap4all_nn, geom = "beeswarm")
plot(shap4all_nn, geom = "profile", variable = "nodes")
plot(shap4all_nn, geom = "profile", variable = "grade")

# (抽样)汇总计算shap
sumshap_nn <- sumshap(
  exper_nn, 
  traindatax, 
  catvars, 
  convars, 
  sampleN=200
)
sumshap_nn$shapvipplot
sumshap_nn$shapplotd
sumshap_nn$shapplotc
