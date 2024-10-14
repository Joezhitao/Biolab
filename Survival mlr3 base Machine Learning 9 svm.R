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

# svm模型
# https://mlr3extralearners.mlr-org.com/reference/mlr_learners_surv.svm.html

# 设定模型
learner_svm <- as_learner(
  po("encode", method = "treatment") %>>%
    ppl(
      "distrcompositor",
      learner = lrn(
        "surv.svm",
        type = "hybrid",
        diff.meth = "makediff3"
      ),
      estimator = "kaplan",
      form = "ph"
    )
)
learner_svm$id <- "svm"
learner_svm

# 超参数寻优空间
learner_svm$param_set$set_values(
  surv.svm.gamma.mu = c(0.5, 0.5)
)
search_space_svm <- ps(
  mu1 = p_dbl(lower = 0, upper = 1),
  mu2 = p_dbl(lower = 0, upper = 1)
)
search_space_svm$extra_trafo <- function(x, param_set) {
  x$surv.svm.gamma.mu = c(x$mu1, x$mu2)
  x$mu1 = NULL
  x$mu2 = NULL
  return(x)
}
search_space_svm

# 多核并行
future::plan("multisession")

# 超参数调优设定
set.seed(42)
tune_svm <- tune(
  tuner = tnr("random_search"),
  task = task_train,
  search_space = search_space_svm,
  learner = learner_svm,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("evals", n_evals = 10)
)
# tune_svm
tune_svm$archive$data %>%
  as.data.frame() %>%
  select(1:3) %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~surv.cindex,
                colorscale = 'Jet',
                showscale = T),
    dimensions = list(
      list(label = 'mu1', values = ~mu1),
      list(label = 'mu2', values = ~mu2)
    )
  ) %>%
  plotly::layout(title = "SVM HPO Guided by C-Index")


# 训练
learner_svm$param_set$values <-
  tune_svm$result_learner_param_vals
set.seed(42)
learner_svm$train(task_train)

# 模型概况
learner_svm

###################################################

# 预测训练集
predtrain_svm <- learner_svm$predict(task_train)
predprobtrain_svm <- predprob(
  pred = predtrain_svm, 
  preddata = traindata, 
  etime = "rfstime",
  estatus = "status",
  model = "svm", 
  dataset = "train", 
  timepoints =itps
)

# 性能指标
predtrain_svm$score(measure_sa)
cindex_bootci(learner_svm, traindata)

evaltrain_svm <- eval4sa(
  predprob = predprobtrain_svm,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "svm",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "nne",  # quantile
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltrain_svm$auc
evaltrain_svm$rocplot
evaltrain_svm$brierscore
evaltrain_svm$brierscoretest
evaltrain_svm$calibrationplot
evaltrain_svm$riskplot

sadca(
  predprob = predprobtrain_svm,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "svm",
  dataset = "train",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 预测测试集
predtest_svm <- learner_svm$predict(task_test)
predprobtest_svm <- predprob(
  pred = predtest_svm, 
  preddata = testdata, 
  etime = "rfstime",
  estatus = "status",
  model = "svm", 
  dataset = "test", 
  timepoints =itps
)
# 性能指标
predtest_svm$score(measure_sa)
cindex_bootci(learner_svm, testdata)

evaltest_svm <- eval4sa(
  predprob = predprobtest_svm,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "svm",
  dataset = "test",
  timepoints = itps,
  plotcalimethod = "nne",  # quantile
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltest_svm$auc
evaltest_svm$rocplot
evaltest_svm$brierscore
evaltest_svm$brierscoretest
evaltest_svm$calibrationplot
evaltest_svm$riskplot

sadca(
  predprob = predprobtest_svm,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "svm",
  dataset = "test",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 保存结果用于比较
save(predprobtrain_svm,
     evaltrain_svm,
     predprobtest_svm,
     evaltest_svm,
     file = "F:/Mach_learn_data/mlr_model/svm.RData")

# 保存结果用于shiny网页计算器
traindata4sadata <- traindata
learner_svm4sadata <- learner_svm
save(traindata4sadata,
     learner_svm4sadata,
     file = "F:/Mach_learn_data/shiny/svm.RData")

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
prednew_svm <- learner_svm$predict_newdata(newdata)
predprobnew_svm <- predprob(
  pred = prednew_svm, 
  preddata = newdata, 
  etime = "rfstime",
  estatus = "status",
  model = "svm", 
  dataset = "new", 
  timepoints =itps
)

# 性能指标
prednew_svm$score(measure_sa)
cindex_bootci(learner_svm, newdata)
evalnew_svm <- eval4sa(
  predprob = predprobnew_svm,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "svm",
  dataset = "new",
  timepoints = itps,
  plotcalimethod = "nne",  # quantile
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evalnew_svm$auc
evalnew_svm$rocplot
evalnew_svm$brierscore
evalnew_svm$brierscoretest
evalnew_svm$calibrationplot
evalnew_svm$riskplot
sadca(
  predprob = predprobnew_svm,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "svm",
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
exper_svm <- survex::explain_survival(
  learner_svm, 
  data = traindatax,
  y = traindatay,
  predict_function = risk_pred4,
  predict_survival_function = surv_pred4,
  predict_cumulative_hazard_function = chf_pred4,
  label = "svm",
  times = itps
)

# 变量重要性
viplot(exper_svm, output_type = "risk")
viplot(exper_svm, output_type = "survival")

# 偏依赖图
pdplot(exper_svm, convars, output_type = "survival")
pdplot(exper_svm, convars, output_type = "chf")
pdplot(exper_svm, convars, output_type = "risk")
pdplot(exper_svm, catvars, output_type = "survival")
pdplot(exper_svm, catvars, output_type = "chf")
pdplot(exper_svm, catvars, output_type = "risk")


# 单样本预测分解
shap4one(exper_svm, traindatax[1,], output_type = "survival")
shap4one(exper_svm, traindatax[1,], output_type = "risk")

# 全局shap
shap4all_svm <- survex::model_survshap(
  exper_svm,
  new_observation = traindatax,
  y_true = traindatay,
  N = 100,
  calculation_method = "kernelshap",
  aggregation_method = "integral",
  output_type = "risk"
)
plot(shap4all_svm)
plot(shap4all_svm, geom = "beeswarm")
plot(shap4all_svm, geom = "profile", variable = "nodes")
plot(shap4all_svm, geom = "profile", variable = "grade")

# (抽样)汇总计算shap
sumshap_svm <- sumshap(
  exper_svm, 
  traindatax, 
  catvars, 
  convars, 
  sampleN=200
)
sumshap_svm$shapvipplot
sumshap_svm$shapplotd
sumshap_svm$shapplotc
