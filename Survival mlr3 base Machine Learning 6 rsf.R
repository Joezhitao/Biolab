##Joe_2024_10_13
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

# 随机森林
# https://mlr3extralearners.mlr-org.com/reference/mlr_learners_surv.rfsrc.html

# 设定模型
learner_rsf <- lrn(
  "surv.rfsrc",
  ntree = to_tune(200, 500),
  mtry = to_tune(3, 5),
  nodesize = to_tune(15, 21)
)
learner_rsf$id <- "rsf"
learner_rsf

# 多核并行
future::plan("multisession")

# 超参数调优设定
set.seed(42)
tune_rsf <- tune(
  tuner = tnr(
    "grid_search", 
    resolution = 3
  ),
  task = task_train,
  learner = learner_rsf,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)
tune_rsf
as.data.table(tune_rsf$archive) %>%
  as.data.frame() %>%
  select(1:4) %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~surv.cindex, 
                colorscale = 'Jet', 
                showscale = T),
    dimensions = list(
      list(label = 'ntree', values = ~ntree),
      list(label = 'mtry', values = ~mtry),
      list(label = 'nodesize', values = ~nodesize)
    )
  ) %>%
  plotly::layout(title = "RSF HPO Guided by C-Index",
                 font = list(family = "serif"))


# 训练最终模型
learner_rsf$param_set$values <-
  tune_rsf$result_learner_param_vals
set.seed(42)
learner_rsf$train(task_train)
learner_rsf

# 模型概况
learner_rsf
learner_rsf$model

###################################################

# 预测训练集
predtrain_rsf <- learner_rsf$predict(task_train)
predprobtrain_rsf <- predprob(
  pred = predtrain_rsf, 
  preddata = traindata, 
  etime = "rfstime",
  estatus = "status",
  model = "rsf", 
  dataset = "train", 
  timepoints =itps
)

# 性能指标
predtrain_rsf$score(measure_sa)
cindex_bootci(learner_rsf, traindata)

evaltrain_rsf <- eval4sa(
  predprob = predprobtrain_rsf,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "rsf",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltrain_rsf$auc
evaltrain_rsf$rocplot
evaltrain_rsf$brierscore
evaltrain_rsf$brierscoretest
evaltrain_rsf$calibrationplot
evaltrain_rsf$riskplot

sadca(
  predprob = predprobtrain_rsf,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "rsf",
  dataset = "train",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 预测测试集
predtest_rsf <- learner_rsf$predict(task_test)
predprobtest_rsf <- predprob(
  pred = predtest_rsf, 
  preddata = testdata, 
  etime = "rfstime",
  estatus = "status",
  model = "rsf", 
  dataset = "test", 
  timepoints =itps
)
# 性能指标
predtest_rsf$score(measure_sa)
cindex_bootci(learner_rsf, testdata)

evaltest_rsf <- eval4sa(
  predprob = predprobtest_rsf,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "rsf",
  dataset = "test",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltest_rsf$auc
evaltest_rsf$rocplot
evaltest_rsf$brierscore
evaltest_rsf$brierscoretest
evaltest_rsf$calibrationplot
evaltest_rsf$riskplot

sadca(
  predprob = predprobtest_rsf,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "rsf",
  dataset = "test",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 保存结果用于比较
save(predprobtrain_rsf,
     evaltrain_rsf,
     predprobtest_rsf,
     evaltest_rsf,
     file = "F:/Mach_learn_data/mlr_model/rsf.RData")

# 保存结果用于shiny网页计算器
traindata4sadata <- traindata
learner_rsf4sadata <- learner_rsf
save(traindata4sadata,
     learner_rsf4sadata,
     file = "F:/Mach_learn_data/shiny/rsf.RData")

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
prednew_rsf <- learner_rsf$predict_newdata(newdata)
predprobnew_rsf <- predprob(
  pred = prednew_rsf, 
  preddata = newdata, 
  etime = "rfstime",
  estatus = "status",
  model = "rsf", 
  dataset = "new", 
  timepoints =itps
)

# 性能指标
prednew_rsf$score(measure_sa)
cindex_bootci(learner_rsf, newdata)
evalnew_rsf <- eval4sa(
  predprob = predprobnew_rsf,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "rsf",
  dataset = "new",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evalnew_rsf$auc
evalnew_rsf$rocplot
evalnew_rsf$brierscore
evalnew_rsf$brierscoretest
evalnew_rsf$calibrationplot
evalnew_rsf$riskplot
sadca(
  predprob = predprobnew_rsf,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "rsf",
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
exper_rsf <- survex::explain(
  learner_rsf$model, 
  data = traindatax,
  y = traindatay,
  times = itps
)

# 变量重要性
viplot(exper_rsf, output_type = "risk")
viplot(exper_rsf, output_type = "survival")

# 偏依赖图
pdplot(exper_rsf, convars, output_type = "survival")
pdplot(exper_rsf, convars, output_type = "chf")
pdplot(exper_rsf, convars, output_type = "risk")
pdplot(exper_rsf, catvars, output_type = "survival")
pdplot(exper_rsf, catvars, output_type = "chf")
pdplot(exper_rsf, catvars, output_type = "risk")


# 单样本预测分解
shap4one(exper_rsf, traindatax[1,], output_type = "survival")
shap4one(exper_rsf, traindatax[1,], output_type = "risk")

# 全局shap
shap4all_rsf <- survex::model_survshap(
  exper_rsf,
  new_observation = traindatax,
  y_true = traindatay,
  N = 100,
  calculation_method = "kernelshap",
  aggregation_method = "integral",
  output_type = "risk"
)
plot(shap4all_rsf)
plot(shap4all_rsf, geom = "beeswarm")
plot(shap4all_rsf, geom = "profile", variable = "nodes")
plot(shap4all_rsf, geom = "profile", variable = "grade")

# (抽样)汇总计算shap
sumshap_rsf <- sumshap(
  exper_rsf, 
  traindatax, 
  catvars, 
  convars, 
  sampleN=200
)
sumshap_rsf$shapvipplot
sumshap_rsf$shapplotd
sumshap_rsf$shapplotc
