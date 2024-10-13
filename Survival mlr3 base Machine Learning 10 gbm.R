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

# gbm模型
# https://mlr3extralearners.mlr-org.com/reference/mlr_learners_surv.gbm.html

# 设定模型
learner_gbm <- as_learner(
  po("encode", method = "treatment") %>>%
    ppl(
      "distrcompositor",
      learner = lrn(
        "surv.gbm",
        n.trees = to_tune(100, 500),
        interaction.depth = to_tune(1, 5),
        n.minobsinnode = to_tune(5, 21),
        shrinkage = to_tune(0.001, 0.1)
      ),
      estimator = "kaplan",
      form = "ph"
    )
)
learner_gbm$id <- "gbm" 
learner_gbm

# 多核并行
future::plan("multisession")

# 超参数调优设定
set.seed(42)
tune_gbm <- tune(
  tuner = tnr("grid_search", resolution = 3),
  task = task_train,
  learner = learner_gbm,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)
tune_gbm
as.data.table(tune_gbm$archive) %>%
  as.data.frame() %>%
  select(1:5) %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~surv.cindex, 
                colorscale = 'Jet', 
                showscale = T),
    dimensions = list(
      list(label = 'ntree', values = ~surv.gbm.n.trees),
      list(label = 'interaction.depth', values = ~surv.gbm.interaction.depth),
      list(label = 'n.minobsinnode', values = ~surv.gbm.n.minobsinnode),
      list(label = 'shrinkage', values = ~surv.gbm.shrinkage)
    )
  ) %>%
  plotly::layout(title = "GBM HPO Guided by C-Index",
                 font = list(family = "serif"))

# 训练
learner_gbm$param_set$values <-
  tune_gbm$result_learner_param_vals
set.seed(42)
learner_gbm$train(task_train)

# 模型概况
learner_gbm
learner_gbm$model$surv.gbm$model


###################################################

# 预测训练集
predtrain_gbm <- learner_gbm$predict(task_train)
predprobtrain_gbm <- predprob(
  pred = predtrain_gbm, 
  preddata = traindata, 
  etime = "rfstime",
  estatus = "status",
  model = "gbm", 
  dataset = "train", 
  timepoints =itps
)

# 性能指标
predtrain_gbm$score(measure_sa)
cindex_bootci(learner_gbm, traindata)

evaltrain_gbm <- eval4sa(
  predprob = predprobtrain_gbm,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "gbm",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltrain_gbm$auc
evaltrain_gbm$rocplot
evaltrain_gbm$brierscore
evaltrain_gbm$brierscoretest
evaltrain_gbm$calibrationplot
evaltrain_gbm$riskplot

sadca(
  predprob = predprobtrain_gbm,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "gbm",
  dataset = "train",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 预测测试集
predtest_gbm <- learner_gbm$predict(task_test)
predprobtest_gbm <- predprob(
  pred = predtest_gbm, 
  preddata = testdata, 
  etime = "rfstime",
  estatus = "status",
  model = "gbm", 
  dataset = "test", 
  timepoints =itps
)
# 性能指标
predtest_gbm$score(measure_sa)
cindex_bootci(learner_gbm, testdata)

evaltest_gbm <- eval4sa(
  predprob = predprobtest_gbm,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "gbm",
  dataset = "test",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltest_gbm$auc
evaltest_gbm$rocplot
evaltest_gbm$brierscore
evaltest_gbm$brierscoretest
evaltest_gbm$calibrationplot
evaltest_gbm$riskplot

sadca(
  predprob = predprobtest_gbm,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "gbm",
  dataset = "test",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 保存结果用于比较
save(predprobtrain_gbm,
     evaltrain_gbm,
     predprobtest_gbm,
     evaltest_gbm,
     file = "F:/Mach_learn_data/gbm.RData")

# 保存结果用于shiny网页计算器
traindata4sadata <- traindata
learner_gbm4sadata <- learner_gbm
save(traindata4sadata,
     learner_gbm4sadata,
     file = "F:/Mach_learn_data/mlr_model/shiny_gbm.RData")

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
prednew_gbm <- learner_gbm$predict_newdata(newdata)
predprobnew_gbm <- predprob(
  pred = prednew_gbm, 
  preddata = newdata, 
  etime = "rfstime",
  estatus = "status",
  model = "gbm", 
  dataset = "new", 
  timepoints =itps
)

# 性能指标
prednew_gbm$score(measure_sa)
cindex_bootci(learner_gbm, newdata)
evalnew_gbm <- eval4sa(
  predprob = predprobnew_gbm,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "gbm",
  dataset = "new",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evalnew_gbm$auc
evalnew_gbm$rocplot
evalnew_gbm$brierscore
evalnew_gbm$brierscoretest
evalnew_gbm$calibrationplot
evalnew_gbm$riskplot
sadca(
  predprob = predprobnew_gbm,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "gbm",
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
exper_gbm <- survex::explain_survival(
  learner_gbm, 
  data = traindatax,
  y = traindatay,
  predict_function = risk_pred,
  predict_survival_function = surv_pred,
  predict_cumulative_hazard_function = chf_pred,
  label = "gbm",
  times = itps
)

# 变量重要性
viplot(exper_gbm, output_type = "risk")
viplot(exper_gbm, output_type = "survival")

# 偏依赖图
pdplot(exper_gbm, convars, output_type = "survival")
pdplot(exper_gbm, convars, output_type = "chf")
pdplot(exper_gbm, convars, output_type = "risk")
pdplot(exper_gbm, catvars, output_type = "survival")
pdplot(exper_gbm, catvars, output_type = "chf")
pdplot(exper_gbm, catvars, output_type = "risk")


# 单样本预测分解
shap4one(exper_gbm, traindatax[1,], output_type = "survival")
shap4one(exper_gbm, traindatax[1,], output_type = "risk")

# 全局shap
shap4all_gbm <- survex::model_survshap(
  exper_gbm,
  new_observation = traindatax,
  y_true = traindatay,
  N = 100,
  calculation_method = "kernelshap",
  aggregation_method = "integral",
  output_type = "risk"
)
plot(shap4all_gbm)
plot(shap4all_gbm, geom = "beeswarm")
plot(shap4all_gbm, geom = "profile", variable = "nodes")
plot(shap4all_gbm, geom = "profile", variable = "grade")

# (抽样)汇总计算shap
sumshap_gbm <- sumshap(
  exper_gbm, 
  traindatax, 
  catvars, 
  convars, 
  sampleN=200
)
sumshap_gbm$shapvipplot
sumshap_gbm$shapplotd
sumshap_gbm$shapplotc
