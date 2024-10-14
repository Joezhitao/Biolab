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

# lasso模型
learner_lasso <- as_learner(
  po("encode", method = "treatment") %>>%
    lrn(
      "surv.cv_glmnet",
      type.measure = "C",
      alpha = 1
    )
)
learner_lasso

########################

# 决策树
learner_ctree <- auto_tuner(
  tuner = tnr(
    "grid_search", 
    resolution = 4,
    batch_size = 4
  ),
  learner = lrn(
    "surv.ctree",
    alpha = to_tune(c(0.01, 0.05, 0.1)),
    minbucket = to_tune(5, 25)
  ),
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)
learner_ctree

########################

# 随机森林
learner_rsf <- auto_tuner(
  tuner = tnr(
    "grid_search", 
    resolution = 3,
    batch_size = 3
  ),
  learner = lrn(
    "surv.rfsrc",
    ntree = to_tune(200, 500),
    mtry = to_tune(3, 5),
    nodesize = to_tune(15, 21)
  ),
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)
learner_rsf

###################################################

# stacking集成
# 设定
learner_stack <- as_learner(
  gunion(
    list(
      po("learner_cv", learner_ctree),
      po("learner_cv", learner_rsf),
      po("nop")
    )
  ) %>>%
    po("featureunion") %>>%
    po(
      "select", 
      selector = selector_name(c(
        task_train$feature_names, 
        "surv.rfsrc.tuned.crank",
        "surv.ctree.tuned.crank"
      ))
    ) %>>%
    po("encode", method = "treatment") %>>%
    learner_lasso
)
learner_stack

# 训练
set.seed(42)
learner_stack$train(task_train)
learner_stack

###################################################

# 预测训练集
predtrain_stack <- learner_stack$predict(task_train)
predprobtrain_stack <- predprob(
  pred = predtrain_stack, 
  preddata = traindata, 
  etime = "rfstime",
  estatus = "status",
  model = "stack", 
  dataset = "train", 
  timepoints =itps
)

# 性能指标
predtrain_stack$score(measure_sa)
cindex_bootci(learner_stack, traindata)

evaltrain_stack <- eval4sa(
  predprob = predprobtrain_stack,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "stack",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltrain_stack$auc
evaltrain_stack$rocplot
evaltrain_stack$brierscore
evaltrain_stack$brierscoretest
evaltrain_stack$calibrationplot
evaltrain_stack$riskplot

sadca(
  predprob = predprobtrain_stack,
  preddata = traindata,
  etime = "rfstime",
  estatus = "status",
  model = "stack",
  dataset = "train",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 预测测试集
predtest_stack <- learner_stack$predict(task_test)
predprobtest_stack <- predprob(
  pred = predtest_stack, 
  preddata = testdata, 
  etime = "rfstime",
  estatus = "status",
  model = "stack", 
  dataset = "test", 
  timepoints =itps
)
# 性能指标
predtest_stack$score(measure_sa)
cindex_bootci(learner_stack, testdata)

evaltest_stack <- eval4sa(
  predprob = predprobtest_stack,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "stack",
  dataset = "test",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evaltest_stack$auc
evaltest_stack$rocplot
evaltest_stack$brierscore
evaltest_stack$brierscoretest
evaltest_stack$calibrationplot
evaltest_stack$riskplot

sadca(
  predprob = predprobtest_stack,
  preddata = testdata,
  etime = "rfstime",
  estatus = "status",
  model = "stack",
  dataset = "test",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)

# 保存结果用于比较
save(predprobtrain_stack,
     evaltrain_stack,
     predprobtest_stack,
     evaltest_stack,
     file = "F:/Mach_learn_data/mlr_model/stack.RData")

# 保存结果用于shiny网页计算器
traindata4sadata <- traindata
learner_stack4sadata <- learner_stack
save(traindata4sadata,
     learner_stack4sadata,
     file = "F:/Mach_learn_data/shiny/stack.RData")

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
prednew_stack <- learner_stack$predict_newdata(newdata)
predprobnew_stack <- predprob(
  pred = prednew_stack, 
  preddata = newdata, 
  etime = "rfstime",
  estatus = "status",
  model = "stack", 
  dataset = "new", 
  timepoints =itps
)

# 性能指标
prednew_stack$score(measure_sa)
cindex_bootci(learner_stack, newdata)

evalnew_stack <- eval4sa(
  predprob = predprobnew_stack,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "stack",
  dataset = "new",
  timepoints = itps,
  plotcalimethod = "quantile",  # nne
  bw4nne = NULL,
  q4quantile = 5,
  cutoff = "median"
)
evalnew_stack$auc
evalnew_stack$rocplot
evalnew_stack$brierscore
evalnew_stack$brierscoretest
evalnew_stack$calibrationplot
evalnew_stack$riskplot
sadca(
  predprob = predprobnew_stack,
  preddata = newdata,
  etime = "rfstime",
  estatus = "status",
  model = "stack",
  dataset = "new",
  timepoints = itps,
  timepoint = 365, 
  xrange = 0:100 / 100
)


