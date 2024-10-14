#Joe_2024_10_12
# 加载包
# library(devtools)
# devtools::install_github("mlr-org/mlr3proba")
# devtools::install_github("mlr-org/mlr3extralearners@*release")
rm(list = ls())
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

# coxph模型设定
learner_coxph <- lrn("surv.coxph")
learner_coxph$id <- "coxph"
learner_coxph

# lasso模型设定
learner_lasso <- as_learner(
  po("encode", method = "treatment") %>>%
    lrn(
      "surv.cv_glmnet",
      type.measure = "C",
      alpha = 1
    )
)
learner_lasso$id <- "lasso"
learner_lasso

# 随机森林模型
learner_rsf_set <- lrn(
  "surv.rfsrc",
  ntree = to_tune(200, 500),
  mtry = to_tune(3, 5),
  nodesize = to_tune(15, 21)
)

# 超参数调优设定
learner_rsf <- auto_tuner(
  tuner = tnr("grid_search", resolution = 3, batch_size = 3),
  learner = learner_rsf_set,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)
learner_rsf$id <- "rsf"
learner_rsf

# 比较设定
bm_design <- benchmark_grid(
  tasks = list(task_train, task_test),
  learners = list(learner_coxph, learner_lasso, learner_rsf),
  resamplings = rsmps("cv", folds = 5)
)
bm_design

# 实施比较
set.seed(42)
bm_result <- benchmark(bm_design)
measures <- msrs(c("surv.cindex"))
bm_result$aggregate(measures)
autoplot(bm_result)

bm_agg <- mlr3benchmark::as_benchmark_aggr(
  bm_result,
  measure = msr("surv.cindex")
)
bm_agg$rank_data(minimize = F)
bm_agg$friedman_test()
autoplot(bm_agg, type = "mean")

