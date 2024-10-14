# Joe_2024_10_12
# 以下代码来自hdnom包的说明文档
# https://cran.r-project.org/web/packages/hdnom/vignettes/hdnom.html
# 部分地方做了修改，仅供参考
# hdnom包支持一系列带惩罚的生存分析模型，具体可参看网页链接

# 加载包
library(survival)
library(survminer)
library(hdnom)

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
  dplyr::filter(rfstime > 0) # 剔除时间非正的样本

# 数据概况
skimr::skim(sadata)

# 感兴趣的时间节点
range(unique(sadata$rfstime))
itps <- c(365 * c(1, 3, 5))
itps
table(cut(sadata$rfstime, c(0, itps, Inf)))

###################################################

# 数据拆分构建任务对象
set.seed(42)
datasplit <- rsample::initial_split(
  sadata, prop = 0.8, strata = rfstime, breaks = 10
)
traindata <- rsample::training(datasplit) %>%
  dplyr::sample_n(nrow(.))
testdata <- rsample::testing(datasplit) %>%
  dplyr::sample_n(nrow(.))

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

###################################################

# 数据预处理
library(recipes)
datarecipe_lasso <- recipe(rfstime + status ~ ., traindata) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep()
datarecipe_lasso

# 按方处理训练集和测试集
traindata2 <- bake(datarecipe_lasso, new_data = NULL) %>%
  dplyr::select(rfstime, status, everything())
testdata2 <- bake(datarecipe_lasso, new_data = testdata) %>%
  dplyr::select(rfstime, status, everything())

# 拆分自变量部分和因变量部分
# 自变量部分处理为矩阵
trainx <- as.matrix(traindata2[, -c(1, 2)])
testx <- as.matrix(testdata2[, -c(1, 2)])
# 因变量部分处理为生存对象
trainy <- Surv(traindata2$rfstime, traindata2$status)
testy <- Surv(testdata2$rfstime, testdata2$status)

###################################################

# 多核并行
library(doParallel)
cl <- makePSOCKcluster(3)
registerDoParallel(cl)
# 结束多核并行
# stopCluster(cl)

# 自适应弹性网络
fit_aen <- fit_aenet(
  x = trainx, 
  y = trainy, 
  nfolds = 10, 
  rule = "lambda.1se", 
  seed = c(42,43)
)
fit_aen
# 自变量系数
coef(fit_aen$model)  # 从中可以看到模型最终保留的自变量

# 列线图nomogram
nom_aen <- as_nomogram(
  fit_aen,
  x = trainx,
  time = traindata2$rfstime,
  event = traindata2$status,
  pred.at = 365*5,
  funlabel = "5年生存概率"
)
nom_aen
plot(nom_aen)

#################################################

# 内部验证
inval_aen <- validate(
  # 数据
  x = trainx,
  time = traindata2$rfstime,
  event = traindata2$status,
  # 模型
  model.type = "aenet",
  alpha = fit_aen$alpha,
  lambda = fit_aen$lambda,
  pen.factor = fit_aen$pen_factor,
  # 验证设定
  method = "bootstrap",
  boot.times = 15,
  tauc.type = "UNO",
  tauc.time = itps,
  seed = 42,
  trace = T
)
inval_aen
summary(inval_aen)
plot(inval_aen)
# 实线是均值、虚线是中位数

# 外部验证
outval_aen <- validate_external(
  # 模型
  fit_aen,
  x = trainx,
  time = traindata2$rfstime,
  event = traindata2$status,
  # 外部数据
  x_new = testx,
  time_new = testdata2$rfstime,
  event_new = testdata2$status,
  # 验证设定
  tauc.type = "UNO",
  tauc.time = itps
)
outval_aen
summary(outval_aen)
plot(outval_aen)

#################################################

# 内部校准
incal_aen <- calibrate(
  # 数据
  x = trainx,
  time = traindata2$rfstime,
  event = traindata2$status,
  # 模型
  model.type = "aenet",
  alpha = fit_aen$alpha,
  lambda = fit_aen$lambda,
  pen.factor = fit_aen$pen_factor,
  # 校准设定
  method = "bootstrap",
  boot.times = 15,
  pred.at = itps[1],
  ngroup = 5,
  seed = 42,
  trace = T
)
incal_aen
summary(incal_aen)
# 校准曲线
plot(incal_aen, 
     xlim = c(0.9, 1), 
     ylim = c(0.9, 1))
# 不同校准组别的KM曲线
kmplot(incal_aen,
       group.name = as.character(1:5),  # 对应之前的ngroup
       time.at = itps)
# 不同校准组别的差异检验
logrank_incal_aen <- logrank_test(incal_aen)
logrank_incal_aen

# 外部校准
outcal_aen <- calibrate_external(
  # 模型
  fit_aen,
  x = trainx,
  time = traindata2$rfstime,
  event = traindata2$status,
  # 外部数据
  x_new = testx,
  time_new = testdata2$rfstime,
  event_new = testdata2$status,
  # 验证设定
  pred.at = itps[1],
  ngroup = 3
)
outcal_aen
summary(outcal_aen)
# 校准曲线
plot(outcal_aen, 
     xlim = c(0.9, 1), 
     ylim = c(0.9, 1))
# 不同校准组别的KM曲线
kmplot(outcal_aen,
       group.name = as.character(1:3),  # 对应之前的ngroup
       time.at = itps)
# 不同校准组别的差异检验
logrank_outcal_aen <- logrank_test(outcal_aen)
logrank_outcal_aen

#################################################

# 模型比较-验证比较
valcomp <- compare_by_validate(
  # 数据
  x = trainx,
  time = traindata2$rfstime,
  event = traindata2$status,
  # 待比较的模型
  model.type = c("aenet", "enet"),
  # 验证设定
  method = "bootstrap",
  boot.times = 15,
  tauc.type = "UNO",
  tauc.time = itps,
  seed = 42,
  trace = T
)
valcomp
plot(valcomp, interval = T)

#################################################

# 模型比较-校准比较
calcomp <- compare_by_calibrate(
  # 数据
  x = trainx,
  time = traindata2$rfstime,
  event = traindata2$status,
  # 模型
  model.type = c("aenet", "enet"),
  # 校准设定
  method = "bootstrap",
  boot.times = 15,
  pred.at = itps[1],
  ngroup = 5,
  seed = 42,
  trace = T
)
calcomp
plot(calcomp, 
     xlim = c(0.9, 1), 
     ylim = c(0.9, 1))

#################################################

# 预测新样本
predtest_aen <- predict(
  # 模型
  fit_aen,
  x = trainx,
  y = trainy,
  # 外部数据
  newx = testx,
  # 验证设定
  pred.at = itps
)
predtest_aen



