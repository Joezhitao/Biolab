#Joe_2024_10_12
# 加载包
library(tidyverse)
library(survival)
library(survminer)
library(rms)
library(riskRegression)

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

###################################################

# 生存曲线
sc_km <- survfit(Surv(rfstime, status) ~ 1, data = traindata2)
sc_km
# summary(sc_km)
plot(sc_km)
ggsurvplot(sc_km)

###################################################

# 构建公式
coxformula <- as.formula(
  paste0(
    "Surv(rfstime, status) ~ ",
    paste(colnames(traindata2)[-c(1,2)], collapse = " + ")
  )
)
coxformula

# 采用survival包coxph函数构建模型
fit_coxph <- 
  coxph(coxformula, data = as.data.frame(traindata2), x = T)
fit_coxph
summary(fit_coxph)
# 森林图
ggforest(fit_coxph)

# 模型诊断
# 1-比例风险假设-Schoenfeld残差
test.ph <- cox.zph(fit_coxph)
test.ph
ggcoxzph(test.ph)
# 违反比例风险假设的情形可以采用以下方法
# 引入时间和自变量的交互项
# 分层
# 2-强影响点或者离群点-Deviance residual
ggcoxdiagnostics(
  fit_coxph, 
  type = "dfbeta",
  linear.predictions = FALSE, 
  ggtheme = theme_bw()
)
ggcoxdiagnostics(
  fit_coxph, 
  type = "deviance",
  linear.predictions = FALSE, 
  ggtheme = theme_bw()
)
# 3-对数风险和自变量是否存在非线性关系-Martingale残差-针对连续型自变量
ggcoxfunctional(
  coxformula, 
  data = as.data.frame(traindata2)
)

predictCox(fit_coxph, times = itps, newdata = traindata2)

# 采用rms包cph函数构建coxph模型
options(datadist = datadist(traindata2))
fit_cph <- cph(
  coxformula, 
  data = as.data.frame(traindata2),
  model = T,
  x = T, 
  y = T, 
  surv = T
)
fit_cph
# 生存函数
surv <- Survival(fit_cph)
# 列线图-nomogram
plot(
  nomogram(
    fit_cph, 
    lp = F,
    fun = list(function(x) surv(365, x),
               function(x) surv(730, x)),
    fun.at = list(c(0.3, 0.5, 0.9),
                  c(0.1, 0.5, 0.9)),
    funlabel = c("1-year Survival Prob",
                 "2-year Survival Prob")
  )
)


