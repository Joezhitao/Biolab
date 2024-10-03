# R语言机器学习-特征筛选-无监督筛选+LASSO筛选

# 读取数据
simdata <- readr::read_csv("F:/Mach_learn_data/simdata.csv")
colnames(simdata) 
# 修正变量类型
# 将分类变量转换为factor
for(i in c(7, 8, 20, 27, 28)){ 
  simdata[[i]] <- factor(simdata[[i]])
}
skimr::skim(simdata)
DataExplorer::plot_correlation(
  simdata,
  cor_args = list(use = "pairwise.complete.obs"),
  title = "所有变量的相关性",
  ggtheme = ggplot2::theme_minimal()
)

# 数据拆分
set.seed(42)
datasplit <- 
  rsample::initial_split(simdata, prop = 0.75, strata = Cl)
traindata <- rsample::training(datasplit)
testdata <- rsample::testing(datasplit)

#################################################################

# 无监督筛选
library(recipes)
# 配方
datarecipe <- recipe(Y+Cl+Cl3 ~ ., traindata) %>%
  # 主观剔除，样本id之类的变量或者主观分析要剔除的变量
  # step_rm(id) %>%
  # 按缺失比例剔除变量，threshold设定缺失比例界限
  step_filter_missing(all_predictors(), threshold = 0) %>%
  # 剔除零方差变量，所有样本取值一样的变量
  step_zv(all_predictors()) %>%
  # 分类独热编码，后续部分算法要求所有自变量为数值
  step_dummy(all_nominal_predictors()) %>%
  # 剔除近零方差变量，取值水平少，个别取值水平占绝大多数的变量
  step_nzv(all_predictors(), freq_cut = 95/5, unique_cut = 10) %>%
  # 剔除高度相关的变量
  step_corr(all_predictors()) %>%
  prep()
# 以上所述每一步处理了哪些变量
datarecipe

# 处理
trainx <- bake(datarecipe, new_data = traindata) %>%
  select(-Y, -Cl, -Cl3)
trainy <- traindata$Y
traincl <- traindata$Cl
traincl3 <- traindata$Cl3

testx <- bake(datarecipe, new_data = testdata) %>%
  select(-Y, -Cl, -Cl3)
testy <- testdata$Y
testcl <- testdata$Cl
testcl3 <- testdata$Cl3

#################################################################

# 有监督筛选
# 选其一或者多种结合均可

#################################


# lasso模型
# 自变量因变量分别输入
# 自变量部分为数值型矩阵

#################################  回归

# 交叉验证寻找最优的lambda值
library(glmnet)
set.seed(42)
cvlasso_reg <- cv.glmnet(
  x = as.matrix(trainx), 
  y = trainy, 
  family = "gaussian", 
  alpha = 1,
  type.measure = "mse"
)
# 模型概要
cvlasso_reg
# lambda参数与评价指标之间的关系图示
plot(cvlasso_reg, xlab = "Log Lambda")

# 系数轨迹
set.seed(42)
lasso_reg <- glmnet(
  x = as.matrix(trainx), 
  y = trainy,
  family = "gaussian",
  alpha = 1 
)
# 模型概要
lasso_reg
# lambda参数与各个变量的系数之间的关系图示
plot(lasso_reg, xvar = "lambda")


# lambda.1se对应的非零变量
library(tidyverse)
predict(cvlasso_reg, type = "coef", s = "lambda.1se") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("X") %>%
  filter(X != "(Intercept)") %>%
  filter(lambda.1se != 0) %>%
  pull(X)

# lambda.min对应的非零变量
predict(cvlasso_reg, type = "coef", s = "lambda.min") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("X") %>%
  filter(X != "(Intercept)") %>%
  filter(lambda.min != 0) %>%
  pull(X)

#################################  二分类

# 交叉验证寻找最优的lambda值
library(glmnet)
set.seed(42)
cvlasso_cls2 <- cv.glmnet(
  x = as.matrix(trainx), 
  y = traincl, 
  family = "binomial", 
  alpha = 1,
  type.measure = "auc"
)
# 模型概要
cvlasso_cls2
# lambda参数与评价指标之间的关系图示
plot(cvlasso_cls2, xlab = "Log Lambda")

# 系数轨迹
set.seed(42)
lasso_cls2 <- glmnet(
  x = as.matrix(trainx), 
  y = traincl,
  family = "binomial",
  alpha = 1 
)
# 模型概要
lasso_cls2
# lambda参数与各个变量的系数之间的关系图示
plot(lasso_cls2, xvar = "lambda")


# lambda.1se对应的非零变量
library(tidyverse)
predict(cvlasso_cls2, type = "coef", s = "lambda.1se") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("X") %>%
  filter(X != "(Intercept)") %>%
  filter(lambda.1se != 0) %>%
  pull(X)

# lambda.min对应的非零变量
predict(cvlasso_cls2, type = "coef", s = "lambda.min") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("X") %>%
  filter(X != "(Intercept)") %>%
  filter(lambda.min != 0) %>%
  pull(X)


#################################  多分类

# 交叉验证寻找最优的lambda值
library(glmnet)
set.seed(42)
cvlasso_cls3 <- cv.glmnet(
  x = as.matrix(trainx), 
  y = traincl3, 
  family = "multinomial", 
  type.multinomial = "grouped", 
  alpha = 1,
  type.measure = "class",
  maxit = 1e+06  # 最大迭代次数，1e+06是科学计数法，即10的6次方
)
# 模型概要
cvlasso_cls3
# lambda参数与评价指标之间的关系图示
plot(cvlasso_cls3, xlab = "Log Lambda")

# 系数轨迹
set.seed(42)
lasso_cls3 <- glmnet(
  x = as.matrix(trainx), 
  y = traincl3,
  family = "multinomial",
  type.multinomial = "grouped",
  alpha = 1 ,
  maxit = 1e+06  # 最大迭代次数，1e+06是科学计数法，即1乘以10的6次方
)
# 模型概要
lasso_cls3
# lambda参数与各个变量的系数之间的关系图示
plot(lasso_cls3, xvar = "lambda")


# lambda.1se对应的非零变量
library(tidyverse)
coefs1se <- predict(cvlasso_cls3, type = "coef", s = "lambda.1se") %>%
  lapply(function(x) as.data.frame(as.matrix(x)))
coefs1se[[1]] %>%
  rownames_to_column("X") %>%
  filter(X != "(Intercept)") %>%
  filter(`1` != 0) %>%
  pull(X) -> coefs1
coefs1se[[2]] %>%
  rownames_to_column("X") %>%
  filter(X != "(Intercept)") %>%
  filter(`1` != 0) %>%
  pull(X) -> coefs2
coefs3 <- coefs1se[[3]] %>%
  rownames_to_column("X") %>%
  filter(X != "(Intercept)") %>%
  filter(`1` != 0) %>%
  pull(X)
# 以上三者可以取并集unique(c(...))

selxs_lasso <- unique(c(coefs1, coefs2, coefs3))
trainx_lasso <- trainx[, selxs_lasso]

# lambda.min对应的非零变量
coefsmin <- predict(cvlasso_cls3, type = "coef", s = "lambda.min") %>%
  lapply(function(x) as.data.frame(as.matrix(x)))
coefsmin[[1]] %>%
  rownames_to_column("X") %>%
  filter(X != "(Intercept)") %>%
  filter(`1` != 0) %>%
  pull(X)
coefsmin[[2]] %>%
  rownames_to_column("X") %>%
  filter(X != "(Intercept)") %>%
  filter(`1` != 0) %>%
  pull(X)
coefsmin[[3]] %>%
  rownames_to_column("X") %>%
  filter(X != "(Intercept)") %>%
  filter(`1` != 0) %>%
  pull(X)
# 以上三者可以取并集unique(c(...))

