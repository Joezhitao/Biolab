# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---LASSO岭回归弹性网络
# https://parsnip.tidymodels.org/reference/details_logistic_reg_glmnet.html

##############################################################

# install.packages("tidymodels")
library(tidymodels)
source("H:/Biolab/Biolab/tidyfuncs4cls2_v18.R")

# 多核并行
library(doParallel)
registerDoParallel(
  makePSOCKcluster(
    max(1, (parallel::detectCores(logical = F))-1)
  )
)

# 读取数据
# file.choose()
Heart <- readr::read_csv("F:/Mach_learn_data/data-Heart4cls2.csv")
colnames(Heart) 
# 修正变量类型
# 将分类变量转换为factor
for(i in c(3,4,7,8,10,12,14,15)){ 
  Heart[[i]] <- factor(Heart[[i]])
}

# 删除无关变量在此处进行
Heart$Id <- NULL
# 删除含有缺失值的样本在此处进行，填充缺失值在后面
Heart <- na.omit(Heart)
# Heart <- Heart %>%
#   drop_na(Thal)

# 数据概况
skimr::skim(Heart)    

# 设定阳性类别和阴性类别
yourpositivelevel <- "Yes"
yournegativelevel <- "No"
# 转换因变量的因子水平，将阳性类别设定为第二个水平
levels(Heart$AHD)
table(Heart$AHD)
Heart$AHD <- factor(
  Heart$AHD,
  levels = c(yournegativelevel, yourpositivelevel)
)
levels(Heart$AHD)
table(Heart$AHD)

##############################################################

# 数据拆分
set.seed(42)
datasplit <- initial_split(Heart, prop = 0.75, strata = AHD)
traindata <- training(datasplit) %>%
  sample_n(nrow(.))
testdata <- testing(datasplit) %>%
  sample_n(nrow(.))

# 重抽样设定-5折交叉验证
set.seed(42)
folds <- vfold_cv(traindata, v = 5, strata = AHD)
folds

# 数据预处理配方
datarecipe_lasso <- recipe(formula = AHD ~ ., traindata) %>%
  step_dummy(all_nominal_predictors(), naming = new_dummy_names) %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_predictors())
datarecipe_lasso


# 设定模型
model_lasso <- logistic_reg(
  mode = "classification",
  engine = "glmnet",
  mixture = 1,   # LASSO
  penalty = tune()
)
model_lasso

# workflow
wk_lasso <- 
  workflow() %>%
  add_recipe(datarecipe_lasso) %>%
  add_model(model_lasso)
wk_lasso

##############################################################
#########################  超参数寻优2选1-贝叶斯优化

# 贝叶斯优化超参数
set.seed(42)
tune_lasso <- wk_lasso %>%
  tune_bayes(
    resamples = folds,
    initial = 10,
    iter = 50,
    metrics = metricset_cls2,
    control = control_bayes(save_pred = T, 
                            verbose = T,
                            no_improve = 10,
                            uncertain = 5,
                            event_level = "second",
                            parallel_over = "everything",
                            save_workflow = T)
  )

########################  超参数寻优结束


# 交叉验证结果
eval_tune_lasso <- tune_lasso %>%
  collect_metrics()
eval_tune_lasso

# 图示
autoplot(tune_lasso)

# 经过交叉验证得到的最优超参数
hpbest_lasso <- tune_lasso %>%
  select_by_one_std_err(metric = "roc_auc", desc(penalty))
hpbest_lasso

# 采用最优超参数组合训练最终模型
set.seed(42)
final_lasso <- wk_lasso %>%
  finalize_workflow(hpbest_lasso) %>%
  fit(traindata)
final_lasso

##################################################################

# 训练集预测评估
predtrain_lasso <- eval4cls2(
  model = final_lasso, 
  dataset = traindata, 
  yname = "AHD", 
  modelname = "LASSO", 
  datasetname = "traindata",
  cutoff = "yueden",
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtrain_lasso$prediction
predtrain_lasso$predprobplot
predtrain_lasso$rocplot
predtrain_lasso$prplot
predtrain_lasso$caliplot
predtrain_lasso$cmplot
predtrain_lasso$metrics
predtrain_lasso$diycutoff
predtrain_lasso$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtrain_lasso$proc)
pROC::ci.auc(predtrain_lasso$proc)

# 预测评估测试集预测评估
predtest_lasso <- eval4cls2(
  model = final_lasso, 
  dataset = testdata, 
  yname = "AHD", 
  modelname = "LASSO", 
  datasetname = "testdata",
  cutoff = predtrain_lasso$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtest_lasso$prediction
predtest_lasso$predprobplot
predtest_lasso$rocplot
predtest_lasso$prplot
predtest_lasso$caliplot
predtest_lasso$cmplot
predtest_lasso$metrics
predtest_lasso$diycutoff
predtest_lasso$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtest_lasso$proc)
pROC::ci.auc(predtest_lasso$proc)

# ROC比较检验
pROC::roc.test(predtrain_lasso$proc, predtest_lasso$proc)


# 合并训练集和测试集上ROC曲线
predtrain_lasso$rocresult %>%
  bind_rows(predtest_lasso$rocresult) %>%
  mutate(dataAUC = paste(data, curvelab),
         dataAUC = forcats::as_factor(dataAUC)) %>%
  ggplot(aes(x = 1-specificity,
             y = sensitivity, 
             color = dataAUC)) +
  geom_path(linewidth = 1) +
  geom_abline(linetype = "dashed") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = seq(0, 1, by = 0.2)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = seq(0, 1, by = 0.2)) +
  labs(color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        text = element_text(family = "serif"))

# 合并训练集和测试集上PR曲线
predtrain_lasso$prresult %>%
  bind_rows(predtest_lasso$prresult) %>%
  mutate(dataAUC = paste(data, curvelab),
         dataAUC = forcats::as_factor(dataAUC)) %>%
  ggplot(aes(x = recall,
             y = precision, 
             color = dataAUC)) +
  geom_path(linewidth = 1) +
  geom_abline(linetype = "dashed", slope = -1, intercept = 1) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = seq(0, 1, by = 0.2)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = seq(0, 1, by = 0.2)) +
  labs(color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        text = element_text(family = "serif"))

# 合并训练集和测试集上校准曲线
predtrain_lasso$caliresult %>%
  bind_rows(predtest_lasso$caliresult) %>%
  mutate(data = forcats::as_factor(data)) %>%
  ggplot(aes(x = predprobgroup,
             y = Fraction, 
             color = data)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, pch = 15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.1),
                     labels = seq(0, 1, by = 0.1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.1),
                     labels = seq(0, 1, by = 0.1)) +
  labs(x = "Bin Midpoint", y = "Event Rate", color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        text = element_text(family = "serif"))

# 合并训练集和测试集上性能指标
predtrain_lasso$metrics %>%
  bind_rows(predtest_lasso$metrics) %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)


# 最优超参数交叉验证的结果
evalcv_lasso <- bestcv4cls2(
  wkflow = wk_lasso,
  tuneresult = tune_lasso,
  hpbest = hpbest_lasso,
  yname = "AHD",
  modelname = "LASSO",
  v = 5,
  positivelevel = yourpositivelevel
)
evalcv_lasso$cvroc
evalcv_lasso$cvpr
evalcv_lasso$evalcv


###################################################################

# 预测新样本
# file.choose()
newHeart <- readr::read_csv("F:/Mach_learn_data/newdata-Heart4cls2.csv")
# 修正变量类型-将分类变量转换为factor
for(i in c(3,4,7,8,10,12,14,15)){ 
  newHeart[[i]] <- factor(newHeart[[i]])
}
# 删除无关变量在此处进行
newHeart$Id <- NULL
# 删除含有缺失值的样本在此处进行，填充缺失值在后面
newHeart <- na.omit(newHeart)
# newHeart <- newHeart %>%
#   drop_na(Thal)
# 转换因变量的因子水平，将阳性类别设定为第二个水平
newHeart$AHD <- factor(
  newHeart$AHD,
  levels = c(yournegativelevel, yourpositivelevel)
)
# 预测
predresult <- newHeart %>%
  bind_cols(predict(final_lasso, new_data = newHeart, type = "prob"))%>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Yes >= predtrain_lasso$diycutoff, 
             yourpositivelevel, 
             yournegativelevel)
    )
  )
# readr::write_excel_csv(predresult, "弹性网络二分类预测结果.csv")
# 评估指标
prednew_lasso <- eval4cls2(
  model = final_lasso, 
  dataset = newHeart, 
  yname = "AHD", 
  modelname = "LASSO", 
  datasetname = "newHeart",
  cutoff = predtrain_lasso$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
prednew_lasso$prediction
prednew_lasso$predprobplot
prednew_lasso$rocplot
prednew_lasso$prplot
prednew_lasso$caliplot
prednew_lasso$cmplot
prednew_lasso$metrics
prednew_lasso$diycutoff
prednew_lasso$ksplot

###################################################################

# 自变量数据集
colnames(traindata)
traindatax <- traindata %>%
  dplyr::select(-AHD)
colnames(traindatax)

# 分类型、连续型自变量名称
catvars <- getcategory(traindatax)
convars <- getcontinuous(traindatax)

# 提取最终的算法模型
final_lasso2 <- final_lasso %>%
  extract_fit_engine()
final_lasso2

# 非零系数自变量
tidy(final_lasso) %>%
  filter(term != "(Intercept)", estimate != 0) %>%
  pull(term)

######################## DALEX解释对象

explainer_lasso <- DALEXtra::explain_tidymodels(
  final_lasso, 
  data = traindatax,
  y = ifelse(traindata$AHD == yourpositivelevel, 1, 0),
  type = "classification",
  label = "LASSO"
)
# 变量重要性
vip_lasso <- viplot(explainer_lasso)
vipdata_lasso <- vip_lasso$data
vip_lasso$plot

# 变量偏依赖图
pdplot(explainer_lasso, convars)
pdplot(explainer_lasso, "Age")
pdplot(explainer_lasso, catvars)
pdplot(explainer_lasso, "Thal")

###################################### iml解释对象

predictor_lasso <- iml::Predictor$new(
  final_lasso, 
  data = traindatax,
  y = traindata$AHD,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>%
      rename_with(~gsub(".pred_", "", .x))
  },
  type = "prob"
)
# 交互作用
interact_lasso <- iml::Interaction$new(predictor_lasso)
plot(interact_lasso) +
  theme_minimal()

interact_lasso_1vo <- 
  iml::Interaction$new(predictor_lasso, feature = "Age")
plot(interact_lasso_1vo) +
  theme_minimal()

interact_lasso_1v1 <- iml::FeatureEffect$new(
  predictor_lasso, 
  feature = c("Age", "Oldpeak"),
  method = "pdp"
)
plot(interact_lasso_1v1) +
  scale_fill_viridis_c() +
  labs(fill = "") +
  theme_minimal()

###################################### lime单样本预测分解

explainer_lasso <- lime::lime(
  traindatax,
  lime::as_classifier(final_lasso, c(yournegativelevel, yourpositivelevel))
)
explanation_lasso <- lime::explain(
  traindatax[1,],  # 训练集第1个样本
  explainer_lasso, 
  n_labels = 2, 
  n_features = ncol(traindatax)
)
lime::plot_features(explanation_lasso)

######################## fastshap包

shapresult_lasso <- shap4cls2(
  finalmodel = final_lasso,
  predfunc = function(model, newdata) {
    predict(model, newdata, type = "prob") %>%
      dplyr::select(ends_with(yourpositivelevel)) %>%
      pull()
  },
  datax = traindatax,
  datay = traindata$AHD,
  yname = "AHD",
  flname = catvars,
  lxname = convars
)

# 基于shap的变量重要性
shapresult_lasso$shapvipplot

# 单样本预测分解
shap41 <- shapviz::shapviz(
  shapresult_lasso$shapley,
  X = traindatax
)
shapviz::sv_force(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
shapviz::sv_waterfall(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))

# 所有分类变量的shap图示
shapresult_lasso$shapplotd_facet
shapresult_lasso$shapplotd_one
# 所有连续变量的shap图示
shapresult_lasso$shapplotc_facet
shapresult_lasso$shapplotc_one
shapresult_lasso$shapplotc_one2
# 单变量shap图示
sdplot(shapresult_lasso, "Thal", "AHD")
sdplot(shapresult_lasso, "Age", "AHD")

# 所有变量一张图
# shap变量重要性
shapresult_lasso$shapvipplot_unity
# shap依赖图
shapresult_lasso$shapplot_unity


#################################################################

# 学习曲线
lcN <- 
  floor(seq(nrow(traindata)%/%2, nrow(traindata), length = 10))
lcresult_lasso <- data.frame()
for (i in lcN) {
  
  set.seed(i)
  traindatai <- traindata[sample(nrow(traindata), i), ]
  
  i_lasso <-  wk_lasso %>%
    finalize_workflow(hpbest_lasso) %>%
    fit(traindatai)
  
  predtrain_i_lasso <- traindatai %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_lasso, new_data = traindatai, type = "prob")) %>%
    yardstick::roc_auc(AHD,
                       paste0(".pred_", yourpositivelevel), 
                       event_level = "second") %>%
    mutate(datasetname = paste("traindata", i, sep = "-"))
  
  predtest_i_lasso <- testdata %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_lasso, new_data = testdata, type = "prob")) %>%
    yardstick::roc_auc(AHD, 
                       paste0(".pred_", yourpositivelevel),
                       event_level = "second") %>%
    mutate(datasetname = paste("testdata", i, sep = "-"))
  
  predi <- bind_rows(predtrain_i_lasso, predtest_i_lasso)
  lcresult_lasso <- rbind(lcresult_lasso, predi)
  print(i)
}
# 图示
lcresult_lasso %>%
  separate(datasetname, into = c("dataset", "N"), sep = "-") %>%
  mutate(N = as.numeric(N),
         dataset = forcats::as_factor(dataset)) %>%
  ggplot(aes(x = N, y = .estimate, color = dataset)) +
  geom_point() +
  geom_smooth(se = F, method = 'loess', formula = 'y ~ x') +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Samples in traindata", y = "ROCAUC", color = "") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = "serif"))

######################################################################

# 保存评估结果
save(datarecipe_lasso,
     model_lasso,
     wk_lasso,
     tune_lasso,
     predtrain_lasso,
     predtest_lasso,
     evalcv_lasso,
     vipdata_lasso,
     file = "F:/Mach_learn_data/Decision Tree/evalresult_lasso.RData")

# 保存模型结果供shiny部署之用，本课程不包括shiny内容
final_lasso_heart <- final_lasso
traindata_heart <- traindata
save(final_lasso_heart,
     traindata_heart,
     file = "F:/Mach_learn_data/Decision Tree/shiny_lasso_heart.RData")
