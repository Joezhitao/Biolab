# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---ridge岭回归弹性网络
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
datarecipe_ridge <- recipe(formula = AHD ~ ., traindata) %>%
  step_dummy(all_nominal_predictors(), naming = new_dummy_names) %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_predictors())
datarecipe_ridge


# 设定模型
model_ridge <- logistic_reg(
  mode = "classification",
  engine = "glmnet",
  mixture = 0,  # 岭回归
  penalty = tune()
)
model_ridge

# workflow
wk_ridge <- 
  workflow() %>%
  add_recipe(datarecipe_ridge) %>%
  add_model(model_ridge)
wk_ridge

##############################################################
#########################  超参数寻优2选1-贝叶斯优化

# 贝叶斯优化超参数
set.seed(42)
tune_ridge <- wk_ridge %>%
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
eval_tune_ridge <- tune_ridge %>%
  collect_metrics()
eval_tune_ridge

# 图示
autoplot(tune_ridge)

# 经过交叉验证得到的最优超参数
hpbest_ridge <- tune_ridge %>%
  select_by_one_std_err(metric = "roc_auc", desc(penalty))
hpbest_ridge

# 采用最优超参数组合训练最终模型
set.seed(42)
final_ridge <- wk_ridge %>%
  finalize_workflow(hpbest_ridge) %>%
  fit(traindata)
final_ridge

##################################################################

# 训练集预测评估
predtrain_ridge <- eval4cls2(
  model = final_ridge, 
  dataset = traindata, 
  yname = "AHD", 
  modelname = "Ridge", 
  datasetname = "traindata",
  cutoff = "yueden",
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtrain_ridge$prediction
predtrain_ridge$predprobplot
predtrain_ridge$rocplot
predtrain_ridge$prplot
predtrain_ridge$caliplot
predtrain_ridge$cmplot
predtrain_ridge$metrics
predtrain_ridge$diycutoff
predtrain_ridge$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtrain_ridge$proc)
pROC::ci.auc(predtrain_ridge$proc)

# 预测评估测试集预测评估
predtest_ridge <- eval4cls2(
  model = final_ridge, 
  dataset = testdata, 
  yname = "AHD", 
  modelname = "Ridge", 
  datasetname = "testdata",
  cutoff = predtrain_ridge$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtest_ridge$prediction
predtest_ridge$predprobplot
predtest_ridge$rocplot
predtest_ridge$prplot
predtest_ridge$caliplot
predtest_ridge$cmplot
predtest_ridge$metrics
predtest_ridge$diycutoff
predtest_ridge$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtest_ridge$proc)
pROC::ci.auc(predtest_ridge$proc)

# ROC比较检验
pROC::roc.test(predtrain_ridge$proc, predtest_ridge$proc)


# 合并训练集和测试集上ROC曲线
predtrain_ridge$rocresult %>%
  bind_rows(predtest_ridge$rocresult) %>%
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
predtrain_ridge$prresult %>%
  bind_rows(predtest_ridge$prresult) %>%
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
predtrain_ridge$caliresult %>%
  bind_rows(predtest_ridge$caliresult) %>%
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
predtrain_ridge$metrics %>%
  bind_rows(predtest_ridge$metrics) %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)


# 最优超参数交叉验证的结果
evalcv_ridge <- bestcv4cls2(
  wkflow = wk_ridge,
  tuneresult = tune_ridge,
  hpbest = hpbest_ridge,
  yname = "AHD",
  modelname = "Ridge",
  v = 5,
  positivelevel = yourpositivelevel
)
evalcv_ridge$cvroc
evalcv_ridge$cvpr
evalcv_ridge$evalcv


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
  bind_cols(predict(final_ridge, new_data = newHeart, type = "prob"))%>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Yes >= predtrain_ridge$diycutoff, 
             yourpositivelevel, 
             yournegativelevel)
    )
  )
# readr::write_excel_csv(predresult, "弹性网络二分类预测结果.csv")
# 评估指标
prednew_ridge <- eval4cls2(
  model = final_ridge, 
  dataset = newHeart, 
  yname = "AHD", 
  modelname = "Ridge", 
  datasetname = "newHeart",
  cutoff = predtrain_ridge$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
prednew_ridge$prediction
prednew_ridge$predprobplot
prednew_ridge$rocplot
prednew_ridge$prplot
prednew_ridge$caliplot
prednew_ridge$cmplot
prednew_ridge$metrics
prednew_ridge$diycutoff
prednew_ridge$ksplot

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
final_ridge2 <- final_ridge %>%
  extract_fit_engine()
final_ridge2

# 非零系数自变量
tidy(final_ridge) %>%
  filter(term != "(Intercept)", estimate != 0) %>%
  pull(term)

######################## DALEX解释对象

explainer_ridge <- DALEXtra::explain_tidymodels(
  final_ridge, 
  data = traindatax,
  y = ifelse(traindata$AHD == yourpositivelevel, 1, 0),
  type = "classification",
  label = "Ridge"
)
# 变量重要性
vip_ridge <- viplot(explainer_ridge)
vipdata_ridge <- vip_ridge$data
vip_ridge$plot

# 变量偏依赖图
pdplot(explainer_ridge, convars)
pdplot(explainer_ridge, "Age")
pdplot(explainer_ridge, catvars)
pdplot(explainer_ridge, "Thal")

###################################### iml解释对象

predictor_ridge <- iml::Predictor$new(
  final_ridge, 
  data = traindatax,
  y = traindata$AHD,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>%
      rename_with(~gsub(".pred_", "", .x))
  },
  type = "prob"
)
# 交互作用
interact_ridge <- iml::Interaction$new(predictor_ridge)
plot(interact_ridge) +
  theme_minimal()

interact_ridge_1vo <- 
  iml::Interaction$new(predictor_ridge, feature = "Age")
plot(interact_ridge_1vo) +
  theme_minimal()

interact_ridge_1v1 <- iml::FeatureEffect$new(
  predictor_ridge, 
  feature = c("Age", "Oldpeak"),
  method = "pdp"
)
plot(interact_ridge_1v1) +
  scale_fill_viridis_c() +
  labs(fill = "") +
  theme_minimal()

###################################### lime单样本预测分解

explainer_ridge <- lime::lime(
  traindatax,
  lime::as_classifier(final_ridge, c(yournegativelevel, yourpositivelevel))
)
explanation_ridge <- lime::explain(
  traindatax[1,],  # 训练集第1个样本
  explainer_ridge, 
  n_labels = 2, 
  n_features = ncol(traindatax)
)
lime::plot_features(explanation_ridge)

######################## fastshap包

shapresult_ridge <- shap4cls2(
  finalmodel = final_ridge,
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
shapresult_ridge$shapvipplot

# 单样本预测分解
shap41 <- shapviz::shapviz(
  shapresult_ridge$shapley,
  X = traindatax
)
shapviz::sv_force(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
shapviz::sv_waterfall(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))

# 所有分类变量的shap图示
shapresult_ridge$shapplotd_facet
shapresult_ridge$shapplotd_one
# 所有连续变量的shap图示
shapresult_ridge$shapplotc_facet
shapresult_ridge$shapplotc_one
shapresult_ridge$shapplotc_one2
# 单变量shap图示
sdplot(shapresult_ridge, "Thal", "AHD")
sdplot(shapresult_ridge, "Age", "AHD")

# 所有变量一张图
# shap变量重要性
shapresult_ridge$shapvipplot_unity
# shap依赖图
shapresult_ridge$shapplot_unity


#################################################################

# 学习曲线
lcN <- 
  floor(seq(nrow(traindata)%/%2, nrow(traindata), length = 10))
lcresult_ridge <- data.frame()
for (i in lcN) {
  
  set.seed(i)
  traindatai <- traindata[sample(nrow(traindata), i), ]
  
  i_ridge <-  wk_ridge %>%
    finalize_workflow(hpbest_ridge) %>%
    fit(traindatai)
  
  predtrain_i_ridge <- traindatai %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_ridge, new_data = traindatai, type = "prob")) %>%
    yardstick::roc_auc(AHD,
                       paste0(".pred_", yourpositivelevel), 
                       event_level = "second") %>%
    mutate(datasetname = paste("traindata", i, sep = "-"))
  
  predtest_i_ridge <- testdata %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_ridge, new_data = testdata, type = "prob")) %>%
    yardstick::roc_auc(AHD, 
                       paste0(".pred_", yourpositivelevel),
                       event_level = "second") %>%
    mutate(datasetname = paste("testdata", i, sep = "-"))
  
  predi <- bind_rows(predtrain_i_ridge, predtest_i_ridge)
  lcresult_ridge <- rbind(lcresult_ridge, predi)
  print(i)
}
# 图示
lcresult_ridge %>%
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
save(datarecipe_ridge,
     model_ridge,
     wk_ridge,
     tune_ridge,
     predtrain_ridge,
     predtest_ridge,
     evalcv_ridge,
     vipdata_ridge,
     file = "F:/Mach_learn_data/Decision Tree/evalresult_ridge.RData")

# 保存模型结果供shiny部署之用，本课程不包括shiny内容
final_ridge_heart <- final_ridge
traindata_heart <- traindata
save(final_ridge_heart,
     traindata_heart,
     file = "F:/Mach_learn_data/Decision Tree/shiny_ridge_heart.RData")

