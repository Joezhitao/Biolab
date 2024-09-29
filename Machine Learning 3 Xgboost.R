# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---xgboost
# https://parsnip.tidymodels.org/reference/details_boost_tree_xgboost.html

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
datarecipe_xgboost <- recipe(formula = AHD ~ ., traindata) %>%
  step_dummy(all_nominal_predictors(), naming = new_dummy_names) %>%
  step_nzv(all_predictors())
datarecipe_xgboost


# 设定模型
model_xgboost <- boost_tree(
  mode = "classification",
  engine = "xgboost",
  mtry = tune(),
  trees = 1000,
  min_n = tune(),
  tree_depth = tune(),
  learn_rate = tune(),
  loss_reduction = tune(),
  sample_size = tune(),
  stop_iter = 25
) %>%
  set_args(validation = 0.2,
           event_level = "second")
model_xgboost

# workflow
wk_xgboost <- 
  workflow() %>%
  add_recipe(datarecipe_xgboost) %>%
  add_model(model_xgboost)
wk_xgboost

#############################################################
#########################  超参数寻优2选1-贝叶斯优化

# 更新超参数范围
param_xgboost <- model_xgboost %>%
  extract_parameter_set_dials() %>%
  update(mtry = mtry(c(2, 10)))

# 贝叶斯优化超参数
set.seed(42)
tune_xgboost <- wk_xgboost %>%
  tune_bayes(
    resamples = folds,
    initial = 10,
    iter = 50,
    param_info = param_xgboost,
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
eval_tune_xgboost <- tune_xgboost %>%
  collect_metrics()
eval_tune_xgboost

# 图示
# autoplot(tune_xgboost)
eval_tune_xgboost %>% 
  filter(.metric == "roc_auc") %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~mean, colorscale = 'Jet', showscale = T),
    dimensions = list(
      list(label = 'mtry', values = ~mtry),
      list(label = 'min_n', values = ~min_n),
      list(label = 'tree_depth', values = ~tree_depth),
      list(label = 'learn_rate', values = ~learn_rate),
      list(label = 'loss_reduction', values = ~loss_reduction),
      list(label = 'sample_size', values = ~sample_size)
    )
  ) %>%
  plotly::layout(title = "xgboost HPO Guided by AUCROC",
                 font = list(family = "serif"))

# 经过交叉验证得到的最优超参数
hpbest_xgboost <- tune_xgboost %>%
  select_by_one_std_err(metric = "roc_auc", desc(min_n))
hpbest_xgboost

# 采用最优超参数组合训练最终模型
set.seed(42)
final_xgboost <- wk_xgboost %>%
  finalize_workflow(hpbest_xgboost) %>%
  fit(traindata)
final_xgboost

##################################################################

# 训练集预测评估
predtrain_xgboost <- eval4cls2(
  model = final_xgboost, 
  dataset = traindata, 
  yname = "AHD", 
  modelname = "Xgboost", 
  datasetname = "traindata",
  cutoff = "yueden",
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtrain_xgboost$prediction
predtrain_xgboost$predprobplot
predtrain_xgboost$rocplot
predtrain_xgboost$prplot
predtrain_xgboost$caliplot
predtrain_xgboost$cmplot
predtrain_xgboost$metrics
predtrain_xgboost$diycutoff
predtrain_xgboost$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtrain_xgboost$proc)
pROC::ci.auc(predtrain_xgboost$proc)

# 预测评估测试集预测评估
predtest_xgboost <- eval4cls2(
  model = final_xgboost, 
  dataset = testdata, 
  yname = "AHD", 
  modelname = "Xgboost", 
  datasetname = "testdata",
  cutoff = predtrain_xgboost$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtest_xgboost$prediction
predtest_xgboost$predprobplot
predtest_xgboost$rocplot
predtest_xgboost$prplot
predtest_xgboost$caliplot
predtest_xgboost$cmplot
predtest_xgboost$metrics
predtest_xgboost$diycutoff
predtest_xgboost$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtest_xgboost$proc)
pROC::ci.auc(predtest_xgboost$proc)

# ROC比较检验
pROC::roc.test(predtrain_xgboost$proc, predtest_xgboost$proc)


# 合并训练集和测试集上ROC曲线
predtrain_xgboost$rocresult %>%
  bind_rows(predtest_xgboost$rocresult) %>%
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
predtrain_xgboost$prresult %>%
  bind_rows(predtest_xgboost$prresult) %>%
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
predtrain_xgboost$caliresult %>%
  bind_rows(predtest_xgboost$caliresult) %>%
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
predtrain_xgboost$metrics %>%
  bind_rows(predtest_xgboost$metrics) %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)


# 最优超参数交叉验证的结果
evalcv_xgboost <- bestcv4cls2(
  wkflow = wk_xgboost,
  tuneresult = tune_xgboost,
  hpbest = hpbest_xgboost,
  yname = "AHD",
  modelname = "Xgboost",
  v = 5,
  positivelevel = yourpositivelevel
)
evalcv_xgboost$cvroc
evalcv_xgboost$cvpr
evalcv_xgboost$evalcv


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
  bind_cols(predict(final_xgboost, new_data = newHeart, type = "prob"))%>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Yes >= predtrain_xgboost$diycutoff, 
             yourpositivelevel, 
             yournegativelevel)
    )
  )
# readr::write_excel_csv(predresult, "Xgboost二分类预测结果.csv")
# 评估指标
prednew_xgboost <- eval4cls2(
  model = final_xgboost, 
  dataset = newHeart, 
  yname = "AHD", 
  modelname = "Xgboost", 
  datasetname = "newHeart",
  cutoff = predtrain_xgboost$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
prednew_xgboost$prediction
prednew_xgboost$predprobplot
prednew_xgboost$rocplot
prednew_xgboost$prplot
prednew_xgboost$caliplot
prednew_xgboost$cmplot
prednew_xgboost$metrics
prednew_xgboost$diycutoff
prednew_xgboost$ksplot

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
final_xgboost2 <- final_xgboost %>%
  extract_fit_engine()
final_xgboost2

# 变量重要性
importance_matrix <- 
  xgboost::xgb.importance(model = final_xgboost2)
print(importance_matrix)
xgboost::xgb.plot.importance(
  importance_matrix = importance_matrix,
  measure = "Cover",
  col = "skyblue",
  family = "serif"
)
# SHAP
xgboost::xgb.plot.shap(
  data = as.matrix(final_xgboost %>%
                     extract_recipe() %>%
                     bake(new_data = traindata) %>%
                     dplyr::select(-AHD)), 
  model = final_xgboost2,
  top_n = 5,
  las = 1,
  family = "serif"
)

######################## DALEX解释对象

explainer_xgboost <- DALEXtra::explain_tidymodels(
  final_xgboost, 
  data = traindatax,
  y = ifelse(traindata$AHD == yourpositivelevel, 1, 0),
  type = "classification",
  label = "Xgboost"
)
# 变量重要性
vip_xgboost <- viplot(explainer_xgboost)
vipdata_xgboost <- vip_xgboost$data
vip_xgboost$plot

# 变量偏依赖图
pdplot(explainer_xgboost, convars)
pdplot(explainer_xgboost, "Age")
pdplot(explainer_xgboost, catvars)
pdplot(explainer_xgboost, "Thal")

###################################### iml解释对象

predictor_xgboost <- iml::Predictor$new(
  final_xgboost, 
  data = traindatax,
  y = traindata$AHD,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>%
      rename_with(~gsub(".pred_", "", .x))
  },
  type = "prob"
)
# 交互作用
interact_xgboost <- iml::Interaction$new(predictor_xgboost)
plot(interact_xgboost) +
  theme_minimal()

interact_xgboost_1vo <- 
  iml::Interaction$new(predictor_xgboost, feature = "Age")
plot(interact_xgboost_1vo) +
  theme_minimal()

interact_xgboost_1v1 <- iml::FeatureEffect$new(
  predictor_xgboost, 
  feature = c("Age", "Oldpeak"),
  method = "pdp"
)
plot(interact_xgboost_1v1) +
  scale_fill_viridis_c() +
  labs(fill = "") +
  theme_minimal()

###################################### lime单样本预测分解

explainer_xgboost <- lime::lime(
  traindatax,
  lime::as_classifier(final_xgboost, c(yournegativelevel, yourpositivelevel))
)
explanation_xgboost <- lime::explain(
  traindatax[1,],  # 训练集第1个样本
  explainer_xgboost, 
  n_labels = 2, 
  n_features = ncol(traindatax)
)
lime::plot_features(explanation_xgboost)

######################## fastshap包

shapresult_xgboost <- shap4cls2(
  finalmodel = final_xgboost,
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
shapresult_xgboost$shapvipplot

# 单样本预测分解
shap41 <- shapviz::shapviz(
  shapresult_xgboost$shapley,
  X = traindatax
)
shapviz::sv_force(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
shapviz::sv_waterfall(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))

# 所有分类变量的shap图示
shapresult_xgboost$shapplotd_facet
shapresult_xgboost$shapplotd_one
# 所有连续变量的shap图示
shapresult_xgboost$shapplotc_facet
shapresult_xgboost$shapplotc_one
shapresult_xgboost$shapplotc_one2
# 单变量shap图示
sdplot(shapresult_xgboost, "Thal", "AHD")
sdplot(shapresult_xgboost, "Age", "AHD")

# 所有变量一张图
# shap变量重要性
shapresult_xgboost$shapvipplot_unity
# shap依赖图
shapresult_xgboost$shapplot_unity

#################################################################

# shap交互作用
traindatax2 <- final_xgboost %>%
  extract_recipe() %>%
  bake(new_data = traindata) %>%
  dplyr::select(-AHD)
colnames(traindatax2)
shapley_xgboost <- shapviz::shapviz(
  final_xgboost2,
  X_pred = as.matrix(traindatax2),
  X = traindatax2,
  interactions = T
)
# 所有变量，横轴对应shap
shapviz::sv_interaction(
  shapley_xgboost,
  max_display = ncol(traindatax)
)
# 指定变量，y轴对应shap interaction
shapviz::sv_dependence(
  shapley_xgboost, 
  v = c("Thal_normal", "Ca"),
  color_var = c("Slope_2", "Age"),
  interactions = T
)
# 指定变量，颜色对应shap interaction
shapviz::sv_dependence2D(
  shapley_xgboost, 
  x = c("Thal_normal", "Ca"),
  y = c("Slope_2", "Age"),
  interactions = T
)


#################################################################

# 学习曲线
lcN <- 
  floor(seq(nrow(traindata)%/%2, nrow(traindata), length = 10))
lcresult_xgboost <- data.frame()
for (i in lcN) {
  
  set.seed(i)
  traindatai <- traindata[sample(nrow(traindata), i), ]
  
  i_xgboost <-  wk_xgboost %>%
    finalize_workflow(hpbest_xgboost) %>%
    fit(traindatai)
  
  predtrain_i_xgboost <- traindatai %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_xgboost, new_data = traindatai, type = "prob")) %>%
    yardstick::roc_auc(AHD,
                       paste0(".pred_", yourpositivelevel), 
                       event_level = "second") %>%
    mutate(datasetname = paste("traindata", i, sep = "-"))
  
  predtest_i_xgboost <- testdata %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_xgboost, new_data = testdata, type = "prob")) %>%
    yardstick::roc_auc(AHD, 
                       paste0(".pred_", yourpositivelevel),
                       event_level = "second") %>%
    mutate(datasetname = paste("testdata", i, sep = "-"))
  
  predi <- bind_rows(predtrain_i_xgboost, predtest_i_xgboost)
  lcresult_xgboost <- rbind(lcresult_xgboost, predi)
  print(i)
}
# 图示
lcresult_xgboost %>%
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
save(datarecipe_xgboost,
     model_xgboost,
     wk_xgboost,
     hpgrid_xgboost,  # 如果采用贝叶斯优化则替换为 param_xgboost
     tune_xgboost,
     predtrain_xgboost,
     predtest_xgboost,
     evalcv_xgboost,
     vipdata_xgboost,
     file = "F:/Mach_learn_data/Decision Tree/evalresult_xgboost.RData")

# 保存模型结果供shiny部署之用，本课程不包括shiny内容
final_xgboost_heart <- final_xgboost
traindata_heart <- traindata
save(final_xgboost_heart,
     traindata_heart,
     file = "F:/Mach_learn_data/Decision Tree/shiny_xgboost_heart.RData")
