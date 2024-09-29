# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---KNN
# https://parsnip.tidymodels.org/reference/details_nearest_neighbor_kknn.html

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
datarecipe_knn <- recipe(formula = AHD ~ ., traindata) %>%
  step_dummy(all_nominal_predictors(), naming = new_dummy_names) %>%
  step_nzv(all_predictors()) %>% 
  step_normalize(all_predictors())
datarecipe_knn


# 设定模型
model_knn <- nearest_neighbor(
  mode = "classification",
  engine = "kknn",
  
  neighbors = tune(),
  weight_func = tune(),
  dist_power = 2
)
model_knn

# workflow
wk_knn <- 
  workflow() %>%
  add_recipe(datarecipe_knn) %>%
  add_model(model_knn)
wk_knn

##############################################################
############################  超参数寻优2选1-网格搜索

# 超参数寻优网格
set.seed(42)
hpgrid_knn <- parameters(
  neighbors(range = c(3, 11)),
  weight_func()
) %>%
  # grid_regular(levels = c(5)) # 常规网格
  grid_random(size = 20) # 随机网格
# grid_latin_hypercube(size = 10) # 拉丁方网格
# grid_max_entropy(size = 10) # 最大熵网格
hpgrid_knn
# 网格也可以自己手动生成expand.grid()
# 交叉验证网格搜索过程
set.seed(42)
tune_knn <- wk_knn %>%
  tune_grid(
    resamples = folds,
    grid = hpgrid_knn,
    metrics = metricset_cls2,
    control = control_grid(save_pred = T, 
                           verbose = T,
                           event_level = "second",
                           parallel_over = "everything",
                           save_workflow = T)
  )

#########################  超参数寻优2选1-贝叶斯优化

# 更新超参数范围
param_knn <- model_knn %>%
  extract_parameter_set_dials() %>%
  update(neighbors = neighbors(c(5, 35)),
         weight_func = weight_func(c("rectangular",  "triangular")))

# 贝叶斯优化超参数
set.seed(42)
tune_knn <- wk_knn %>%
  tune_bayes(
    resamples = folds,
    initial = 20,
    iter = 50,
    param_info = param_knn,
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
eval_tune_knn <- tune_knn %>%
  collect_metrics()
eval_tune_knn

# 图示
# autoplot(tune_knn)
eval_tune_knn %>% 
  filter(.metric == "roc_auc") %>%
  mutate(weight_func2 = as.numeric(as.factor(weight_func))) %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~mean, colorscale = 'Jet', showscale = T),
    dimensions = list(
      list(label = 'neighbors', values = ~neighbors),
      list(label = 'weight_func', values = ~weight_func2,
           range = c(1,length(unique(eval_tune_knn$weight_func))), 
           tickvals = 1:length(unique(eval_tune_knn$weight_func)),
           ticktext = sort(unique(eval_tune_knn$weight_func)))
    )
  ) %>%
  plotly::layout(title = "KNN HPO Guided by AUCROC",
                 font = list(family = "serif"))

# 经过交叉验证得到的最优超参数
hpbest_knn <- tune_knn %>%
  select_by_one_std_err(metric = "roc_auc", desc(neighbors))
hpbest_knn

# 采用最优超参数组合训练最终模型
set.seed(42)
final_knn <- wk_knn %>%
  finalize_workflow(hpbest_knn) %>%
  fit(traindata)
final_knn

##################################################################

# 训练集预测评估
predtrain_knn <- eval4cls2(
  model = final_knn, 
  dataset = traindata, 
  yname = "AHD", 
  modelname = "KNN", 
  datasetname = "traindata",
  cutoff = "yueden",
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtrain_knn$prediction
predtrain_knn$predprobplot
predtrain_knn$rocplot
predtrain_knn$prplot
predtrain_knn$caliplot
predtrain_knn$cmplot
predtrain_knn$metrics
predtrain_knn$diycutoff
predtrain_knn$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtrain_knn$proc)
pROC::ci.auc(predtrain_knn$proc)

# 预测评估测试集预测评估
predtest_knn <- eval4cls2(
  model = final_knn, 
  dataset = testdata, 
  yname = "AHD", 
  modelname = "KNN", 
  datasetname = "testdata",
  cutoff = predtrain_knn$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtest_knn$prediction
predtest_knn$predprobplot
predtest_knn$rocplot
predtest_knn$prplot
predtest_knn$caliplot
predtest_knn$cmplot
predtest_knn$metrics
predtest_knn$diycutoff
predtest_knn$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtest_knn$proc)
pROC::ci.auc(predtest_knn$proc)

# ROC比较检验
pROC::roc.test(predtrain_knn$proc, predtest_knn$proc)


# 合并训练集和测试集上ROC曲线
predtrain_knn$rocresult %>%
  bind_rows(predtest_knn$rocresult) %>%
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
predtrain_knn$prresult %>%
  bind_rows(predtest_knn$prresult) %>%
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
predtrain_knn$caliresult %>%
  bind_rows(predtest_knn$caliresult) %>%
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
predtrain_knn$metrics %>%
  bind_rows(predtest_knn$metrics) %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)


# 最优超参数交叉验证的结果
evalcv_knn <- bestcv4cls2(
  wkflow = wk_knn,
  tuneresult = tune_knn,
  hpbest = hpbest_knn,
  yname = "AHD",
  modelname = "KNN",
  v = 5,
  positivelevel = yourpositivelevel
)
evalcv_knn$cvroc
evalcv_knn$cvpr
evalcv_knn$evalcv


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
  bind_cols(predict(final_knn, new_data = newHeart, type = "prob"))%>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Yes >= predtrain_knn$diycutoff, 
             yourpositivelevel, 
             yournegativelevel)
    )
  )
# readr::write_excel_csv(predresult, "KNN二分类预测结果.csv")
# 评估指标
prednew_knn <- eval4cls2(
  model = final_knn, 
  dataset = newHeart, 
  yname = "AHD", 
  modelname = "KNN", 
  datasetname = "newHeart",
  cutoff = predtrain_knn$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
prednew_knn$prediction
prednew_knn$predprobplot
prednew_knn$rocplot
prednew_knn$prplot
prednew_knn$caliplot
prednew_knn$cmplot
prednew_knn$metrics
prednew_knn$diycutoff
prednew_knn$ksplot

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
final_knn2 <- final_knn %>%
  extract_fit_engine()
final_knn2


######################## DALEX解释对象

explainer_knn <- DALEXtra::explain_tidymodels(
  final_knn, 
  data = traindatax,
  y = ifelse(traindata$AHD == yourpositivelevel, 1, 0),
  type = "classification",
  label = "KNN"
)
# 变量重要性
vip_knn <- viplot(explainer_knn)
vipdata_knn <- vip_knn$data
vip_knn$plot

# 变量偏依赖图
pdplot(explainer_knn, convars)
pdplot(explainer_knn, "Age")
pdplot(explainer_knn, catvars)
pdplot(explainer_knn, "Thal")

###################################### iml解释对象

predictor_knn <- iml::Predictor$new(
  final_knn, 
  data = traindatax,
  y = traindata$AHD,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>%
      rename_with(~gsub(".pred_", "", .x))
  },
  type = "prob"
)
# 交互作用
interact_knn <- iml::Interaction$new(predictor_knn)
plot(interact_knn) +
  theme_minimal()

interact_knn_1vo <- 
  iml::Interaction$new(predictor_knn, feature = "Age")
plot(interact_knn_1vo) +
  theme_minimal()

interact_knn_1v1 <- iml::FeatureEffect$new(
  predictor_knn, 
  feature = c("Age", "Oldpeak"),
  method = "pdp"
)
plot(interact_knn_1v1) +
  scale_fill_viridis_c() +
  labs(fill = "") +
  theme_minimal()

###################################### lime单样本预测分解

explainer_knn <- lime::lime(
  traindatax,
  lime::as_classifier(final_knn, c(yournegativelevel, yourpositivelevel))
)
explanation_knn <- lime::explain(
  traindatax[1,],  # 训练集第1个样本
  explainer_knn, 
  n_labels = 2, 
  n_features = ncol(traindatax)
)
lime::plot_features(explanation_knn)

######################## fastshap包

shapresult_knn <- shap4cls2(
  finalmodel = final_knn,
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
shapresult_knn$shapvipplot

# 单样本预测分解
shap41 <- shapviz::shapviz(
  shapresult_knn$shapley,
  X = traindatax
)
shapviz::sv_force(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
shapviz::sv_waterfall(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))

# 所有分类变量的shap图示
shapresult_knn$shapplotd_facet
shapresult_knn$shapplotd_one
# 所有连续变量的shap图示
shapresult_knn$shapplotc_facet
shapresult_knn$shapplotc_one
shapresult_knn$shapplotc_one2
# 单变量shap图示
sdplot(shapresult_knn, "Thal", "AHD")
sdplot(shapresult_knn, "Age", "AHD")

# 所有变量一张图
# shap变量重要性
shapresult_knn$shapvipplot_unity
# shap依赖图
shapresult_knn$shapplot_unity


#################################################################

# 学习曲线
lcN <- 
  floor(seq(nrow(traindata)%/%2, nrow(traindata), length = 10))
lcresult_knn <- data.frame()
for (i in lcN) {
  
  set.seed(i)
  traindatai <- traindata[sample(nrow(traindata), i), ]
  
  i_knn <-  wk_knn %>%
    finalize_workflow(hpbest_knn) %>%
    fit(traindatai)
  
  predtrain_i_knn <- traindatai %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_knn, new_data = traindatai, type = "prob")) %>%
    yardstick::roc_auc(AHD,
                       paste0(".pred_", yourpositivelevel), 
                       event_level = "second") %>%
    mutate(datasetname = paste("traindata", i, sep = "-"))
  
  predtest_i_knn <- testdata %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_knn, new_data = testdata, type = "prob")) %>%
    yardstick::roc_auc(AHD, 
                       paste0(".pred_", yourpositivelevel),
                       event_level = "second") %>%
    mutate(datasetname = paste("testdata", i, sep = "-"))
  
  predi <- bind_rows(predtrain_i_knn, predtest_i_knn)
  lcresult_knn <- rbind(lcresult_knn, predi)
  print(i)
}
# 图示
lcresult_knn %>%
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
save(datarecipe_knn,
     model_knn,
     wk_knn,
     tune_knn,
     predtrain_knn,
     predtest_knn,
     evalcv_knn,
     vipdata_knn,
     file = "F:/Mach_learn_data/Decision Tree/evalresult_knn.RData")

# 保存模型结果供shiny部署之用，本课程不包括shiny内容
final_knn_heart <- final_knn
traindata_heart <- traindata
save(final_knn_heart,
     traindata_heart,
     file = "F:/Mach_learn_data/Decision Tree/shiny_knn_heart.RData")
