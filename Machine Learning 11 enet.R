# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---enet岭回归弹性网络
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
datarecipe_enet <- recipe(formula = AHD ~ ., traindata) %>%
  step_dummy(all_nominal_predictors(), naming = new_dummy_names) %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_predictors())
datarecipe_enet


# 设定模型
model_enet <- logistic_reg(
  mode = "classification",
  engine = "glmnet",
  mixture = tune(),
  penalty = tune()
)
model_enet

# workflow
wk_enet <- 
  workflow() %>%
  add_recipe(datarecipe_enet) %>%
  add_model(model_enet)
wk_enet

##############################################################
#########################  超参数寻优2选1-贝叶斯优化

# 贝叶斯优化超参数
set.seed(42)
tune_enet <- wk_enet %>%
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
eval_tune_enet <- tune_enet %>%
  collect_metrics()
eval_tune_enet

# 图示
# autoplot(tune_enet)
eval_tune_enet %>% 
  filter(.metric == "roc_auc") %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~mean, colorscale = 'Jet', showscale = T),
    dimensions = list(
      list(label = 'mixture', values = ~mixture),
      list(label = 'penalty', values = ~penalty)
    )
  ) %>%
  plotly::layout(title = "ENet HPO Guided by AUCROC",
                 font = list(family = "serif"))

# 经过交叉验证得到的最优超参数
hpbest_enet <- tune_enet %>%
  select_by_one_std_err(metric = "roc_auc", desc(penalty))
hpbest_enet

# 采用最优超参数组合训练最终模型
set.seed(42)
final_enet <- wk_enet %>%
  finalize_workflow(hpbest_enet) %>%
  fit(traindata)
final_enet

##################################################################

# 训练集预测评估
predtrain_enet <- eval4cls2(
  model = final_enet, 
  dataset = traindata, 
  yname = "AHD", 
  modelname = "ENet", 
  datasetname = "traindata",
  cutoff = "yueden",
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtrain_enet$prediction
predtrain_enet$predprobplot
predtrain_enet$rocplot
predtrain_enet$prplot
predtrain_enet$caliplot
predtrain_enet$cmplot
predtrain_enet$metrics
predtrain_enet$diycutoff
predtrain_enet$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtrain_enet$proc)
pROC::ci.auc(predtrain_enet$proc)

# 预测评估测试集预测评估
predtest_enet <- eval4cls2(
  model = final_enet, 
  dataset = testdata, 
  yname = "AHD", 
  modelname = "ENet", 
  datasetname = "testdata",
  cutoff = predtrain_enet$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtest_enet$prediction
predtest_enet$predprobplot
predtest_enet$rocplot
predtest_enet$prplot
predtest_enet$caliplot
predtest_enet$cmplot
predtest_enet$metrics
predtest_enet$diycutoff
predtest_enet$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtest_enet$proc)
pROC::ci.auc(predtest_enet$proc)

# ROC比较检验
pROC::roc.test(predtrain_enet$proc, predtest_enet$proc)


# 合并训练集和测试集上ROC曲线
predtrain_enet$rocresult %>%
  bind_rows(predtest_enet$rocresult) %>%
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
predtrain_enet$prresult %>%
  bind_rows(predtest_enet$prresult) %>%
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
predtrain_enet$caliresult %>%
  bind_rows(predtest_enet$caliresult) %>%
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
predtrain_enet$metrics %>%
  bind_rows(predtest_enet$metrics) %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)


# 最优超参数交叉验证的结果
evalcv_enet <- bestcv4cls2(
  wkflow = wk_enet,
  tuneresult = tune_enet,
  hpbest = hpbest_enet,
  yname = "AHD",
  modelname = "ENet",
  v = 5,
  positivelevel = yourpositivelevel
)
evalcv_enet$cvroc
evalcv_enet$cvpr
evalcv_enet$evalcv


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
  bind_cols(predict(final_enet, new_data = newHeart, type = "prob"))%>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Yes >= predtrain_enet$diycutoff, 
             yourpositivelevel, 
             yournegativelevel)
    )
  )
# readr::write_excel_csv(predresult, "弹性网络二分类预测结果.csv")
# 评估指标
prednew_enet <- eval4cls2(
  model = final_enet, 
  dataset = newHeart, 
  yname = "AHD", 
  modelname = "ENet", 
  datasetname = "newHeart",
  cutoff = predtrain_enet$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
prednew_enet$prediction
prednew_enet$predprobplot
prednew_enet$rocplot
prednew_enet$prplot
prednew_enet$caliplot
prednew_enet$cmplot
prednew_enet$metrics
prednew_enet$diycutoff
prednew_enet$ksplot

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
final_enet2 <- final_enet %>%
  extract_fit_engine()
final_enet2

# 非零系数自变量
tidy(final_enet) %>%
  filter(term != "(Intercept)", estimate != 0) %>%
  pull(term)

######################## DALEX解释对象

explainer_enet <- DALEXtra::explain_tidymodels(
  final_enet, 
  data = traindatax,
  y = ifelse(traindata$AHD == yourpositivelevel, 1, 0),
  type = "classification",
  label = "ENet"
)
# 变量重要性
vip_enet <- viplot(explainer_enet)
vipdata_enet <- vip_enet$data
vip_enet$plot

# 变量偏依赖图
pdplot(explainer_enet, convars)
pdplot(explainer_enet, "Age")
pdplot(explainer_enet, catvars)
pdplot(explainer_enet, "Thal")

###################################### iml解释对象

predictor_enet <- iml::Predictor$new(
  final_enet, 
  data = traindatax,
  y = traindata$AHD,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>%
      rename_with(~gsub(".pred_", "", .x))
  },
  type = "prob"
)
# 交互作用
interact_enet <- iml::Interaction$new(predictor_enet)
plot(interact_enet) +
  theme_minimal()

interact_enet_1vo <- 
  iml::Interaction$new(predictor_enet, feature = "Age")
plot(interact_enet_1vo) +
  theme_minimal()

interact_enet_1v1 <- iml::FeatureEffect$new(
  predictor_enet, 
  feature = c("Age", "Oldpeak"),
  method = "pdp"
)
plot(interact_enet_1v1) +
  scale_fill_viridis_c() +
  labs(fill = "") +
  theme_minimal()

###################################### lime单样本预测分解

explainer_enet <- lime::lime(
  traindatax,
  lime::as_classifier(final_enet, c(yournegativelevel, yourpositivelevel))
)
explanation_enet <- lime::explain(
  traindatax[1,],  # 训练集第1个样本
  explainer_enet, 
  n_labels = 2, 
  n_features = ncol(traindatax)
)
lime::plot_features(explanation_enet)

######################## fastshap包

shapresult_enet <- shap4cls2(
  finalmodel = final_enet,
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
shapresult_enet$shapvipplot

# 单样本预测分解
shap41 <- shapviz::shapviz(
  shapresult_enet$shapley,
  X = traindatax
)
shapviz::sv_force(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
shapviz::sv_waterfall(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))

# 所有分类变量的shap图示
shapresult_enet$shapplotd_facet
shapresult_enet$shapplotd_one
# 所有连续变量的shap图示
shapresult_enet$shapplotc_facet
shapresult_enet$shapplotc_one
shapresult_enet$shapplotc_one2
# 单变量shap图示
sdplot(shapresult_enet, "Thal", "AHD")
sdplot(shapresult_enet, "Age", "AHD")

# 所有变量一张图
# shap变量重要性
shapresult_enet$shapvipplot_unity
# shap依赖图
shapresult_enet$shapplot_unity


#################################################################

# 学习曲线
lcN <- 
  floor(seq(nrow(traindata)%/%2, nrow(traindata), length = 10))
lcresult_enet <- data.frame()
for (i in lcN) {
  
  set.seed(i)
  traindatai <- traindata[sample(nrow(traindata), i), ]
  
  i_enet <-  wk_enet %>%
    finalize_workflow(hpbest_enet) %>%
    fit(traindatai)
  
  predtrain_i_enet <- traindatai %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_enet, new_data = traindatai, type = "prob")) %>%
    yardstick::roc_auc(AHD,
                       paste0(".pred_", yourpositivelevel), 
                       event_level = "second") %>%
    mutate(datasetname = paste("traindata", i, sep = "-"))
  
  predtest_i_enet <- testdata %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_enet, new_data = testdata, type = "prob")) %>%
    yardstick::roc_auc(AHD, 
                       paste0(".pred_", yourpositivelevel),
                       event_level = "second") %>%
    mutate(datasetname = paste("testdata", i, sep = "-"))
  
  predi <- bind_rows(predtrain_i_enet, predtest_i_enet)
  lcresult_enet <- rbind(lcresult_enet, predi)
  print(i)
}
# 图示
lcresult_enet %>%
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
save(datarecipe_enet,
     model_enet,
     wk_enet,
     tune_enet,
     predtrain_enet,
     predtest_enet,
     evalcv_enet,
     vipdata_enet,
     file = "F:/Mach_learn_data/Decision Tree/evalresult_enet.RData")

# 保存模型结果供shiny部署之用，本课程不包括shiny内容
final_enet_heart <- final_enet
traindata_heart <- traindata
save(final_enet_heart,
     traindata_heart,
     file = "F:/Mach_learn_data/Decision Tree/shiny_enet_heart.RData")
