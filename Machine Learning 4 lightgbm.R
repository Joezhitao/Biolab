# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---lightgbm
# https://parsnip.tidymodels.org/reference/details_boost_tree_lightgbm.html

##############################################################

# install.packages("tidymodels")
library(tidymodels)
library(bonsai)
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
datarecipe_lightgbm <- recipe(formula = AHD ~ ., traindata) %>%
  step_dummy(all_nominal_predictors(), naming = new_dummy_names) %>%
  step_nzv(all_predictors())
datarecipe_lightgbm


# 设定模型
model_lightgbm <- boost_tree(
  mode = "classification",
  engine = "lightgbm",
  tree_depth = tune(),
  trees = tune(),
  learn_rate = tune(),
  mtry = tune(),
  min_n = tune(),
  loss_reduction = tune()
)
model_lightgbm

# workflow
wk_lightgbm <- 
  workflow() %>%
  add_recipe(datarecipe_lightgbm) %>%
  add_model(model_lightgbm)
wk_lightgbm

##############################################################
############################  超参数寻优2选1-网格搜索

# 超参数寻优网格
set.seed(42)
hpgrid_lightgbm <- parameters(
  tree_depth(range = c(1, 3)),
  trees(range = c(100, 500)),
  learn_rate(range = c(-3, -1)),
  mtry(range = c(2, 8)),
  min_n(range = c(5, 10)),
  loss_reduction(range = c(-3, 0))
) %>%
  # grid_regular(levels = c(3, 2, 2, 3, 2, 2)) # 常规网格
  grid_random(size = 20) # 随机网格
# grid_latin_hypercube(size = 10) # 拉丁方网格
# grid_max_entropy(size = 10) # 最大熵网格
hpgrid_lightgbm
# 网格也可以自己手动生成expand.grid()

# 交叉验证网格搜索过程
set.seed(42)
tune_lightgbm <- wk_lightgbm %>%
  tune_grid(
    resamples = folds,
    grid = hpgrid_lightgbm,
    metrics = metricset_cls2,
    control = control_grid(save_pred = T, 
                           verbose = T,
                           event_level = "second",
                           parallel_over = "everything",
                           save_workflow = T)
  )

#########################  超参数寻优2选1-贝叶斯优化

# 更新超参数范围
param_lightgbm <- model_lightgbm %>%
  extract_parameter_set_dials() %>%
  update(mtry = mtry(c(2, 10)),
         min_n = min_n(c(15, 55)))

# 贝叶斯优化超参数
set.seed(42)
tune_lightgbm <- wk_lightgbm %>%
  tune_bayes(
    resamples = folds,
    initial = 10,
    iter = 50,
    param_info = param_lightgbm,
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
eval_tune_lightgbm <- tune_lightgbm %>%
  collect_metrics()
eval_tune_lightgbm

# 图示
# autoplot(tune_lightgbm)
eval_tune_lightgbm %>% 
  filter(.metric == "roc_auc") %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~mean, colorscale = 'Jet', showscale = T),
    dimensions = list(
      list(label = 'mtry', values = ~mtry),
      list(label = 'trees', values = ~trees),
      list(label = 'min_n', values = ~min_n),
      list(label = 'tree_depth', values = ~tree_depth),
      list(label = 'learn_rate', values = ~learn_rate),
      list(label = 'loss_reduction', values = ~loss_reduction)
    )
  ) %>%
  plotly::layout(title = "lightgbm HPO Guided by AUCROC",
                 font = list(family = "serif"))

# 经过交叉验证得到的最优超参数
hpbest_lightgbm <- tune_lightgbm %>%
  select_by_one_std_err(metric = "roc_auc", desc(min_n))
hpbest_lightgbm

# 采用最优超参数组合训练最终模型
set.seed(42)
final_lightgbm <- wk_lightgbm %>%
  finalize_workflow(hpbest_lightgbm) %>%
  fit(traindata)
final_lightgbm

##################################################################

# 训练集预测评估
predtrain_lightgbm <- eval4cls2(
  model = final_lightgbm, 
  dataset = traindata, 
  yname = "AHD", 
  modelname = "Lightgbm", 
  datasetname = "traindata",
  cutoff = "yueden",
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtrain_lightgbm$prediction
predtrain_lightgbm$predprobplot
predtrain_lightgbm$rocplot
predtrain_lightgbm$prplot
predtrain_lightgbm$caliplot
predtrain_lightgbm$cmplot
predtrain_lightgbm$metrics
predtrain_lightgbm$diycutoff
predtrain_lightgbm$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtrain_lightgbm$proc)
pROC::ci.auc(predtrain_lightgbm$proc)

# 预测评估测试集预测评估
predtest_lightgbm <- eval4cls2(
  model = final_lightgbm, 
  dataset = testdata, 
  yname = "AHD", 
  modelname = "Lightgbm", 
  datasetname = "testdata",
  cutoff = predtrain_lightgbm$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtest_lightgbm$prediction
predtest_lightgbm$predprobplot
predtest_lightgbm$rocplot
predtest_lightgbm$prplot
predtest_lightgbm$caliplot
predtest_lightgbm$cmplot
predtest_lightgbm$metrics
predtest_lightgbm$diycutoff
predtest_lightgbm$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtest_lightgbm$proc)
pROC::ci.auc(predtest_lightgbm$proc)

# ROC比较检验
pROC::roc.test(predtrain_lightgbm$proc, predtest_lightgbm$proc)


# 合并训练集和测试集上ROC曲线
predtrain_lightgbm$rocresult %>%
  bind_rows(predtest_lightgbm$rocresult) %>%
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
predtrain_lightgbm$prresult %>%
  bind_rows(predtest_lightgbm$prresult) %>%
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
predtrain_lightgbm$caliresult %>%
  bind_rows(predtest_lightgbm$caliresult) %>%
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
predtrain_lightgbm$metrics %>%
  bind_rows(predtest_lightgbm$metrics) %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)


# 最优超参数交叉验证的结果
evalcv_lightgbm <- bestcv4cls2(
  wkflow = wk_lightgbm,
  tuneresult = tune_lightgbm,
  hpbest = hpbest_lightgbm,
  yname = "AHD",
  modelname = "Lightgbm",
  v = 5,
  positivelevel = yourpositivelevel
)
evalcv_lightgbm$cvroc
evalcv_lightgbm$cvpr
evalcv_lightgbm$evalcv


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
  bind_cols(predict(final_lightgbm, new_data = newHeart, type = "prob"))%>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Yes >= predtrain_lightgbm$diycutoff, 
             yourpositivelevel, 
             yournegativelevel)
    )
  )
# readr::write_excel_csv(predresult, "Lightgbm二分类预测结果.csv")
# 评估指标
prednew_lightgbm <- eval4cls2(
  model = final_lightgbm, 
  dataset = newHeart, 
  yname = "AHD", 
  modelname = "Lightgbm", 
  datasetname = "newHeart",
  cutoff = predtrain_lightgbm$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
prednew_lightgbm$prediction
prednew_lightgbm$predprobplot
prednew_lightgbm$rocplot
prednew_lightgbm$prplot
prednew_lightgbm$caliplot
prednew_lightgbm$cmplot
prednew_lightgbm$metrics
prednew_lightgbm$diycutoff
prednew_lightgbm$ksplot

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
final_lightgbm2 <- final_lightgbm %>%
  extract_fit_engine()
final_lightgbm2

# 保存lightgbm模型比较特殊
model_file <- 
  tempfile(pattern = "Lightgbm", tmpdir = ".", fileext = ".txt")
lightgbm::lgb.save(final_lightgbm2, model_file)

# # 加载也需要自己的函数
# load_booster <- lightgbm::lgb.load(file.choose())

# 变量重要性
lightgbm::lgb.importance(final_lightgbm2, percentage = T)
lightgbm::lgb.plot.importance(
  lightgbm::lgb.importance(final_lightgbm2, percentage = T)
)

# 变量对预测的贡献
lightgbm::lgb.interprete(
  final_lightgbm2, 
  as.matrix(final_lightgbm %>%
              extract_recipe() %>%
              bake(new_data = traindata) %>%
              dplyr::select(-AHD)), 
  1:2
)
lightgbm::lgb.plot.interpretation(
  lightgbm::lgb.interprete(
    final_lightgbm2, 
    as.matrix(final_lightgbm %>%
                extract_recipe() %>%
                bake(new_data = traindata) %>%
                dplyr::select(-AHD)),
    1:2
  )[[1]]
)

######################## DALEX解释对象

explainer_lightgbm <- DALEXtra::explain_tidymodels(
  final_lightgbm, 
  data = traindatax,
  y = ifelse(traindata$AHD == yourpositivelevel, 1, 0),
  type = "classification",
  label = "Lightgbm"
)
# 变量重要性
vip_lightgbm <- viplot(explainer_lightgbm)
vipdata_lightgbm <- vip_lightgbm$data
vip_lightgbm$plot

# 变量偏依赖图
pdplot(explainer_lightgbm, convars)
pdplot(explainer_lightgbm, "Age")
pdplot(explainer_lightgbm, catvars)
pdplot(explainer_lightgbm, "Thal")

###################################### iml解释对象

predictor_lightgbm <- iml::Predictor$new(
  final_lightgbm, 
  data = traindatax,
  y = traindata$AHD,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>%
      rename_with(~gsub(".pred_", "", .x))
  },
  type = "prob"
)
# 交互作用
interact_lightgbm <- iml::Interaction$new(predictor_lightgbm)
plot(interact_lightgbm) +
  theme_minimal()

interact_lightgbm_1vo <- 
  iml::Interaction$new(predictor_lightgbm, feature = "Age")
plot(interact_lightgbm_1vo) +
  theme_minimal()

interact_lightgbm_1v1 <- iml::FeatureEffect$new(
  predictor_lightgbm, 
  feature = c("Age", "Oldpeak"),
  method = "pdp"
)
plot(interact_lightgbm_1v1) +
  scale_fill_viridis_c() +
  labs(fill = "") +
  theme_minimal()

###################################### lime单样本预测分解

explainer_lightgbm <- lime::lime(
  traindatax,
  lime::as_classifier(final_lightgbm, c(yournegativelevel, yourpositivelevel))
)
explanation_lightgbm <- lime::explain(
  traindatax[1,],  # 训练集第1个样本
  explainer_lightgbm, 
  n_labels = 2, 
  n_features = ncol(traindatax)
)
lime::plot_features(explanation_lightgbm)

######################## fastshap包

shapresult_lightgbm <- shap4cls2(
  finalmodel = final_lightgbm,
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
shapresult_lightgbm$shapvipplot

# 单样本预测分解
shap41 <- shapviz::shapviz(
  shapresult_lightgbm$shapley,
  X = traindatax
)
shapviz::sv_force(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
shapviz::sv_waterfall(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))

# 所有分类变量的shap图示
shapresult_lightgbm$shapplotd_facet
shapresult_lightgbm$shapplotd_one
# 所有连续变量的shap图示
shapresult_lightgbm$shapplotc_facet
shapresult_lightgbm$shapplotc_one
shapresult_lightgbm$shapplotc_one2
# 单变量shap图示
sdplot(shapresult_lightgbm, "Thal", "AHD")
sdplot(shapresult_lightgbm, "Age", "AHD")

# 所有变量一张图
# shap变量重要性
shapresult_lightgbm$shapvipplot_unity
# shap依赖图
shapresult_lightgbm$shapplot_unity

#################################################################

# shap交互作用
traindata2 <- final_lightgbm %>%
  extract_recipe() %>%
  bake(new_data = traindata)
colnames(traindata2)
unify_lightgbm <- treeshap::lightgbm.unify(
  final_lightgbm2,
  traindata2
)
set.seed(42)
treeshap_lightgbm <- treeshap::treeshap(
  unify_lightgbm, 
  traindata2 %>%
    dplyr::select(-AHD),
  interactions = T
)
shapley_lightgbm <- shapviz::shapviz(
  treeshap_lightgbm,
  X = traindata2 %>%
    dplyr::select(-AHD)
)
# 所有变量，横轴对应shap
shapviz::sv_interaction(
  shapley_lightgbm,
  max_display = ncol(traindatax)
) +
  theme_bw()
# 指定变量，y轴对应shap interaction
shapviz::sv_dependence(
  shapley_lightgbm, 
  v = c("Thal_normal", "Ca"),
  color_var = c("Slope_2", "Age"),
  interactions = T
)
# 指定变量，颜色对应shap interaction
shapviz::sv_dependence2D(
  shapley_lightgbm, 
  x = c("Thal_normal", "Ca"),
  y = c("Slope_2", "Age"),
  interactions = T
)

#################################################################

# 学习曲线
lcN <- 
  floor(seq(nrow(traindata)%/%2, nrow(traindata), length = 10))
lcresult_lightgbm <- data.frame()
for (i in lcN) {
  
  set.seed(i)
  traindatai <- traindata[sample(nrow(traindata), i), ]
  
  i_lightgbm <-  wk_lightgbm %>%
    finalize_workflow(hpbest_lightgbm) %>%
    fit(traindatai)
  
  predtrain_i_lightgbm <- traindatai %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_lightgbm, new_data = traindatai, type = "prob")) %>%
    yardstick::roc_auc(AHD,
                       paste0(".pred_", yourpositivelevel), 
                       event_level = "second") %>%
    mutate(datasetname = paste("traindata", i, sep = "-"))
  
  predtest_i_lightgbm <- testdata %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_lightgbm, new_data = testdata, type = "prob")) %>%
    yardstick::roc_auc(AHD, 
                       paste0(".pred_", yourpositivelevel),
                       event_level = "second") %>%
    mutate(datasetname = paste("testdata", i, sep = "-"))
  
  predi <- bind_rows(predtrain_i_lightgbm, predtest_i_lightgbm)
  lcresult_lightgbm <- rbind(lcresult_lightgbm, predi)
  print(i)
}
# 图示
lcresult_lightgbm %>%
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
save(datarecipe_lightgbm,
     model_lightgbm,
     wk_lightgbm,
     param_lightgbm, # 如果采用贝叶斯优化则替换为 param_lightgbm
     tune_lightgbm,
     predtrain_lightgbm,
     predtest_lightgbm,
     evalcv_lightgbm,
     vipdata_lightgbm,
     file = "F:/Mach_learn_data/Decision Tree/evalresult_lightgbm.RData")

# 保存模型结果供shiny部署之用，本课程不包括shiny内容
final_lightgbm_heart <- final_lightgbm
traindata_heart <- traindata
save(final_lightgbm_heart,
     traindata_heart,
     file = "F:/Mach_learn_data/Decision Tree/shiny_lightgbm_heart.RData")
