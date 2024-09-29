# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---随机森林
# https://parsnip.tidymodels.org/reference/details_rand_forest_randomForest.html

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
datarecipe_rf <- recipe(formula = AHD ~ ., traindata)
datarecipe_rf


# 设定模型
model_rf <- rand_forest(
  mode = "classification",
  engine = "randomForest", # ranger
  mtry = tune(),
  trees = tune(),
  min_n = tune()
) %>%
  set_args(importance = T)
model_rf

# workflow
wk_rf <- 
  workflow() %>%
  add_recipe(datarecipe_rf) %>%
  add_model(model_rf)
wk_rf

##############################################################

#########################  超参数寻优2选1-贝叶斯优化

# 更新超参数范围
param_rf <- model_rf %>%
  extract_parameter_set_dials() %>%
  update(mtry = mtry(c(2, 10)),
         trees = trees(c(100, 1000)),
         min_n = min_n(c(7, 55)))

# 贝叶斯优化超参数
set.seed(42)
tune_rf <- wk_rf %>%
  tune_bayes(
    resamples = folds,
    initial = 10,
    iter = 50,
    param_info = param_rf,
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
eval_tune_rf <- tune_rf %>%
  collect_metrics()
eval_tune_rf

# 图示
# autoplot(tune_rf)
eval_tune_rf %>% 
  filter(.metric == "roc_auc") %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~mean, colorscale = 'Jet', showscale = T),
    dimensions = list(
      list(label = 'mtry', values = ~mtry),
      list(label = 'trees', values = ~trees),
      list(label = 'min_n', values = ~min_n)
    )
  ) %>%
  plotly::layout(title = "RF HPO Guided by AUCROC",
                 font = list(family = "serif"))

# 经过交叉验证得到的最优超参数
hpbest_rf <- tune_rf %>%
  select_by_one_std_err(metric = "roc_auc", desc(min_n))
hpbest_rf

# 采用最优超参数组合训练最终模型
set.seed(42)
final_rf <- wk_rf %>%
  finalize_workflow(hpbest_rf) %>%
  fit(traindata)
final_rf

##################################################################

# 训练集预测评估
predtrain_rf <- eval4cls2(
  model = final_rf, 
  dataset = traindata, 
  yname = "AHD", 
  modelname = "RF", 
  datasetname = "traindata",
  cutoff = "yueden",
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtrain_rf$prediction
predtrain_rf$predprobplot
predtrain_rf$rocplot
predtrain_rf$prplot
predtrain_rf$caliplot
predtrain_rf$cmplot
predtrain_rf$metrics
predtrain_rf$diycutoff
predtrain_rf$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtrain_rf$proc)
pROC::ci.auc(predtrain_rf$proc)

# 预测评估测试集预测评估
predtest_rf <- eval4cls2(
  model = final_rf, 
  dataset = testdata, 
  yname = "AHD", 
  modelname = "RF", 
  datasetname = "testdata",
  cutoff = predtrain_rf$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtest_rf$prediction
predtest_rf$predprobplot
predtest_rf$rocplot
predtest_rf$prplot
predtest_rf$caliplot
predtest_rf$cmplot
predtest_rf$metrics
predtest_rf$diycutoff
predtest_rf$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtest_rf$proc)
pROC::ci.auc(predtest_rf$proc)

# ROC比较检验
pROC::roc.test(predtrain_rf$proc, predtest_rf$proc)


# 合并训练集和测试集上ROC曲线
predtrain_rf$rocresult %>%
  bind_rows(predtest_rf$rocresult) %>%
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
predtrain_rf$prresult %>%
  bind_rows(predtest_rf$prresult) %>%
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
predtrain_rf$caliresult %>%
  bind_rows(predtest_rf$caliresult) %>%
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
predtrain_rf$metrics %>%
  bind_rows(predtest_rf$metrics) %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)


# 最优超参数交叉验证的结果
evalcv_rf <- bestcv4cls2(
  wkflow = wk_rf,
  tuneresult = tune_rf,
  hpbest = hpbest_rf,
  yname = "AHD",
  modelname = "RF",
  v = 5,
  positivelevel = yourpositivelevel
)
evalcv_rf$cvroc
evalcv_rf$cvpr
evalcv_rf$evalcv


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
  bind_cols(predict(final_rf, new_data = newHeart, type = "prob"))%>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Yes >= predtrain_rf$diycutoff, 
             yourpositivelevel, 
             yournegativelevel)
    )
  )
# readr::write_excel_csv(predresult, "随机森林二分类预测结果.csv")
# 评估指标
prednew_rf <- eval4cls2(
  model = final_rf, 
  dataset = newHeart, 
  yname = "AHD", 
  modelname = "RF", 
  datasetname = "newHeart",
  cutoff = predtrain_rf$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
prednew_rf$prediction
prednew_rf$predprobplot
prednew_rf$rocplot
prednew_rf$prplot
prednew_rf$caliplot
prednew_rf$cmplot
prednew_rf$metrics
prednew_rf$diycutoff
prednew_rf$ksplot

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
final_rf2 <- final_rf %>%
  extract_fit_engine()
final_rf2

# 随机森林树的棵树与误差演变
plot(final_rf2, 
     main = "随机森林树的棵树与误差演变", 
     ylim = c(0, 1), 
     las = 1,
     lwd = 2,
     family = "serif")
legend("top", 
       legend = colnames(final_rf2$err.rate),
       lty = 1:3,
       lwd = 2,
       col = 1:3,
       horiz = T)

# 变量重要性
randomForest::importance(final_rf2)
randomForest::varImpPlot(
  final_rf2, 
  main = "变量重要性", 
  family = "serif"
)
# 偏依赖图
randomForest::partialPlot(
  final_rf2, 
  pred.data = as.data.frame(traindatax), 
  x.var = ChestPain,
  which.class = yourpositivelevel,
  las = 1,
  family = "serif"
)

######################## DALEX解释对象

explainer_rf <- DALEXtra::explain_tidymodels(
  final_rf, 
  data = traindatax,
  y = ifelse(traindata$AHD == yourpositivelevel, 1, 0),
  type = "classification",
  label = "RF"
)
# 变量重要性
vip_rf <- viplot(explainer_rf)
vipdata_rf <- vip_rf$data
vip_rf$plot

# 变量偏依赖图
pdplot(explainer_rf, convars)
pdplot(explainer_rf, "Age")
pdplot(explainer_rf, catvars)
pdplot(explainer_rf, "Thal")

###################################### iml解释对象

predictor_rf <- iml::Predictor$new(
  final_rf, 
  data = traindatax,
  y = traindata$AHD,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>%
      rename_with(~gsub(".pred_", "", .x))
  },
  type = "prob"
)
# 交互作用
interact_rf <- iml::Interaction$new(predictor_rf)
plot(interact_rf) +
  theme_minimal()

interact_rf_1vo <- 
  iml::Interaction$new(predictor_rf, feature = "Age")
plot(interact_rf_1vo) +
  theme_minimal()

interact_rf_1v1 <- iml::FeatureEffect$new(
  predictor_rf, 
  feature = c("Age", "Oldpeak"),
  method = "pdp"
)
plot(interact_rf_1v1) +
  scale_fill_viridis_c() +
  labs(fill = "") +
  theme_minimal()

###################################### lime单样本预测分解

explainer_rf <- lime::lime(
  traindatax,
  lime::as_classifier(final_rf, c(yournegativelevel, yourpositivelevel))
)
explanation_rf <- lime::explain(
  traindatax[1,],  # 训练集第1个样本
  explainer_rf, 
  n_labels = 2, 
  n_features = ncol(traindatax)
)
lime::plot_features(explanation_rf)

######################## fastshap包

shapresult_rf <- shap4cls2(
  finalmodel = final_rf,
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
shapresult_rf$shapvipplot

# 单样本预测分解
shap41 <- shapviz::shapviz(
  shapresult_rf$shapley,
  X = traindatax
)
shapviz::sv_force(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
shapviz::sv_waterfall(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))

# 所有分类变量的shap图示
shapresult_rf$shapplotd_facet
shapresult_rf$shapplotd_one
# 所有连续变量的shap图示
shapresult_rf$shapplotc_facet
shapresult_rf$shapplotc_one
shapresult_rf$shapplotc_one2
# 单变量shap图示
sdplot(shapresult_rf, "Thal", "AHD")
sdplot(shapresult_rf, "Age", "AHD")

# 所有变量一张图
# shap变量重要性
shapresult_rf$shapvipplot_unity
# shap依赖图
shapresult_rf$shapplot_unity

#################################################################

# shap交互作用
traindata3 <-traindata %>%
  mutate(AHD = factor(ifelse(AHD == yourpositivelevel, 1, 0)))
set.seed(42)
final_rf3 <- wk_rf %>%
  finalize_workflow(hpbest_rf) %>%
  fit(traindata3) %>%
  extract_fit_engine()


unify_rf <- treeshap::randomForest.unify(
  final_rf3,
  traindata3
)

set.seed(42)
treeshap_rf <- treeshap::treeshap(
  unify_rf, 
  traindatax,
  interactions = T
)
shapley_rf <- shapviz::shapviz(
  treeshap_rf,
  X = traindatax
)
# 所有变量，横轴对应shap
shapviz::sv_interaction(
  shapley_rf,
  max_display = ncol(traindatax)
)
# 指定变量，y轴对应shap interaction
shapviz::sv_dependence(
  shapley_rf, 
  v = c("Thal", "Ca"),
  color_var = c("Slope", "Age"),
  interactions = T
)
# 指定变量，颜色对应shap interaction
shapviz::sv_dependence2D(
  shapley_rf, 
  x = c("Thal", "Ca"),
  y = c("Slope", "Age"),
  interactions = T
)

#################################################################

# 学习曲线
lcN <- 
  floor(seq(nrow(traindata)%/%2, nrow(traindata), length = 10))
lcresult_rf <- data.frame()
for (i in lcN) {
  
  set.seed(i)
  traindatai <- traindata[sample(nrow(traindata), i), ]
  
  i_rf <-  wk_rf %>%
    finalize_workflow(hpbest_rf) %>%
    fit(traindatai)
  
  predtrain_i_rf <- traindatai %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_rf, new_data = traindatai, type = "prob")) %>%
    yardstick::roc_auc(AHD,
                       paste0(".pred_", yourpositivelevel), 
                       event_level = "second") %>%
    mutate(datasetname = paste("traindata", i, sep = "-"))
  
  predtest_i_rf <- testdata %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_rf, new_data = testdata, type = "prob")) %>%
    yardstick::roc_auc(AHD, 
                       paste0(".pred_", yourpositivelevel),
                       event_level = "second") %>%
    mutate(datasetname = paste("testdata", i, sep = "-"))
  
  predi <- bind_rows(predtrain_i_rf, predtest_i_rf)
  lcresult_rf <- rbind(lcresult_rf, predi)
  print(i)
}
# 图示
lcresult_rf %>%
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
save(datarecipe_rf,
     model_rf,
     wk_rf,
     param_rf,   # 如果采用贝叶斯优化则替换为 param_rf
     tune_rf,
     predtrain_rf,
     predtest_rf,
     evalcv_rf,
     vipdata_rf,
     file = "F:/Mach_learn_data/Decision Tree/evalresult_rf.RData")

# 保存模型结果供shiny部署之用，本课程不包括shiny内容
final_rf_heart <- final_rf
traindata_heart <- traindata
save(final_rf_heart,
     traindata_heart,
     file = "F:/Mach_learn_data/Decision Tree/shiny_rf_heart.RData")
