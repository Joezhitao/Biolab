# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---SVM
# https://parsnip.tidymodels.org/reference/details_svm_rbf_kernlab.html

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
datarecipe_svm <- recipe(formula = AHD ~ ., traindata) %>%
  step_dummy(all_nominal_predictors(), naming = new_dummy_names) %>%
  step_nzv(all_predictors()) %>% 
  step_normalize(all_predictors())
datarecipe_svm


# 设定模型
model_svm <- svm_rbf(
  mode = "classification",
  engine = "kernlab",
  cost = tune(),
  rbf_sigma = tune()
) %>%
  set_args(class.weights = c("No" = 1, "Yes" = 2)) 
# 此处NoYes是因变量的取值水平，换数据之后要相应更改
# 后面的数字12表示分类错误时的成本权重，无需设定时都等于1即可
model_svm

# workflow
wk_svm <- 
  workflow() %>%
  add_recipe(datarecipe_svm) %>%
  add_model(model_svm)
wk_svm

##############################################################
#########################  超参数寻优2选1-贝叶斯优化

# 贝叶斯优化超参数
set.seed(42)
tune_svm <- wk_svm %>%
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
eval_tune_svm <- tune_svm %>%
  collect_metrics()
eval_tune_svm

# 图示
# autoplot(tune_svm)
eval_tune_svm %>% 
  filter(.metric == "roc_auc") %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~mean, colorscale = 'Jet', showscale = T),
    dimensions = list(
      list(label = 'cost', values = ~cost),
      list(label = 'rbf_sigma', values = ~rbf_sigma,
           font = list(family = "serif"))
    )
  ) %>%
  plotly::layout(title = "SVM HPO Guided by AUCROC",
                 font = list(family = "serif"))

# 经过交叉验证得到的最优超参数
hpbest_svm <- tune_svm %>%
  select_best(metric = "roc_auc")
hpbest_svm

# 采用最优超参数组合训练最终模型
set.seed(42)
final_svm <- wk_svm %>%
  finalize_workflow(hpbest_svm) %>%
  fit(traindata)
final_svm

##################################################################

# 训练集预测评估
predtrain_svm <- eval4cls2(
  model = final_svm, 
  dataset = traindata, 
  yname = "AHD", 
  modelname = "SVM", 
  datasetname = "traindata",
  cutoff = "yueden",
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtrain_svm$prediction
predtrain_svm$predprobplot
predtrain_svm$rocplot
predtrain_svm$prplot
predtrain_svm$caliplot
predtrain_svm$cmplot
predtrain_svm$metrics
predtrain_svm$diycutoff
predtrain_svm$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtrain_svm$proc)
pROC::ci.auc(predtrain_svm$proc)

# 预测评估测试集预测评估
predtest_svm <- eval4cls2(
  model = final_svm, 
  dataset = testdata, 
  yname = "AHD", 
  modelname = "SVM", 
  datasetname = "testdata",
  cutoff = predtrain_svm$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtest_svm$prediction
predtest_svm$predprobplot
predtest_svm$rocplot
predtest_svm$prplot
predtest_svm$caliplot
predtest_svm$cmplot
predtest_svm$metrics
predtest_svm$diycutoff
predtest_svm$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtest_svm$proc)
pROC::ci.auc(predtest_svm$proc)

# ROC比较检验
pROC::roc.test(predtrain_svm$proc, predtest_svm$proc)


# 合并训练集和测试集上ROC曲线
predtrain_svm$rocresult %>%
  bind_rows(predtest_svm$rocresult) %>%
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
predtrain_svm$prresult %>%
  bind_rows(predtest_svm$prresult) %>%
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
predtrain_svm$caliresult %>%
  bind_rows(predtest_svm$caliresult) %>%
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
predtrain_svm$metrics %>%
  bind_rows(predtest_svm$metrics) %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)


# 最优超参数交叉验证的结果
evalcv_svm <- bestcv4cls2(
  wkflow = wk_svm,
  tuneresult = tune_svm,
  hpbest = hpbest_svm,
  yname = "AHD",
  modelname = "SVM",
  v = 5,
  positivelevel = yourpositivelevel
)
evalcv_svm$cvroc
evalcv_svm$cvpr
evalcv_svm$evalcv


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
  bind_cols(predict(final_svm, new_data = newHeart, type = "prob"))%>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Yes >= predtrain_svm$diycutoff, 
             yourpositivelevel, 
             yournegativelevel)
    )
  )
# readr::write_excel_csv(predresult, "SVM二分类预测结果.csv")
# 评估指标
prednew_svm <- eval4cls2(
  model = final_svm, 
  dataset = newHeart, 
  yname = "AHD", 
  modelname = "SVM", 
  datasetname = "newHeart",
  cutoff = predtrain_svm$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
prednew_svm$prediction
prednew_svm$predprobplot
prednew_svm$rocplot
prednew_svm$prplot
prednew_svm$caliplot
prednew_svm$cmplot
prednew_svm$metrics
prednew_svm$diycutoff
prednew_svm$ksplot

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
final_svm2 <- final_svm %>%
  extract_fit_engine()
final_svm2


######################## DALEX解释对象

explainer_svm <- DALEXtra::explain_tidymodels(
  final_svm, 
  data = traindatax,
  y = ifelse(traindata$AHD == yourpositivelevel, 1, 0),
  type = "classification",
  label = "SVM"
)
# 变量重要性
vip_svm <- viplot(explainer_svm)
vipdata_svm <- vip_svm$data
vip_svm$plot

# 变量偏依赖图
pdplot(explainer_svm, convars)
pdplot(explainer_svm, "Age")
pdplot(explainer_svm, catvars)
pdplot(explainer_svm, "Thal")

###################################### iml解释对象

predictor_svm <- iml::Predictor$new(
  final_svm, 
  data = traindatax,
  y = traindata$AHD,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>%
      rename_with(~gsub(".pred_", "", .x))
  },
  type = "prob"
)
# 交互作用
interact_svm <- iml::Interaction$new(predictor_svm)
plot(interact_svm) +
  theme_minimal()

interact_svm_1vo <- 
  iml::Interaction$new(predictor_svm, feature = "Age")
plot(interact_svm_1vo) +
  theme_minimal()

interact_svm_1v1 <- iml::FeatureEffect$new(
  predictor_svm, 
  feature = c("Age", "Oldpeak"),
  method = "pdp"
)
plot(interact_svm_1v1) +
  scale_fill_viridis_c() +
  labs(fill = "") +
  theme_minimal()

###################################### lime单样本预测分解

explainer_svm <- lime::lime(
  traindatax,
  lime::as_classifier(final_svm, c(yournegativelevel, yourpositivelevel))
)
explanation_svm <- lime::explain(
  traindatax[1,],  # 训练集第1个样本
  explainer_svm, 
  n_labels = 2, 
  n_features = ncol(traindatax)
)
lime::plot_features(explanation_svm)

######################## fastshap包

shapresult_svm <- shap4cls2(
  finalmodel = final_svm,
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
shapresult_svm$shapvipplot

# 单样本预测分解
shap41 <- shapviz::shapviz(
  shapresult_svm$shapley,
  X = traindatax
)
shapviz::sv_force(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
shapviz::sv_waterfall(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))

# 所有分类变量的shap图示
shapresult_svm$shapplotd_facet
shapresult_svm$shapplotd_one
# 所有连续变量的shap图示
shapresult_svm$shapplotc_facet
shapresult_svm$shapplotc_one
shapresult_svm$shapplotc_one2
# 单变量shap图示
sdplot(shapresult_svm, "Thal", "AHD")
sdplot(shapresult_svm, "Age", "AHD")

# 所有变量一张图
# shap变量重要性
shapresult_svm$shapvipplot_unity
# shap依赖图
shapresult_svm$shapplot_unity

#################################################################

# 学习曲线
lcN <- 
  floor(seq(nrow(traindata)%/%2, nrow(traindata), length = 10))
lcresult_svm <- data.frame()
for (i in lcN) {
  
  set.seed(i)
  traindatai <- traindata[sample(nrow(traindata), i), ]
  
  i_svm <-  wk_svm %>%
    finalize_workflow(hpbest_svm) %>%
    fit(traindatai)
  
  predtrain_i_svm <- traindatai %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_svm, new_data = traindatai, type = "prob")) %>%
    yardstick::roc_auc(AHD,
                       paste0(".pred_", yourpositivelevel), 
                       event_level = "second") %>%
    mutate(datasetname = paste("traindata", i, sep = "-"))
  
  predtest_i_svm <- testdata %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_svm, new_data = testdata, type = "prob")) %>%
    yardstick::roc_auc(AHD, 
                       paste0(".pred_", yourpositivelevel),
                       event_level = "second") %>%
    mutate(datasetname = paste("testdata", i, sep = "-"))
  
  predi <- bind_rows(predtrain_i_svm, predtest_i_svm)
  lcresult_svm <- rbind(lcresult_svm, predi)
  print(i)
}
# 图示
lcresult_svm %>%
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
save(datarecipe_svm,
     model_svm,
     wk_svm,# 如果采用贝叶斯优化则删除这一行
     tune_svm,
     predtrain_svm,
     predtest_svm,
     evalcv_svm,
     vipdata_svm,
     file = "F:/Mach_learn_data/Decision Tree/evalresult_svm.RData")

# 保存模型结果供shiny部署之用，本课程不包括shiny内容
final_svm_heart <- final_svm
traindata_heart <- traindata
save(final_svm_heart,
     traindata_heart,
     file = "F:/Mach_learn_data/Decision Tree/shiny_svm_heart.RData")
