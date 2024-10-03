# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---logistic回归
# https://parsnip.tidymodels.org/reference/details_logistic_reg_glm.html

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
datarecipe_logistic <- recipe(AHD ~ ., traindata)
datarecipe_logistic


# 设定模型
model_logistic <- logistic_reg(
  mode = "classification",
  engine = "glm"
)
model_logistic

# workflow
wk_logistic <- 
  workflow() %>%
  add_recipe(datarecipe_logistic) %>%
  add_model(model_logistic)
wk_logistic

# 训练模型
set.seed(42)
final_logistic <- wk_logistic %>%
  fit(traindata)
final_logistic

##################################################################

# 训练集预测评估
predtrain_logistic <- eval4cls2(
  model = final_logistic, 
  dataset = traindata, 
  yname = "AHD", 
  modelname = "Logistic", 
  datasetname = "traindata",
  cutoff = "yueden",
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtrain_logistic$prediction
predtrain_logistic$predprobplot
predtrain_logistic$rocplot
predtrain_logistic$prplot
predtrain_logistic$caliplot
predtrain_logistic$cmplot
predtrain_logistic$metrics
predtrain_logistic$diycutoff
predtrain_logistic$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtrain_logistic$proc)
pROC::ci.auc(predtrain_logistic$proc)

# 预测评估测试集预测评估
predtest_logistic <- eval4cls2(
  model = final_logistic, 
  dataset = testdata, 
  yname = "AHD", 
  modelname = "Logistic", 
  datasetname = "testdata",
  cutoff = predtrain_logistic$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtest_logistic$prediction
predtest_logistic$predprobplot
predtest_logistic$rocplot
predtest_logistic$prplot
predtest_logistic$caliplot
predtest_logistic$cmplot
predtest_logistic$metrics
predtest_logistic$diycutoff
predtest_logistic$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtest_logistic$proc)
pROC::ci.auc(predtest_logistic$proc)

# ROC比较检验
pROC::roc.test(predtrain_logistic$proc, predtest_logistic$proc)


# 合并训练集和测试集上ROC曲线
predtrain_logistic$rocresult %>%
  bind_rows(predtest_logistic$rocresult) %>%
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
predtrain_logistic$prresult %>%
  bind_rows(predtest_logistic$prresult) %>%
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
predtrain_logistic$caliresult %>%
  bind_rows(predtest_logistic$caliresult) %>%
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
predtrain_logistic$metrics %>%
  bind_rows(predtest_logistic$metrics) %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)

##################################################################

# 交叉验证
set.seed(42)
cv_logistic <- 
  wk_logistic %>%
  fit_resamples(
    folds,
    metrics = metricset_cls2,
    control = control_resamples(save_pred = T,
                                verbose = T,
                                event_level = "second",
                                parallel_over = "everything",
                                save_workflow = T)
  )
cv_logistic

# 交叉验证指标结果
evalcv_logistic <- list()
# 评估指标设定
metrictemp <- metric_set(yardstick::roc_auc, yardstick::pr_auc)
evalcv_logistic$evalcv <- 
  collect_predictions(cv_logistic) %>%
  group_by(id) %>%
  metrictemp(AHD, .pred_Yes, event_level = "second") %>%
  group_by(.metric) %>%
  mutate(model = "logistic",
         mean = mean(.estimate),
         sd = sd(.estimate)/sqrt(length(folds$splits)))
evalcv_logistic$evalcv

# 交叉验证预测结果图示
# ROC
evalcv_logistic$cvroc <- 
  collect_predictions(cv_logistic) %>%
  group_by(id) %>%
  roc_curve(AHD, .pred_Yes, event_level = "second") %>%
  ungroup() %>%
  left_join(evalcv_logistic$evalcv %>% filter(.metric == "roc_auc"), 
            by = "id") %>%
  mutate(idAUC = paste(id, " ROCAUC:", round(.estimate, 4)),
         idAUC = forcats::as_factor(idAUC)) %>%
  ggplot(aes(x = 1-specificity, y = sensitivity, color = idAUC)) +
  geom_path(linewidth = 1) +
  geom_abline(linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  labs(color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        text = element_text(family = "serif"))
evalcv_logistic$cvroc

# PR
evalcv_logistic$cvpr <- 
  collect_predictions(cv_logistic) %>%
  group_by(id) %>%
  pr_curve(AHD, .pred_Yes, event_level = "second") %>%
  ungroup() %>%
  left_join(evalcv_logistic$evalcv %>% filter(.metric == "pr_auc"), 
            by = "id") %>%
  mutate(idAUC = paste(id, " PRAUC:", round(.estimate, 4)),
         idAUC = forcats::as_factor(idAUC)) %>%
  ggplot(aes(x = recall, y = precision, color = idAUC)) +
  geom_path(linewidth = 1) +
  geom_abline(linetype = "dashed", intercept = 1, slope = -1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  labs(color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        text = element_text(family = "serif"))
evalcv_logistic$cvpr


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
  bind_cols(predict(final_logistic, new_data = newHeart, type = "prob")) %>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Yes >= predtrain_logistic$diycutoff, 
             yourpositivelevel, 
             yournegativelevel)
    )
  )
# readr::write_excel_csv(predresult, "Logistic二分类预测结果.csv")
# 评估指标
prednew_logistic <- eval4cls2(
  model = final_logistic, 
  dataset = newHeart, 
  yname = "AHD", 
  modelname = "Logistic", 
  datasetname = "newHeart",
  cutoff = predtrain_logistic$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
prednew_logistic$prediction
prednew_logistic$predprobplot
prednew_logistic$rocplot
prednew_logistic$prplot
prednew_logistic$caliplot
prednew_logistic$cmplot
prednew_logistic$metrics
prednew_logistic$diycutoff
prednew_logistic$ksplot

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
final_logistic2 <- final_logistic %>%
  extract_fit_engine()
final_logistic2


######################## DALEX解释对象

explainer_logistic <- DALEXtra::explain_tidymodels(
  final_logistic, 
  data = traindatax,
  y = ifelse(traindata$AHD == yourpositivelevel, 1, 0),
  type = "classification",
  label = "Logistic"
)
# 变量重要性
vip_logistic <- viplot(explainer_logistic)
vipdata_logistic <- vip_logistic$data
vip_logistic$plot

# 变量偏依赖图
pdplot(explainer_logistic, convars)
pdplot(explainer_logistic, "Age")
pdplot(explainer_logistic, catvars)
pdplot(explainer_logistic, "Thal")

###################################### iml解释对象

predictor_logistic <- iml::Predictor$new(
  final_logistic, 
  data = traindatax,
  y = traindata$AHD,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>%
      rename_with(~gsub(".pred_", "", .x))
  },
  type = "prob"
)
# 交互作用
interact_logistic <- iml::Interaction$new(predictor_logistic)
plot(interact_logistic) +
  theme_minimal()

interact_logistic_1vo <- 
  iml::Interaction$new(predictor_logistic, feature = "Age")
plot(interact_logistic_1vo) +
  theme_minimal()

interact_logistic_1v1 <- iml::FeatureEffect$new(
  predictor_logistic, 
  feature = c("Age", "Oldpeak"),
  method = "pdp"
)
plot(interact_logistic_1v1) +
  scale_fill_viridis_c() +
  labs(fill = "") +
  theme_minimal()

###################################### lime单样本预测分解

explainer_logistic <- lime::lime(
  traindatax,
  lime::as_classifier(final_logistic, c(yournegativelevel, yourpositivelevel))
)
explanation_logistic <- lime::explain(
  traindatax[1,],  # 训练集第1个样本
  explainer_logistic, 
  n_labels = 2, 
  n_features = ncol(traindatax)
)
lime::plot_features(explanation_logistic)

######################## fastshap包

shapresult_logistic <- shap4cls2(
  finalmodel = final_logistic,
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
shapresult_logistic$shapvipplot

# 单样本预测分解
shap41 <- shapviz::shapviz(
  shapresult_logistic$shapley,
  X = traindatax
)
shapviz::sv_force(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
shapviz::sv_waterfall(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
# 所有分类变量的shap图示
shapresult_logistic$shapplotd_facet
shapresult_logistic$shapplotd_one
# 所有连续变量的shap图示
shapresult_logistic$shapplotc_facet
shapresult_logistic$shapplotc_one
shapresult_logistic$shapplotc_one2
# 单变量shap图示
sdplot(shapresult_logistic, "Thal", "AHD")
sdplot(shapresult_logistic, "Age", "AHD")

# 所有变量一张图
# shap变量重要性
shapresult_logistic$shapvipplot_unity
# shap依赖图
shapresult_logistic$shapplot_unity


#################################################################

# 学习曲线
lcN <- 
  floor(seq(nrow(traindata)%/%2, nrow(traindata), length = 10))
lcresult_logistic <- data.frame()
for (i in lcN) {
  
  set.seed(i)
  traindatai <- traindata[sample(nrow(traindata), i), ]
  
  i_logistic <-  wk_logistic %>%
    fit(traindatai)
  
  predtrain_i_logistic <- traindatai %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_logistic, new_data = traindatai, type = "prob")) %>%
    yardstick::roc_auc(AHD,
                       paste0(".pred_", yourpositivelevel), 
                       event_level = "second") %>%
    mutate(datasetname = paste("traindata", i, sep = "-"))
  
  predtest_i_logistic <- testdata %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_logistic, new_data = testdata, type = "prob")) %>%
    yardstick::roc_auc(AHD, 
                       paste0(".pred_", yourpositivelevel),
                       event_level = "second") %>%
    mutate(datasetname = paste("testdata", i, sep = "-"))
  
  predi <- bind_rows(predtrain_i_logistic, predtest_i_logistic)
  lcresult_logistic <- rbind(lcresult_logistic, predi)
  print(i)
}
# 图示
lcresult_logistic %>%
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
save(datarecipe_logistic,
     model_logistic,
     wk_logistic,
     cv_logistic,
     predtrain_logistic,
     predtest_logistic,
     evalcv_logistic,
     vipdata_logistic,
     file = "F:/Mach_learn_data/Decision Tree/evalresult_logistic.RData")


# 保存模型结果供shiny部署之用，本课程不包括shiny内容
final_logistic_heart <- final_logistic
traindata_heart <- traindata
save(final_logistic_heart,
     traindata_heart,
     file = "F:/Mach_learn_data/Decision Tree/shiny_logistic_heart.RData")
