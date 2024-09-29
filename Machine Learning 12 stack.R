# Joe---R语言tidymodels包机器学习分类与回归模型---二分类
# stacking
# 注意事项
# stacks包只能采用lasso作为第二层元模型，且第一层基础模型不能使用lightgbm
# stacks包的结果不具备可重复性，可能是因为多核并行，请注意保存结果
# 

##############################################################

# install.packages("tidymodels")
library(tidymodels)
library(bonsai)
source("H:/Biolab/Biolab/tidyfuncs4cls2_v18.R")

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

################################################################

# stacking模型中可能包括不同类型的模型，也可能包括相同类型的模型的不同配置

# 用于构建stacking模型的样本自变量值是候选基础模型的样本外预测结果

library(stacks)

##############################

# 也可以是之前单个模型建模的结果
load("F:/Mach_learn_data/Decision Tree/evalresult_knn.RData")
load("F:/Mach_learn_data/Decision Tree/evalresult_rf.RData")
load("F:/Mach_learn_data/Decision Tree/evalresult_logistic.RData")
models_stack <- 
  stacks() %>% 
  add_candidates(tune_knn) %>%
  add_candidates(tune_rf) %>%
  add_candidates(cv_logistic)
models_stack

##############################

# 拟合stacking元模型——stack
set.seed(42)
meta_stack <- blend_predictions(
  models_stack, 
  penalty = 10^seq(-2, -0.5, length = 20),
  control = control_grid(save_pred = T, 
                         verbose = T,
                         event_level = "second",
                         parallel_over = "everything",
                         save_workflow = T)
)
meta_stack
autoplot(meta_stack) +
  theme_bw() +
  theme(text = element_text(family = "serif"))

# 拟合选定的基础模型
set.seed(42)
final_stack <- fit_members(meta_stack)
final_stack

######################################################

# 应用stacking模型预测并评估
# 训练集
predtrain_stack <- eval4cls2(
  model = final_stack, 
  dataset = traindata, 
  yname = "AHD", 
  modelname = "Stacking", 
  datasetname = "traindata",
  cutoff = "yueden",
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtrain_stack$prediction
predtrain_stack$predprobplot
predtrain_stack$rocplot
predtrain_stack$prplot
predtrain_stack$caliplot
predtrain_stack$cmplot
predtrain_stack$metrics
predtrain_stack$diycutoff
predtrain_stack$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtrain_stack$proc)
pROC::ci.auc(predtrain_stack$proc)

# 测试集
predtest_stack <- eval4cls2(
  model = final_stack, 
  dataset = testdata, 
  yname = "AHD", 
  modelname = "Stacking", 
  datasetname = "testdata",
  cutoff = predtrain_stack$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtest_stack$prediction
predtest_stack$predprobplot
predtest_stack$rocplot
predtest_stack$prplot
predtest_stack$caliplot
predtest_stack$cmplot
predtest_stack$metrics
predtest_stack$diycutoff
predtest_stack$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtest_stack$proc)
pROC::ci.auc(predtest_stack$proc)

# ROC比较检验
pROC::roc.test(predtrain_stack$proc, predtest_stack$proc)


# 合并训练集和测试集上ROC曲线
predtrain_stack$rocresult %>%
  bind_rows(predtest_stack$rocresult) %>%
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
predtrain_stack$prresult %>%
  bind_rows(predtest_stack$prresult) %>%
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
predtrain_stack$caliresult %>%
  bind_rows(predtest_stack$caliresult) %>%
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
predtrain_stack$metrics %>%
  bind_rows(predtest_stack$metrics) %>%
  select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)

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
  bind_cols(predict(final_stack, new_data = newHeart, type = "prob"))%>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Yes >= predtrain_stack$diycutoff, 
             yourpositivelevel, 
             yournegativelevel)
    )
  )
# readr::write_excel_csv(predresult, "stacking二分类预测结果.csv")
# 评估指标
prednew_stack <- eval4cls2(
  model = final_stack, 
  dataset = newHeart, 
  yname = "AHD", 
  modelname = "Stacking", 
  datasetname = "newHeart",
  cutoff = predtrain_stack$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
prednew_stack$prediction
prednew_stack$predprobplot
prednew_stack$rocplot
prednew_stack$prplot
prednew_stack$caliplot
prednew_stack$cmplot
prednew_stack$metrics
prednew_stack$diycutoff
prednew_stack$ksplot

###################################################################

# 自变量数据集
colnames(traindata)
traindatax <- traindata %>%
  dplyr::select(-AHD)
colnames(traindatax)

# 分类型、连续型自变量名称
catvars <- getcategory(traindatax)
convars <- getcontinuous(traindatax)

######################## DALEX解释对象

explainer_stack <- DALEXtra::explain_tidymodels(
  final_stack, 
  data = traindatax,
  y = ifelse(traindata$AHD == yourpositivelevel, 1, 0),
  type = "classification",
  label = "Stacking"
)
# 变量重要性
vip_stack <- viplot(explainer_stack)
vipdata_stack <- vip_stack$data
vip_stack$plot

# 变量偏依赖图
pdplot(explainer_stack, convars)
pdplot(explainer_stack, "Age")
pdplot(explainer_stack, catvars)
pdplot(explainer_stack, "Thal")

###################################### iml解释对象

predictor_stack <- iml::Predictor$new(
  final_stack, 
  data = traindatax,
  y = traindata$AHD,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>%
      rename_with(~gsub(".pred_", "", .x))
  },
  type = "prob"
)
# 交互作用
interact_stack <- iml::Interaction$new(predictor_stack)
plot(interact_stack) +
  theme_minimal()

interact_stack_1vo <- 
  iml::Interaction$new(predictor_stack, feature = "Age")
plot(interact_stack_1vo) +
  theme_minimal()

interact_stack_1v1 <- iml::FeatureEffect$new(
  predictor_stack, 
  feature = c("Age", "Oldpeak"),
  method = "pdp"
)
plot(interact_stack_1v1) +
  scale_fill_viridis_c() +
  labs(fill = "") +
  theme_minimal()

###################################### lime单样本预测分解

explainer_stack <- lime::lime(
  traindatax,
  lime::as_classifier(final_stack, c(yournegativelevel, yourpositivelevel))
)
explanation_stack <- lime::explain(
  traindatax[1,],  # 训练集第1个样本
  explainer_stack, 
  n_labels = 2, 
  n_features = ncol(traindatax)
)
lime::plot_features(explanation_stack)

######################## fastshap包

shapresult_stack <- shap4cls2(
  finalmodel = final_stack,
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
shapresult_stack$shapvipplot

# 单样本预测分解
shap41 <- shapviz::shapviz(
  shapresult_stack$shapley,
  X = traindatax
)
shapviz::sv_force(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
shapviz::sv_waterfall(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))

# 所有分类变量的shap图示
shapresult_stack$shapplotd_facet
shapresult_stack$shapplotd_one
# 所有连续变量的shap图示
shapresult_stack$shapplotc_facet
shapresult_stack$shapplotc_one
shapresult_stack$shapplotc_one2
# 单变量shap图示
sdplot(shapresult_stack, "Thal", "AHD")
sdplot(shapresult_stack, "Age", "AHD")

# 所有变量一张图
# shap变量重要性
shapresult_stack$shapvipplot_unity
# shap依赖图
shapresult_stack$shapplot_unity

######################################################################

# 保存评估结果
save(predtrain_stack,
     predtest_stack,
     vipdata_stack,
     file = "F:/Mach_learn_data/Decision Tree/evalresult_stack.RData")

# 保存模型结果供shiny部署之用，本课程不包括shiny内容
final_stack_heart <- final_stack
traindata_heart <- traindata
save(final_stack_heart,
     traindata_heart,
     file = "F:/Mach_learn_data/Decision Tree/shiny_stack_heart.RData")
