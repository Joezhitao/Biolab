# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---决策树
# https://parsnip.tidymodels.org/reference/details_decision_tree_rpart.html
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
datarecipe_dt <- recipe(formula = AHD ~ ., traindata)
datarecipe_dt


# 设定模型
model_dt <- decision_tree(
  mode = "classification",
  engine = "rpart",
  tree_depth = tune(),
  min_n = tune(),
  cost_complexity = tune()
) %>%
  set_args(model=T)
model_dt

# workflow
wk_dt <- 
  workflow() %>%
  add_recipe(datarecipe_dt) %>%
  add_model(model_dt)
wk_dt

##############################################################
############################  超参数寻优2选1-网格搜索

# 超参数寻优网格
set.seed(42)
hpgrid_dt <- parameters(
  tree_depth(range = c(3, 7)),
  min_n(range = c(5, 10)),
  cost_complexity(range = c(-6, -3))
) %>%
  # grid_regular(levels = c(3, 2, 4)) # 常规网格
  grid_random(size = 20) # 随机网格
# grid_latin_hypercube(size = 10) # 拉丁方网格
# grid_max_entropy(size = 10) # 最大熵网格
hpgrid_dt
log10(hpgrid_dt$cost_complexity)
# 网格也可以自己手动生成expand.grid()
# hpgrid_dt <- expand.grid(
#   tree_depth = c(2:5),
#   min_n = c(5, 11),
#   cost_complexity = 10^(-5:-1)
# )

# 交叉验证网格搜索过程
set.seed(42)
tune_dt <- wk_dt %>%
  tune_grid(
    resamples = folds,
    grid = hpgrid_dt,
    metrics = metricset_cls2,
    control = control_grid(save_pred = T, 
                           verbose = T,
                           event_level = "second",
                           parallel_over = "everything",
                           save_workflow = T)
  )

#########################  超参数寻优2选1-贝叶斯优化

# 贝叶斯优化超参数
set.seed(42)
tune_dt <- wk_dt %>%
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
eval_tune_dt <- tune_dt %>%
  collect_metrics()
eval_tune_dt

# 图示
# autoplot(tune_dt)
eval_tune_dt %>% 
  filter(.metric == "roc_auc") %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~mean, colorscale = 'Jet', showscale = T),
    dimensions = list(
      list(label = 'cost_complexity', values = ~cost_complexity),
      list(label = 'tree_depth', values = ~tree_depth),
      list(label = 'min_n', values = ~min_n)
    )
  ) %>%
  plotly::layout(title = "DT HPO Guided by AUCROC",
                 font = list(family = "serif"))

# 经过交叉验证得到的最优超参数
hpbest_dt <- tune_dt %>%
  select_by_one_std_err(metric = "roc_auc", desc(cost_complexity))
hpbest_dt

# 采用最优超参数组合训练最终模型
set.seed(42)
final_dt <- wk_dt %>%
  finalize_workflow(hpbest_dt) %>%
  fit(traindata)
final_dt

##################################################################

# 训练集预测评估
predtrain_dt <- eval4cls2(
  model = final_dt, 
  dataset = traindata, 
  yname = "AHD", 
  modelname = "DT", 
  datasetname = "traindata",
  cutoff = "yueden",
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtrain_dt$prediction
predtrain_dt$predprobplot
predtrain_dt$rocplot
predtrain_dt$prplot
predtrain_dt$caliplot
predtrain_dt$cmplot
predtrain_dt$metrics
predtrain_dt$diycutoff
predtrain_dt$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtrain_dt$proc)
pROC::ci.auc(predtrain_dt$proc)

# 预测评估测试集预测评估
predtest_dt <- eval4cls2(
  model = final_dt, 
  dataset = testdata, 
  yname = "AHD", 
  modelname = "DT", 
  datasetname = "testdata",
  cutoff = predtrain_dt$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtest_dt$prediction
predtest_dt$predprobplot
predtest_dt$rocplot
predtest_dt$prplot
predtest_dt$caliplot
predtest_dt$cmplot
predtest_dt$metrics
predtest_dt$diycutoff
predtest_dt$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtest_dt$proc)
pROC::ci.auc(predtest_dt$proc)

# ROC比较检验
pROC::roc.test(predtrain_dt$proc, predtest_dt$proc)


# 合并训练集和测试集上ROC曲线
predtrain_dt$rocresult %>%
  bind_rows(predtest_dt$rocresult) %>%
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
predtrain_dt$prresult %>%
  bind_rows(predtest_dt$prresult) %>%
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
predtrain_dt$caliresult %>%
  bind_rows(predtest_dt$caliresult) %>%
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
predtrain_dt$metrics %>%
  bind_rows(predtest_dt$metrics) %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)


# 最优超参数交叉验证的结果
evalcv_dt <- bestcv4cls2(
  wkflow = wk_dt,
  tuneresult = tune_dt,
  hpbest = hpbest_dt,
  yname = "AHD",
  modelname = "DT",
  v = 5,
  positivelevel = yourpositivelevel
)
evalcv_dt$cvroc
evalcv_dt$cvpr
evalcv_dt$evalcv


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
  bind_cols(predict(final_dt, new_data = newHeart, type = "prob"))%>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Yes >= predtrain_dt$diycutoff, 
             yourpositivelevel, 
             yournegativelevel)
    )
  )
# readr::write_excel_csv(predresult, "决策树二分类预测结果.csv")
# 评估指标
prednew_dt <- eval4cls2(
  model = final_dt, 
  dataset = newHeart, 
  yname = "AHD", 
  modelname = "DT", 
  datasetname = "newHeart",
  cutoff = predtrain_dt$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
prednew_dt$prediction
prednew_dt$predprobplot
prednew_dt$rocplot
prednew_dt$prplot
prednew_dt$caliplot
prednew_dt$cmplot
prednew_dt$metrics
prednew_dt$diycutoff
prednew_dt$ksplot

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
final_dt2 <- final_dt %>%
  extract_fit_engine()
final_dt2

# 树形图
rpart.plot::rpart.plot(final_dt2, family = "serif")

# 变量重要性
final_dt2$variable.importance
# par(mar = c(10, 3, 1, 1))
barplot(final_dt2$variable.importance, las = 2, family = "serif")

######################## DALEX解释对象

explainer_dt <- DALEXtra::explain_tidymodels(
  final_dt, 
  data = traindatax,
  y = ifelse(traindata$AHD == yourpositivelevel, 1, 0),
  type = "classification",
  label = "DT"
)
# 变量重要性
vip_dt <- viplot(explainer_dt)
vipdata_dt <- vip_dt$data
vip_dt$plot

# 变量偏依赖图
pdplot(explainer_dt, convars)
pdplot(explainer_dt, "Age")
pdplot(explainer_dt, catvars)
pdplot(explainer_dt, "Thal")

###################################### iml解释对象

predictor_dt <- iml::Predictor$new(
  final_dt, 
  data = traindatax,
  y = traindata$AHD,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>%
      rename_with(~gsub(".pred_", "", .x))
  },
  type = "prob"
)
# 交互作用
interact_dt <- iml::Interaction$new(predictor_dt)
plot(interact_dt) +
  theme_minimal()

interact_dt_1vo <- 
  iml::Interaction$new(predictor_dt, feature = "Age")
plot(interact_dt_1vo) +
  theme_minimal()

interact_dt_1v1 <- iml::FeatureEffect$new(
  predictor_dt, 
  feature = c("Age", "Oldpeak"),
  method = "pdp"
)
plot(interact_dt_1v1) +
  scale_fill_viridis_c() +
  labs(fill = "") +
  theme_minimal()

###################################### lime单样本预测分解

explainer_dt <- lime::lime(
  traindatax,
  lime::as_classifier(final_dt, c(yournegativelevel, yourpositivelevel))
)
explanation_dt <- lime::explain(
  traindatax[1,],  # 训练集第1个样本
  explainer_dt, 
  n_labels = 2, 
  n_features = ncol(traindatax)
)
lime::plot_features(explanation_dt)

######################## fastshap包

shapresult_dt <- shap4cls2(
  finalmodel = final_dt,
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
shapresult_dt$shapvipplot

# 单样本预测分解
shap41 <- shapviz::shapviz(
  shapresult_dt$shapley,
  X = traindatax
)
shapviz::sv_force(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
shapviz::sv_waterfall(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))

# 所有分类变量的shap图示
shapresult_dt$shapplotd_facet
shapresult_dt$shapplotd_one
# 所有连续变量的shap图示
shapresult_dt$shapplotc_facet
shapresult_dt$shapplotc_one
shapresult_dt$shapplotc_one2
# 单变量shap图示
sdplot(shapresult_dt, "Thal", "AHD")
sdplot(shapresult_dt, "Age", "AHD")

# 所有变量一张图
# shap变量重要性
shapresult_dt$shapvipplot_unity
# shap依赖图
shapresult_dt$shapplot_unity


#################################################################

# 学习曲线
lcN <- 
  floor(seq(nrow(traindata)%/%2, nrow(traindata), length = 10))
lcresult_dt <- data.frame()
for (i in lcN) {
  
  set.seed(i)
  traindatai <- traindata[sample(nrow(traindata), i), ]
  
  i_dt <-  wk_dt %>%
    finalize_workflow(hpbest_dt) %>%
    fit(traindatai)
  
  predtrain_i_dt <- traindatai %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_dt, new_data = traindatai, type = "prob")) %>%
    yardstick::roc_auc(AHD,
                       paste0(".pred_", yourpositivelevel), 
                       event_level = "second") %>%
    mutate(datasetname = paste("traindata", i, sep = "-"))
  
  predtest_i_dt <- testdata %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_dt, new_data = testdata, type = "prob")) %>%
    yardstick::roc_auc(AHD, 
                       paste0(".pred_", yourpositivelevel),
                       event_level = "second") %>%
    mutate(datasetname = paste("testdata", i, sep = "-"))
  
  predi <- bind_rows(predtrain_i_dt, predtest_i_dt)
  lcresult_dt <- rbind(lcresult_dt, predi)
  print(i)
}
# 图示
lcresult_dt %>%
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
save(datarecipe_dt,
     model_dt,
     wk_dt,
     #hpgrid_dt, # 如果采用贝叶斯优化则删掉这一行
     tune_dt,
     predtrain_dt,
     predtest_dt,
     evalcv_dt,
     vipdata_dt,
     file = "F:/Mach_learn_data/Decision Tree/evalresult_dt.RData")

# 保存模型结果供shiny部署之用，本课程不包括shiny内容
final_dt_heart <- final_dt
traindata_heart <- traindata
save(final_dt_heart,
     traindata_heart,
     file = "F:/Mach_learn_data/Decision Tree/shiny_dt_heart.RData")
