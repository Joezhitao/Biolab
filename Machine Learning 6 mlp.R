# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---单隐藏层神经网络
# https://parsnip.tidymodels.org/reference/details_mlp_nnet.html

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
datarecipe_mlp <- recipe(formula = AHD ~ ., traindata) %>%
  step_dummy(all_nominal_predictors(), naming = new_dummy_names) %>%
  step_nzv(all_predictors()) %>% 
  step_range(all_predictors())
datarecipe_mlp


# 设定模型
model_mlp <- mlp(
  mode = "classification",
  engine = "nnet",
  hidden_units = tune(),
  penalty = tune(),
  epochs = tune()
) %>%
  set_args(MaxNWts = 10000)
model_mlp

# workflow
wk_mlp <- 
  workflow() %>%
  add_recipe(datarecipe_mlp) %>%
  add_model(model_mlp)
wk_mlp

##############################################################
############################  超参数寻优2选1-网格搜索

# 超参数寻优网格
set.seed(42)
hpgrid_mlp <- parameters(
  hidden_units(range = c(15, 24)),
  penalty(range = c(-3, 0)),
  epochs(range = c(50, 150))
) %>%
  grid_regular(levels = 3) # 常规网格
# grid_random(size = 5) # 随机网格
# grid_latin_hypercube(size = 10) # 拉丁方网格
# grid_max_entropy(size = 10) # 最大熵网格
hpgrid_mlp
# 网格也可以自己手动生成expand.grid()

# 交叉验证网格搜索过程
set.seed(42)
tune_mlp <- wk_mlp %>%
  tune_grid(
    resamples = folds,
    grid = hpgrid_mlp,
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
tune_mlp <- wk_mlp %>%
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
eval_tune_mlp <- tune_mlp %>%
  collect_metrics()
eval_tune_mlp

# 图示
# autoplot(tune_mlp)
eval_tune_mlp %>% 
  filter(.metric == "roc_auc") %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~mean, colorscale = 'Jet', showscale = T),
    dimensions = list(
      list(label = 'hidden_units', values = ~hidden_units),
      list(label = 'penalty', values = ~penalty),
      list(label = 'epochs', values = ~epochs)
    )
  ) %>%
  plotly::layout(title = "MLP HPO Guided by AUCROC",
                 font = list(family = "serif"))

# 经过交叉验证得到的最优超参数
hpbest_mlp <- tune_mlp %>%
  select_by_one_std_err(metric = "roc_auc", desc(penalty))
hpbest_mlp

# 采用最优超参数组合训练最终模型
set.seed(42)
final_mlp <- wk_mlp %>%
  finalize_workflow(hpbest_mlp) %>%
  fit(traindata)
final_mlp

##################################################################

# 训练集预测评估
predtrain_mlp <- eval4cls2(
  model = final_mlp, 
  dataset = traindata, 
  yname = "AHD", 
  modelname = "MLP", 
  datasetname = "traindata",
  cutoff = "yueden",
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtrain_mlp$prediction
predtrain_mlp$predprobplot
predtrain_mlp$rocplot
predtrain_mlp$prplot
predtrain_mlp$caliplot
predtrain_mlp$cmplot
predtrain_mlp$metrics
predtrain_mlp$diycutoff
predtrain_mlp$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtrain_mlp$proc)
pROC::ci.auc(predtrain_mlp$proc)

# 预测评估测试集预测评估
predtest_mlp <- eval4cls2(
  model = final_mlp, 
  dataset = testdata, 
  yname = "AHD", 
  modelname = "MLP", 
  datasetname = "testdata",
  cutoff = predtrain_mlp$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
predtest_mlp$prediction
predtest_mlp$predprobplot
predtest_mlp$rocplot
predtest_mlp$prplot
predtest_mlp$caliplot
predtest_mlp$cmplot
predtest_mlp$metrics
predtest_mlp$diycutoff
predtest_mlp$ksplot

# pROC包auc值及其置信区间
pROC::auc(predtest_mlp$proc)
pROC::ci.auc(predtest_mlp$proc)

# ROC比较检验
pROC::roc.test(predtrain_mlp$proc, predtest_mlp$proc)


# 合并训练集和测试集上ROC曲线
predtrain_mlp$rocresult %>%
  bind_rows(predtest_mlp$rocresult) %>%
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
predtrain_mlp$prresult %>%
  bind_rows(predtest_mlp$prresult) %>%
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
predtrain_mlp$caliresult %>%
  bind_rows(predtest_mlp$caliresult) %>%
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
predtrain_mlp$metrics %>%
  bind_rows(predtest_mlp$metrics) %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)


# 最优超参数交叉验证的结果
evalcv_mlp <- bestcv4cls2(
  wkflow = wk_mlp,
  tuneresult = tune_mlp,
  hpbest = hpbest_mlp,
  yname = "AHD",
  modelname = "MLP",
  v = 5,
  positivelevel = yourpositivelevel
)
evalcv_mlp$cvroc
evalcv_mlp$cvpr
evalcv_mlp$evalcv


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
  bind_cols(predict(final_mlp, new_data = newHeart, type = "prob"))%>%
  mutate(
    .pred_class = factor(
      ifelse(.pred_Yes >= predtrain_mlp$diycutoff, 
             yourpositivelevel, 
             yournegativelevel)
    )
  )
# readr::write_excel_csv(predresult, "神经网络二分类预测结果.csv")
# 评估指标
prednew_mlp <- eval4cls2(
  model = final_mlp, 
  dataset = newHeart, 
  yname = "AHD", 
  modelname = "MLP", 
  datasetname = "newHeart",
  cutoff = predtrain_mlp$diycutoff,
  positivelevel = yourpositivelevel,
  negativelevel = yournegativelevel
)
prednew_mlp$prediction
prednew_mlp$predprobplot
prednew_mlp$rocplot
prednew_mlp$prplot
prednew_mlp$caliplot
prednew_mlp$cmplot
prednew_mlp$metrics
prednew_mlp$diycutoff
prednew_mlp$ksplot

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
final_mlp2 <- final_mlp %>%
  extract_fit_engine()
final_mlp2

# 神经网络拓扑结构图
NeuralNetTools::plotnet(final_mlp2, pad_x = 0.5)

# 变量重要性
NeuralNetTools::garson(final_mlp2) +
  coord_flip() +
  theme(text = element_text(family = "serif"))
NeuralNetTools::olden(final_mlp2) +
  coord_flip() +
  theme(text = element_text(family = "serif"))

######################## DALEX解释对象

explainer_mlp <- DALEXtra::explain_tidymodels(
  final_mlp, 
  data = traindatax,
  y = ifelse(traindata$AHD == yourpositivelevel, 1, 0),
  type = "classification",
  label = "MLP"
)
# 变量重要性
vip_mlp <- viplot(explainer_mlp)
vipdata_mlp <- vip_mlp$data
vip_mlp$plot

# 变量偏依赖图
pdplot(explainer_mlp, convars)
pdplot(explainer_mlp, "Age")
pdplot(explainer_mlp, catvars)
pdplot(explainer_mlp, "Thal")

###################################### iml解释对象

predictor_mlp <- iml::Predictor$new(
  final_mlp, 
  data = traindatax,
  y = traindata$AHD,
  predict.function = function(model, newdata){
    predict(model, newdata, type = "prob") %>%
      rename_with(~gsub(".pred_", "", .x))
  },
  type = "prob"
)
# 交互作用
interact_mlp <- iml::Interaction$new(predictor_mlp)
plot(interact_mlp) +
  theme_minimal()

interact_mlp_1vo <- 
  iml::Interaction$new(predictor_mlp, feature = "Age")
plot(interact_mlp_1vo) +
  theme_minimal()

interact_mlp_1v1 <- iml::FeatureEffect$new(
  predictor_mlp, 
  feature = c("Age", "Oldpeak"),
  method = "pdp"
)
plot(interact_mlp_1v1) +
  scale_fill_viridis_c() +
  labs(fill = "") +
  theme_minimal()

###################################### lime单样本预测分解

explainer_mlp <- lime::lime(
  traindatax,
  lime::as_classifier(final_mlp, c(yournegativelevel, yourpositivelevel))
)
explanation_mlp <- lime::explain(
  traindatax[1,],  # 训练集第1个样本
  explainer_mlp, 
  n_labels = 2, 
  n_features = ncol(traindatax)
)
lime::plot_features(explanation_mlp)

######################## fastshap包

shapresult_mlp <- shap4cls2(
  finalmodel = final_mlp,
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
shapresult_mlp$shapvipplot

# 单样本预测分解
shap41 <- shapviz::shapviz(
  shapresult_mlp$shapley,
  X = traindatax
)
shapviz::sv_force(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))
shapviz::sv_waterfall(shap41, row_id = 1)  +  # 训练集第1个样本
  theme(text = element_text(family = "serif"))

# 所有分类变量的shap图示
shapresult_mlp$shapplotd_facet
shapresult_mlp$shapplotd_one
# 所有连续变量的shap图示
shapresult_mlp$shapplotc_facet
shapresult_mlp$shapplotc_one
shapresult_mlp$shapplotc_one2
# 单变量shap图示
sdplot(shapresult_mlp, "Thal", "AHD")
sdplot(shapresult_mlp, "Age", "AHD")

# 所有变量一张图
# shap变量重要性
shapresult_mlp$shapvipplot_unity
# shap依赖图
shapresult_mlp$shapplot_unity


#################################################################

# 学习曲线
lcN <- 
  floor(seq(nrow(traindata)%/%2, nrow(traindata), length = 10))
lcresult_mlp <- data.frame()
for (i in lcN) {
  
  set.seed(i)
  traindatai <- traindata[sample(nrow(traindata), i), ]
  
  i_mlp <-  wk_mlp %>%
    finalize_workflow(hpbest_mlp) %>%
    fit(traindatai)
  
  predtrain_i_mlp <- traindatai %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_mlp, new_data = traindatai, type = "prob")) %>%
    yardstick::roc_auc(AHD,
                       paste0(".pred_", yourpositivelevel), 
                       event_level = "second") %>%
    mutate(datasetname = paste("traindata", i, sep = "-"))
  
  predtest_i_mlp <- testdata %>%
    dplyr::select(AHD) %>%
    bind_cols(predict(i_mlp, new_data = testdata, type = "prob")) %>%
    yardstick::roc_auc(AHD, 
                       paste0(".pred_", yourpositivelevel),
                       event_level = "second") %>%
    mutate(datasetname = paste("testdata", i, sep = "-"))
  
  predi <- bind_rows(predtrain_i_mlp, predtest_i_mlp)
  lcresult_mlp <- rbind(lcresult_mlp, predi)
  print(i)
}
# 图示
lcresult_mlp %>%
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
save(datarecipe_mlp,
     model_mlp,
     wk_mlp,
     tune_mlp,
     predtrain_mlp,
     predtest_mlp,
     evalcv_mlp,
     vipdata_mlp,
     file = "F:/Mach_learn_data/Decision Tree/evalresult_mlp.RData")

# 保存模型结果供shiny部署之用，本课程不包括shiny内容
final_mlp_heart <- final_mlp
traindata_heart <- traindata
save(final_mlp_heart,
     traindata_heart,
     file = "F:/Mach_learn_data/Decision Tree/shiny_mlp_heart.RData")
