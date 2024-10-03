# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---嵌套式重采样，用于模型评估和比较

#############################################################

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

# 变量类型修正后数据概况
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

#############################################################

# 数据预处理配方
datarecipe_dt <- recipe(formula = AHD ~ ., Heart)
datarecipe_dt

# 设定模型
model_dt <- decision_tree(
  mode = "classification",
  engine = "rpart",
  tree_depth = tune(),
  min_n = tune(),
  cost_complexity = tune()
) %>%
  set_args(model=TRUE)
model_dt

# workflow
wk_dt <- 
  workflow() %>%
  add_recipe(datarecipe_dt) %>%
  add_model(model_dt)
wk_dt

# 超参数寻优网格
set.seed(42)
hpgrid_dt <- parameters(
  tree_depth(range = c(3, 7)),
  min_n(range = c(5, 10)),
  cost_complexity(range = c(-6, -1))
) %>%
  grid_random(size = 10) # 随机网格
hpgrid_dt


#############################################################

# 嵌套式重采样设定
oifolds <- nested_cv(
  Heart,
  # 外层3折交叉验证，不建议使用bootstrap
  outside = vfold_cv(v = 3, strata = AHD),
  # 内层5折交叉验证，可以使用bootstrap
  inside = vfold_cv(v = 5, strata = AHD)
)
oifolds
oifolds$splits[[1]]
oifolds$inner_resamples[[1]]
oifolds$splits[[2]]
oifolds$inner_resamples[[2]]
oifolds$splits[[3]]
oifolds$inner_resamples[[3]]

#############################################################

nestedmetric <- list()
nestedpred <- list()
for (i in seq_along(oifolds$splits)) {
  
  # 外层训练集上采用内层重采样方式做超参数调优
  ############################  超参数寻优2选1-网格搜索
  # i = 1
  # 交叉验证网格搜索过程
  set.seed(42)
  tune_outi <- wk_dt %>%
    tune_grid(
      resamples = oifolds$inner_resamples[[i]],
      grid = hpgrid_dt,
      metrics = metricset_cls2,
      control = control_grid(save_pred = T, 
                             verbose = T,
                             event_level = "second",
                             parallel_over = "everything",
                             save_workflow = T)
    )
  
  # #########################  超参数寻优2选1-贝叶斯优化
  # 
  # # 贝叶斯优化超参数
  # set.seed(42)
  # tune_outi <- wk_dt %>%
  #   tune_bayes(
  #     resamples = oifolds$inner_resamples[[i]],
  #     initial = 10,
  #     iter = 50,
  #     metrics = metricset_cls2,
  #     control = control_bayes(save_pred = T,
  #                            verbose = T,
  #                            no_improve = 10,
  #                            event_level = "second",
  #                            parallel_over = "everything",
  #                            save_workflow = T)
  #   )
  
  ########################  超参数寻优结束
  
  # 经过交叉验证得到的最优超参数
  hpbest_outi <- tune_outi %>%
    select_best(metric = "roc_auc")
  # hpbest_outi
  
  # 采用最优超参数在外层训练集上训练模型并预测外层测试集
  set.seed(42)
  model_outi <- wk_dt %>%
    finalize_workflow(hpbest_outi) %>%
    last_fit(
      oifolds$splits[[i]], 
      metrics = metricset_cls2,
      control = control_last_fit(event_level = "second")
    )
  
  # 外层测试集上的性能指标
  nestedmetric[[i]] <- model_outi$.metrics %>%
    bind_rows() %>%
    mutate(.out = i,
           .model = "DT") %>%
    select(.model, .out, .metric, .estimate)
  # 外层测试集上的预测结果
  nestedpred[[i]] <- model_outi$.predictions %>%
    bind_rows() %>%
    mutate(.out = i,
           .model = "DT") %>%
    select(.model, .out, .row, .pred_No, .pred_Yes, .obs = AHD)
  # model_outi$.workflow
}

# 外层测试集性能指标
nestedmetric_all <- bind_rows(nestedmetric)
nestedmetric_all %>%
  group_by(.model, .metric) %>%
  summarise(mean = mean(.estimate),
            sd = sd(.estimate)/sqrt(length(oifolds$splits))) %>%
  ungroup()

# 外层测试集预测结果汇总
nestedpred_all <- bind_rows(nestedpred)
nestedpred_all %>%
  group_by(.out) %>%
  roc_curve(.obs, 
            ends_with(yourpositivelevel), 
            event_level = "second") %>%
  ungroup() %>%
  left_join(filter(nestedmetric_all, .metric == "roc_auc")) %>%
  mutate(idAUC = paste("Fold", .out, " ROCAUC:", round(.estimate, 4)),
         idAUC = forcats::as_factor(idAUC)) %>%
  ggplot(aes(x = 1-specificity, y = sensitivity, 
             color = idAUC)) +
  geom_path(linewidth = 1) +
  geom_abline(linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  labs(color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        text = element_text(family = "serif"))

#############################################################

# 最后用于外部预测验证和解释的模型基于所用数据做超参数调优和拟合

