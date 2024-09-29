# Joe---R语言tidymodels包机器学习分类与回归模型---二分类---模型比较

#############################################################
# remotes::install_github("tidymodels/probably")
library(tidymodels)

# 加载各个模型的评估结果
evalfiles <- list.files("F:/Mach_learn_data/Decision Tree/", full.names = T)
lapply(evalfiles, load, .GlobalEnv)

# 横向比较的模型个数
nmodels <- 12
cols4model <- rainbow(nmodels)  # 模型统一配色

#############################################################

# 各个模型在训练集上的性能指标
evaltrain <- bind_rows(
  lapply(list(predtrain_logistic, predtrain_dt, 
              predtrain_lasso, predtrain_ridge, predtrain_enet,
              predtrain_knn, predtrain_lightgbm, predtrain_rf,
              predtrain_xgboost, predtrain_svm, predtrain_mlp,
              predtrain_stack), 
         "[[", 
         "metrics")
) %>%
  mutate(model = forcats::as_factor(model))
evaltrain
# 平行线图
evaltrain_max <-   evaltrain %>% 
  filter(!(.metric %in% c("detection_prevalence"))) %>%
  group_by(.metric) %>%
  slice_max(.estimate)
evaltrain_min <-   evaltrain %>% 
  filter(!(.metric %in% c("detection_prevalence"))) %>%
  group_by(.metric) %>%
  slice_min(.estimate)

evaltrain %>%
  filter(!(.metric %in% c("detection_prevalence"))) %>%
  ggplot(aes(x = .metric, y = .estimate, color = model)) +
  geom_point() +
  geom_line(aes(group = model)) +
  # ggrepel::geom_text_repel(evaltrain_max,
  #                          mapping = aes(label = model),
  #                          nudge_y = 0.05,
  #                          angle = 90,
  #                          show.legend = F) +
  # ggrepel::geom_text_repel(evaltrain_min,
  #                          mapping = aes(label = model),
  #                          nudge_y = -0.05,
  #                          angle = 90,
  #                          show.legend = F) +
  scale_color_manual(values = cols4model) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "", y = "", color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(family = "serif"))

# 指标热图
evaltrain %>%
  filter(!(.metric %in% c("detection_prevalence"))) %>%
  dplyr::select(model, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate) %>%
  mutate(model = reorder(model, roc_auc)) %>%
  pivot_longer(cols = -1) %>%
  group_by(name) %>%
  mutate(valuescale = (value-min(value)) / (max(value)-min(value))) %>%
  ungroup() %>%
  ggplot(aes(x = name, y = model, fill = valuescale)) +
  geom_tile(color = "white", show.legend = F) +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient(low = "green", high = "red") +
  labs(x = "", y = "", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        text = element_text(family = "serif"))

# 各个模型在训练集上的性能指标表格
evaltrain2 <- evaltrain %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)
evaltrain2

# 各个模型在训练集上的性能指标图示
# ROCAUC
evaltrain2 %>%
  ggplot(aes(x = model, y = roc_auc, fill = model)) +
  geom_col(width = 0.5, show.legend = F) +
  geom_text(aes(label = round(roc_auc, 2)), 
            nudge_y = -0.03) +
  scale_fill_manual(values = cols4model) +
  labs(x = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "serif"))

#############################

# 各个模型在训练集上的ROC
roctrain <- bind_rows(
  lapply(list(predtrain_logistic, predtrain_dt, 
              predtrain_lasso, predtrain_ridge, predtrain_enet,
              predtrain_knn, predtrain_lightgbm, predtrain_rf,
              predtrain_xgboost, predtrain_svm, predtrain_mlp,
              predtrain_stack), 
         "[[", 
         "rocresult")
) %>%
  mutate(model = forcats::as_factor(model))
roctrain

roctrain %>%
  mutate(modelauc = paste(model,  curvelab),
         modelauc = forcats::as_factor(modelauc)) %>%
  ggplot(aes(x = 1-specificity, y = sensitivity, color = modelauc)) +
  geom_path(linewidth = 1) +
  geom_abline(linetype = "dashed") +
  scale_color_manual(values = cols4model) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  labs(color = "", title = paste0("ROCs on traindata")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "serif"))

# 各个模型在训练集上的PR
prtrain <- bind_rows(
  lapply(list(predtrain_logistic, predtrain_dt, 
              predtrain_lasso, predtrain_ridge, predtrain_enet,
              predtrain_knn, predtrain_lightgbm, predtrain_rf,
              predtrain_xgboost, predtrain_svm, predtrain_mlp,
              predtrain_stack), 
         "[[", 
         "prresult")
) %>%
  mutate(model = forcats::as_factor(model))
prtrain

prtrain %>%
  mutate(modelauc = paste(model, curvelab),
         modelauc = forcats::as_factor(modelauc)) %>%
  ggplot(aes(x = recall, y = precision, color = modelauc)) +
  geom_path(linewidth = 1) +
  geom_abline(linetype = "dashed", slope = -1, intercept = 1) +
  scale_color_manual(values = cols4model) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  labs(color = "", title = paste0("PRs on traindata")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "serif"))

############################

# 各个模型在训练集上的预测概率
predtrain <- bind_rows(
  lapply(list(predtrain_logistic, predtrain_dt, 
              predtrain_lasso, predtrain_ridge, predtrain_enet,
              predtrain_knn, predtrain_lightgbm, predtrain_rf,
              predtrain_xgboost, predtrain_svm, predtrain_mlp,
              predtrain_stack), 
         "[[", 
         "prediction")
) %>%
  mutate(model = forcats::as_factor(model))
predtrain

# 各个模型在训练集上的预测概率---宽数据
predtrain2 <- predtrain %>%
  dplyr::select(-.pred_No) %>%
  mutate(id = rep(1:nrow(predtrain_logistic$prediction), 
                  length(unique(predtrain$model)))) %>%
  pivot_wider(id_cols = c(id, .obs), 
              names_from = model, 
              values_from = .pred_Yes) %>%
  dplyr::select(id, .obs, sort(unique(predtrain$model)))
predtrain2

############################


# 各个模型在训练集上的校准曲线
# 校准曲线附加置信区间
predtrain %>%
  probably::cal_plot_breaks(.obs, 
                            .pred_Yes, 
                            event_level = "second", 
                            num_breaks = 5,  # 可以改大改小
                            .by = model) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_color_manual(values = cols4model) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        text = element_text(family = "serif"))

# 各个模型在训练集上的校准曲线
calitrain <- bind_rows(
  lapply(list(predtrain_logistic, predtrain_dt, 
              predtrain_lasso, predtrain_ridge, predtrain_enet,
              predtrain_knn, predtrain_lightgbm, predtrain_rf,
              predtrain_xgboost, predtrain_svm, predtrain_mlp,
              predtrain_stack), 
         "[[", 
         "caliresult")
) %>%
  mutate(model = forcats::as_factor(model))
calitrain

calitrain %>%
  mutate(model = forcats::as_factor(model)) %>%
  ggplot(aes(x = predprobgroup, y = Fraction, color = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, pch = 15) +
  geom_abline(linetype = "dashed") +
  scale_color_manual(values = cols4model) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  labs(color = "", x = "Bin Midpoint", y = "Event Rate",
       title = paste0("calibration on traindata")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "serif"))


############################

# 各个模型在训练集上的DCA
traindca_obj <- dcurves::dca(as.formula(
  paste0(".obs ~ ", 
         paste(colnames(predtrain2)[3:ncol(predtrain2)], 
               collapse = " + "))
),
data = predtrain2,
thresholds = seq(0, 1, by = 0.01)
)
plot(traindca_obj, smooth = T, span = 0.5) +
  scale_color_manual(values = c("black", "grey", cols4model)) +
  labs(title = "DCA on traindata") +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        text = element_text(family = "serif"))

#############################################################

# 各个模型在测试集上的性能指标
evaltest <- bind_rows(
  lapply(list(predtest_logistic, predtest_dt, 
              predtest_lasso, predtest_ridge, predtest_enet,
              predtest_knn, predtest_lightgbm, predtest_rf,
              predtest_xgboost, predtest_svm, predtest_mlp,
              predtest_stack), 
         "[[", 
         "metrics")
) %>%
  mutate(model = forcats::as_factor(model))
evaltest
# 平行线图
evaltest_max <-   evaltest %>% 
  filter(!(.metric %in% c("detection_prevalence"))) %>%
  group_by(.metric) %>%
  slice_max(.estimate)
evaltest_min <-   evaltest %>% 
  filter(!(.metric %in% c("detection_prevalence"))) %>%
  group_by(.metric) %>%
  slice_min(.estimate)

evaltest %>%
  filter(!(.metric %in% c("detection_prevalence"))) %>%
  ggplot(aes(x = .metric, y = .estimate, color = model)) +
  geom_point() +
  geom_line(aes(group = model)) +
  # ggrepel::geom_text_repel(evaltest_max, 
  #                          mapping = aes(label = model), 
  #                          nudge_y = 0.05,
  #                          angle = 90,
  #                          show.legend = F) +
  # ggrepel::geom_text_repel(evaltest_min, 
  #                          mapping = aes(label = model), 
  #                          nudge_y = -0.05,
  #                          angle = 90,
  #                          show.legend = F) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = cols4model) +
  labs(x = "", y = "", color = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        text = element_text(family = "serif"))

# 指标热图
evaltest %>%
  filter(!(.metric %in% c("detection_prevalence"))) %>%
  dplyr::select(model, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate) %>%
  mutate(model = reorder(model, roc_auc)) %>%
  pivot_longer(cols = -1) %>%
  group_by(name) %>%
  mutate(valuescale = (value-min(value)) / (max(value)-min(value))) %>%
  ungroup() %>%
  ggplot(aes(x = name, y = model, fill = valuescale)) +
  geom_tile(color = "white", show.legend = F) +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient(low = "green", high = "red") +
  labs(x = "", y = "", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        text = element_text(family = "serif"))

# 各个模型在测试集上的性能指标表格
evaltest2 <- evaltest %>%
  dplyr::select(-.estimator) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)
evaltest2

# 各个模型在测试集上的性能指标图示
# ROCAUC
evaltest2 %>%
  ggplot(aes(x = model, y = roc_auc, fill = model)) +
  geom_col(width = 0.5, show.legend = F) +
  geom_text(aes(label = round(roc_auc, 2)), 
            nudge_y = -0.03) +
  scale_fill_manual(values = cols4model) +
  labs(x = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "serif"))

#############################

# 各个模型在测试集上的ROC
roctest <- bind_rows(
  lapply(list(predtest_logistic, predtest_dt,
              predtest_lasso, predtest_ridge, predtest_enet,
              predtest_knn, predtest_lightgbm, predtest_rf,
              predtest_xgboost, predtest_svm, predtest_mlp,
              predtest_stack), 
         "[[", 
         "rocresult")
) %>%
  mutate(model = forcats::as_factor(model))
roctest

roctest %>%
  mutate(modelauc = paste(model,  curvelab),
         modelauc = forcats::as_factor(modelauc)) %>%
  ggplot(aes(x = 1-specificity, y = sensitivity, color = modelauc)) +
  geom_path(linewidth = 1) +
  geom_abline(linetype = "dashed") +
  scale_color_manual(values = cols4model) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  labs(color = "", title = paste0("ROCs on testdata")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "serif"))

# 各个模型在测试集上的PR
prtest <- bind_rows(
  lapply(list(predtest_logistic, predtest_dt,
              predtest_lasso, predtest_ridge, predtest_enet,
              predtest_knn, predtest_lightgbm, predtest_rf,
              predtest_xgboost, predtest_svm, predtest_mlp,
              predtest_stack), 
         "[[", 
         "prresult")
) %>%
  mutate(model = forcats::as_factor(model))
prtest

prtest %>%
  mutate(modelauc = paste(model, curvelab),
         modelauc = forcats::as_factor(modelauc)) %>%
  ggplot(aes(x = recall, y = precision, color = modelauc)) +
  geom_path(linewidth = 1) +
  geom_abline(linetype = "dashed", slope = -1, intercept = 1) +
  scale_color_manual(values = cols4model) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  labs(color = "", title = paste0("PRs on testdata")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "serif"))

############################

# 各个模型在测试集上的预测概率
predtest <- bind_rows(
  lapply(list(predtest_logistic, predtest_dt,
              predtest_lasso, predtest_ridge, predtest_enet,
              predtest_knn, predtest_lightgbm, predtest_rf,
              predtest_xgboost, predtest_svm, predtest_mlp,
              predtest_stack), 
         "[[", 
         "prediction")
) %>%
  mutate(model = forcats::as_factor(model))
predtest

# 各个模型在测试集上的预测概率---宽数据
predtest2 <- predtest %>%
  dplyr::select(-.pred_No) %>%
  mutate(id = rep(1:nrow(predtest_logistic$prediction), 
                  length(unique(predtest$model)))) %>%
  pivot_wider(id_cols = c(id, .obs), 
              names_from = model, 
              values_from = .pred_Yes) %>%
  dplyr::select(id, .obs, sort(unique(predtest$model)))
predtest2

############################


# 各个模型在测试集上的校准曲线
# 校准曲线附加置信区间
predtest %>%
  probably::cal_plot_breaks(.obs, 
                            .pred_Yes, 
                            event_level = "second", 
                            num_breaks = 5,  # 可以改大改小
                            .by = model) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_color_manual(values = cols4model) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        text = element_text(family = "serif"))

# 各个模型在测试集上的校准曲线
calitest <- bind_rows(
  lapply(list(predtest_logistic, predtest_dt, 
              predtest_lasso, predtest_ridge, predtest_enet,
              predtest_knn, predtest_lightgbm, predtest_rf,
              predtest_xgboost, predtest_svm, predtest_mlp,
              predtest_stack), 
         "[[", 
         "caliresult")
) %>%
  mutate(model = forcats::as_factor(model))
calitest

calitest %>%
  mutate(model = forcats::as_factor(model)) %>%
  ggplot(aes(x = predprobgroup, y = Fraction, color = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, pch = 15) +
  geom_abline(linetype = "dashed") +
  scale_color_manual(values = cols4model) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), 
                     breaks = seq(0, 1, by = 0.2),
                     labels = c(0, seq(0.2, 0.8, by = 0.2), 1)) +
  labs(color = "", x = "Bin Midpoint", y = "Event Rate",
       title = paste0("calibration on testdata")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(family = "serif"))


############################

# 各个模型在测试集上的DCA
testdca_obj <- dcurves::dca(as.formula(
  paste0(".obs ~ ", 
         paste(colnames(predtest2)[3:ncol(predtest2)], 
               collapse = " + "))
),
data = predtest2,
thresholds = seq(0, 1, by = 0.01)
)
plot(testdca_obj, smooth = T, span = 0.5) +
  scale_color_manual(values = c("black", "grey", cols4model)) +
  labs(title = "DCA on testdata") +
  theme(panel.grid = element_blank(),
        legend.position = "inside",
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        text = element_text(family = "serif"))

#############################################################

# 各个模型交叉验证的各折指标点线图
evalcv <- bind_rows(
  lapply(list(evalcv_logistic, evalcv_dt,
              evalcv_lasso, evalcv_ridge, evalcv_enet,
              evalcv_knn, evalcv_lightgbm, evalcv_rf,
              evalcv_xgboost, evalcv_svm, evalcv_mlp), 
         "[[", 
         "evalcv")
) %>%
  mutate(
    model = forcats::as_factor(model),
    modelperf = paste0(model, "(", round(mean, 2),"±",
                       round(sd,2), ")")
  )
evalcv

# ROC
evalcvroc_max <-   evalcv %>% 
  filter(.metric == "roc_auc") %>%
  group_by(id) %>%
  slice_max(.estimate)
evalcvroc_min <-   evalcv %>% 
  filter(.metric == "roc_auc") %>%
  group_by(id) %>%
  slice_min(.estimate)
evalcv %>%
  filter(.metric == "roc_auc") %>%
  ggplot(aes(x = id, y = .estimate, 
             group = modelperf, color = modelperf)) +
  geom_point() +
  geom_line() +
  ggrepel::geom_text_repel(evalcvroc_max, 
                           mapping = aes(label = model), 
                           nudge_y = 0.01,
                           show.legend = F) +
  ggrepel::geom_text_repel(evalcvroc_min, 
                           mapping = aes(label = model), 
                           nudge_y = -0.01,
                           show.legend = F) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = cols4model) +
  labs(x = "", y = "ROCAUC", color = "") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        text = element_text(family = "serif"))

# PR
evalcvpr_max <-   evalcv %>% 
  filter(.metric == "pr_auc") %>%
  group_by(id) %>%
  slice_max(.estimate)
evalcvpr_min <-   evalcv %>% 
  filter(.metric == "pr_auc") %>%
  group_by(id) %>%
  slice_min(.estimate)
evalcv %>%
  filter(.metric == "pr_auc") %>%
  ggplot(aes(x = id, y = .estimate, 
             group = modelperf, color = modelperf)) +
  geom_point() +
  geom_line() +
  ggrepel::geom_text_repel(evalcvpr_max, 
                           mapping = aes(label = model), 
                           nudge_y = 0.01,
                           show.legend = F) +
  ggrepel::geom_text_repel(evalcvpr_min, 
                           mapping = aes(label = model), 
                           nudge_y = -0.01,
                           show.legend = F) +
  scale_y_continuous(limits = c(0.5, 1)) +
  scale_color_manual(values = cols4model) +
  labs(x = "", y = "PRAUC", color = "") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        text = element_text(family = "serif"))

# 各个模型交叉验证的指标平均值图(带上下限)
# ROC
evalcv %>%
  filter(.metric == "roc_auc") %>%
  group_by(model) %>%
  sample_n(size = 1) %>%
  ungroup() %>%
  ggplot(aes(x = model, y = mean, color = model)) +
  geom_point(size = 2, show.legend = F) +
  # geom_line(group = 1) +
  geom_errorbar(aes(ymin = mean-sd, 
                    ymax = mean+sd),
                width = 0.1, 
                linewidth = 1.2,
                show.legend = F) +
  scale_y_continuous(limits = c(0.5, 1)) +
  scale_color_manual(values = cols4model) +
  labs(x = "", y = "cv roc_auc") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        text = element_text(family = "serif"))

# PR
evalcv %>%
  filter(.metric == "pr_auc") %>%
  group_by(model) %>%
  sample_n(size = 1) %>%
  ungroup() %>%
  ggplot(aes(x = model, y = mean, color = model)) +
  geom_point(size = 2, show.legend = F) +
  # geom_line(group = 1) +
  geom_errorbar(aes(ymin = mean-sd, 
                    ymax = mean+sd),
                width = 0.1, 
                linewidth = 1.2,
                show.legend = F) +
  scale_y_continuous(limits = c(0.5, 1)) +
  scale_color_manual(values = cols4model) +
  labs(x = "", y = "cv pr_auc") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        text = element_text(family = "serif"))

