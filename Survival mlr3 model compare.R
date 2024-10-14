#Joe_2024_10_12
rm(list = ls())
# 加载包
library(tidyverse)

# 加载各个模型的评估结果
evalfiles <- list.files("F:/Mach_learn_data/mlr_model/", full.names = T)
lapply(evalfiles, load, .GlobalEnv)

# 横向比较的模型个数
nmodels <- 11
cols4model <- rainbow(nmodels)  # 模型统一配色

#############################################################

# 各个模型在训练集上的评估结果
# 各个模型对训练集各个时点的AUC
trainauc_all <- bind_rows(
  lapply(list(evaltrain_coxph, evaltrain_ctree, 
              evaltrain_gbm, evaltrain_lasso, 
              evaltrain_ridge, evaltrain_enet, 
              evaltrain_nn, evaltrain_rsf,
              evaltrain_svm, evaltrain_xgboost, evaltrain_stack), 
         "[[", 
         "auc")
)
trainauc_all

trainauc_min <- trainauc_all %>%
  group_by(times) %>%
  slice_min(AUC, n = 1)
trainauc_max <- trainauc_all %>%
  group_by(times) %>%
  slice_max(AUC, n = 1)

trainauc_all %>%
  ggplot(aes(x = times, y = AUC, color = model)) +
  geom_point() +
  geom_line(aes(group = model)) +
  ggrepel::geom_text_repel(trainauc_max,
                           mapping = aes(label = model),
                           nudge_y = 0.05,
                           show.legend = F) +
  ggrepel::geom_text_repel(trainauc_min,
                           mapping = aes(label = model),
                           nudge_y = -0.05,
                           show.legend = F) +
  scale_x_continuous(breaks = unique(trainauc_all$times)) +
  scale_color_manual(values = c(cols4model)) +
  labs(x = "Time", color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "serif"))

#############################

# 各个模型对训练集各个时点的BS
trainbs_all <- bind_rows(
  lapply(list(evaltrain_coxph, evaltrain_ctree, 
              evaltrain_gbm, evaltrain_lasso, 
              evaltrain_ridge, evaltrain_enet, 
              evaltrain_nn, evaltrain_rsf,
              evaltrain_svm, evaltrain_xgboost, evaltrain_stack), 
         "[[", 
         "brierscore")
)
trainbs_all

trainbs_min <- trainbs_all %>%
  filter(model != "Null model") %>%
  group_by(times) %>%
  slice_min(Brier, n = 1)
trainbs_max <- trainbs_all %>%
  filter(model != "Null model") %>%
  group_by(times) %>%
  slice_max(Brier, n = 1)

trainbs_all %>%
  filter(model != "Null model") %>%
  ggplot(aes(x = times, y = Brier, color = model)) +
  geom_point() +
  geom_line(aes(group = model)) +
  ggrepel::geom_text_repel(trainbs_max,
                           mapping = aes(label = model),
                           nudge_y = 0.05,
                           show.legend = F) +
  ggrepel::geom_text_repel(trainbs_min,
                           mapping = aes(label = model),
                           nudge_y = -0.05,
                           show.legend = F) +
  scale_x_continuous(breaks = unique(trainbs_all$times)) +
  scale_color_manual(values = c(cols4model)) +
  labs(x = "Time") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "serif"))


##################################

# 各个模型在训练集上的ROC
rocsall <- bind_rows(
  lapply(list(evaltrain_coxph, evaltrain_ctree, 
              evaltrain_gbm, evaltrain_lasso, 
              evaltrain_ridge, evaltrain_enet, 
              evaltrain_nn, evaltrain_rsf,
              evaltrain_svm, evaltrain_xgboost, evaltrain_stack), 
         "[[", 
         "roc")
)
rocsall %>%
  left_join(trainauc_all, by = c("model", "times")) %>%
  filter(times == 1825) %>% # 时间点可以替换
  mutate(mtAUC = paste0(model, ", T=", times, 
                        ", AUC=", round(AUC, 4))) %>%
  ggplot(aes(x = FPR,
             y = TPR, 
             color = forcats::as_factor(mtAUC))) +
  geom_path(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = c(cols4model)) +
  labs(color = "", x = "1 - Specificity", y = "Sensitivity") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = "serif"))

#################################

# 各个模型在训练集上的DCA
predprobtrain_all <- bind_rows(
  predprobtrain_coxph, 
  predprobtrain_ctree, 
  predprobtrain_gbm, 
  predprobtrain_lasso, 
  predprobtrain_ridge, 
  predprobtrain_enet, 
  predprobtrain_nn, 
  predprobtrain_rsf,
  predprobtrain_svm,
  predprobtrain_xgboost,
  predprobtrain_stack
)
predprobtrain_all

predprobtrain_all %>%
  select(`1825`, time, status, model) %>%  # 时间点1825可以替换
  mutate(id = rep(1:nrow(predprobtrain_coxph), nmodels),
         `1825` = 1 - `1825`) %>%
  pivot_wider(names_from = model, values_from = `1825`) %>%
  select(-id) %>%
  dcurves::dca(
    survival::Surv(time, status) ~ ., 
    data = .,
    time = 1825, # 60跟上面对应
    thresholds = 0:100 / 100  # 范围可以改
  ) %>%
  plot(smooth = T, span = 0.5) +
  scale_color_manual(values = c("black", "grey", cols4model)) +
  labs(title = "DCA on traindata") +
  theme(panel.grid = element_blank(), 
        legend.position = "inside",
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        text = element_text(family = "serif"))



#############################################################

# 各个模型在测试集上的评估结果
names(evaltest_svm)
# 各个模型对测试集各个时点的AUC
testauc_all <- bind_rows(
  lapply(list(evaltest_coxph, evaltest_ctree, 
              evaltest_gbm, evaltest_lasso, 
              evaltest_ridge, evaltest_enet, 
              evaltest_nn, evaltest_rsf,
              evaltest_svm, evaltest_xgboost, evaltest_stack), 
         "[[", 
         "auc")
)
testauc_all

testauc_min <- testauc_all %>%
  group_by(times) %>%
  slice_min(AUC, n = 1)
testauc_max <- testauc_all %>%
  group_by(times) %>%
  slice_max(AUC, n = 1)

testauc_all %>%
  ggplot(aes(x = times, y = AUC, color = model)) +
  geom_point() +
  geom_line(aes(group = model)) +
  ggrepel::geom_text_repel(testauc_max,
                           mapping = aes(label = model),
                           nudge_y = 0.05,
                           show.legend = F) +
  ggrepel::geom_text_repel(testauc_min,
                           mapping = aes(label = model),
                           nudge_y = -0.05,
                           show.legend = F) +
  scale_x_continuous(breaks = unique(testauc_all$times)) +
  scale_color_manual(values = c(cols4model)) +
  labs(x = "Time", color = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "serif"))

#############################

# 各个模型对测试集各个时点的BS
testbs_all <- bind_rows(
  lapply(list(evaltest_coxph, evaltest_ctree, 
              evaltest_gbm, evaltest_lasso, 
              evaltest_ridge, evaltest_enet,
              evaltest_nn, evaltest_rsf,
              evaltest_svm, evaltest_xgboost, evaltest_stack), 
         "[[", 
         "brierscore")
)
testbs_all

testbs_min <- testbs_all %>%
  filter(model != "Null model") %>%
  group_by(times) %>%
  slice_min(Brier, n = 1)
testbs_max <- testbs_all %>%
  filter(model != "Null model") %>%
  group_by(times) %>%
  slice_max(Brier, n = 1)

testbs_all %>%
  filter(model != "Null model") %>%
  ggplot(aes(x = times, y = Brier, color = model)) +
  geom_point() +
  geom_line(aes(group = model)) +
  ggrepel::geom_text_repel(testbs_max,
                           mapping = aes(label = model),
                           nudge_y = 0.05,
                           show.legend = F) +
  ggrepel::geom_text_repel(testbs_min,
                           mapping = aes(label = model),
                           nudge_y = -0.05,
                           show.legend = F) +
  scale_x_continuous(breaks = unique(testbs_all$times)) +
  scale_color_manual(values = c(cols4model)) +
  labs(x = "Time") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(family = "serif"))


##################################

# 各个模型在测试集上的ROC
rocsall <- bind_rows(
  lapply(list(evaltest_coxph, evaltest_ctree, 
              evaltest_gbm, evaltest_lasso, 
              evaltest_ridge, evaltest_enet,
              evaltest_nn, evaltest_rsf,
              evaltest_svm, evaltest_xgboost, evaltest_stack), 
         "[[", 
         "roc")
)
rocsall %>%
  left_join(testauc_all, by = c("model", "times")) %>%
  filter(times == 1825) %>% # 时间点可以替换
  mutate(mtAUC = paste0(model, ", T=", times, 
                        ", AUC=", round(AUC, 4))) %>%
  ggplot(aes(x = FPR,
             y = TPR, 
             color = forcats::as_factor(mtAUC))) +
  geom_path(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = c(cols4model)) +
  labs(color = "", x = "1 - Specificity", y = "Sensitivity") +
  theme_bw() +
  theme(legend.position = "inside",
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = "serif"))

#################################

# 各个模型在测试集上的DCA
predprobtest_all <- bind_rows(
  predprobtest_coxph, 
  predprobtest_ctree, 
  predprobtest_gbm, 
  predprobtest_lasso, 
  predprobtest_ridge, 
  predprobtest_enet, 
  predprobtest_nn, 
  predprobtest_rsf,
  predprobtest_svm,
  predprobtest_xgboost,
  predprobtest_stack
)
predprobtest_all

predprobtest_all %>%
  select(`1825`, time, status, model) %>%  # 时间点1825可以替换
  mutate(id = rep(1:nrow(predprobtest_coxph), nmodels),
         `1825` = 1 - `1825`) %>%
  pivot_wider(names_from = model, values_from = `1825`) %>%
  select(-id) %>%
  dcurves::dca(
    survival::Surv(time, status) ~ ., 
    data = .,
    time = 1825, # 60跟上面对应
    thresholds = 0:100 / 100  # 范围可以改
  ) %>%
  plot(smooth = T, span = 0.5) +
  scale_color_manual(values = c("black", "grey", cols4model)) +
  labs(title = "DCA on testdata") +
  theme(panel.grid = element_blank(), 
        legend.position = "inside",
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.key = element_blank(), 
        text = element_text(family = "serif"))

