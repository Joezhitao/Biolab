predprob <- function(pred, preddata, etime, estatus, model, dataset, timepoints){
  resulti <- pred$distr[
    1:nrow(preddata)
  ]$survival(timepoints) %>%
    t() %>%
    as.data.frame() %>%
    mutate(crank = pred$crank,
           time = preddata[[etime]],
           status = preddata[[estatus]],
           model = model,
           dataset = dataset)
  return(resulti)
}

cindex_bootci <- function(model, data4pred){
  cindexs <- c()
  for (i in 1:500) {
    set.seed(i)
    samplei <- sample(
      1:nrow(data4pred), nrow(data4pred), replace = T
    )
    datai <- data4pred[samplei, ]
    cindexi <- as.numeric(
      model$predict_newdata(datai)$score(msr("surv.cindex"))
    )
    cindexs <- c(cindexs, cindexi)
  }
  return(quantile(cindexs, c(0.025, 0.975)))
}
  
sadca <- function(predprob, preddata, etime, estatus, model, dataset, timepoints, timepoint, xrange = 0:100 / 100){
  tpat <- which(timepoints == timepoint)
  data.frame(time = preddata[[etime]],
             status = preddata[[estatus]],
             Model = 1-predprob[[tpat]]) %>% 
    dcurves::dca(
      Surv(time, status) ~ Model, 
      data = .,
      time = timepoint, 
      label = list(Model = model), 
      thresholds = xrange
    ) %>%
    plot(smooth = T) +
    scale_color_manual(values = c("black", "grey", "red")) +
    theme(panel.grid = element_blank(), 
          legend.position = "inside",
          legend.justification = c(1,1),
          legend.background = element_blank(),
          legend.key = element_blank(), 
          text = element_text(family = "serif"))
}

eval4sa <- function(
    predprob, preddata, etime, estatus, model, dataset, timepoints,
    plotcalimethod, bw4nne, q4quantile, cutoff
) {
  
  if(cutoff == "median"){
    predprob2 <- predprob %>%
      dplyr::mutate(risk = ifelse(crank > median(crank), "High", "Low")) %>%
      dplyr::mutate(risk = as.factor(risk)) %>%
      as_tibble()
  } else {
    predprob2 <- predprob %>%
      dplyr::mutate(risk = ifelse(crank > cutoff, "High", "Low")) %>%
      dplyr::mutate(risk = as.factor(risk)) %>%
      as_tibble()
  }
  
  riskfit <-  survfit(Surv(time, status) ~ risk, data = predprob2)
  riskplot <- survminer::ggsurvplot(
      riskfit, 
      data = predprob2,
      pval=TRUE,
      pval.coord = c(0.1, 0.5),
      risk.table=TRUE,
      ggtheme = survminer::theme_survminer() +
        theme(text = element_text(family = "serif")),
      font.family = "serif"
    )
  
  survobj <- as.formula(
    paste0('Surv(', etime, ', ', estatus, ') ~ 1')
  )
  
  aucdf <- list()
  rocdf <- list()
  auc2df <- list()
  roc2df <- list()
  bsdf <- list()
  bscdf <- list()
  caldf <- list()
  
  for (i in seq_along(timepoints)) {
    datai <- data.frame(vi = 1-predprob[[i]])
    colnames(datai) <- model
    
    set.seed(1220)
    score_obj <- riskRegression::Score(
      as.list(datai),
      formula = survobj, 
      data = preddata,
      
      metrics=c("auc", "brier"),
      summary = c("IPA"),
      plots = c("roc", "calibrate"),
      
      times = timepoints[i],
      conf.int = T
    )
    aucdf[[i]] <- score_obj$AUC$score %>%
      mutate(dataset = dataset)
    rocdf[[i]] <- score_obj$ROC$plotframe %>%
      mutate(dataset = dataset)
    
    bsdf[[i]] <- score_obj$Brier$score %>%
      mutate(dataset = dataset)
    bscdf[[i]] <- score_obj$Brier$contrasts %>%
      mutate(dataset = dataset)
    plotcali <- riskRegression::plotCalibration(
      score_obj,
      cens.method = "local",
      method = plotcalimethod,
      bandwidth = bw4nne,
      q = q4quantile,
      plot = F
    )
    caldf[[i]] <- plotcali$plotFrames[[1]] %>%
      mutate(dataset = dataset,
             times = timepoints[i],
             model = model)
    
  }
  
  rocdf_plus1 <- data.frame(
    model = model,
    times = timepoints,
    risk = -Inf,
    TPR = 0,
    FPR = 0,
    dataset = dataset
  )

  dataauc <- bind_rows(aucdf)
  dataroc <- bind_rows(rocdf, rocdf_plus1) %>%
    arrange(TPR)
  plotroc <- dataroc %>%
    left_join(dataauc, by = c("model", "times", "dataset")) %>%
    mutate(tauc = paste0("T=", times, ", AUC=", round(AUC, 3),
                         "(", round(lower, 3), "~", 
                         round(upper, 3), ")"),
           tauc = forcats::as_factor(tauc)) %>%
    ggplot(aes(x = FPR, y = TPR, group = tauc, color = tauc)) +
    geom_line(linewidth = 1) +
    geom_abline(color = "grey") +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.2)) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.2)) +
    labs(x = "1-Specificity", y = "Sensitivity", color = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "inside",
          legend.justification = c(1,0),
          legend.background = element_blank(),
          legend.key = element_blank(), 
          text = element_text(family = "serif"))
  
  databs <- bind_rows(bsdf)
  databsc <- bind_rows(bscdf)
  datacal <- bind_rows(caldf)
  rownames(datacal) <- NULL
  plotcal <- datacal %>%
    mutate(t = paste0("T=", times),
           t = forcats::as_factor(t)) %>%
    ggplot(aes(x = Pred, y = Obs, group = t, color = t)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3, 
               alpha = ifelse(plotcalimethod == "nne", 0, 1)) +
    geom_abline() +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.2)) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.2)) +
    labs(x = "Predicted risk", 
         y = "Estimated actual risk", 
         color = "") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "inside",
          legend.justification = c(1,0),
          legend.background = element_blank(),
          legend.key = element_blank(), 
          text = element_text(family = "serif"))
  
  return(list(
    riskplot = riskplot,
    auc = dataauc,
    roc = dataroc,
    rocplot = plotroc,
    brierscore = databs,
    brierscoretest = databsc,
    calibration = datacal,
    calibrationplot = plotcal
  ))
}

viplot <- function(explainer, output_type = "risk"){
  if (output_type == "survival") {
    set.seed(1220)
    vipdata <- survex::model_parts(
      explainer,
      type = "ratio",
      N = 100,
      output_type = "survival"
    )
    vipdata$result %>%
      filter(`_permutation_` == 0) %>%
      rename("Time" = "_times_") %>%
      select(Time, all_of(colnames(explainer$data))) %>%
      pivot_longer(cols = -1) %>%
      mutate(name = tidytext::reorder_within(name, value, Time)) %>%
      ggplot(aes(x = value-1, y = name)) +
      geom_col(fill = "skyblue", color = "grey", width = 0.5) +
      geom_vline(xintercept = 0, color = "black") +
      tidytext::scale_y_reordered() +
      labs(x = "Importance", y = "") +
      facet_wrap(~Time, scales = "free_y", labeller = "label_both") +
      theme_minimal() +
      theme(text = element_text(family = "serif"),
            axis.text.y = element_text(size = 13),
            panel.grid.minor.x = element_blank())
  } else if (output_type == "risk") {
    set.seed(1220)
    vipdata <- survex::model_parts(
      explainer,
      type = "ratio",
      N = 100,
      loss_function = survex::loss_one_minus_c_index,
      output_type = "risk"
    )
    plot(vipdata, 
         show_boxplots = F,
         max_vars = ncol(explainer$data)+1) +
      labs(subtitle = NULL) +
      theme(text = element_text(family = "serif"),
            axis.text.y = element_text(size = 13))
  }
}

pdplot <- function(explainer, vars, output_type = "risk") {
  set.seed(1220)
  partialdata <- survex::model_profile(
    explainer,
    variables = vars,
    type = "partial",
    output_type = output_type
  )
  plot(partialdata, 
       facet_ncol = 2) +
    labs(subtitle = NULL) +
    theme(text = element_text(family = "serif"))
}

shap4one <- function(explainer, dataone, output_type) {
  if (output_type == "survival") {
    set.seed(1220)
    shap41data <- survex::predict_parts(
      explainer, 
      new_observation = dataone,
      N = 100,
      output_type = "survival",
      type = "survshap"
    )
    plot(shap41data, 
         max_vars = ncol(dataone)) +
      labs(subtitle = NULL) +
      scale_color_discrete(guide = guide_legend(title = NULL, ncol = 1)) +
      theme(legend.position = "right",
            legend.text = element_text(size = 13),
            text = element_text(family = "serif"))
  } else if (output_type == "risk") {
    set.seed(1220)
    shap41data <- survex::predict_parts(
      explainer, 
      new_observation = dataone,
      N = 100,
      output_type = "risk",
      type = "shap"
    )
    plot(shap41data, 
         show_boxplots = F,
         max_vars = ncol(dataone)) +
      theme(legend.position = "none",
            text = element_text(family = "serif"))
  }
}



sumshap <- function(explainer, datax, catvars, convars, sampleN){
  set.seed(1220)
  sampleidx <- sample(1:nrow(datax), sampleN)
  samplex <- datax[sampleidx, ]
  sampleshapis <- list()
  jdt <- progress::progress_bar$new(total = sampleN)
  for (i in 1:sampleN) {
    set.seed(1220)
    sampleshapis[[i]] <- survex::predict_parts(
      explainer, 
      new_observation = samplex[i, ],
      N = 100,
      output_type = "risk",
      type = "shap"
    ) %>%
      as.data.frame() %>%
      group_by(variable_name, variable_value) %>%
      summarise(shap = mean(contribution), .groups = "drop") %>%
      mutate(id = i)
    jdt$tick()
  }
  sampleshap <- bind_rows(sampleshapis)
  
  # shap重要性
  shapvipdata <- sampleshap %>%
    group_by(variable_name) %>%
    summarise(meanabsshap = mean(abs(shap)), .groups = "drop") %>%
    mutate(variable_name = reorder(variable_name, meanabsshap))
  shapvipplot <- shapvipdata %>%
    ggplot(aes(x = meanabsshap, y = variable_name)) +
    geom_col(aes(fill = meanabsshap), show.legend = F, color = "grey") +
    labs(y = "", x = "Average |Aggregated SurvSHAP(t)|") +
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    theme_minimal() +
    theme(text = element_text(family = "serif"),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_text(size = 13))
  
  # 分类变量图
  library(ggh4x)
  shapplotd <- sampleshap %>%
    filter(variable_name %in% catvars) %>%
    mutate(variable_name = factor(
      variable_name, 
      levels = rev(levels(shapvipdata$variable_name))
    )) %>%
    ggplot(aes(x = interaction(variable_value, variable_name),
               y = shap)) +
    geom_boxplot(aes(fill = variable_name),
                 show.legend = F) +
    geom_hline(yintercept = 0, color = "grey10") +
    scale_x_discrete(NULL, guide = "axis_nested") +
    scale_colour_viridis_c() +
    labs(x = "", y = "Aggregated SurvSHAP(t)") + 
    theme_bw() +
    theme(text = element_text(family = "serif"),
          axis.text.x = element_text(size = 13))
  
  # 连续变量shap图
  shapplotc <- sampleshap %>%
    filter(variable_name %in% convars) %>%
    mutate(variable_value = as.numeric(variable_value)) %>%
    dplyr::group_by(variable_name) %>%
    dplyr::mutate(
      variable_value = (variable_value - min(variable_value)) / (max(variable_value) - min(variable_value)),
      variable_name = factor(variable_name, levels = levels(shapvipdata$variable_name))
    ) %>%
    dplyr::arrange(variable_value) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x = shap, y = variable_name, color = variable_value)) +
    ggbeeswarm::geom_quasirandom(width = 0.1) +
    geom_vline(xintercept = 0) +
    scale_color_gradient2(
      low = "red", 
      mid = "darkmagenta",
      high = "dodgerblue", 
      midpoint = 0.5,
      breaks = c(0, 1), 
      labels = c("Low", "High"), 
      guide = guide_colorbar(barwidth = 0.5,
                             barheight = length(convars)*2, 
                             ticks = F,
                             title.position = "right",
                             title.hjust = 0.5)
    ) +
    labs(x = "Aggregated SurvSHAP(t)",
         y = "", 
         color = "Feature value") +
    theme_minimal() +
    theme(legend.title = element_text(angle = -90),
          text = element_text(family = "serif"),
          axis.text.y = element_text(size = 13))
  
  return(list(shapley = sampleshap,
              shapvipdata = shapvipdata,
              shapvipplot = shapvipplot,
              shapplotd = shapplotd,
              shapplotc = shapplotc))
}

##############################################################

risk_pred <- function(model, newdata){
  model$predict_newdata(newdata)$crank
} 
surv_pred <- function(model, newdata, times){
  t(model$predict_newdata(newdata)$distr$survival(times))
} 
chf_pred <- function(model, newdata, times){
  t(model$predict_newdata(newdata)$distr$cumHazard(times))
} 

########################################

risk_pred2 <- function(model, newdata){
  model$predict_newdata(newdata)$crank
} 
surv_pred2 <- function(model, newdata, times){
  ers <- (exp(model$predict_newdata(newdata)$crank))
  km0 <- survfit(traindatay ~ 1)
  km01 <- 
    sapply(times, function(x){which.max(km0$time[km0$time<x])})
  km02 <- km0$surv[km01]
  si <- list()
  for (i in seq_along(km02)) {
    si[[i]] <- km02[i] ^ ers
  }
  names(si) <- times
  return(bind_cols(si))
} 

chf_pred2 <- function(model, newdata, times) {
  survex::survival_to_cumulative_hazard(
    surv_pred2(model, newdata, times)
  )
}

########################################

risk_pred3 <- function(model, newdata){
  ri <- predict(
    model,
    newdata %>% select(sort(colnames(newdata))),
    type = "risk"
  )
  return(ri)
}

surv_pred3 <- function(model, newdata, times){
  si <- predict(
    model,
    newdata %>% select(sort(colnames(newdata)))
  )
  sitimes <- as.numeric(colnames(si))
  itpswhere <- 
    sapply(times, function(x){which.max(sitimes[sitimes<x])})
  si2 <- matrix(si[,itpswhere], nrow = nrow(newdata))
  colnames(si2) <- times
  return(si2)
} 

chf_pred3 <- function(model, newdata, times) {
  survex::survival_to_cumulative_hazard(
    surv_pred3(model, newdata, times)
  )
}

#################################################

risk_pred4 <- function(model, newdata){
  if (nrow(newdata)==1) {
    newdata2 <- rbind(newdata, newdata)
    crank <- model$predict_newdata(newdata2)$crank[1]
  } else {
    crank <- model$predict_newdata(newdata)$crank
  }
  return(crank)
} 
surv_pred4 <- function(model, newdata, times){
  ers <- (exp(risk_pred4(model, newdata)))
  km0 <- survfit(traindatay ~ 1)
  km01 <- 
    sapply(times, function(x){which.max(km0$time[km0$time<x])})
  km02 <- km0$surv[km01]
  si <- list()
  for (i in seq_along(km02)) {
    si[[i]] <- km02[i] ^ ers
  }
  names(si) <- times
  return(bind_cols(si))
} 

chf_pred4 <- function(model, newdata, times) {
  survex::survival_to_cumulative_hazard(
    surv_pred4(model, newdata, times)
  )
}
