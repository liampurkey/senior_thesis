#Readme: This R file provides utils necessary to implement the synthetic difference in differences 
#estimator (Arkhangelsky et al 2020) in R

#Demean data for weighted two-way fixed effects estimation using lm
sdd.preparedata = function(outcome, treatment_var, unit_var, time_var, time_pre, data) {
  
  prepared_data = data %>%
    dplyr::mutate(treatment_id = !!as.name(treatment_var)) %>%
    dplyr::mutate(!!as.name(treatment_var) := case_when(!!as.name(time_var) <= time_pre ~ 0,
                                                       !!as.name(time_var) > time_pre ~ !!as.name(treatment_var))) %>%
    dplyr::group_by(!!as.name(unit_var)) %>%
    dplyr::mutate(unit_outcome_mean = mean(!!as.name(outcome))) %>%
    dplyr::mutate(unit_treatment_mean = mean(!!as.name(treatment_var))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!as.name(time_var)) %>%
    dplyr::mutate(time_outcome_mean = mean(!!as.name(outcome))) %>%
    dplyr::mutate(time_treatment_mean = mean(!!as.name(treatment_var))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(outcome_mean = mean(!!as.name(outcome))) %>%
    dplyr::mutate(treatment_mean = mean(!!as.name(treatment_var))) %>%
    dplyr::mutate(outcome_twfe = !!as.name(outcome) - unit_outcome_mean - time_outcome_mean + outcome_mean) %>%
    dplyr::mutate(treatment_twfe = !!as.name(treatment_var) - unit_treatment_mean - time_treatment_mean + treatment_mean) %>%
    dplyr::select(-c(unit_outcome_mean, unit_treatment_mean, time_outcome_mean, time_treatment_mean, outcome_mean, treatment_mean))
  
  return(prepared_data)
  
}

#Estimate Jackknife variance
sdd.jackknife = function(outcome, treatment_var, unit_var, time_var, cluster_var = NULL, tau_hat, time_pre, cluster = FALSE, data) {

  if (isFALSE(cluster)) {
    
    units = data %>%
      select(unit_var) %>%
      as.matrix() %>%
      as.vector() %>%
      unique()
    
    jackknife_taus = c()
    
    for (i in 1:length(units)) {
      
      jackknife_data = data %>%
        dplyr::filter(!!as.name(unit_var) != units[i]) %>%
        sdd.preparedata(outcome = outcome, treatment_var = treatment_var, unit_var = unit_var, time_var = time_var, time_pre = time_pre, data = .)
      
      jackknife_tau = lm(outcome_twfe ~ treatment_twfe, data = jackknife_data, weights = reg_weight)$coefficients[2]
      
      jackknife_taus = append(jackknife_taus, jackknife_tau)
      
    }
    
  } else if (isTRUE(cluster)) {
    
    units = data %>%
      select(cluster_var) %>%
      as.matrix() %>%
      as.vector() %>%
      unique()
    
    jackknife_taus = c()
    
    for (i in 1:length(units)) {
      
      jackknife_data = data %>%
        dplyr::filter(!!as.name(cluster_var) != units[i]) %>%
        sdd.preparedata(outcome = outcome, treatment_var = treatment_var, unit_var = unit_var, time_var = time_var, time_pre = time_pre, data = .)
      
      jackknife_tau = lm(outcome_twfe ~ treatment_twfe, data = jackknife_data, weights = reg_weight)$coefficients[2]
      
      jackknife_taus = append(jackknife_taus, jackknife_tau)
      
    }
    
  }
  
  tau_vector = rep(tau_hat, length(units))
  jackknife_variance = ((length(units) - 1)/length(units)) * sum((jackknife_taus - tau_vector)^2)
  
  return(jackknife_variance)
  
}

#Produce statistics 
sdd.statistics = function(tau_hat, tau_variance, unit_var, time_var, cluster_var = NULL, cluster = FALSE, data) {
  
  t_statistic = tau_hat/sqrt(tau_variance)
  
  if (isFALSE(cluster)) {
    t_df = nrow(data) - (nrow(unique(data %>% dplyr::select(all_of(unit_var)) %>% as.matrix)) - 1) - 
      (nrow(unique(data %>% dplyr::select(all_of(time_var)) %>% as.matrix)) - 1) - 2
  } else if (isTRUE(cluster)) {
    t_df = nrow(unique(data %>% dplyr::select(all_of(cluster_var)) %>% as.matrix)) - 1
  }
  
  p_value = stats::pt(-abs(t_statistic), df = t_df, lower.tail = TRUE) + stats::pt(abs(t_statistic), df = t_df, lower.tail = FALSE)
  
  tau_ci = c(tau_hat - 1.96*sqrt(tau_variance), tau_hat + 1.96*sqrt(tau_variance))
  
  return(list("t_statistic" = t_statistic, "t_df" = t_df, "p_value" = p_value, "tau_ci" = tau_ci))
}

#Calculate Normalized Root Mean Squared Error and Covariate Normalized Differences
sdd.predictionerror = function(outcome, treatment_var, unit_var, time_var, time_pre, covariates = NULL, data) {
  
  omega_nrmse = data %>% 
    dplyr::filter(!!as.name(time_var) <= time_pre) %>%
    dplyr::select(all_of(c(outcome, "unit_weight", "unit_shift", treatment_var, unit_var, time_var))) %>%
    dplyr::mutate(!!as.name(treatment_var) := case_when(!!as.name(treatment_var) == 1 ~ "treated",
                                                        !!as.name(treatment_var) == 0 ~ "control")) %>%
    dplyr::group_by(!!as.name(treatment_var), !!as.name(time_var)) %>%
    dplyr::summarize(outcome_mean = sum((!!as.name(outcome) + unit_shift)*unit_weight)) %>%
    tidyr::spread(key = !!as.name(treatment_var), value = outcome_mean) %>%
    dplyr::summarize(output = sqrt(mean((treated - control)^2)) / mean(treated)) %>%
    dplyr::mutate(rownames = "omega") %>%
    tibble::column_to_rownames("rownames") %>%
    as.matrix() 
  
  lambda_nrmse = data %>% 
    dplyr::filter(!!as.name(treatment_var) == 0) %>%
    dplyr::select(all_of(c(outcome, "time_weight", "time_shift", unit_var, time_var))) %>%
    dplyr::mutate(time_post = case_when(!!as.name(time_var) > time_pre ~ "post",
                                        !!as.name(time_var) <= time_pre ~ "pre")) %>%
    dplyr::group_by(time_post, !!as.name(unit_var)) %>%
    dplyr::summarize(outcome_mean = sum((!!as.name(outcome) + time_shift)*time_weight)) %>%
    tidyr::spread(key = time_post, value = outcome_mean) %>%
    dplyr::summarize(output = sqrt(mean((post - pre)^2)) / mean(post)) %>%
    dplyr::mutate(rownames = "lambda") %>%
    tibble::column_to_rownames("rownames") %>%
    as.matrix() 
  
  if (!is.null(covariates)) {
    
    covariates_difference = data %>%
      dplyr::select(all_of(c(unit_var, treatment_var, "unit_weight", covariates))) %>%
      dplyr::distinct() %>%
      dplyr::group_by(!!as.name(treatment_var)) %>%
      dplyr::summarize(across(all_of(covariates), .fns = ~sum(.x * unit_weight))) %>%
      dplyr::ungroup() %>%
      dplyr::summarize(across(all_of(covariates), .fns = ~(lag(.x, order_by = !!as.name(treatment_var)) - .x)/.x)) %>%
      dplyr::filter(!if_any(all_of(covariates), is.na)) %>%
      as.matrix() %>%
      t()
    
    differences = rbind(omega_nrmse, lambda_nrmse, covariates_difference)
    
  } else if (is.null(covariates)) {
    
    differences = rbind(omega_nrmse, lambda_nrmse)
    
  }
    
  return(differences)
  
}

#Format Tex Output
sdd.format = function(tau_hat, tau_t, tau_p, differences, treatment_var, unit_var, output_name, data) {
  
  if (tau_p <= 0.01) {
    tau_hat = paste("$", round(tau_hat, 2), "^{***}$", sep = "")
  } else if (tau_p <= 0.05) {
    tau_hat = paste("$", round(tau_hat, 2), "^{**}$", sep = "")
  } else if (tau_p <= 0.1) {
    tau_hat = paste("$", round(tau_hat, 2), "^{*}$", sep = "")
  } else if (tau_p > 0.1) {
    tau_hat = paste("$", round(tau_hat, 2), "$", sep = "")
  }
  
  tau_hat = as.matrix(tau_hat)
  
  tau_t = as.matrix(paste("(", round(tau_t, 2), ")", sep = ""))
  
  differences = round(differences, 4)
  
  units = as.matrix(nrow(data %>% select(unit_var) %>% distinct()))
  obs = nrow(data)
  
  weight_norm = data %>%
    dplyr::filter(!!as.name(treatment_var) == 0) %>%
    dplyr::select(all_of(c(unit_var, "unit_weight"))) %>%
    dplyr::distinct() %>%
    dplyr::select("unit_weight") %>%
    base::norm("2") %>%
    round(4) %>%
    as.matrix()
  
  table_out = rbind(tau_hat, tau_t, differences, weight_norm, units, obs)
  rownames(table_out) = c(treatment_var, "", rownames(differences), "Weight Norm", "Units", "Observations")
  
  table_out = as.data.frame(table_out) %>%
    tibble::rownames_to_column("Name") %>%
    rename(!!as.name(output_name) := output)
  
  return(table_out)
  
}

#Calculate Cluster-Level Weights
sdd.clusters = function(treatment_var, unit_var, cluster_var, data) {
  
  cluster_data = data %>%
    dplyr::filter(!!as.name(treatment_var) == 0) %>%
    dplyr::select(all_of(c(cluster_var, unit_var, "unit_weight"))) %>%
    dplyr::distinct() %>%
    dplyr::group_by(!!as.name(cluster_var)) %>% 
    dplyr::summarize(cluster_weight = sum(unit_weight)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster_rank = dplyr::dense_rank(desc(cluster_weight))) %>%
    dplyr::filter(cluster_rank <= 10) %>%
    dplyr::select(-cluster_rank) %>%
    dplyr::arrange(desc(cluster_weight)) %>%
    janitor::adorn_totals("row") %>%
    as.data.frame()
  
  return(cluster_data)
  
}
