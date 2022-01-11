#Readme: This script implements the synthetic difference in differences estimator 
#in Arkhangelsky et al (2020)

sdd.estimate = function(outcome, treatment_var, unit_var, time_var, cluster_var = NULL,
                        time_pre, time_init, covariates = NULL, cluster = FALSE, output_name, data) {
  
  #Calculate Regularization Parameter
  zeta = sdd.regularization(outcome = outcome, treatment_var = treatment_var, unit_var = unit_var, time_var = time_var,
                            time_pre = time_pre, covariates = covariates, data = data)
  
  #Calculate Unit Weights
  unit_weight_data = sdd.unitweights(outcome = outcome, treatment_var = treatment_var, unit_var = unit_var, time_var = time_var,
                                 time_pre = time_pre, time_init = time_init, zeta = zeta, covariates = covariates, data = data)$omega_out
  
  #Calculate Unit Weights
  time_weight_data = sdd.timeweights(outcome = outcome, treatment_var = treatment_var, unit_var = unit_var, time_var = time_var,
                                 time_pre = time_pre, data = data)$lambda_out
  
  #Merge Weights with Input Data
  weight_data = data %>%
    inner_join(unit_weight_data, by = "GEOID") %>%
    inner_join(time_weight_data, by = "year") %>%
    mutate(reg_weight = unit_weight*time_weight) 
  
  estimation_data = weight_data %>%
    sdd.preparedata(outcome = outcome, treatment_var = treatment_var, unit_var = unit_var, time_var = time_var, time_pre = time_pre, data = .)
  
  #Estimate Model 
  tau_hat = lm(outcome_twfe ~ treatment_twfe, weights = reg_weight, data = estimation_data)$coefficients[2]
   
  #Estimate Variance
  tau_variance = sdd.jackknife(outcome = outcome, treatment_var = treatment_var, unit_var = unit_var, time_var = time_var, cluster_var = cluster_var, 
                               tau_hat = tau_hat, time_pre = time_pre, cluster = cluster, data = weight_data)
  
  #Calculate Statistics
  output_statistics = sdd.statistics(tau_hat = tau_hat, tau_variance = tau_variance, unit_var = unit_var, time_var = time_var, cluster_var = cluster_var, 
                                     cluster = cluster, data = data)
  
  #Calculate Treatment - Control Balance
  balances = sdd.predictionerror(outcome = outcome, treatment_var = treatment_var, unit_var = unit_var, time_var = time_var, 
                                 time_pre = time_pre, covariates = covariates, data = weight_data)
  
  #Create Plot Data
  plot_data = weight_data %>%
    dplyr::select(all_of(c(outcome, treatment_var, time_var, "unit_weight", "unit_shift", "time_weight"))) %>%
    dplyr::group_by(!!as.name(treatment_var), !!as.name(time_var)) %>%
    dplyr::summarize(outcome_mean = sum((!!as.name(outcome) + unit_shift)*unit_weight)) %>%
    dplyr::mutate(!!as.name(treatment_var) := case_when(!!as.name(treatment_var) == 1 ~ "treated",
                                                        !!as.name(treatment_var) == 0 ~ "control")) %>%
    tidyr::spread(key = !!as.name(treatment_var), value = outcome_mean) %>%
    dplyr::inner_join(weight_data %>% dplyr::select(all_of(c(time_var, "time_weight"))) %>% dplyr::distinct(), by = time_var)
  
  #Create Tex Output
  tex_table = sdd.format(tau_hat = tau_hat, tau_t = output_statistics$t_statistic, tau_p = output_statistics$p_value, differences = balances, 
                         treatment_var = treatment_var, unit_var = unit_var, output_name = output_name, data = weight_data)
  
  #Create Cluster Importance Data
  if (!is.null(cluster_var)) {
    
    cluster_data = sdd.clusters(treatment_var = treatment_var, unit_var = unit_var, cluster_var = cluster_var, data = weight_data)
    
    return(list("tau_hat" = tau_hat, "statistics" = output_statistics, "plot_data" = plot_data, "tex_table" = tex_table, "reg_data" = estimation_data, 
                "cluster_data" = cluster_data))
    
  }
  
  return(list("tau_hat" = tau_hat, "statistics" = output_statistics, "plot_data" = plot_data, "tex_table" = tex_table, "reg_data" = estimation_data))
}