#Readme: This R file provides a tidy function to calculate the regularization parameter in
#Arkhangelsky et al (2020)

sdd.regularization = function(outcome, treatment_var, unit_var, time_var, time_pre, covariates = NULL, data) {
  
  if (is.null(covariates)) {
    
    zeta = data %>%
      dplyr::filter(!!as.name(time_var) <= time_pre) %>% #Keep pre-treatment periods
      dplyr::select(tidyselect::all_of(c(unit_var, time_var, treatment_var, outcome))) %>% #Select variables to keep
      tidyr::gather(key = "variable", value = "value", -c(!!as.name(unit_var), !!as.name(time_var), !!as.name(treatment_var))) %>% #Transform data to long format
      dplyr::select(-variable) %>% #Drop variable name
      dplyr::group_by(!!as.name(treatment_var), !!as.name(time_var)) %>% #Group by treatment status and year
      dplyr::summarise(mean = mean(value), var = var(value)) %>% #Calculate group means and variances
      dplyr::ungroup() %>%
      dplyr::group_by(!!as.name(time_var)) %>% #Group by year
      dplyr::summarize(normalized_difference = abs((mean - dplyr::lag(mean, order_by = !!as.name(treatment_var)))/
                                                     sqrt((var + dplyr::lag(var, order_by = !!as.name(treatment_var)))/2))) %>% #Calculate normalized difference between treatment groups
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(normalized_difference)) %>% #Drop NAs
      dplyr::summarize(output = 1/max(normalized_difference)) %>% #Set regularization parameter to inverse of maximum normalized difference
      as.matrix() %>%
      as.vector()
    
  } else if (!is.null(covariates)) {
    
    zeta = data %>%
      dplyr::filter(!!as.name(time_var) <= time_pre) %>% #Keep pre-treatment periods
      dplyr::select(tidyselect::all_of(c(unit_var, time_var, treatment_var, outcome, covariates))) %>% #Select variables to keep
      tidyr::gather(key = "variable", value = "value", -c(!!as.name(unit_var), !!as.name(time_var), !!as.name(treatment_var))) %>% #Transform data to long format
      dplyr::mutate(rownames = dplyr::case_when(variable == outcome ~ as.character(year), #Create rowname variable equal to outcome year or covariate name
                                                variable != outcome ~ variable)) %>%
      dplyr::select(-!!as.name(time_var), -variable) %>% 
      dplyr::distinct() %>% #Drop duplicates of time-invariant covariates
      dplyr::group_by(!!as.name(treatment_var), rownames) %>% #Group by treatment status and variable name
      dplyr::summarise(mean = mean(value), var = var(value)) %>% #Calculate group means and variances
      dplyr::ungroup() %>%
      dplyr::group_by(rownames) %>% #Group by variable name
      dplyr::summarize(normalized_difference = abs((mean - dplyr::lag(mean, order_by = !!as.name(treatment_var)))/
                                                     sqrt((var + dplyr::lag(var, order_by = !!as.name(treatment_var)))/2))) %>% #Calculate normalized difference between treatment groups
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(normalized_difference)) %>% #Drop NAs 
      summarize(output = 1/max(normalized_difference)) %>% #Set regularization parameter to inverse of maximum normalized difference
      as.matrix() %>%
      as.vector()
    
  }
   
  return(zeta)
  
}

