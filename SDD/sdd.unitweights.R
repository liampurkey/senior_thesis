#Readme: This R file provides a tidy function to calculate the unit weights in Arkhangelsky et al (2020)
#with the option to include covariates

sdd.unitweights = function(outcome, treatment_var, unit_var, time_var, time_pre, time_init, zeta, covariates = NULL, data) {
  
  #The sdd.unitweights takes in a balanced panel and regularization
  #parameter and computes the SDD unit weights by passing
  #a (t + c) x (i + 1) matrix of outcomes and a (t + c) x 1 vector of average outcomes 
  #to R's CVXR solve
  
  #Check that there are no missing observations
  if (is.null(covariates)) {
    
    inputs = data %>%
      dplyr::select(tidyselect::all_of(c(outcome, treatment_var, unit_var, time_var))) %>%
      as.matrix()
    
    missing_data = any(is.na(inputs))

  } else if (!is.null(covariates)) {
    
    inputs = data %>%
      dplyr::select(tidyselect::all_of(c(outcome, treatment_var, unit_var, time_var, covariates))) %>%
      as.matrix()
    
    missing_data = any(is.na(inputs))
    
  }

  if (isTRUE(missing_data)) {
    print("The input data contains missing observations")
    stop()
  }
  
  #Prepare CVXR inputs without covariates
  if (is.null(covariates)) {
    
    #Create target vector
    target_vec = data %>%
      dplyr::filter(!!as.name(treatment_var) == 1) %>% #Keep treated units
      dplyr::filter(!!as.name(time_var) <= time_pre) %>% #Keep pre-treatment periods
      dplyr::select(tidyselect::all_of(c(time_var, outcome))) %>% #Select time variable and outcome
      dplyr::group_by(!!as.name(time_var)) %>% 
      dplyr::summarize(mean_outcome = mean(!!as.name(outcome))) %>% #Calculate average pre-treatment outcome for each pre-treatment period
      dplyr::ungroup() %>%
      dplyr::arrange(!!as.name(time_var)) %>% #Order rows from initial pre-treatment period to final pre-treatment period
      tibble::column_to_rownames(time_var) %>% #Set rownames to time variable
      as.matrix()
    
    #Create input matrix
    input_matrix = data %>%
      dplyr::filter(!!as.name(treatment_var) == 0) %>% #Keep control units
      dplyr::filter(!!as.name(time_var) <= time_pre) %>% #Keep pre-treatment periods
      dplyr::select(tidyselect::all_of(c(unit_var, time_var, outcome))) %>% #Select unit variable, time variable and outcome
      tidyr::spread(key = !!as.name(unit_var), value = !!as.name(outcome)) %>% #Transform data to T x N_co data frame
      dplyr::mutate(constant = 1) %>% #Add constant to data
      dplyr::relocate(constant, .after = !!as.name(time_var)) %>% #Set location of constant for input into solver
      dplyr::arrange(!!as.name(time_var)) %>% #Order rows from initial pre-treatment period to final pre-treatment period 
      tibble::column_to_rownames(time_var) %>% #Set rownames to time variable
      as.matrix()
    
    #Check that rows (time variable) are in the same order for both inputs
    if (isFALSE(identical(rownames(input_matrix), rownames(target_vec)))) {
      print("Row names do not match")
      stop()
    }
    
  } else if (!is.null(covariates)) {
    
    #Create target vector
    target_vec = data %>% 
      dplyr::filter(!!as.name(treatment_var) == 1) %>% #Keep treated units
      dplyr::filter(!!as.name(time_var) <= time_pre) %>% #Keep pre-treatment periods
      dplyr::select(tidyselect::all_of(c(time_var, outcome, covariates))) %>% #Select time variable, outcome, and covariates
      dplyr::group_by(!!as.name(time_var)) %>%
      dplyr::summarise(dplyr::across(tidyselect::all_of(c(outcome, covariates)), .fns = mean)) %>% #Summarize average pre-treatment outcome and covariates for each pre-treatment period
      dplyr::ungroup() %>%
      tidyr::gather(key = "variable", value = "value", -!!as.name(time_var)) %>% #Collapse data to long format while keeping time indicator variable 
      dplyr::mutate(rownames = dplyr::case_when(variable == outcome ~ as.character(!!as.name(time_var)),
                                                variable != outcome ~ variable)) %>% #Set rownames variable to outcome year or covariate name
      dplyr::select(-!!as.name(time_var), -variable) %>%
      dplyr::distinct() %>% #Drop duplicates of time-invariant covariates
      dplyr::arrange(rownames) %>% #Order rows from initial pre-treatment period to final pre-treatment period and then by covariate name
      tibble::column_to_rownames("rownames") %>% #Set rownames to rownames variable 
      as.matrix()
    
    #Create input matrix
    input_matrix = data %>%
      dplyr::filter(!!as.name(treatment_var) == 0) %>% #Keep treated units
      dplyr::filter(!!as.name(time_var) <= time_pre) %>% #Keep pre-treatment periods
      dplyr::select(tidyselect::all_of(c(unit_var, time_var, outcome, covariates))) %>% #Select unit variable, time variable, outcome, and covariates
      tidyr::gather(key = "variable", value = "value", -c(!!as.name(unit_var), !!as.name(time_var))) %>% #Transform to long format while keeping unit and time indicators
      dplyr::mutate(rownames = dplyr::case_when(variable == outcome ~ as.character(year), 
                                                variable != outcome ~ variable)) %>% #Set rownames variable to outcome year or covariate name
      dplyr::select(-!!as.name(time_var), -variable) %>%
      dplyr::distinct() %>% #Drop duplicates of time-invariant covariates
      tidyr::spread(key = !!as.name(unit_var), value = value) %>% #Transform data to (t + c) x N_co matrix
      dplyr::mutate(constant = ifelse(rownames %in% covariates, 0, 1)) %>% #Add constant to data frame
      dplyr::relocate(constant, .after = rownames) %>% #Set location of constant for input into solver
      dplyr::arrange(rownames) %>% #Order rows from initial pre-treatment period to final pre-treatment period and then by covariate name
      tibble::column_to_rownames("rownames") %>% #Set rownames to rownames variable 
      as.matrix()
      
    #Check that rows (time variable and covariates) are in the same order for both inputs
    if (isFALSE(identical(rownames(input_matrix), rownames(target_vec)))) {
      print("Row names do not match")
      stop()
    }
    
    #Create inverse-variance weight vector
    inverse_variance_weights = as.data.frame(cbind(rownames(input_matrix), input_matrix[,-1]))%>% 
      dplyr::rename(rownames = V1) %>%
      tidyr::gather(key = unit_var, value = "value", -rownames) %>%
      dplyr::group_by(rownames) %>%
      dplyr::summarize(inverse_variance = 1/var(value)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(rownames) %>% #Order rows from initial pre-treatment period to final pre-treatment period and then by covariate name
      tibble::column_to_rownames(var = "rownames") %>% 
      as.matrix()
    
    #Check that rows (time variable and covariates) are in the same order for input matrix and inverse variance weights
    if (isFALSE(identical(rownames(inverse_variance_weights), rownames(input_matrix)))) {
      print("Row names do not match")
      stop()
    }
    
  }
  
  #Create matrix that maps solution vector to unit weight vector
  map_matrix = cbind(rep(0, ncol(input_matrix) - 1), diag(ncol(input_matrix) - 1))
  
  #Use CVXR to produce unit weigts
  omega = CVXR::Variable(ncol(input_matrix))
  if (is.null(covariates)) {
    objective = sum((input_matrix %*% omega - target_vec)^2) + zeta*(time_pre - time_init)*sum((map_matrix %*% omega)^2)
    constraints = list(sum(map_matrix %*% omega) == 1, map_matrix %*% omega >= 0)
    problem = CVXR::Problem(Minimize(objective), constraints)
    result = CVXR::solve(problem)
    if (result$status != "optimal") {
      print("An optimal Solution was not obtained")
      stop()
    }
    omega_hat = result$getValue(omega)
  } else if (!is.null(covariates)) {
    objective = sum(inverse_variance_weights*(input_matrix %*% omega - target_vec)^2) + zeta*(time_pre - time_init)*sum((map_matrix %*% omega)^2)
    constraints = list(sum(map_matrix %*% omega) == 1, map_matrix %*% omega >= 0)
    problem = CVXR::Problem(Minimize(objective), constraints)
    result = CVXR::solve(problem)
    if (result$status != "optimal") {
      print("An optimal Solution was not obtained")
      stop()
    }
    omega_hat = result$getValue(omega)
  }
  
  #Check that weights sum to one
  if (!(sum(omega_hat[-1,]) >= 0.99 & sum(omega_hat[-1,]) <= 1.01)) {
    print("Unit weights do not sum to one")
    stop()
  }
    
  #Format Output
  omega_out = cbind(colnames(input_matrix[,-1]), omega_hat[-1,], rep(omega_hat[1,], length(omega_hat[-1,]))) %>%
    as.data.frame() %>%
    dplyr::mutate(V1 = as.double(V1)) %>%
    dplyr::mutate(V2 = case_when(as.double(V2) < 0 ~ 0,
                          as.double(V2) >= 0 ~ as.double(V2))) %>%
    dplyr::mutate(V3 = as.double(V3)) %>%
    dplyr::rename(!!as.name(unit_var) := V1) %>%
    dplyr::rename(unit_weight = V2) %>%
    dplyr::rename(unit_shift = V3) %>%
    dplyr::full_join(data %>% 
                       dplyr::filter(!!as.name(treatment_var) == 1) %>% 
                       dplyr::select(!!as.name(unit_var)) %>%
                       dplyr::distinct(), by = unit_var) %>%
    dplyr::mutate(unit_weight = case_when(is.na(unit_weight) ~ 1/nrow(data %>% 
                                                                          dplyr::filter(!!as.name(treatment_var) == 1) %>% 
                                                                          dplyr::select(!!as.name(unit_var)) %>%
                                                                          dplyr::distinct()),
                                          !is.na(unit_weight) ~ unit_weight)) %>%
    dplyr::mutate(unit_shift = case_when(is.na(unit_shift) ~ 0,
                                         !is.na(unit_shift) ~ unit_shift))
  
  return(list("omega_out" = omega_out, "target_vec" = target_vec, "input_matrix" = input_matrix))
}