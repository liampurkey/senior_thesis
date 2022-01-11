#Readme: This R file provides a tidy function to calculate the time weights in Arkhangelsky et al (2020)

sdd.timeweights = function(outcome, treatment_var, unit_var, time_var, time_pre, data) {
  
  #Check that there are no missing observations
  inputs = data %>%
    dplyr::select(tidyselect::all_of(c(outcome, treatment_var, unit_var, time_var))) %>%
    as.matrix()
  
  missing_data = any(is.na(inputs))
  
  if (isTRUE(missing_data)) {
    print("The input data contains missing observations")
    stop()
  }
  
  #Calclulate inputs into convex solver
  target_vec = data %>%
    dplyr::filter(!!as.name(treatment_var) == 0) %>% #Keep comparison units
    dplyr::filter(!!as.name(time_var) > time_pre) %>% #Keep post-treatment periods
    dplyr::select(tidyselect::all_of(c(unit_var, outcome))) %>%
    dplyr::group_by(!!as.name(unit_var)) %>%
    dplyr::summarise(mean_outcome = mean(!!as.name(outcome))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(!!as.name(unit_var)) %>%
    tibble::column_to_rownames(unit_var) %>%
    as.matrix()
    
  input_matrix = data %>%
    dplyr::filter(!!as.name(treatment_var) == 0) %>% #Keep comparison units
    dplyr::filter(!!as.name(time_var) <= time_pre) %>% #Keep post-treatment periods
    dplyr::select(tidyselect::all_of(c(unit_var, time_var, outcome))) %>%
    tidyr::spread(key = !!as.name(time_var), value = !!as.name(outcome)) %>%
    dplyr::mutate(constant = 1) %>%
    dplyr::relocate(constant, .after = !!as.name(unit_var)) %>%
    dplyr::arrange(!!as.name(unit_var)) %>%
    tibble::column_to_rownames(unit_var) %>%
    as.matrix()
  
  #Check that rows (time variable) are in the same order for both inputs
  if (isFALSE(identical(rownames(input_matrix), rownames(target_vec)))) {
    print("Row names do not match")
    stop()
  }
  
  #Create matrix that maps solution vector to time weight vector
  map_matrix = cbind(rep(0, ncol(input_matrix) - 1), diag(ncol(input_matrix) - 1))
  
  #Use CVXR to produce unit weigts
  lambda = CVXR::Variable(ncol(input_matrix))
  objective = CVXR::cvxr_norm(input_matrix %*% lambda - target_vec, 2) 
  constraints = list(sum(map_matrix %*% lambda) == 1, map_matrix %*% lambda >= 0)
  problem = CVXR::Problem(Minimize(objective), constraints)
  result = CVXR::solve(problem)
  if (result$status != "optimal") {
    print("An optimal Solution was not obtained")
    stop()
  }
  lambda_hat = result$getValue(lambda)
  
  #Check that weights sum to one
  if (!(sum(lambda_hat[-1,]) >= 0.99 & sum(lambda_hat[-1,]) <= 1.01)) {
    print("Unit weights do not sum to one")
    stop()
  }
  
  #Format output  
  lambda_out = cbind(colnames(input_matrix[,-1]), lambda_hat[-1,], rep(lambda_hat[1,], length(lambda_hat[-1,]))) %>%
    as.data.frame() %>%
    dplyr::mutate(V1 = as.double(V1)) %>%
    dplyr::mutate(V2 = case_when(as.double(V2) < 0 ~ 0,
                                 as.double(V2) >= 0 ~ as.double(V2))) %>%
    dplyr::mutate(V3 = as.double(V3)) %>%
    dplyr::rename(!!as.name(time_var) := V1) %>%
    dplyr::rename(time_weight = V2) %>%
    dplyr::rename(time_shift = V3) %>%
    dplyr::full_join(data %>% 
                       dplyr::filter(!!as.name(time_var) > time_pre) %>% 
                       dplyr::select(!!as.name(time_var)) %>%
                       dplyr::distinct(), by = time_var) %>%
    dplyr::mutate(time_weight = case_when(is.na(time_weight) ~ 1/nrow(data %>% 
                                                                        dplyr::filter(!!as.name(time_var) > time_pre) %>% 
                                                                        dplyr::select(!!as.name(time_var)) %>%
                                                                        dplyr::distinct()),
                                          !is.na(time_weight) ~ time_weight)) %>%
    dplyr::mutate(time_shift = case_when(is.na(time_shift) ~ 0,
                                         !is.na(time_shift) ~ time_shift))
  
  return(list("lambda_out" = lambda_out, "target_vec" = target_vec, "input_matrix" = input_matrix))
    
}