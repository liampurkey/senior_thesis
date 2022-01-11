#The ps.estimate function is my implementation of the propenity score algorithm
#presented in Imbens (2014) (Matching Methods in Practice: Three Examples).


#Function to Determine Best Single Additional Variable Model from Set of Covariates
ps.best.single.term = function(outcome, baseline.covariates, covariates, data, clin = 1, cqua = 2.7, linear) {
  
  #Calculate Baseline Model
  baseline.formula = paste(outcome, paste(baseline.covariates, collapse = " + "), sep = " ~ ")
  baseline.model = glm(baseline.formula, family="binomial"(link = "logit"), data = data)
  
  #Initialize Null Condition
  best.variable = NULL
  
  if (isTRUE(linear)) {
    
    #Select Linear Terms to Include
    for(i in 1:length(covariates)) {
      
      #Estimate Test Model
      test.model = glm(paste(baseline.formula, covariates[i], sep = " + "), family="binomial"(link = "logit"), data = data)
      
      #Calculate Likelihood Ratio
      lr = lrtest(baseline.model, test.model)
      
      #Determine if Covariate Meets Matching Cutoff 
      if (is.null(best.variable) & lr$Chisq[2] < clin) {
        break
      } else if (is.null(best.variable) & lr$Chisq[2] >= clin) {
        best.variable = covariates[i]
        lrstat = lr$Chisq[2]
      } else if (!(is.null(best.variable)) & lr$Chisq[2] > lrstat) {
        best.variable = covariates[i]
        lrstat = lr$Chisq[2]
      }
    }
  } else if (isFALSE(linear)) {
    
    #Select Linear Terms to Include
    for(i in 1:length(covariates)) {
      
      #Estimate Test Model
      test.model = glm(paste(baseline.formula, covariates[i], sep = " + "), family="binomial"(link = "logit"), data = data)
      
      #Calculate Likelihood Ratio
      lr = lrtest(baseline.model, test.model)
      
      #Determine if Covariate Meets Matching Cutoff 
      if (is.null(best.variable) & lr$Chisq[2] < cqua) {
        break
      } else if (is.null(best.variable) & lr$Chisq[2] >= cqua) {
        best.variable = covariates[i]
        lrstat = lr$Chisq[2]
      } else if (!(is.null(best.variable)) & lr$Chisq[2] > lrstat) {
        best.variable = covariates[i]
        lrstat = lr$Chisq[2]
      }
    }
  }
  
  
  return(best.variable)
  
}

#Recursively determine covariates to include by comparing log likelihoods 
ps.best.terms = function(outcome, baseline.covariates, covariates, data, clin = 1, cqua = 2.7, linear) {
  
  term = ps.best.single.term(outcome = outcome, baseline.covariates = baseline.covariates,
                             covariates = covariates, data = data, clin = clin, cqua = cqua, linear = linear)
  
  if(is.null(term)) {
    
    return(c(baseline.covariates, term))
    
  } else if (!is.null(term)) {
    
    baseline.covariates = c(baseline.covariates, term)
    covariates = setdiff(covariates, term)
    
    return(ps.best.terms(outcome = outcome, baseline.covariates = baseline.covariates, 
                         covariates = covariates, data = data, clin = clin, cqua = cqua, linear = linear))
  }
  
}

#Calculate Propensity Score 
ps.estimate = function(outcome, baseline.covariates = NULL, covariates, data, clin = 1, cqua = 2.7) {
  
  if (is.null(baseline.covariates)) {
    baseline.covariates = c("1")
  } 
  
  #Calculate Optimal Linear Terms 
  linear.terms = ps.best.terms(outcome = outcome, baseline.covariates = baseline.covariates, 
                               covariates = covariates, data = data, clin = clin, cqua = cqua, linear = TRUE)
  
  #Calculate Possible Interactions and Squared Terms 
  interactions = combn(linear.terms, 2)
  interactions = paste(interactions[1,], interactions[2,], sep = "*")
  squared = paste(linear.terms, linear.terms, sep = "*")
  quadratic.choices = c(interactions, squared)
  
  #Calculate Optimal Quadratic Terms
  quadratic.terms = ps.best.terms(outcome = outcome, baseline.covariates = linear.terms, 
                                  covariates = quadratic.choices, data = data, clin = clin, cqua = cqua, linear = FALSE)
  
  #Join Terms
  pscore.terms = unique(c(linear.terms, quadratic.terms))
  
  #Estimate Propensity Score
  pscore.model = glm(paste(outcome, paste(pscore.terms, collapse = " + "), sep = " ~ "), family="binomial"(link = "logit"), data = data)
  
  #Return Propensity Score
  return(pscore.model)
}