#Readme: This script contains functions used in Purkey (2021). This file is modefied to construct 
#radii around gentrifying neighborhoods but not around non-gentrifying neighborhoods. The original
#functions are found in the Depricated folder.

###Adjust for Inflation

to_2010_dollars = function(dataset, variables) {
  for (i in 1:length(variables)) {
    varname = variables[i]
    dataset = dataset %>%
      mutate(!!varname := !!as.name(varname) * CPI_ratio)
  }
  return(dataset)
}

### Normalize Variables

normalize_drop = function(dataset, variables, total, newname = NULL) {
  #inputs must be strings
  if (is.null(newname)) {
    for (i in 1:length(variables)) {
      var = variables[i]
      varname = paste("p", var, sep = "_")
      dataset = dataset %>% 
        mutate(!!varname := !!as.name(var)/!!as.name(total)) %>%
        select(-!!as.name(var))
    }
  } else {
    for (i in 1:length(variables)) {
      var = variables[i]
      varname = paste(newname, var, sep = "_")
      dataset = dataset %>% 
        mutate(!!varname := !!as.name(var)/!!as.name(total)) %>%
        select(-!!as.name(var))
    }
  }
  return(dataset)
}

normalize = function(dataset, variables, total, newname = NULL) {
  #inputs must be strings
  if (is.null(newname)) {
    for (i in 1:length(variables)) {
      var = variables[i]
      varname = paste("p", var, sep = "_")
      dataset = dataset %>% 
        mutate(!!varname := !!as.name(var)/!!as.name(total)) 
    }
  } else {
    for (i in 1:length(variables)) {
      var = variables[i]
      varname = paste(newname, var, sep = "_")
      dataset = dataset %>% 
        mutate(!!varname := !!as.name(var)/!!as.name(total))
    }
  }
  return(dataset)
}

### Sum Multiple Columns

sumcols_drop = function(dataset, variables, varname) {
  dataset = dataset %>%
    mutate(!!as.name(varname) := rowSums(dataset[,variables])) 
  return(dataset)
}

sumcols = function(dataset, variables, varname) {
  dataset = dataset %>%
    mutate(!!as.name(varname) := rowSums(across(variables)))
  return(dataset)
}

### Create Percentage Changes

pc_changes = function(dataset, variables) {
  dataset = dataset %>%
    group_by(GEOID)
    
  for (i in 1:length(variables)) {
    var = variables[i]
    pc_varname = paste("pc", var, sep = "_")
    dataset = dataset %>%
      mutate(!!pc_varname := (!!as.name(var) - lag(!!as.name(var)))/lag(!!as.name(var)))
  }
  
  dataset = dataset %>%
    ungroup()
  
  return(dataset)
  
}

create_pc_changes = function(dataset, variables, newname = NULL, start_year, end_year) {
  
  start_data = dataset %>%
    filter(year == start_year) %>%
    select(-c(year, tract_area)) %>%
    setNames(paste0(names(.), "_2010")) %>%
    rename(GEOID = GEOID_2010) %>%
    rename(City = City_2010) %>%
    rename(State = State_2010)
  
  end_data = dataset %>%
    filter(year == end_year) %>%
    select(-c(year, tract_area)) %>%
    setNames(paste0(names(.), "_2018")) %>%
    rename(GEOID = GEOID_2018) %>%
    rename(City = City_2018) %>%
    rename(State = State_2018) 
  
  pc_changes_data = pc_changes(dataset, variables)
  pc_changes_vars = names(pc_changes_data)
  
  changes_vars = grep("pc_", pc_changes_vars, value = TRUE)
  
  keep = c("GEOID", "City", "State", "tract_area", changes_vars)
  
  pc_changes_data = pc_changes_data %>% 
    filter(year == end_year) %>%
    select(all_of(keep))
  
  data = inner_join(pc_changes_data, start_data, by = c("GEOID", "City", "State"))
  data = inner_join(data, end_data, by = c("GEOID", "City", "State"))
  
  return(data)
}

### Clean Units

clean_units <- function(x){
  attr(x,"units") <- NULL
  class(x) <- setdiff(class(x),"units")
  return(x)
}

### Make Maps 

map_gen = function(dataset, cities, radius) {
  
  for (i in 1:length(cities)) {
    
    city = cities[i]
    
    map_data = dataset %>%
      filter(City == city)
    
    gen_data = cities_gen %>%
      filter(City == city)
    
    map = tm_shape(map_data) + tm_polygons("edugentrify", style = "cat", 
                                           palette = c("#DEEBF7", "#FFECEB", "#FDBEB9"), 
                                           labels = c("Not in Sample", "Not Gentrifying", "Gentrifying"), 
                                           title = "Neighborhoods") + 
      tm_layout(legend.outside = TRUE, legend.text.size = .8, legend.title.size = 1) +
      tm_shape(gen_data) + tm_polygons("edugentrify", style = "cat", 
                                       palette = c("red4"), 
                                       labels = c("Gentrifying Tract"), 
                                       title = "Tracts")
    
    outputdir = paste("/Users/liampurkey/Desktop/Honors/Results/Maps", city, sep = "/")
    
    mapname = paste(city, "gen", radius, sep = "_")
    mapname = paste(mapname, "png", sep = ".")
      
    mappath = paste(outputdir, mapname, sep = "/") 
      
    tmap_save(tm = map, filename = mappath)
    
  }
  
}

#Calculate Differences and T-statistics 

test.diff = function(data, variables) {
  
  test_vec = c()
  
  data_1 = data %>%
    filter(edugentrify == 1)
  
  data_0 = data %>%
    filter(edugentrify == 0)
  
  for (i in 1:length(variables)) {
    
    vec_1 = data_1 %>%
      select(variables[i]) %>%
      filter(!is.na(!!as.name(variables[i]))) %>%
      filter(!is.nan(!!as.name(variables[i]))) %>%
      filter(!is.infinite(!!as.name(variables[i]))) %>%
      as.matrix() %>%
      as.vector() 
    
    vec_0 = data_0 %>%
      select(variables[i]) %>%
      filter(!is.na(!!as.name(variables[i]))) %>%
      filter(!is.nan(!!as.name(variables[i]))) %>%
      filter(!is.infinite(!!as.name(variables[i]))) %>%
      as.matrix() %>%
      as.vector()
    
    test = t.test(vec_1, vec_0)
    diff = round(test$estimate[1] - test$estimate[2], 2)
    t = paste("(", round(test$statistic, 2), ")", sep = "")
    
    if (test$p.value <= 0.01) {
      diff = paste("$", diff, "^{***}$", sep = "")
    } else if (test$p.value <= 0.05) {
      diff = paste("$", diff, "^{**}$", sep = "")
    } else if (test$p.value <= 0.1) {
      diff = paste("$", diff, "^{*}$", sep = "")
    }
    
    out = c(diff, t)
    names(out) = c(variables[i], variables[i])
      
    test_vec = append(test_vec, out)
    
  }
  return(test_vec)
}

spread.variables = function(variable, years, unit_var, time_var, data) {
  
  for (i in 1:length(years)) {
    
    data = data %>%
      dplyr::group_by(!!as.name(unit_var)) %>%
      mutate(!!as.name(paste(variable, years[i], sep = "_")) := ifelse(!!as.name(time_var) == years[i], !!as.name(variable), NA)) %>%
      fill(!!as.name(paste(variable, years[i], sep = "_")), .direction = "downup") %>%
      dplyr::ungroup()
  }
  
  return(data)
}
