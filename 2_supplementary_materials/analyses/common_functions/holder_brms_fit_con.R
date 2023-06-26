holder_brms_fit_con <- function(data = NA, 
                                brms_fit = NA, 
                                lod1 = NA, 
                                lod2 = NA, 
                                mod = c("GoodnessOfFit","Posterior_grandmean")){
  
  # browser()
  
  if (mod == "GoodnessOfFit") {
    
    ids <- unique(data$ID)
    
    # do Parallel
    cl <- makeCluster(4)
    registerDoParallel(cl)
    
    temp <- foreach(i = ids) %dopar% {
      
      library(dplyr)
      library(tidyr)
      library(brms)
      
      temp_data <- data %>% filter(ID == i)
      op <- posterior_epred(object = brms_fit,
                            newdata = temp_data,
                            re_formula = NULL) %>% 
        brms::posterior_summary(probs = c(0.025, 0.975)) %>% 
        as_tibble() %>% 
        bind_cols(temp_data) %>% 
        mutate(pos_sday = if_else(Q2.5 > lod1, 1, 0), # when lower bound higher than lod1
               pos_eday = if_else(Q97.5 < lod1, 1, 0)) # when higher bound lower than lod1
      
      gc()
      
      return(op)
      
    }
  }
  
  if (mod == "Posterior_grandmean") {
    
    num_dose <- unique(data$dose_cat)
    
    cl <- makeCluster(4)
    registerDoParallel(cl)
    
    temp <- foreach(i = num_dose) %dopar% {
      
      library(dplyr)
      library(tidyr)
      library(brms)
      
      temp_data <- data %>% filter(dose_cat == i)
      op <- posterior_epred(object = brms_fit,
                            newdata = temp_data,
                            re_formula = NA) %>% 
        brms::posterior_summary(probs = c(0.025, 0.975)) %>% 
        as_tibble() %>% 
        bind_cols(temp_data) %>% 
        mutate(pos_sday = if_else(Q2.5 > lod1, 1, 0), # when lower bound higher than lod1
               pos_eday = if_else(Q97.5 < lod1, 1, 0)) # when higher bound lower than lod1
      
      gc()
      
      return(op)
      
    }
  }
  
  stopCluster(cl)
  return(temp)
}
# end