holder_brms_fit_op_con <- function(data = NA, brms_fit = NA, dose_log = NA){
  
  # browser()
  
  min_days <- min(data$days)
  madays <- max(data$days)
  my_days <- seq(min_days,madays,by = 1/24)
  
  group_data <- list()
  group_curves <- list()
  group_auc <- list()
  group_tpkv <- list()
  group_oneset <- list()
  group_duration <- list()
  
  # extract 5000 samples of population level
  # posterior samples for 2nd outcome summary
  mysamples = 5000
  
  for (i in 1:nrow(dose_log)) {
    # data 
    condition = dose_log$dose_cat[i]
    group_data[[paste0("dose",i)]] <- crossing(days = my_days, dose_log) %>% 
      filter(dose_cat == condition)
    
    # group fitted values and curves
    set.seed(123456789)
    
    # it has mysamples rows, and my_days columns, log(y)
    group_temp <- posterior_epred(object = brms_fit,
                                  ndraws = mysamples,
                                  newdata = group_data[[paste0("dose",i)]],
                                  re_formula = NA) 
    
    # it has my_days rows, and mysamples columns, y
    temp <- as_tibble( t( exp(group_temp) ))
    colnames(temp) <- paste0("nsamples_",seq(1,mysamples))
    
    # total virus shedding
    group_auc[[paste0("dose",i)]] <- apply(temp, 2, trapz, x = my_days)  %>% 
      posterior_summary(probs = c(0.025, 0.975)) %>% 
      as_tibble() %>% 
      mutate(dose = dose_log$dose[i],
             dose_cat = dose_log$dose_cat[i],
             dose_log = dose_log$dose_log[i])
    
    # time of peak shedding
    group_tpkv[[paste0("dose",i)]] <- apply(temp, 2, function(x) my_days[x == max(x)]) %>% 
      posterior_summary(probs = c(0.025, 0.975)) %>% 
      as_tibble() %>% 
      mutate(dose = dose_log$dose[i],
             dose_cat = dose_log$dose_cat[i],
             dose_log = dose_log$dose_log[i])
    
    pcr_times <- apply(temp, 2, function(x) if_else(x >= 225000, TRUE, FALSE))
    
    # oneset time of POS
    group_oneset[[paste0("dose",i)]] <- apply(pcr_times, 2, function(x) min(my_days[x])) %>% 
      posterior_summary(probs = c(0.025, 0.975)) %>% 
      as_tibble() %>% 
      mutate(dose = dose_log$dose[i],
             dose_cat = dose_log$dose_cat[i],
             dose_log = dose_log$dose_log[i])
    
    # duration of POS
    group_duration[[paste0("dose",i)]] <- apply(pcr_times, 2, function(x) max(my_days[x]) - min(my_days[x])) %>% 
      posterior_summary(probs = c(0.025, 0.975)) %>% 
      as_tibble() %>% 
      mutate(dose = dose_log$dose[i],
             dose_cat = dose_log$dose_cat[i],
             dose_log = dose_log$dose_log[i])
    
    gc()
  }
  
  # row bind predictions
  # ===============================================
  ts_onset <- bind_rows(group_oneset)
  ts_duration <- bind_rows(group_duration)
  ts_tpkv <- bind_rows(group_tpkv)
  ts_auc <- bind_rows(group_auc)
  
  return(list(ts_onset = ts_onset, ts_duration = ts_duration, ts_tpkv = ts_tpkv, ts_auc = ts_auc))
  
}