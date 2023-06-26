holder_brms_fit_op_cat <- function(data = NA, brms_fit = NA, dose_cat = NA){
  
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
  
  # to predict 2nd outcome, so we used a new individual (e.g. ID0), who
  # theoretically similar to individuals in the trial but has not been enrolled
  # previously, be given to a certain amount of virus
  
  for (i in dose_cat) {
    # data 
    group_data[[paste0("dose",i)]] <- crossing(days = my_days, dose_cat = i) 
    
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
    
    # for each mysamples (columns), calculating 2nd outcomes
    # total virus shedding
    group_auc[[paste0("dose",i)]] <- apply(temp, 2, trapz, x = my_days)  %>% 
      posterior_summary(probs = c(0.025, 0.975)) %>% 
      as_tibble() %>% 
      mutate(dose = i)
    
    # time of peak shedding
    group_tpkv[[paste0("dose",i)]] <- apply(temp, 2, function(x) my_days[x == max(x)]) %>% 
      posterior_summary(probs = c(0.025, 0.975)) %>% 
      as_tibble() %>% 
      mutate(dose = i)
    
    pcr_times <- apply(temp, 2, function(x) if_else(x >= 225000, TRUE, FALSE))
    
    # oneset time of POS
    group_oneset[[paste0("dose",i)]] <- apply(pcr_times, 2, function(x) min(my_days[x])) %>% 
      posterior_summary(probs = c(0.025, 0.975)) %>% 
      as_tibble() %>% 
      mutate(dose = i)
    
    # duration of POS
    group_duration[[paste0("dose",i)]] <- apply(pcr_times, 2, function(x) max(my_days[x]) - min(my_days[x])) %>% 
      posterior_summary(probs = c(0.025, 0.975)) %>% 
      as_tibble() %>% 
      mutate(dose = i)
    
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
# end