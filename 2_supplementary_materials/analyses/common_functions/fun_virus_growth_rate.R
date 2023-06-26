fun_virus_growth_rate <- function(x){
  
  # because most peak after 1st day, so we choose 0.5th and 0.625th days as the two point of the
  # growth linear endpoints
  day0.625 <- x %>% 
    filter(days < 0.626, days > 0.624) 
  
  day0.5 <- x %>% 
    filter(days < 0.51, days > 0.49) 
  
  rate = (day0.625$Estimate - day0.5$Estimate) / (day0.625$days - day0.5$days)

  growth_rate = data.frame(ID = unique(x$ID), 
                           dose_cat  = unique(x$dose_cat ),
                           growth_rate = rate)

  return(growth_rate)
}