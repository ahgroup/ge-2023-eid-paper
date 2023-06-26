mycoef <- function(x){
  # browser()
  
  if ("dose_log" %in% rownames(fixef(x))) {
    op <- fixef(x)[rownames(fixef(x)) == "dose_log"] %>% t() %>% as.data.frame()
    colnames(op)= colnames(fixef(x))
    op$variable = rownames(fixef(x))[rownames(fixef(x)) == "dose_log"]
    op$op = paste0(round(op$Estimate,2), " (95%CI, ", round(op$Q2.5,2), " to ", round(op$Q97.5,2), ")")
  }
  
  if ("dose_catMedium" %in% rownames(fixef(x))) {
    op <- fixef(x) %>% as.data.frame()
    colnames(op) = colnames(fixef(x))
    op$variable = rownames(fixef(x))
    op$op = paste0(round(op$Estimate,2), " (95%CI, ", round(op$Q2.5,2), " to ", round(op$Q97.5,2), ")")
  }
  
  if ("growth_rate" %in% rownames(fixef(x))) {
    op <- fixef(x)[rownames(fixef(x)) == "growth_rate"] %>% t() %>% as.data.frame()
    colnames(op) = colnames(fixef(x))
    op$variable = rownames(fixef(x))[rownames(fixef(x)) == "growth_rate"]
    op$op = paste0(round(op$Estimate,2), " (95%CI, ", round(op$Q2.5,2), " to ", round(op$Q97.5,2), ")")
  }
  
  if ("outcome_log" %in% rownames(fixef(x))) {
    op <- fixef(x)[rownames(fixef(x)) == "outcome_log"] %>% t() %>% as.data.frame()
    colnames(op) = colnames(fixef(x))
    op$variable = rownames(fixef(x))[rownames(fixef(x)) == "outcome_log"]
    op$op = paste0(round(op$Estimate,2), " (95%CI, ", round(op$Q2.5,2), " to ", round(op$Q97.5,2), ")")
  }
  
  return(op)
}