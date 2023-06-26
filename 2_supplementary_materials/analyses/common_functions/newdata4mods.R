newdata_cat <- data.frame(dose_cat = c("Low_L","Low","Medium","High"))

newdata_cnt19 <- data.frame(dose = c(4.8,12,48,75,300,1200,4800),
                            dose_cat = as.character(c(4.8,12,48,75,300,1200,4800)),
                            dose_log = log(c(4.8,12,48,75,300,1200,4800)) - log(48))

newdata_cnt20 <- data.frame(dose = c(0.48,1.2,4.8,12,48,75,300,1200,4800),
                            dose_cat = as.character(c(0.48,1.2,4.8,12,48,75,300,1200,4800)),
                            dose_log = log(c(0.48,1.2,4.8,12,48,75,300,1200,4800)) - log(48))
