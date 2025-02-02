##########################################
# GSE71014
##########################################

e1 <- new.env()
load('data/RData/GSE71014.RData', envir = e1)

e1$times <- as.double(gsub('.*:', '', e1$GSE71014.p$characteristics_ch1.2))
e1$event <- event_by_contains(e1$GSE71014.p$characteristics_ch1.3, "event (1: dead, 0:alive): 0", 
                              "event (1: dead, 0:alive): 1", fixed = TRUE)


GSE71014 <- list(X = e1$GSE71014.expr,
                  Y = data.frame(time = e1$times,
                                 death = e1$event),
                  NAME = 'GSE71014')

saveRDS(GSE71014, file = 'LEUKEMIA_LIB/data/GSE71014.rds')
