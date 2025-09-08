source('functions.R')
load('data/RData/GSE37642_1.RData')
times <- as.double(gsub('.*:', '', GSE37642_1.p$characteristics_ch1.5))  / 365*12
death <- event_by_contains(GSE37642_1.p$characteristics_ch1.6, "alive", "dead")
X <- GSE37642_1.expr
risk <- get_favorable("37642_570", colnames(X))$ELN2017
risk[risk%in% 'MDS'] <- NA
GSE37642_1 <- list(X = X, Y = data.frame(risk = risk, time = times, death = death),
                   NAME = 'GSE37642 (GPL570)', p = GSE37642_1.p)
saveRDS(GSE37642_1, file = 'LEUKEMIA_LIB/data/GSE37642_1.rds')


rm2()
source('functions.R')
load('data/RData/GSE37642_2.RData')
times <- as.double(gsub('.*:', '', GSE37642_2.p$characteristics_ch1.5))  / 365*12
death <- event_by_contains(GSE37642_2.p$characteristics_ch1.6, "alive", "dead")
X <- GSE37642_2.expr
risk <- get_favorable("37642_96", colnames(X))$ELN2017
risk[risk%in% 'MDS'] <- NA
GSE37642_2 <- list(X = X, Y = data.frame(risk = risk, time = times, death = death),
                   NAME = 'GSE37642 (GPL96)',
                   p = GSE37642_2.p)
saveRDS(GSE37642_2, file = 'LEUKEMIA_LIB/data/GSE37642_2.rds')
