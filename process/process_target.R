source('functions.R')

load('data/RData/target_aml.RData')

stopifnot(all(colnames(target_logCPM) == target_aml_clin$PATIENT_ID))

X <- target_logCPM
risk <- target_aml_clin$RISK_GROUP
risk[risk == "Unknown"] <- NA

time <- target_aml_clin$OS_DAYS / 365*12
death <- event_by_contains(target_aml_clin$OS_STATUS, "0:LIVING", "1:DECEASED")

library(ensembldb)
library(EnsDb.Hsapiens.v86)
keys <- keys(EnsDb.Hsapiens.v86)
s <- ensembldb::select(EnsDb.Hsapiens.v86, keys = keys, 
            columns=c("GENEID","SYMBOL"),keytype="GENEID")

s <- dplyr::rename(s, ID = GENEID)

rownames(X) <- gsub('\\.\\d*$', '', rownames(X))

TARGET <- list(X = X, Y = data.frame(risk, time, death, gender = tolower(target_aml_clin$SEX)),
               NAME = 'TARGET', p = target_aml_clin)

saveRDS(TARGET, file = 'LEUKEMIA_LIB/data/target.rds')
saveRDS(s, file = 'LEUKEMIA_LIB/data/ensemble.rds')
