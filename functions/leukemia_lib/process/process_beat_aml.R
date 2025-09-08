##########################################
# BEAT AML
##########################################
source('functions.R')

library(AMLbeatR)
e1 <- new.env()
load('data/RData/beatAML.RData', envir = e1)
rm(av, rpkm, v, envir = e1)

e1$d=tidyClin(e1$clin) #672 rows in clin, one for each measurement; 562 rows in d, one for each patient
e1$beat_aml_expr <- e1$cpm[,-1:-2]
rownames(e1$beat_aml_expr) <- e1$cpm$Symbol 
e1$beat_aml_expr <- log2(e1$beat_aml_expr + 1)

e1$fav_all <- get_favorable('BEAT', e1$clin$LabId)
e1$clin <- cbind(e1$fav_all, e1$clin)

e1$myclin_bma <- e1$clin %>% dplyr::filter(specimenGroups == "Initial Acute Leukemia Diagnosis",
                                           specimenType == "Bone Marrow Aspirate")

e1$myclin_pb <- e1$clin %>% dplyr::filter(specimenGroups == "Initial Acute Leukemia Diagnosis",
                                          specimenType == "Peripheral Blood")


## any overlap??
ids <- c(e1$myclin_bma$PatientId, e1$myclin_pb$PatientId) %>% unique()
#vc <- vennCounts(cbind(bma = ids%in% e1$myclin_bma$PatientId, pb = ids%in% e1$myclin_pb$PatientId))
#vennDiagram(vc)

# remove duplicates
cat("remove duplicates from BEAT BMA/PB")
dups <- intersect(e1$myclin_bma$PatientId, e1$myclin_pb$PatientId)
e1$myclin_bma <- e1$myclin_bma %>% dplyr::filter(!PatientId %in% dups)
e1$myclin_pb <- e1$myclin_pb %>% dplyr::filter(!PatientId %in% dups)

e1$x <- e1$beat_aml_expr

get_common <- function(x, f) {
  common <- intersect(f$ID, colnames(x))
  m1 <- match(common, f$ID)
  m2 <- match(common, colnames(x))
  f <- f[m1,]
  x <- x[,m2]
  list(x = x, f = f)
}

e1$CC_bma <- get_common(e1$x, e1$myclin_bma)
e1$CC_pb <- get_common(e1$x, e1$myclin_pb)

risk_bma <- e1$CC_bma$f$Favorability
risk_pb <- e1$CC_pb$f$Favorability


e1$pb_status <- event_by_contains(e1$CC_pb$f$vitalStatus, "Alive", "Dead")
e1$bma_status <- event_by_contains(e1$CC_bma$f$vitalStatus, "Alive", "Dead")

e1$CC_bma$f$overallSurvival <- e1$CC_bma$f$overallSurvival / 365*12
e1$CC_pb$f$overallSurvival <- e1$CC_pb$f$overallSurvival / 365*12

clean_risk <- function(x) {
  x[x == 'FavorableOrIntermediate'] <- NA
  x[x == 'IntermediateOrAdverse'] <- NA
  x
}

risk_bma <- clean_risk(risk_bma)
risk_pb <- clean_risk(risk_pb)

BMA <- list(X = e1$CC_bma$x, 
            Y = data.frame(risk = risk_bma,
                           time = e1$CC_bma$f$overallSurvival,
                           death = e1$bma_status),
            NAME = 'BEAT AML (BMA)', p = e1$CC_bma$f
)
        
PB <- list(X = e1$CC_pb$x, 
            Y = data.frame(risk = risk_pb,
                           time = e1$CC_pb$f$overallSurvival,
                           death = e1$pb_status),
           NAME = 'BEAT AML (PB)', p = e1$CC_pb$f
)

saveRDS(BMA, file = 'LEUKEMIA_LIB/data/BEAT_BMA.rds')
saveRDS(PB, file = 'LEUKEMIA_LIB/data/BEAT_PB.rds')

