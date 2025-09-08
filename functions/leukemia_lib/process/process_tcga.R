
##############################################################################
# TCGA - read in and process raw data (download from http://firebrowse.org/)
# for LEUKEMIA_LIB
##############################################################################

library(stringr)
library(edgeR)
library(FirebrowseR)
library(readr)
source('functions.R')

aml_expr <- read_delim("data/LAML.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt",
                       "\t", escape_double = FALSE, trim_ws = TRUE)

g <- grep('raw_count', aml_expr[1,])
aml_expr <- aml_expr[,c(1,g)]
aml_expr <- aml_expr[-1,]

genes <- aml_expr$`Hybridization REF`
aml_expr <- aml_expr[,-1]
aml_expr <- apply(aml_expr, 2, as.double)

rownames(aml_expr) <- genes

## process data
dge <- DGEList(counts=aml_expr)
keep <- filterByExpr(aml_expr)
dge <- dge[keep,keep.lib.sizes = TRUE]
dge <- calcNormFactors(dge)

logCPM <- cpm(dge, log=TRUE, prior.count=1)

samples <- colnames(logCPM)
colnames(logCPM) <- str_match(samples, 'TCGA-[:alnum:]{2}-[:alnum:]{4}')

favorable <- get_favorable('TCGA', colnames(logCPM))

# get patient data
aml_p = Samples.Clinical(cohort = 'LAML', format="tsv", 
                         page_size = 500)

# select samples
m <- match(colnames(logCPM), aml_p$tcga_participant_barcode)
aml_p <- aml_p[m,]

# km curve
time <- aml_p$days_to_last_followup
death<- event_by_contains(aml_p$vital_status, 'alive', 'dead')
time[death == 1] <- aml_p$days_to_death[death==1]

time <- time / 365*12

X <- logCPM
Y <- data.frame(risk = aml_p$acute_myeloid_leukemia_calgb_cytogenetics_risk_category, time = time, death = death,
                gender = aml_p$gender)
NOTES <- list(death = 'os')

TCGA <- list(X = X, Y = Y, NOTES = NOTES, NAME = 'TCGA',
             p = aml_p)
saveRDS(TCGA, file = 'LEUKEMIA_LIB/data/TCGA.rds')


