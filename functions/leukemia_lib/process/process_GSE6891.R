source('functions.R')

e1 <- new.env()

# cohort 1
load('data/RData/GSE6891_training.RData', envir = e1)
#e1$x1 <- get_gene(e1$GSE6891_training.expr, e1$GPL570, 'ALDH1A1', combine = FALSE)
e1$risk <- get_favorable("6891", colnames(e1$GSE6891_training.expr))

risk <- e1$risk$Favorability

stopifnot(all(rownames(e1$GSE6891_training.p) == colnames(e1$GSE6891_training.expr)))

GSE6891_1 <- list(X = e1$GSE6891_training.expr,
                Y = data.frame(risk = risk, gender = e1$GSE6891_training.p$`gender:ch1`),
                NAME = 'GSE6891 (Cohort #1)',
                p = e1$GSE6891_training.p)

e1 <- new.env()
# cohort 2
load('data/RData/GSE6891_validation.RData',envir = e1)
e1$favorable <- get_favorable("6891", colnames(e1$GSE6891_validation.expr))
risk <- e1$favorable$Favorability

stopifnot(all(rownames(e1$GSE6891_validation.p) == colnames(e1$GSE6891_validation.expr)))

GSE6891_2 <- list(X = e1$GSE6891_validation.expr,
                  Y = data.frame(risk = risk, gender = e1$GSE6891_validation.p$`gender:ch1`),
                  NAME = 'GSE6891 (Cohort #2)',
                  p = e1$GSE6891_validation.p)

saveRDS(GSE6891_1, file = 'LEUKEMIA_LIB/data/GSE6891_1.rds')
saveRDS(GSE6891_2, file = 'LEUKEMIA_LIB/data/GSE6891_2.rds')
