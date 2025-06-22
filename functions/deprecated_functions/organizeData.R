library(dplyr)
library(survival)
library(GGally)
library(ggplot2)
library(readxl)
library(ggsignif)
library(tidyverse)
library(MASS)
library(mongolite)

connect_mongo1 <- function(url = "mongodb://localhost:27017/") {
  url2 <- url
  function() {
    m <- mongo(url2)
  }
}


# if your mongoDB port varies, add connection URL to connect_mongo1 param
connect_mongo <- connect_mongo1()

get_TCGA_data <- function() {
  # read .RDS into scope
  TCGA <- readRDS("C:/Users/nateg/Desktop/R Files/AML-BET-main/data/TCGA.rds")
  # remove miscellaneous
  TCGA$X <- TCGA$X[-(1:8),]
  
  for (x in 1:length(rownames(TCGA$X))) {
    query_name <- rownames(TCGA$X)[x]
    name_to_replace <- sub("\\|.*", "", rownames(TCGA$X)[x])
    if (name_to_replace %in% rownames(TCGA$X)) {
      avg_original <- rowMeans(TCGA$X[name_to_replace, ,drop = FALSE])
      avg_new <- rowMeans(TCGA$X[query_name, ,drop = FALSE])
      if (avg_original > avg_new) {
        name_to_replace <- query_name
      } else {rownames(TCGA$X)[rownames(TCGA$X) == name_to_replace] <- query_name}
    }
    rownames(TCGA$X)[x] <- name_to_replace
  }
  
  connect_mongo()
  TCGA_expr <- mongo("TCGA_expr")
  
  for (x in 1:length(rownames(TCGA$X))) {
    if ("|" %in% rownames(TCGA$X)[x]) {next}
    expr <- as.vector(TCGA$X[x,])
    expr <- paste0(expr, collapse = ",")
    
    str <- paste0('{"gene" :"', rownames(TCGA$X)[x],'", "expr" : [', expr, ']}')
    TCGA_expr$insert(str)
  }
  
}

get_TCGA_clinical <- function() {
  TCGA <- readRDS("C:/Users/nateg/Desktop/R Files/AML-BET-main/data/TCGA.rds")
  rownames(TCGA$Y) <- TCGA$p[,1]
  
  TCGA$Y$risk[TCGA$Y$risk == "intermediate/normal"] <- "intermediate"
  
  connect_mongo()
  TCGA_clinical <- mongo("TCGA_clinical")
  
  for (x in 1:length(rownames(TCGA$Y))) {
    str <- paste0('{ "id" : "', rownames(TCGA$Y)[x], '",')
    if (!is.na(TCGA$Y[x,1])) {
      str <- paste0(str, '"risk" : "', TCGA$Y[x,1], '",')
    }
    if (!is.na(TCGA$Y[x,2])) {
      str <- paste0(str, '"time" : ', TCGA$Y[x,2], ',')
    }
    if (!is.na(TCGA$Y[x,3])) {
      str <- paste0(str, '"death" : ', TCGA$Y[x,3], ',')
    }
    if (!is.na(TCGA$Y[x,4])) {
      str <- paste0(str, '"gender" : "', TCGA$Y[x,4], '"')
    }
    str <- paste0(str, ' }')
    
    TCGA_clinical$insert(str)
  }
}

get_BEAT_BMA_data <- function() {
  BEAT_BMA <- readRDS("C:/Users/nateg/Desktop/R Files/AML-BET-main/data/BEAT_BMA.rds")
  
  connect_mongo()
  BEAT_BMA_expr <- mongo("BEAT_BMA_expr")
  
  for (x in 1:length(rownames(BEAT_BMA$X))) {
    expr <- as.vector(BEAT_BMA$X[x,])
    expr <- paste0(expr, collapse = ",")
    
    str <- paste0('{"gene" :"', rownames(BEAT_BMA$X)[x],'", "expr" : [', expr, ']}')
    BEAT_BMA_expr$insert(str)
  }
}

get_BEAT_BMA_clinical <- function() {
  BEAT_BMA <- readRDS("C:/Users/nateg/Desktop/R Files/AML-BET-main/data/BEAT_BMA.rds")
  rownames(BEAT_BMA$Y) <- rownames(BEAT_BMA$p)
  
  BEAT_BMA$Y$risk[BEAT_BMA$Y$risk == "Intermediate"] <- "intermediate"
  BEAT_BMA$Y$risk[BEAT_BMA$Y$risk == "Favorable"] <- "favorable"
  BEAT_BMA$Y$risk[BEAT_BMA$Y$risk == "Adverse"] <- "poor"
  
  connect_mongo()
  BEAT_BMA_clinical <- mongo("BEAT_BMA_clinical")
  
  for (x in 1:length(rownames(BEAT_BMA$Y))) {
    str <- paste0('{ "id" : "', rownames(BEAT_BMA$Y)[x], '",')
    if (!is.na(BEAT_BMA$Y[x,1])) {
      str <- paste0(str, '"risk" : "', BEAT_BMA$Y[x,1], '",')
    }
    if (!is.na(BEAT_BMA$Y[x,2])) {
      str <- paste0(str, '"time" : ', BEAT_BMA$Y[x,2], ',')
    }
    if (!is.na(BEAT_BMA$Y[x,3])) {
      str <- paste0(str, '"death" : ', BEAT_BMA$Y[x,3])
    } else {str <- substr(str,1,nchar(str)-1)}
    str <- paste0(str, ' }')
    
    BEAT_BMA_clinical$insert(str)
  }
}

get_BEAT_PB_data <- function() {
  BEAT_PB <- readRDS("C:/Users/nateg/Desktop/R Files/AML-BET-main/data/BEAT_PB.rds")
  
  connect_mongo()
  BEAT_PB_expr <- mongo("BEAT_PB_expr")
  
  for (x in 1:length(rownames(BEAT_PB$X))) {
    expr <- as.vector(BEAT_PB$X[x,])
    expr <- paste0(expr, collapse = ",")
    
    str <- paste0('{"gene" :"', rownames(BEAT_PB$X)[x],'", "expr" : [', expr, ']}')
    BEAT_PB_expr$insert(str)
  }
}

get_BEAT_PB_clinical <- function() {
  BEAT_PB <- readRDS("C:/Users/nateg/Desktop/R Files/AML-BET-main/data/BEAT_PB.rds")
  rownames(BEAT_PB$Y) <- rownames(BEAT_PB$p)
  
  BEAT_PB$Y$risk[BEAT_PB$Y$risk == "Intermediate"] <- "intermediate"
  BEAT_PB$Y$risk[BEAT_PB$Y$risk == "Favorable"] <- "favorable"
  BEAT_PB$Y$risk[BEAT_PB$Y$risk == "Adverse"] <- "poor"
  
  connect_mongo()
  BEAT_PB_clinical <- mongo("BEAT_PB_clinical")
  
  for (x in 1:length(rownames(BEAT_PB$Y))) {
    str <- paste0('{ "id" : "', rownames(BEAT_PB$Y)[x], '",')
    if (!is.na(BEAT_PB$Y[x,1])) {
      str <- paste0(str, '"risk" : "', BEAT_PB$Y[x,1], '",')
    }
    if (!is.na(BEAT_PB$Y[x,2])) {
      str <- paste0(str, '"time" : ', BEAT_PB$Y[x,2], ',')
    }
    if (!is.na(BEAT_PB$Y[x,3])) {
      str <- paste0(str, '"death" : ', BEAT_PB$Y[x,3])
    } else {str <- substr(str,1,nchar(str)-1)}
    str <- paste0(str, ' }')
    
    BEAT_PB_clinical$insert(str)
  }
}

get_target_data <- function() {
  ensemble <- readRDS("C:/Users/nateg/Desktop/R Files/AML-BET-main/data/ensemble.rds")
  rownames(ensemble) <- ensemble[,1]
  target <- readRDS("C:/Users/nateg/Desktop/R Files/AML-BET-main/data/target.rds")
  
  connect_mongo()
  target_expr <- mongo("target_expr")
  
  for (x in 1:length(rownames(target$X))) {
    name <- ensemble[rownames(target$X)[x],2]
    
    if (is.na(name)) {next}
    expr <- as.vector(target$X[x,])
    expr <- paste0(expr, collapse = ",")
    
    str <- paste0('{"gene" :"', name,'", "expr" : [', expr, ']}')
    target_expr$insert(str)
  }
}

get_target_clinical <- function() {
  target <- readRDS("C:/Users/nateg/Desktop/R Files/AML-BET-main/data/target.rds")
  rownames(target$Y) <- target$p[,1]
  
  target$Y$risk[target$Y$risk == "Standard"] <- "intermediate"
  target$Y$risk[target$Y$risk == "Low"] <- "favorable"
  target$Y$risk[target$Y$risk == "High"] <- "poor"
  
  connect_mongo()
  target_clinical <- mongo("target_clinical")
  
  for (x in 1:length(rownames(target$Y))) {
    str <- paste0('{ "id" : "', rownames(target$Y)[x], '",')
    if (!is.na(target$Y[x,1])) {
      str <- paste0(str, '"risk" : "', target$Y[x,1], '",')
    }
    if (!is.na(target$Y[x,2])) {
      str <- paste0(str, '"time" : ', target$Y[x,2], ',')
    }
    if (!is.na(target$Y[x,3])) {
      str <- paste0(str, '"death" : ', target$Y[x,3], ',')
    }
    if (!is.na(target$Y[x,4])) {
      str <- paste0(str, '"gender" : "', target$Y[x,4], '"')
    }
    str <- paste0(str, ' }')
    
    target_clinical$insert(str)
  }
}

get_GSE71014_data <- function() {
  load("C:/Users/nateg/Downloads/GitHub/AML-BET/data/RData/GPL10558.RData")
  rownames(GPL10558) <- GPL10558[,1]
  GSE71014 <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/GSE71014.rds")
  GSE71014$X <- na.omit(GSE71014$X)
  
  connect_mongo()
  GSE71014_expr <- mongo("GSE71014_expr")
  
  for (x in 1:length(rownames(GSE71014$X))) {
    name <- GPL10558[rownames(GSE71014$X)[x],2]
    
    if (is.na(name)) {next}
    expr <- as.vector(GSE71014$X[x,])
    expr <- paste0(expr, collapse = ",")
    
    str <- paste0('{"gene" :"', name,'", "expr" : [', expr, ']}')
    GSE71014_expr$insert(str)
  }
}

get_GSE71014_clinical <- function() {
  GSE71014 <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/GSE71014.rds")
  
  connect_mongo()
  GSE71014_clinical <- mongo("GSE71014_clinical")
  
  for (x in 1:length(rownames(GSE71014$Y))) {
    str <- paste0('{ "time" : ', GSE71014$Y[x,1], ', "death" : ', GSE71014$Y[x,2], '}')
    GSE71014_clinical$insert(str)
  }
}

get_GSE6891_1_data <- function() {
  GSE6891_1 <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/GSE6891_1.rds")
  GSE6891_1$X <- data.frame(GSE6891_1$X)
  load("C:/Users/nateg/Downloads/GitHub/AML-BET/data/RData/GPL570.RData")
  
  rownames(GPL570) <- GPL570[,1]
  
  # GSE6891_1_data <- data.frame (
  #   gene <- "builder"
  # )
  # 
  #colnames(GSE6891_1_data)[1] <- "gene"
  
  # for (x in 1:length(colnames(GSE6891_1$X))) {
  #   GSE6891_1_data$name <- 0
  #   colnames(GSE6891_1_data)[x+1] <- colnames(GSE6891_1$X)[x]
  # }

  
  m <- matrix(0,ncol = ncol(GSE6891_1$X), nrow = 1)
  colnames(m) <- colnames(GSE6891_1$X)
  m <- cbind(gene = 0, m)
  GSE6891_1_data <- data.frame(m)


  for (x in 1:length(rownames(GSE6891_1$X))) {
    gene <- GPL570[rownames(GSE6891_1$X)[x],2]
    gene <- unlist(strsplit(genes, ' /// '))
    
    # for (y in 1:length(genes)) {
    #   GSE6891_1_data[nrow(GSE6891_1_data) + 1,] = c(genes[y], GSE6891_1$X[x,])
    # }
    
    rows_to_add <- cbind(gene, GSE6891_1$X[x,] )
    GSE6891_1_data <- rbind(GSE6891_1_data, rows_to_add)
    
    
  }
  
  GSE6891_1_data[1,1] <- NA
  
  known_genes <- unique(GSE6891_1_data[2:length(rownames(GSE6891_1_data)),1])
  
  known_genes <- na.omit(known_genes)
  
  final_data <- data.frame (
    gene <- NA
  )
  
  colnames(final_data)[1] <- "gene"
  
  for (x in 1:length(colnames(GSE6891_1$X))) {
    final_data$name <- 0
    colnames(final_data)[x+1] <- colnames(GSE6891_1$X)[x]
  }
  
  for (x in 1:length(known_genes)) {
    curr_gene <- filter(GSE6891_1_data, gene == known_genes[x])
    mean_expr <- rowMeans(curr_gene[2:length(colnames(curr_gene))])
    expr_to_keep <- which(mean_expr == max(mean_expr))
    final_data[nrow(final_data) + 1,] = curr_gene[expr_to_keep,]
  }
  
  # 1. get a vector of all possible unique genes
  
  # 2. start with an empty expression matrix with nrow = number of unique genes and ncol the number of samples
  
  # 3. For each gene :
        # get the corresponding expression values (possibly the highest mean expression) and assign to
        # expression matrix
  
      # (note: the R way to do this is to create a function that takes gene name and returns the expression value)
  
  connect_mongo()
  GSE6891_1_expr <- mongo("GSE6891_1_expr")
  
  for (x in 2:length(rownames(final_data))) {
    name <- final_data[x,1]
    
    expr <- as.vector(final_data[x,2:length(colnames(final_data))])
    expr <- paste0(expr, collapse = ",")
    
    str <- paste0('{"gene" :"', name,'", "expr" : [', expr, ']}')
    GSE6891_1_expr$insert(str)
  }
}

get_GSE6891_1_clinical <- function() {
  GSE6891_1 <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/GSE6891_1.rds")
  rownames(GSE6891_1$Y) <- rownames(GSE6891_1$p)
  GSE_6891_survival <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/GSE_6891_survival.rds")
  rownames(GSE_6891_survival) <- GSE_6891_survival$volgnummer
  
  GSE6891_1$Y$risk[GSE6891_1$Y$risk == "Intermediate"] <- "intermediate"
  GSE6891_1$Y$risk[GSE6891_1$Y$risk == "Favorable"] <- "favorable"
  GSE6891_1$Y$risk[GSE6891_1$Y$risk == "Adverse"] <- "poor"
  
  GSE6891_1$Y$os <- NA
  GSE6891_1$Y$osi <- NA
  
  for (x in 1:length(rownames(GSE6891_1$Y))) {
    title <- GSE6891_1$p[rownames(GSE6891_1$Y)[x],1]
    title <- substr(title, 5, nchar(title))
    GSE6891_1$Y[x,3] <- GSE_6891_survival[title, 4]
    GSE6891_1$Y[x,4] <- GSE_6891_survival[title, 3]
  }
  
  GSE6891_1$Y$osi[GSE6891_1$Y$osi == "alive"] <- 0
  GSE6891_1$Y$osi[GSE6891_1$Y$osi == "dead"] <- 1
  
  connect_mongo()
  GSE6891_1_clinical <- mongo("GSE6891_1_clinical")
  
  for (x in 1:length(rownames(GSE6891_1$Y))) {
    str <- paste0('{ "id" : "', rownames(GSE6891_1$Y)[x], '",')
    if (!is.na(GSE6891_1$Y[x,1])) {
      str <- paste0(str, '"risk" : "', GSE6891_1$Y[x,1], '",')
    }
    if (!is.na(GSE6891_1$Y[x,4])) {
      str <- paste0(str, '"death" : ', GSE6891_1$Y[x,4], ',')
    }
    if (!is.na(GSE6891_1$Y[x,3])) {
      str <- paste0(str, '"time" : ', GSE6891_1$Y[x,3], ',')
    }
    if (!is.na(GSE6891_1$Y[x,2])) {
      str <- paste0(str, '"gender" : "', GSE6891_1$Y[x,2], '"')
    }
    str <- paste0(str, ' }')
    
    GSE6891_1_clinical$insert(str)
  }
}

get_GSE6891_2_data <- function() {
  GSE6891_2 <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/GSE6891_2.rds")
  GSE6891_2$X <- data.frame(GSE6891_2$X)
  load("C:/Users/nateg/Downloads/GitHub/AML-BET/data/RData/GPL570.RData")
  
  rownames(GPL570) <- GPL570[,1]
  
  GSE6891_2_data <- data.frame (
    gene <- "builder"
  )
  
  colnames(GSE6891_2_data)[1] <- "gene"
  
  for (x in 1:length(colnames(GSE6891_2$X))) {
    GSE6891_2_data$name <- 0
    colnames(GSE6891_2_data)[x+1] <- colnames(GSE6891_2$X)[x]
  }
  
  for (x in 1:length(rownames(GSE6891_2$X))) {
    genes <- GPL570[rownames(GSE6891_2$X)[x],2]
    genes <- unlist(strsplit(genes, ' /// '))
    
    for (y in 1:length(genes)) {
      GSE6891_2_data[nrow(GSE6891_2_data) + 1,] = c(genes[y], GSE6891_2$X[x,])
    }
  }
  
  GSE6891_2_data[1,1] <- NA
  
  known_genes <- unique(GSE6891_2_data[2:length(rownames(GSE6891_2_data)),1])
  
  known_genes <- na.omit(known_genes)
  
  final_data <- data.frame (
    gene <- NA
  )
  
  colnames(final_data)[1] <- "gene"
  
  for (x in 1:length(colnames(GSE6891_2$X))) {
    final_data$name <- 0
    colnames(final_data)[x+1] <- colnames(GSE6891_2$X)[x]
  }
  
  for (x in 1:length(known_genes)) {
    curr_gene <- filter(GSE6891_2_data, gene == known_genes[x])
    mean_expr <- rowMeans(curr_gene[2:length(colnames(curr_gene))])
    expr_to_keep <- which(mean_expr == max(mean_expr))
    final_data[nrow(final_data) + 1,] = curr_gene[expr_to_keep,]
  }
  
  connect_mongo()
  GSE6891_2_expr <- mongo("GSE6891_2_expr")
  
  for (x in 2:length(rownames(final_data))) {
    name <- final_data[x,1]
    
    expr <- as.vector(final_data[x,2:length(colnames(final_data))])
    expr <- paste0(expr, collapse = ",")
    
    str <- paste0('{"gene" :"', name,'", "expr" : [', expr, ']}')
    GSE6891_2_expr$insert(str)
  }
}

get_GSE6891_2_clinical <- function() {
  GSE6891_2 <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/GSE6891_2.rds")
  rownames(GSE6891_2$Y) <- rownames(GSE6891_2$p)
  GSE_6891_survival <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/GSE_6891_survival.rds")
  rownames(GSE_6891_survival) <- GSE_6891_survival$volgnummer
  
  GSE6891_2$Y$risk[GSE6891_2$Y$risk == "Intermediate"] <- "intermediate"
  GSE6891_2$Y$risk[GSE6891_2$Y$risk == "Favorable"] <- "favorable"
  GSE6891_2$Y$risk[GSE6891_2$Y$risk == "Adverse"] <- "poor"
  
  GSE6891_2$Y$os <- NA
  GSE6891_2$Y$osi <- NA
  
  for (x in 1:length(rownames(GSE6891_2$Y))) {
    title <- GSE6891_2$p[rownames(GSE6891_2$Y)[x],1]
    title <- substr(title, 5, nchar(title))
    GSE6891_2$Y[x,3] <- GSE_6891_survival[title, 4]
    GSE6891_2$Y[x,4] <- GSE_6891_survival[title, 3]
    # Use gsub and match
  }
  
  GSE6891_2$Y$osi[GSE6891_2$Y$osi == "alive"] <- 0
  GSE6891_2$Y$osi[GSE6891_2$Y$osi == "dead"] <- 1
  
  connect_mongo()
  GSE6891_2_clinical <- mongo("GSE6891_2_clinical")
  
  for (x in 1:length(rownames(GSE6891_2$Y))) {
    str <- paste0('{ "id" : "', rownames(GSE6891_2$Y)[x], '",')
    if (!is.na(GSE6891_2$Y[x,1])) {
      str <- paste0(str, '"risk" : "', GSE6891_2$Y[x,1], '",')
    }
    if (!is.na(GSE6891_2$Y[x,4])) {
      str <- paste0(str, '"death" : ', GSE6891_2$Y[x,4], ',')
    }
    if (!is.na(GSE6891_2$Y[x,3])) {
      str <- paste0(str, '"time" : ', GSE6891_2$Y[x,3], ',')
    }
    if (!is.na(GSE6891_2$Y[x,2])) {
      str <- paste0(str, '"gender" : "', GSE6891_2$Y[x,2], '"')
    }
    str <- paste0(str, ' }')
    
    GSE6891_2_clinical$insert(str)
  }
}

get_GSE37642_1_data <- function() {
  GSE37642_1 <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/GSE37642_1.rds")
  GSE37642_1$X <- data.frame(GSE37642_1$X)
  load("C:/Users/nateg/Downloads/GitHub/AML-BET/data/RData/GPL570.RData")
  
  rownames(GPL570) <- GPL570[,1]
  
  GSE37642_1_data <- data.frame (
    gene <- "builder"
  )
  
  colnames(GSE37642_1_data)[1] <- "gene"
  
  for (x in 1:length(colnames(GSE37642_1$X))) {
    GSE37642_1_data$name <- 0
    colnames(GSE37642_1_data)[x+1] <- colnames(GSE37642_1$X)[x]
  }
  
  for (x in 1:length(rownames(GSE37642_1$X))) {
    genes <- GPL570[rownames(GSE37642_1$X)[x],2]
    genes <- unlist(strsplit(genes, ' /// '))
    
    for (y in 1:length(genes)) {
      GSE37642_1_data[nrow(GSE37642_1_data) + 1,] = c(genes[y], GSE37642_1$X[x,])
    }
  }
  
  GSE37642_1_data[1,1] <- NA
  
  known_genes <- unique(GSE37642_1_data[2:length(rownames(GSE37642_1_data)),1])
  
  known_genes <- na.omit(known_genes)
  
  final_data <- data.frame (
    gene <- NA
  )
  
  colnames(final_data)[1] <- "gene"
  
  for (x in 1:length(colnames(GSE37642_1$X))) {
    final_data$name <- 0
    colnames(final_data)[x+1] <- colnames(GSE37642_1$X)[x]
  }
  
  for (x in 1:length(known_genes)) {
    curr_gene <- filter(GSE37642_1_data, gene == known_genes[x])
    mean_expr <- rowMeans(curr_gene[2:length(colnames(curr_gene))])
    expr_to_keep <- which(mean_expr == max(mean_expr))
    final_data[nrow(final_data) + 1,] = curr_gene[expr_to_keep,]
  }
  
  connect_mongo()
  GSE37642_1_expr <- mongo("GSE37642_1_expr")
  
  for (x in 2:length(rownames(final_data))) {
    name <- final_data[x,1]
    
    expr <- as.vector(final_data[x,2:length(colnames(final_data))])
    expr <- paste0(expr, collapse = ",")
    
    str <- paste0('{"gene" :"', name,'", "expr" : [', expr, ']}')
    GSE37642_1_expr$insert(str)
  }
}

get_GSE37642_1_clinical <- function() {
  GSE37642_1 <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/GSE37642_1.rds")
  rownames(GSE37642_1$Y) <- rownames(GSE37642_1$p)
  
  GSE37642_1$Y$risk[GSE37642_1$Y$risk == "NA"] <- NA
  GSE37642_1$Y$risk[GSE37642_1$Y$risk == "adverse"] <- "poor"
  
  connect_mongo()
  GSE37642_1_clinical <- mongo("GSE37642_1_clinical")
  
  for (x in 1:length(rownames(GSE37642_1$Y))) {
    str <- paste0('{ "id" : "', rownames(GSE37642_1$Y)[x], '",')
    if (!is.na(GSE37642_1$Y[x,1])) {
      str <- paste0(str, '"risk" : "', GSE37642_1$Y[x,1], '",')
    }
    if (!is.na(GSE37642_1$Y[x,2])) {
      str <- paste0(str, '"time" : ', GSE37642_1$Y[x,2], ',')
    }
    if (!is.na(GSE37642_1$Y[x,3])) {
      str <- paste0(str, '"death" : ', GSE37642_1$Y[x,3])
    } else {str <- substr(str,1,nchar(str)-1)}
    str <- paste0(str, ' }')
    
    GSE37642_1_clinical$insert(str)
  }
}

get_GSE37642_2_data <- function() {
  GSE37642_2 <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/GSE37642_2.rds")
  GSE37642_2$X <- data.frame(GSE37642_2$X)
  load("C:/Users/nateg/Downloads/GitHub/AML-BET/data/RData/GPL96.RData")
  
  rownames(GPL96) <- GPL96[,1]
  
  GSE37642_2_data <- data.frame (
    gene <- "builder"
  )
  
  colnames(GSE37642_2_data)[1] <- "gene"
  
  for (x in 1:length(colnames(GSE37642_2$X))) {
    GSE37642_2_data$name <- 0
    colnames(GSE37642_2_data)[x+1] <- colnames(GSE37642_2$X)[x]
  }
  
  for (x in 1:length(rownames(GSE37642_2$X))) {
    genes <- GPL96[rownames(GSE37642_2$X)[x],2]
    genes <- unlist(strsplit(genes, ' /// '))
    
    for (y in 1:length(genes)) {
      GSE37642_2_data[nrow(GSE37642_2_data) + 1,] = c(genes[y], GSE37642_2$X[x,])
    }
  }
  
  GSE37642_2_data[1,1] <- NA
  
  known_genes <- unique(GSE37642_2_data[2:length(rownames(GSE37642_2_data)),1])
  
  known_genes <- na.omit(known_genes)
  
  final_data <- data.frame (
    gene <- NA
  )
  
  colnames(final_data)[1] <- "gene"
  
  for (x in 1:length(colnames(GSE37642_2$X))) {
    final_data$name <- 0
    colnames(final_data)[x+1] <- colnames(GSE37642_2$X)[x]
  }
  
  for (x in 1:length(known_genes)) {
    curr_gene <- filter(GSE37642_2_data, gene == known_genes[x])
    mean_expr <- rowMeans(curr_gene[2:length(colnames(curr_gene))])
    expr_to_keep <- which(mean_expr == max(mean_expr))
    final_data[nrow(final_data) + 1,] = curr_gene[expr_to_keep,]
  }
  
  connect_mongo()
  GSE37642_2_expr <- mongo("GSE37642_2_expr")
  
  for (x in 2:length(rownames(final_data))) {
    name <- final_data[x,1]
    
    expr <- as.vector(final_data[x,2:length(colnames(final_data))])
    expr <- paste0(expr, collapse = ",")
    
    str <- paste0('{"gene" :"', name,'", "expr" : [', expr, ']}')
    GSE37642_2_expr$insert(str)
  }
}

get_GSE37642_2_clinical <- function() {
  GSE37642_2 <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/GSE37642_2.rds")
  rownames(GSE37642_2$Y) <- rownames(GSE37642_2$p)
  
  GSE37642_2$Y$risk[GSE37642_2$Y$risk == "NA"] <- NA
  GSE37642_2$Y$risk[GSE37642_2$Y$risk == "adverse"] <- "poor"
  
  connect_mongo()
  GSE37642_2_clinical <- mongo("GSE37642_2_clinical")
  
  for (x in 1:length(rownames(GSE37642_2$Y))) {
    str <- paste0('{ "id" : "', rownames(GSE37642_2$Y)[x], '",')
    if (!is.na(GSE37642_2$Y[x,1])) {
      str <- paste0(str, '"risk" : "', GSE37642_2$Y[x,1], '",')
    }
    if (!is.na(GSE37642_2$Y[x,2])) {
      str <- paste0(str, '"time" : ', GSE37642_2$Y[x,2], ',')
    }
    if (!is.na(GSE37642_2$Y[x,3])) {
      str <- paste0(str, '"death" : ', GSE37642_2$Y[x,3])
    } else {str <- substr(str,1,nchar(str)-1)}
    str <- paste0(str, ' }')
    
    GSE37642_2_clinical$insert(str)
  }
}



stopifnot(all(
      rownames(GSE6891_1$p) == colnames(GSE6891_1$X)
))




