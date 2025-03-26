library(dplyr)
library(survival)
library(GGally)
library(ggplot2)
library(readxl)
library(ggsignif)
library(tidyverse)
library(MASS)
library(mongolite)

# GPL96 is only for GSE37642_2
# GPL10558 is only for GSE71014
# ensemble is only for target

upload_GSE_expr_data <- function(GSE_dataset, GPL_dataset, target = FALSE) {
  GSE <- readRDS(paste0("C:/Users/nateg/Downloads/GitHub/AML-BET/data/", GSE_dataset, ".rds"))
  if (target == FALSE) {
    GPL <- load(paste0("C:/Users/nateg/Downloads/GitHub/AML-BET/data/RData/", GPL_dataset,".RData"))
    if (length(GPL) > 1) {stop("GPL data provided containts more than one object.")}
    GPL <- get(GPL)
  }
  else {
    GPL <- readRDS(paste0("C:/Users/nateg/Downloads/GitHub/AML-BET/data/", GPL_dataset, ".rds"))
  }
  
  # 1. vector of all possible unique genes
  all_genes <- unique(unlist(strsplit(GPL[,2], ' /// ')))
  all_genes <- na.omit(all_genes)
  
  # 2. Calling return_expr_GSE for every entry in all_genes
  expression_data <- t(sapply(all_genes, return_expr_GSE, dataset = GSE, GPL = GPL))
  
  # 3. Upload data to MongoDB
  upload_expr_mongo(GSE_dataset, expression_data)
}

# 3. Helper Function for 3
return_expr_GSE <- function(gene_name, dataset, GPL) {
  all_probes <- GPL[grepl(gene_name, GPL[,2], fixed = TRUE),1]
  probe_validity_check <- match(all_probes, rownames(dataset$X))
  if (all(is.na(probe_validity_check)) == TRUE) {
    return(rep(0, length(colnames(dataset$X))))
  }
  probe_validity_index <- which(!is.na(probe_validity_check))
  all_probes <- all_probes[probe_validity_index]
  mean_expr <- rowMeans(dataset$X[all_probes, , drop = FALSE])
  max_probe <- all_probes[which.max(mean_expr)]
  expr_data_to_return <- dataset$X[max_probe, , drop = FALSE]
}

upload_BEAT_expr_data <- function(dataset_name) {
  BEAT_dataset <- readRDS(paste0("C:/Users/nateg/Downloads/GitHub/AML-BET/data/", dataset_name, ".rds"))
  upload_expr_mongo(dataset_name, BEAT_dataset$X)
}

upload_TCGA_expr_data <- function() {
  TCGA <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/TCGA.rds")
  TCGA$X <- TCGA$X[-(1:8),]
  genes <- sub("\\|.*", "", rownames(TCGA$X))
  expression_data <- t(sapply(genes, return_expr_TCGA))
  upload_expr_mongo("TCGA", expression_data)
}

# Helper function for returning TCGA data
return_expr_TCGA <- function(gene_name) {
  all_entries <- grep(paste0(gene_name, "|"), rownames(TCGA$X), fixed = TRUE, value = TRUE)
  mean_expr <- rowMeans(TCGA$X[all_entries, , drop = FALSE])
  max_probe <- all_entries[which.max(mean_expr)]
  expr_data_to_return <- TCGA$X[max_probe, , drop = FALSE]
}

# Helper function for uploading final expression data
upload_expr_mongo <- function(dataset_name, expression_data) {
  connect_mongo <- connect_mongo1()
  connect_mongo()
  connection <- mongo(paste0(dataset_name, "_expr"))
  
  for (x in 1:length(rownames(expression_data))) {
    if (mean(expression_data[x,]) == 0) {next()}
    expr <- as.vector(expression_data[x,])
    expr <- paste0(expr, collapse = ",")
    
    str <- paste0('{"gene" :"', rownames(expression_data)[x],'", "expr" : [', expr, ']}')
    connection$insert(str)
  }
  print(paste0("Expression data from ", dataset_name, " successfully uploaded to MongoDB,"))
}

upload_clinical_data <- function(dataset_name) {
  dataset <- readRDS(paste0("C:/Users/nateg/Downloads/GitHub/AML-BET/data/", dataset_name, ".rds"))
  rownames(dataset$Y) <- colnames(dataset$X)
  dataset <- dataset$Y
  factors <- colnames(dataset)
  dataset <- standardize_clinical_data(dataset, factors)
  
  connect_mongo <- connect_mongo1()
  connect_mongo()
  connection <- mongo(paste0(dataset_name, "_clinical"))
  
  for (x in 1:length(rownames(dataset))) {
    json_text <- paste0('{ "id" : "', rownames(dataset)[x], '",')
    for (y in 1:length(factors)) {
      if(is.na(dataset[x,y])) {next()}
      if (is.numeric(dataset[x,y])) {
        json_text <- paste0(json_text, ' "', factors[y], '" : ', dataset[x,y], ',')
      }
      else {
        json_text <- paste0(json_text, ' "', factors[y], '" : "', dataset[x,y], '",')
      }
    }
    json_text <- substr(json_text, 1, nchar(json_text)-1)
    json_text <- paste0(json_text, ' }')
    connection$insert(json_text)
  }
}

standardize_clinical_data <- function(dataset, factors) {
  if ("risk" %in% factors) {
    dataset$risk[dataset$risk == "intermediate/normal"] <- "intermediate"
    dataset$risk[dataset$risk == "Intermediate"] <- "intermediate"
    dataset$risk[dataset$risk == "Favorable"] <- "favorable"
    dataset$risk[dataset$risk == "Adverse"] <- "poor"
    dataset$risk[dataset$risk == "Standard"] <- "intermediate"
    dataset$risk[dataset$risk == "Low"] <- "favorable"
    dataset$risk[dataset$risk == "High"] <- "poor"
    dataset$risk[dataset$risk == "NA"] <- NA
    dataset$risk[dataset$risk == "adverse"] <- "poor"
  }
  return(dataset)
}

# Connecting to MongoDB
connect_mongo1 <- function(url = "mongodb://localhost:27017/") {
  url2 <- url
  function() {
    m <- mongo(url2)
  }
}

return_data <- function(dataset, variable) {
  
}