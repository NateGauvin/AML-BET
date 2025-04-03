library(dplyr)
library(ggplot2)
library(tidyverse)
library(MASS)
library(mongolite)
library(argparse)

expression_upload <- function(dataset, db) {
  if      (grepl("GSE", dataset)) {upload_GSE_expr_data(dataset, upload = db)}
  else if (dataset == "target") {upload_GSE_expr_data(dataset, GPL_dataset = "ensemble", target = TRUE, upload = db)}
  else if (grepl("BEAT", dataset)) {upload_BEAT_expr_data(dataset, upload = db)}
  else if (grepl("TCGA", dataset)) {upload_TCGA_expr_data(upload = db)}
  else {stop(paste0("unable to process dataset: ", dataset))}
  return(paste0("Successfully processed/uploaded ", dataset, " expression data."))
}



# GPL96 is only for GSE37642_2
# GPL10558 is only for GSE71014
# ensemble is only for target

upload_GSE_expr_data <- function(GSE_dataset, GPL_dataset = "auto", target = FALSE, upload = TRUE) {
  GSE <- readRDS(paste0("data/", GSE_dataset, ".rds"))
  if (GPL_dataset == "auto" && target == FALSE) {
    if (GSE_dataset == "GSE37642_2") {GPL_dataset <- "GPL96"}
    if (GSE_dataset == "GSE71014") {GPL_dataset <- "GPL10558"}
    else {GPL_dataset <- "GPL570"}
  }
  if (target == FALSE) {
    GPL <- load(paste0("data/RData/", GPL_dataset,".RData"))
    if (length(GPL) > 1) {stop("GPL data provided containts more than one object.")}
    GPL <- get(GPL)
  }
  else {
    GPL <- readRDS(paste0("data/additional_data/", GPL_dataset, ".rds"))
  }
  
  # 1. vector of all possible unique genes
  all_genes <- unique(unlist(strsplit(GPL[,2], ' /// ')))
  all_genes <- na.omit(all_genes)
  
  # 2. Calling return_expr_GSE for every entry in all_genes
  expression_data <- t(sapply(all_genes, return_expr_GSE, dataset = GSE, GPL = GPL))
  
  # 3. Upload data to MongoDB
  if (upload == TRUE) {upload_expr_mongo(GSE_dataset, expression_data)}
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

upload_BEAT_expr_data <- function(dataset_name, upload = TRUE) {
  BEAT_dataset <- readRDS(paste0("data/", dataset_name, ".rds"))
  if(upload == TRUE) {upload_expr_mongo(dataset_name, BEAT_dataset$X)}
}

upload_TCGA_expr_data <- function(upload = TRUE) {
  TCGA <- readRDS("data/TCGA.rds")
  TCGA$X <- TCGA$X[-(1:8),]
  genes <- sub("\\|.*", "", rownames(TCGA$X))
  expression_data <- t(sapply(genes, return_expr_TCGA, TCGA))
  if (upload == TRUE) {upload_expr_mongo("TCGA", expression_data)}
}

# Helper function for returning TCGA data
return_expr_TCGA <- function(gene_name, TCGA) {
  all_entries <- grep(paste0(gene_name, "|"), rownames(TCGA$X), fixed = TRUE, value = TRUE)
  mean_expr <- rowMeans(TCGA$X[all_entries, , drop = FALSE])
  max_probe <- all_entries[which.max(mean_expr)]
  expr_data_to_return <- TCGA$X[max_probe, , drop = FALSE]
}

# Helper function for uploading final expression data
upload_expr_mongo <- function(dataset_name, expression_data) {
  
  print("not yet connected")
  
  connection <- connect_mongo(collection_name = paste0(dataset_name, "_expr"))
  
  print("connection/collection made")
  
  expression_data <- data.frame(expression_data)
  
  for (x in 1:length(rownames(expression_data))) {
    if (rowMeans(expression_data[x,]) == 0) {next()}
    expr <- as.vector(expression_data[x,])
    expr <- paste0(expr, collapse = ",")
    
    str <- paste0('{"gene" :"', rownames(expression_data)[x],'", "expr" : [', expr, ']}')
    connection$insert(str)
  }
}

upload_clinical_data <- function(dataset_name, upload = TRUE) {
  dataset <- readRDS(paste0("data/", dataset_name, ".rds"))
  rownames(dataset$Y) <- colnames(dataset$X)
  if (grepl("GSE6891", dataset_name)) {dataset <- add_survival_data_GSE6891(dataset)}
  dataset <- dataset$Y
  factors <- colnames(dataset)
  dataset <- standardize_clinical_data(dataset, factors)
  
  if (upload == FALSE) {
    return(paste0("Successfully processed ", dataset_name, " clinical data."))
    }
  
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
  return(paste0("Successfully uploaded ", dataset_name, " clinical data."))
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

add_survival_data_GSE6891 <- function(dataset) {
  dataset$Y$time <- NA
  dataset$Y$death <- NA
  
  rownames(dataset$Y)
  
  GSE_6891_survival <- readRDS("data/additional_data/GSE_6891_survival.rds")
  rownames(GSE_6891_survival) <- GSE_6891_survival$volgnummer
  
  for (x in 1:length(rownames(dataset$Y))) {
    title <- dataset$p[rownames(dataset$Y)[x],1]
    title <- substr(title, 5, nchar(title))
    dataset$Y[x,3] <- GSE_6891_survival[title, 4]
    dataset$Y[x,4] <- GSE_6891_survival[title, 3]
    # Use gsub and match
  }
  
  dataset$Y$death[dataset$Y$death == "alive"] <- 0
  dataset$Y$death[dataset$Y$death == "dead"] <- 1
  
  return(dataset)
}

# Connecting to MongoDB
connect_mongo <- function(collection_name, user = "root", pass = "password", host = "localhost:27017") {
  uri <- sprintf("mongodb://%s:%s@%s/", user, pass, host)
  return (mongo(url = uri, db = "database", collection = collection_name))
  
}

return_data <- function(dataset, variable) {
  
}

