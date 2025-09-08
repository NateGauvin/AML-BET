source("functions/mongoConnection.R")
library(ggplot2)
library(survival)
library(survminer)

query_mongo <- function(gene_name, variable, datasets = "all", method) {
  connection <- connect_mongo("temp_collection")
  all_collections <- connection$run('{ "listCollections" : 1, "nameOnly" : true}')
  all_collections <- all_collections$cursor$firstBatch[,1]
  expr_check <- lapply("_expr", grepl, all_collections)
  expr_check <- which(expr_check[[1]])
  all_collections <- all_collections[expr_check]
  all_collections <- gsub('.{5}$', '', all_collections)
  
  #check datasets
  if (datasets != "all") {
    datasets <- strsplit(datasets, ",")
    if (any(is.na(match(datasets[[1]], all_collections)))) {
      stop('Invalid datasets entered')
    }
    all_collections <- datasets
  }
  
  #obtain data from mongo
  all_graph_data <- lapply(all_collections, return_data, gene_name, variable)
  all_graph_data <- lapply(1:length(all_graph_data), make_data_list, all_graph_data, all_collections)
  
  #generate boxplots for data (might need work)
  lapply(all_graph_data, generate_boxplots, gene_name, variable, method)
  
}

make_data_list <- function(index, data, names) {
  data[index] <- c(names[index], data[index])
}

return_data <- function(dataset_name, gene_name, variable) {
  data_connection <- connect_mongo(paste0(dataset_name, "_expr"))
  genes_in_dataset <- data_connection$find(fields = '{"gene" : true, "_id" : false}')
  if (!(gene_name %in% genes_in_dataset[[1]])) {
    print("Gene not found in dataset, returning NA")
    return(NA)
    }
  
  expression_data <- data_connection$find(paste0('{ "gene" : "', gene_name, '"}'))
  expression_data <- as.vector(expression_data[2])
  
  data_connection <- connect_mongo(paste0(dataset_name, "_clinical"))
  possible_factors <- data_connection$find(limit = 1)
  possible_factors <- colnames(possible_factors)[2:length(colnames(possible_factors))]
  if (!(variable %in% possible_factors)) {
    print("Requested variable not found in dataset.")
    return(NA)
    }
  clinical_factors <- data_connection$find(fields = paste0('{"', variable, '" : true, "_id" : false}'))
  
  graph_data <- data.frame(
    expr = expression_data,
    factor = clinical_factors
  )
  colnames(graph_data)[1] <- "expr"
  
  return(graph_data)
}

return_data_km <- function(dataset_name, gene_name) {
  data_connection <- connect_mongo(paste0(dataset_name, "_expr"))
  genes_in_dataset <- data_connection$find(fields = '{"gene" : true, "_id" : false}')
  if (!(gene_name %in% genes_in_dataset[[1]])) {
    print("Gene not found in dataset, returning NA")
    return(NA)
  }
  
  expression_data <- data_connection$find(paste0('{ "gene" : "', gene_name, '"}'))
  expression_data <- as.vector(expression_data[2])
  
  data_connection <- connect_mongo(paste0(dataset_name, "_clinical"))
  possible_factors <- data_connection$find(limit = 1)
  possible_factors <- colnames(possible_factors)[2:length(colnames(possible_factors))]
  if (!("time" %in% possible_factors)) {
    print("Survival data not found in dataset.")
    return(NA)
    }
  if (!("death" %in% possible_factors)) {
    print("Survival data not found in dataset.")
    return(NA)
    }
  time <- data_connection$find(fields = paste0('{"time" : true, "_id" : false}'))
  death <- data_connection$find(fields = paste0('{"death" : true, "_id" : false}'))
  
  graph_data <- data.frame(
    expr = expression_data,
    time = time,
    death = death
  )
  colnames(graph_data) <- c("expr", "time", "death")
  
  return(graph_data)
}

generate_boxplots <- function(data, gene_name, variable, method) {
  if (!is.data.frame(data[[2]])) {return(NA)}
  
  dataset <- data[[1]]
  data <- na.omit(data[[2]])
  
  if      (method == "t-test") {test <- t.test}
  else if (method == "wilcox") {test <- wilcox.test}
  else if (method == "anova")  {test <- aov}
  
  if (method != "anova") {
    data[data[,2] == "intermediate",]  <- NA
    data <- na.omit(data)
  }   
  
  test_data <- test(formula = data$expr ~ data$risk,
                        var.equal = FALSE)
    
  if (method == "t-test") {test_label <- paste0("FC: ", trunc(2**(test_data$estimate[2] - test_data$estimate[1]) * 1000)/1000, ", ")}
  else if (method == "wilcox") {test_label <- paste0("AUC: ", trunc(test_data$statistic / 
                                                                      (length(data[data[,2] == "poor",][[1]]) *
                                                                         length(data[data[,2] == "favorable",][[1]]))
                                                                    * 1000)/1000,", ")}
  else {test_label <- ""}
  
  if (method == "anova") {p_label <- paste0("p-value: ", trunc(summary(test_data)[[1]]$`Pr(>F)`[1] * 1000)/1000)}  
  else {p_label <- paste0("p-value: ", trunc(test_data$p.value * 1000)/1000)}

  # Make sure favorable = "lightblue", poor = "darkred"
  
  ggplot(data, aes_string(x = variable, y = "expr", fill = variable)) + 
    geom_boxplot() + theme_classic() + 
    ggtitle(paste0(dataset, ": ", gene_name, 
                   "\n", test_label, p_label)) + 
    labs(y = "Log Counts per Million", x = paste0(variable)) + 
    theme(legend.position = "none")
}

generate_km <- function(data, gene_name) {
  if (!is.data.frame(data[[2]])) {return(NA)}
  
  dataset <- data[[1]]
  data <- na.omit(data[[2]])
  cutoff <- median(data[,1])
  
  data$group <- ifelse(data[,1] <= cutoff, "Low", 'High')
  
  haz_ratio <- coxph(Surv(time = data[,2], event = data[,3]) ~ data$group, data = data)
  hr <- trunc((1 / exp(haz_ratio$coefficients)) * 1000)/1000
  fit <- survfit(Surv(time = data[,2], event = data[,3]) ~ data$group, data = data)
  
  ggsurvplot(fit, data = data, pval = TRUE,
             title = paste0(dataset, ", ", "HR: ", hr),
             legend.title = "Placeholder", xlab = "Time (Days?)",
             legend.labs = c(paste0("High (> ", trunc(cutoff * 1000)/1000, ")"), 
                             paste0("Low (<= ", trunc(cutoff * 1000)/1000,")")))
}

query_mongo_km <- function(gene_name, datasets = "all") {
  connection <- connect_mongo("temp_collection")
  all_collections <- connection$run('{ "listCollections" : 1, "nameOnly" : true}')
  all_collections <- all_collections$cursor$firstBatch[,1]
  expr_check <- lapply("_expr", grepl, all_collections)
  expr_check <- which(expr_check[[1]])
  all_collections <- all_collections[expr_check]
  all_collections <- gsub('.{5}$', '', all_collections)
  
  #check datasets
  if (datasets != "all") {
    datasets <- strsplit(datasets, ",")
    if (any(is.na(match(datasets[[1]], all_collections)))) {
      stop('Invalid datasets entered')
    }
    all_collections <- datasets
  }
  
  #obtain survival data from mongo
  all_graph_data <- lapply(all_collections, return_data_km, gene_name)
  all_graph_data <- lapply(1:length(all_graph_data), make_data_list, all_graph_data, all_collections)
  
  #generate KM curves
  lapply(all_graph_data, generate_km, gene_name)
}
