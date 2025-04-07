source("functions/mongoConnection.R")
library(ggplot2)

query_mongo <- function(gene_name, variable, datasets = "all") {
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
  names(all_graph_data) <- all_collections
  
  #generate boxplots for data (might need work)
  lapply(all_graph_data,generate_boxplots, gene_name, variable)
  
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
  if (!(variable %in% possible_factors)) {stop("Requested variable not found in dataset.")}
  clinical_factors <- data_connection$find(fields = paste0('{"', variable, '" : true, "_id" : false}'))
  
  graph_data <- data.frame(
    expr = expression_data,
    factor = clinical_factors
  )
  colnames(graph_data)[1] <- "expr"
  
  graph_data[3] <- dataset_name
  
  return(graph_data)
}


generate_boxplots <- function(data, gene_name, variable) {
  data <- na.omit(data)
  
  ggplot(data, aes_string(x = variable, y = "expr", fill = variable)) + 
    geom_boxplot() + theme_classic() + ggtitle(paste0(data[1,3], ": ", gene_name, " by ", variable)) + 
    labs(y = "Log Counts per Million", x = paste0(variable))
}
