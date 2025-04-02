#library(ggplot2)
#library(dplyr)
library(shiny)
library(shinylive)
setwd("C:/Users/nateg/Downloads/GitHub/AML-BET/AML-BET_Prototype")

#beat <- readRDS("C:/Users/nateg/Downloads/GitHub/AML-BET/data/BEAT_BMA.rds")


# riskByExpression <- function(dataset, gene) {
#   filterset <- filter(dataset$X, rownames(dataset$X) == gene)
#   filterset <- rbind(filterset, dataset$Y$risk)
#   rownames(filterset) <- c("expression", "risk")
#   filterset <- t(filterset)
#   filterset <- na.omit(filterset)
#   filterset <- data.frame(filterset)
#   filterset$expression <- as.numeric(filterset$expression)
#   
#   ggplot(filterset, aes(x = factor(risk, levels = c("Favorable", "Intermediate", "Adverse")), y = expression, fill = risk)) + 
#     geom_boxplot() + labs(x = "Level of Risk", y = "Expression (Log Counts per Million)") +
#     ggtitle(paste("Expression of Gene", gene, "by Risk")) + theme_classic(base_size = 18)
# }
