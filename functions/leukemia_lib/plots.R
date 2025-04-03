
plot_km_ds1 <- function(DS, continuous = FALSE, colors = c('darkred', 'darkblue'), ...) {
  plot.shiny.km(DS$Y$time, DS$Y$death, as.double(DS$X[1,]), title = DS$NAME,
                col = colors, continuous = continuous, ...)
}

plot_km_ds <- function(DSList, ...) {
  lapply(DSList, plot_km_ds1, ...)
}

plot_groups_from_DS <- function(DS, y, ...) {
  plot_groups(as.double(DS$X[1,]), y, DS$NAME, ...)
}

plot_groups <- function(x,y, title, no_keep = "", 
                           ylab = 'log2 expression',
                           xlab = '',
                           print = TRUE, violin = FALSE,
                           find = NULL, replace = '', return_fit = FALSE) {
  
  if (is.null(print)) {
    print <- function(...){}
    cat <- function(...){}
  }

  if (is.null(y)) {
    return(NULL)
  }
   
  if (!is.null(find)) {
    y <- gsub(find, replace, y)
  } 
  
  
  y <- toupper(y)
  no_keep <- toupper(no_keep)
  
  df <- data.frame(x,y)
  df <- dplyr::filter(df, !y %in% no_keep & !is.na(y))
 
  cat('\n\nplot for', title, 'with N =', df %>% drop_na %>% nrow(), '\n\n')
  
  fit <- lm(x ~y, data = df)

  if (return_fit) {
      return(fit)
  }
  a <- anova(fit)
  
  p <- round(a$`Pr(>F)`[1], 3)
  if (p < 0.01) {
    p <- '(P < 0.01)'
  } else {
    p <- paste0('(P = ', p, ')')
  }
  
 
 ptype <- geom_boxplot
 if (violin) {
     ptype <- geom_violin
 }

  ggplot(df, aes(x = y, y = x)) + ptype(aes(fill = y)) +
    ggtitle(paste(title, p)) + 
    theme_classic() + theme(legend.position = 'none') +
    labs(x = xlab, y = ylab)# + scale_fill_manual(values = mycols)
  
    # cc <- utils::combn(levels(df$y), 2)
    # compare <- lapply(c(1,3,2), function(i,cc)cc[,i], cc )
    # pvalues <- factor(tt$y[,4] < .05, labels = c('NS', '*'))[c(1,3,2)]
  
    # g1 <- g1 + geom_signif(comparisons = compare, step_increase = .1, 
    #                annotations = pvalues)
  
  }
 
library(scales)
km_by_group1 <- function(DS, y, find = NULL, no_keep = '', replace = '') {
  if (!is.null(find)) {
    y <- gsub(find, replace, y)
  }

  y <- toupper(y)
  no_keep <- toupper(no_keep)

  df <- data.frame(FAB.group = factor(y), time = DS$Y$time, death = DS$Y$death)
  df <- dplyr::filter(df, !y %in% no_keep & !is.na(y))

  fit <- coxph(Surv(time, death) ~ FAB.group, data = df)
  s <- survfit(Surv(time, death) ~ FAB.group, data = df)
  pvalue <- round(summary(fit)$logtest[3],3)
  if (pvalue < 0.001) {
    title <- paste0(DS$NAME, ' (P < ',0.001, ')')
  } else {
    title <- paste0(DS$NAME, ' (P = ',pvalue, ')')
  }
  
  colors <- hue_pal()(length(unique(df$FAB.group)))
  colors[2] <- 'black'
  ggsurv(s, surv.col = colors, order.legend = FALSE) + ggtitle(title) + theme_classic()
}


