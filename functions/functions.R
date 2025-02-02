###################################################################################
# plots km curves for 'high' and 'low' values relative to median (default)
# or using an optimal cutoff if optimal.cut is TRUE 
###################################################################################

library(dplyr)
library(survival)
library(GGally)
library(ggplot2)
library(readxl)
library(ggsignif)
library(tidyverse)
library(MASS)


get_multi_genes <- function(X, pl, genes, combine) {
  s <- lapply(genes, function(x) {
    get_gene(X, pl, x, combine)
  })
  
  do.call("rbind", s)
}

get_gene <- function(X, pl, gene, combine = TRUE) {
  # get expression of gene, using the probe with
  # the highest mean expression if combine is TRUE;
  # otherwise, change label to gene (probe)

  if (length(gene) > 1) {
    stop("get_gene takes a single 'gene'")
  }
  
  pl <- dplyr::filter(pl, pl$ID %in% rownames(X))
  
  qry <- paste0('\\b', gene, '\\b')
  g <- grep(qry, pl[,2])
  
  if (length(g) == 0) {
    return(NULL)
  } 
  m <- match(pl[g,1], rownames(X))
  
  cat('number of matches: ', length(m), '\n')
  
  xx <- X[m,,drop = FALSE]
  
  if (nrow(xx) > 1 && combine) {
    w <- 1
    w <- which.max(rowMeans(xx))
    xx <- xx[w, , drop = FALSE]
  }
  
  if (nrow(xx) == 1 && combine) {
    rownames(xx) <- gene  
  } else {
    rownames(xx) <- paste0(gene, " (", rownames(xx), ")")
  }

  xx
}

plot.shiny.km <- function(time, death, x, title = "", ids = NULL, 
                          subset = rep(TRUE, length(time)), 
                          col = NULL,  xlab = NULL, ylab = NULL, hr.inverse = FALSE, 
                          ret.plot = FALSE, 
                          continuous = TRUE,
                          optimal.cut = FALSE, ...) {

  # if continuous is TRUE, we return data frame of median + continuous results
  # if ret.plot is TRUE, return plot object but do not print it
  
  ## filter out missing data ##
  subset = subset & !is.na(time) & !is.na(death) & !is.na(x) & !is.nan(x)
  x = x[subset]; time = time[subset]; death = death[subset]; 
  if (!is.null(ids)) ids = ids[subset]
  if (length(x) ==0 | length(time) == 0 | length(death) == 0 
      | length(unique(death)) < 2) {
    return(invisible(NULL))
  }
  if (is.null(ids)) ids = 1:length(x)
  km.data = data.frame(ID = ids, X = x, Group = NA, time = time, event = death)

  ## settings for median cutoff ##
  cut = median(x)
  p.adj = NULL
  upper = "upper 50%"; lower = "lower 50%"

  p.adj = NA  # for cutp only
  
  ## find optimal cutoff if specified ##
  if (optimal.cut) {
     mod <- coxph(Surv(time, event) ~ X, data = km.data)
     cc = try(survMisc::cutp(mod), silent = TRUE)
     if (class(cc) %in% "try-error") {
	    return(invisible(NULL))
     }
     cut = cc$X$X[1]
     p.adj = cc$X$p[1]
     percentile = round(sum(x>=cut) / length(x) * 100,2)
     upper = paste0("upper ", percentile, "%")
     lower = paste0("lower ", 100-percentile, "%")
  }
    
  ## split into high and low groups using appropriate cutoff ## 
  newX = x
  newX[x >= cut] = upper
  newX[x < cut] = lower
  expression = as.factor(newX)
  km.data$Group = expression
  
  n = length(levels(expression))
  km.group = try(coxph(Surv(time, death) ~ expression), silent = TRUE)
  if (class(km.group) %in% "try-error") return(invisible(NULL))
  p.km = 1 - pchisq(km.group$score,n-1)
  hr = exp(km.group$coefficients)
  n = km.group$n
  
  if (hr.inverse) hr = 1/hr 
  
  if (continuous) {
    km.con = try(coxph(Surv(time, death) ~ x), silent = TRUE)
    if (class(km.con) %in% "try-error") return(invisible(NULL))
    
    p.con = 1 - pchisq(km.con$score,1)
    
    # use adjusted p if available
    if (!is.na(p.adj)) {
      p.km <- p.adj
    } 
    
    hr.con = exp(km.con$coefficients)
    if (hr.inverse) hr.con <- 1/hr.con
    df <- data.frame(hr_median = hr, p_median = p.km,
                     hr_continuous = hr.con, p_continuous = p.con)
    rownames(df) <- title
    return(df)
  }
  
  
  hr.str = paste("HR = ", round(hr,2), ", ")
  p.str = paste0("P = ", round(p.km,4))
  if (!is.null(p.adj) && !is.na(p.adj) ) {
    p.str = paste0(p.str, ", P(adjusted) = ", round(p.adj,4))
  }
  
  if (title=="") { 
    title = paste(hr.str, p.str, sep = "")
  } else {
    title = paste0(title, " (",hr.str, p.str,")")
  }
  
  if (is.null(xlab)) xlab = "Time"
  if (is.null(ylab)) ylab = "Survival"
  
  if (is.null(col)) {
    stop("need to specify col (color)")
  }
  
  ## plot graph ### ggplot2/GGally form
  km.group1 = survfit(Surv(time, death) ~ expression)
  km.group = ggsurv(km.group1, 
                    main = title, 
                    surv.col = col, cens.col = col, xlab = xlab, ylab = ylab, ...) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) + theme_classic()
  #    theme(plot.title = element_text(vjust = 0.5, colour = "black"))
  
  if (ret.plot) {
    return(km.group)
  }
  plot(km.group)
}

plot_all_kms <- function(X, time, outcome) {
  XX <- vector('list', nrow(X))
  for (i in 1:nrow(X)) {
    XX[[i]] <- plot.shiny.km(time, outcome, X[i,], 
                           col = c("darkblue", "darkred"), 
                           title = rownames(X)[i], 
                           ret.plot = TRUE, xlab = 'Time (months)') + 
                           theme(text = element_text(size = 8))
  }
  return(XX)
}

calculate_all_kms <- function(time, outcome, X) {
  s <- lapply(1:nrow(X), function(i) {
      plot.shiny.km(time, outcome, X[i,], title = rownames(X)[i], continuous = TRUE)
  })
  do.call('rbind', s)
}

# displays table with rows highlighted based on
# hr_median, p_median, hr_con, and p_con, and top 
# borders separating different datasets
displayTable <- function(df) {
  
  ds <- unique(df$dataset)
  
  # determine which rows begin a new dataset
  m <- match(ds, df$dataset)
  new_datasets <- rep(FALSE, nrow(df))
  new_datasets[m] <- TRUE
  
  
  myvals <- rep("", nrow(df))
  myvals[df$p_median < 0.05 & df$hr_median > 1 |
           df$p_continuous < 0.05 & df$hr_continuous > 1] <- 'up'
  myvals[df$p_median < 0.05 & df$hr_median < 1 |
           df$p_continuous < 0.05 & df$hr_continuous < 1] <- 'down'
  myvals[df$p_median < 0.05 & df$p_continuous < 0.05 & 
           sign(1-df$hr_median) != sign(1-df$hr_continuous)] <- 'conflict'
           
  df <- data.frame(df, significant = myvals, new_ds = new_datasets)
  
  hideCols <- ncol(df)
  hideCols <- (hideCols-1):hideCols - 1
  
  datatable(df, rownames = FALSE, options = list(dom = 't',
                               ordering = FALSE,
                               pageLength = nrow(df),
                               
                               ## collapse border so we can add to rows
                               initComplete = JS(
                                 "function(settings, json) {",
                                 "$('table').css({'border-collapse': 'collapse'});",
                                 "}"),
                               
                               columnDefs = list(list(visible=FALSE,
                                                      targets=hideCols)))) %>%
    formatStyle('significant', target = 'row', 
                fontWeight = styleEqual(c('up', 'down'), c('bold', 'bold')),
                backgroundColor = styleEqual(c('up', 'down', 'conflict'), 
                                             c('yellow', 'lightblue', 'pink'))) %>%
    formatStyle('new_ds', target = 'row', 
                borderTop = styleEqual(c(FALSE, TRUE), c('none', '2px solid black'))) %>%
    
    formatRound(columns = ~ hr_median + p_median + hr_continuous + p_continuous, digits = 2)
}

get_results <- function(X, PL, GENE.LIST, times, events) {
  G <- get_multi_genes(X, PL, GENE.LIST, combine = TRUE)
  KM <- plot_all_kms(G, times, events)

  G <- get_multi_genes(X, PL, GENE.LIST, combine = FALSE)
  RES <- calculate_all_kms(times, event, G)
  
  list(KM = KM, RES = RES)
  
}

event_by_contains <-function(x, alive, dead, fixed = FALSE) {
  g1 <- grep(alive, x, fixed = fixed)
  g2 <- grep(dead,x, fixed = fixed)
  if (length(g1) == 0) {
    stop ('no alives are found')
  } else if (length(g2) == 0) {
    stop('no deads are found')
  } else if (length(intersect(g1,g2)) > 0) {
    stop("someone is alive AND dead!")
  }
  event <- rep(NA, length(x))
  event[g1] <- 0
  event[g2] <- 1
  event
}

##################################################
### process favorable data
##################################################
get_favorable <- function(sheet, sample_names) {
  r <- read_excel("data/Favorable_Status.xlsx", sheet = sheet, 
                  col_names = TRUE)
  m <- match(sample_names, r$ID)
  r[m,]
}

favorable_ordering <- function(x) {
  x <- toupper(x)
  if ("ADVERSE" %in% x) {
    return(5)
  } else if ("POOR" %in% x) {
    return(4.5)
  } else if ("HIGH" %in% x) {
    return(4.4)
  }else if ("INTERMEDIATE_2" %in% x) {
    return(4)
  } else if ("INTERMEDIATE_1" %in% x) {
    return(3)
  } else if ("INTERMEDIATE" %in% x) {
    return(2)
  } else if ("INTERMEDIATE/NORMAL" %in% x) {
    return(1.5)
  } else if ("STANDARD" %in% x) {
    return(1.3)
  } else if ("FAVORABLE" %in% x) {
    return(1)
  } else if ("LOW" %in% x) {
    return(0.9)
  }
  
  stop('No assignment for: ', x[1])
}


plot_favorable <- function(x,y, title, no_keep = "NA", 
                           ylab = 'log2 expression',
                           dots = FALSE,
                           lower_p = FALSE, return_plot = FALSE, print = TRUE) {
  
  if (is.null(print)) {
    print <- function(...){}
    cat <- function(...){}
  }

  if (is.null(y)) {
    return(NULL)
  }
    
  y <- tolower(y)
  no_keep <- tolower(no_keep)
  
  df <- data.frame(x,y)
  df <- dplyr::filter(df, !y %in% no_keep & !is.na(y))
 
  df$y <- reorder(df$y, df$y, favorable_ordering)
  
  cat('\n\nfavorabilty for ', title, ' N = ', df %>% drop_na %>% nrow(), '\n\n')
  
  fit <- lm(x ~y, data = df)
  a <- anova(fit)
  
  
  tt <- TukeyHSD(aov(x~y, data = df))
  print(tt)
  
  p <- round(a$`Pr(>F)`[1], 3)
  if (p < 0.01) {
    p <- '(P < 0.01)'
  } else {
    p <- paste0('(P = ', p, ')')
  }
  
  if (lower_p) {
    p <- gsub('P','p', p)
  }
  
  mycols <- c(3,4,2)
  
  if (dots) {
    mycols <- c('darkgreen', 'darkblue', 'darkred')
    g1 <- ggplot(df, aes(x=y, y = x)) +
      geom_point(position = position_jitter(h=0,w=.25), 
                 aes(color = y), na.rm = TRUE) +
      stat_summary(fun="mean", geom="crossbar", width = 0.5) +
      theme_classic() + theme(legend.position = 'none') +
      labs(x = '', y = ylab) + scale_color_manual(values = mycols) +
      ggtitle(paste0(title, p, 'N=', nrow(df)))
  } else {
  g1 <- ggplot(df, aes(x = y, y = x)) + geom_boxplot(aes(fill = y)) +
    ggtitle(paste0(title, p)) + 
    theme_classic() + theme(legend.position = 'none') +
    labs(x = '', y = ylab) + scale_fill_manual(values = mycols)
  
    # cc <- utils::combn(levels(df$y), 2)
    # compare <- lapply(c(1,3,2), function(i,cc)cc[,i], cc )
    # pvalues <- factor(tt$y[,4] < .05, labels = c('NS', '*'))[c(1,3,2)]
  
    # g1 <- g1 + geom_signif(comparisons = compare, step_increase = .1, 
    #                annotations = pvalues)
  
  }
  
  
  print(g1)
  return(g1)

}

# groups -- group labels in order ('favorable', 'intermediate', 'adverse')
plot_km_favorability <- function(survival, status, fav, title = 'title',
                                 groups = NULL,
                                 colors = c('darkgreen', 'darkblue', 'darkred'),
                                 xlab = 'Time (Months)') {
  if (is.null(groups)) {
    stop('groups must be specified')
  }

  fav <- tolower(fav); groups <- tolower(groups)
  fav[!fav%in%groups] <- NA
  fav <- factor(fav, levels = groups)
  
  cat('\n\n', title, ', N = ', sum(!is.na(fav) & !is.na(survival) & !is.na(status)))
  
  pp <- coxph(Surv(survival, status) ~ fav)
  s <- summary(pp)
  p<- s$logtest[3]
  p <- round(p,3)
  if (p < 0.01) {
    p <- 'P < 0.01'
  } 
  
  title <- paste0(title,'\n','(',p,')')
  
  risk <- fav
  s <- survfit(Surv(survival, status) ~ risk)

  g1 <- ggsurv(s, size.est = 1, order.legend = FALSE) + theme_classic() +
          scale_color_manual(values=colors) +
          ggtitle(title) + labs(x = xlab, y = 'Survival probability')

  print(g1)
  g1 + scale_linetype_manual(values = c('twodash', 'dotted', 'solid'))
}

# dd <- data.frame(x = rnorm(50), group = sample(c('a','b','c'),50, replace = TRUE))
# 
# ggplot(dd, aes(x=group, y = x)) +
# geom_point(position = position_jitter(h=0,w=.2), aes(colour = group), na.rm = TRUE) +
#   stat_summary(fun="mean", 
#                geom="crossbar", width=0.5)
#   
#   
#    geom_errorbar(stat = "summary", fun.y = "mean.no.na", width=0.8,
#                  aes(ymax=..y..,ymin=..y..)) + theme_classic()


km_ci <- function(time, death, x, title, use_median = TRUE) {
  subset = !is.na(time) & !is.na(death) & !is.na(x) & !is.nan(x)
  x = x[subset]; time = time[subset]; death = death[subset]; 

  if (use_median) {  
    risk <- x >= median(x)
  } else {
    risk <- x
  }
  
  res <- coxph(Surv(time, death) ~ risk)
  
  ci <- confint(res)
  df <- data.frame(Dataset = title, N = length(x),
    HR = exp(res$coefficients), L=exp(ci[1]), U = exp(ci[2]),
    var = res$var[1])
  df
  
}



km_multi <- function(time, death, x, risk) {
  subset = !is.na(time) & !is.na(death) & !is.na(x) & 
    !is.nan(x) & !is.na(risk)
  x = x[subset]; risk = risk[subset]
  time = time[subset]; death = death[subset]; 
  
  res <- coxph(Surv(time, death) ~ risk + x)
  res
  
}

km_ci_multi_from_DS <- function(DS) {

  time <- DS$Y$time
  death <- DS$Y$death
  risk <- DS$Y$risk
  x <- as.double(DS$X[1,])
  title <- DS$NAME

  if (is.null(time) || is.null(death) || is.null(risk)) {
	return(NULL)
  }

  res <- km_multi(time, death, x, risk)
  ci <- confint(res)
  nr <- nrow(ci)
  df <- data.frame(Dataset = title, N = res$n,
    HR = exp(res$coefficients[nr]), L=exp(ci[nr,1]), U = exp(ci[nr,2]),
    var = res$var[1])
  df
}

# get correlations
get_cor1 <- function(x,y) {
  res <- cor.test(x,y)
  res2 <- cor.test(x,y, method = "spearman")
  data.frame(pearson_est = res$estimate, pearson_p = res$p.value,
             spearman_est = res2$estimate, spearman_p = res2$p.value)
}

get_cor <- function(x, X) {
  K <- apply(X, 1, get_cor1, y = x)
  do.call('rbind', K)
}

plot_linear_svm <- function(x1, x2, risk, keep, ds, rev_colors = FALSE) { 
  
  df <- data.frame(x1, x2, y = factor(risk))
  
  df2 <- dplyr::filter(df, y %in% keep) #%>% dplyr::mutate(y = factor(y))
  
  library(e1071)
  svmfit = svm(y ~ ., data = df2, kernel = "linear", cost = 1, scale = FALSE) #, class.weights = 'inverse')
  
  pred <- fitted(svmfit)
  acc1 <- sum(pred[df2$y == keep[1]] == df2$y[df2$y==keep[1]]) / sum(df2$y == keep[1])
  acc2 <- sum(pred[df2$y == keep[2]] == df2$y[df2$y==keep[2]]) / sum(df2$y == keep[2])
  acc <- sum(acc1,acc2)/2 * 100

  # make.grid = function(x, n = 75) {
  #   grange = apply(x, 2, range)
  #   x1 = seq(from = grange[1,1], to = grange[2,1], length = n)
  #   x2 = seq(from = grange[1,2], to = grange[2,2], length = n)
  #   expand.grid(x1 = x1, x2 = x2)
  # }
  # 
  # xgrid <- make.grid(dplyr::select(df, x1,x2), n = 40)
  # ygrid = predict(svmfit, xgrid)
  # df_grid <- data.frame(xgrid, y = ygrid)
  # 

  beta = drop(t(svmfit$coefs)%*%as.matrix(df2[svmfit$index,1:2]))
  beta0 = svmfit$rho
  
  colors <- rev(c('forestgreen', 'darkred'))
  if (rev_colors) {
    colors <- rev(colors)
  }
  
  ggplot()  + #geom_point(data = df_grid, aes(x1,x2,col=y), size = .1) + 
    geom_point(data = df2, aes(x1,x2,col=y), size = 2) + theme_classic() + 
    geom_abline(slope = -beta[1]/beta[2], intercept = beta0/beta[2]) + 
    labs(col = 'risk', x = 'ATF4 expression', y = 'NOXA expression') +
    scale_color_manual(values = colors) + ggtitle(paste0(ds, '\n(accuracy = ', round(acc),'%)'))
}



plot_qda <- function(x1,x2,y,keep, ds, rev_colors = FALSE) {
  
  # construct the model
  df <- data.frame(x1=x1,x2=x2,y=y) %>% filter(y %in% keep)
  
  mdl <- lda(y ~ x1 + x2, data = df)
  
  pred <- predict(mdl)$class
  
  acc1 <- sum(pred[df$y == keep[1]] == df$y[df$y==keep[1]]) / sum(df$y == keep[1])
  acc2 <- sum(pred[df$y == keep[2]] == df$y[df$y==keep[2]]) / sum(df$y == keep[2])
  acc <- sum(acc1,acc2)/2 * 100
  
  np <- 300
  nd.x <- seq(from = min(df$x1), to = max(df$x1), length.out = np)
  nd.y <- seq(from = min(df$x2), to = max(df$x2), length.out = np)
  nd <- expand.grid(x1 = nd.x, x2 = nd.y)
  
  prd <- as.numeric(predict(mdl, newdata = nd)$class)
  df_contour <- data.frame(x1 = nd$x1, x2 = nd$x2, z = prd )
  
  colors <- rev(c('forestgreen', 'darkred'))
  if (rev_colors) {
    colors <- rev(colors)
  }
  
  g1<-ggplot() + geom_point(data = df %>% drop_na(), aes(x1,x2,col=y), size = 2) + 
    theme_classic() + scale_color_manual(values = c('darkred', 'forestgreen')) +
    stat_contour(data = df_contour, aes(x=x1,y=x2,z=z)) +
    labs(col = 'risk', x = 'ATF4 expression', y = 'NOXA expression') +
    scale_color_manual(values = colors) + ggtitle(paste0(ds, '\n(accuracy = ', round(acc),'%)'))
  
  np <- 40
  nd.x <- seq(from = min(df$x1), to = max(df$x1), length.out = np)
  nd.y <- seq(from = min(df$x2), to = max(df$x2), length.out = np)
  nd <- expand.grid(x1 = nd.x, x2 = nd.y)
  
  df2 <- data.frame(nd, pred = predict(mdl, data.frame(x1 = nd$x1, x2 = nd$x2))$class)
  
  g1 + geom_point(data = df2, aes(x=x1,y=x2,col=pred),size = .3)
  
  
}



fc_risk_ci_from_ds <- function(DS) {
  
  tmp <- DS$Y$risk
  
  if (is.null(tmp)) {
    return(NULL)
  }
  
  y <- rep(NA, length(tmp)) 
  y[tolower(tmp)%in% c('favorable', 'low')] <- 0
  y[tolower(tmp)%in% c('adverse', 'high', 'poor')] <- 1
  
  res <- lm(as.double(DS$X[1,]) ~ y)
  
  n <- sum(!is.na(y))
  
  ci <- confint(res)
  var <- vcov(res)[2,2]
  
  df <- data.frame(Dataset = DS$NAME, N = n,
                   LogFC = res$coefficients[2], L=ci[2,1], U = ci[2,2],
                   var = var)
  df
  
}
