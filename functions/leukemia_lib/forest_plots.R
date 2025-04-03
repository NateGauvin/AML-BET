
# Note: need 1.10:
#https://cran.r-project.org/src/contrib/Archive/forestplot/

library(forestplot)

stopifnot(sessionInfo()$otherPkgs$forestplot$Version == "1.10")


km_ci_from_DS1 <- function(DS, use_median = TRUE) {
  kk <- NULL
  if (!is.null(DS$X) && !is.null(DS$Y$time) && !is.null(DS$Y$death)) {
    kk <- km_ci(DS$Y$time, DS$Y$death, as.double(DS$X[1,]), DS$NAME,
                use_median = use_median)
  } 
  kk
}

km_ci_from_DSList <- function(DSList, use_median = TRUE) {
  lapply(DSList, km_ci_from_DS1, use_median = use_median)
}



## optionally pass a list ('kk') of km_ci values
forest_from_list <- function(DSList, km = TRUE, new_page = TRUE, 
                             use_median = TRUE, kk = NULL, ...) {
 
if (is.null(kk)) {
  if (km) {
    kk <- km_ci_from_DSList(DSList, use_median = use_median)
  } else {
    kk <- lapply(DSList, fc_risk_ci_from_ds)
  }
}

  CI <- do.call(rbind, kk)
  
  if (km) {
    xlab <- 'HR'
    xlog <- TRUE
    logHR <- log(CI$HR)
  } else {
    xlab <- 'logFC'
    xlog <- FALSE
    logHR <- CI$LogFC
  }
  
  log_hr_avg <- sum(CI$N *logHR) / sum(CI$N)
  log_hr_avg <- round(log_hr_avg,2)
  
  v <- sum(CI$N**2*CI$var / (sum(CI$N)**2))
  
  log_hr_l <- log_hr_avg - qnorm(.975)*sqrt(v)
  log_hr_u <- log_hr_avg + qnorm(.975)*sqrt(v)
  
  CI <- rbind(CI, NA)
  
  CI[nrow(CI),1] <- 'Weighted Average'

  if (km) {
    CI$HR <- round(CI$HR,2)
    CI[nrow(CI),3:5] <- c(round(exp(log_hr_avg),2),
                        exp(log_hr_l), exp(log_hr_u))
  } else {
    CI$LogFC <- round(CI$LogFC,2)
    CI[nrow(CI),3:5] <- c(round(log_hr_avg,2),
                          log_hr_l, log_hr_u)
  }
  
  
  CI[nrow(CI),2] <- sum(CI$N, na.rm=TRUE)
  
  cis <- CI[,3:5]
  
  ci_str <- paste0('(',round(CI$L,2), ', ', round(CI$U,2),')')
  
  CI <- cbind(CI, '95% CI' = ci_str)
  
  CI <- rbind(colnames(CI), CI)
  cis <- rbind(NA, cis)
  
  is_summary <- rep(FALSE, ncol(CI[,1:3])*nrow(CI[,1:3]))
  is_summary[1] <- TRUE
  is_summary[nrow(CI)] <- TRUE
  
  forestplot(CI[,c(1:3,7)], cis,
             is.summary = is_summary,
             #align = align,
             col = fpColors(box = "royalblue",
                            line = "darkblue",
                            summary = "royalblue"),
             vertices = TRUE,
             xlog = xlog,
             mar = unit(c(0,24,0,24), 'mm'),
             txt_gp = fpTxtGp(cex = .8,
                              ticks = gpar(fontfamily = "", cex = .8),
                              xlab = gpar(cex = .75)),
             xlab = xlab, new_page = new_page,
             fn.ci_norm = function(size, ...) {
               fpDrawNormalCI(size = size * 2, ...)
             },
             graph.pos = 4,
             hrzl_lines = TRUE, ...
  )
  
}


myPushViewport <- function(num, byrow, setup = FALSE) {
  # if setup is TRUE, then num is number of rows or columns
  # otherwise, move to next row/col according to 'num'
  if (setup && byrow) {
    pushViewport(viewport(layout = grid.layout(num,1)))
    pushViewport(viewport(layout.pos.row = 1))
  } else if (setup) {
    pushViewport(viewport(layout = grid.layout(1,num)))
    pushViewport(viewport(layout.pos.col= 1))
  } else if (byrow) {
    pushViewport(viewport(layout.pos.row = num))
  } else {
    pushViewport(viewport(layout.pos.col = num))
  }
}

side_by_side_forest_plots <-function(l1,  l2, l3 = NULL, km = TRUE, titles = NULL, use_median = TRUE,
                                     byrow = FALSE, ...) { 
  grid.newpage()
  ncol = 3
  if (is.null(l3)) {
    ncol = 2
  }
  
  myPushViewport(ncol, byrow, setup = TRUE)
  forest_from_list(l1, km, FALSE, title = titles[1], use_median = use_median, ...)
  popViewport()
  myPushViewport(2,byrow)
  forest_from_list(l2, km, FALSE,title = titles[2], use_median = use_median, ...)
  
  if (!is.null(l3)) {  
    popViewport(1)
    myPushViewport(3,byrow)
    forest_from_list(l3, km, FALSE, title = titles[3], use_median = use_median, ...)
  }
  
  popViewport()

}


