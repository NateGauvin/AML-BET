my_summary1 <- function(xcol, p_valid = 0.1) {
  tt <- table(xcol)
  valid <- sum(tt > p_valid*length(xcol))
  if (sum(valid) >= 2) {
    return (tt)
  } else {
    return(NULL)
  }
}

summarize_table <- function(df, p_valid = 0.1, show_null = FALSE) {
  a <- apply(df, 2, my_summary1, p_valid = p_valid)
  if (!show_null) {
      a[lengths(a)>0]
  }
}

ds_filter_by_gender1 <-function(DS, gender) {
  if (is.null(DS$Y$gender)) {
    return(NULL)
  }
  
  keep <- DS$Y$gender == gender
  DS$Y <- DS$Y[keep,]
  DS$X <- DS$X[,keep,drop = FALSE]
  if (!is.null(DS$p)) {
      DS$p <- DS$p[keep,]
  }

  DS
}

ds_filter_by_gender <- function(DSList, gender) {
  res <- lapply(DSList, ds_filter_by_gender1, gender = gender)
  res[lengths(res) > 0]
}


