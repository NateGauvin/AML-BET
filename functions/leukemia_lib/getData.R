
source('functions/functions.R')

# returns gene1 / (gene2 + gene3) if ratio = TRUE;
# otherwise returns gene1 + gene2
get_all_data <- function(gene1, gene2 = NULL, gene3 = NULL, no_y = FALSE, 
                         RNA_SEQ_only = FALSE,
                         SPLIT_DENOVO = FALSE,
                         NO_RECURRENT = FALSE, ...) {

  if (NO_RECURRENT) SPLIT_DENOVO <- TRUE
  
  if (RNA_SEQ_only) {
    LL <- list(
      TCGA = get_tcga(gene1, gene2, gene3, ...),
      TARGET = get_target(gene1, gene2, gene3, ...),
      BEAT_BMA = get_beat_aml('BMA', gene1, gene2, gene3, ...),
      BEAT_PB = get_beat_aml('PB', gene1, gene2, gene3, ...)  
    )
  } else {
    LL <- list(
      TCGA = get_tcga(gene1, gene2, gene3, ...),
      TARGET = get_target(gene1, gene2, gene3, ...),
      BEAT_BMA = get_beat_aml('BMA', gene1, gene2, gene3, ...),
      BEAT_PB = get_beat_aml('PB', gene1, gene2, gene3, ...),
      GSE37642_1 = get_GSE37642(1, gene1, gene2, gene3, ...),
      GSE37642_2 = get_GSE37642(2,gene1, gene2, gene3, ...),
      GSE6891_1 = get_GSE6891(1, gene1, gene2, gene3, ...),
      GSE6891_2 = get_GSE6891(2, gene1, gene2, gene3, ...),
      GSE71014 = get_GSE71014(gene1, gene2, gene3, ...)
    )
  }
  
  if (SPLIT_DENOVO) {
    
    pb_denovo <- filter_ds(LL$BEAT_PB, LL$BEAT_PB$p$isDenovo == TRUE,
                               name = 'BEAT AML PB, de novo')
    
    pb_recurrent <- filter_ds(LL$BEAT_PB, LL$BEAT_PB$p$isDenovo == FALSE,
                                  name = 'BEAT AML PB, recurrent')
    
    bma_denovo <- filter_ds(LL$BEAT_BMA, LL$BEAT_BMA$p$isDenovo == TRUE,
                                name = 'BEAT AML BMA, de novo')
    
    bma_recurrent <- filter_ds(LL$BEAT_BMA, LL$BEAT_BMA$p$isDenovo == FALSE,
                                   name = 'BEAT AML BMA, recurrent')
    
    LL$BEAT_BMA <- NULL
    LL$BEAT_PB <- NULL
    
    LL$BEAT_BMA_DENOVO <- bma_denovo
    LL$BEAT_PB_DENOVO <- pb_denovo
    
    if (!NO_RECURRENT) {
      LL$BEAT_BMA_RECURRENT <- bma_recurrent
      LL$BEAT_PB_RECURRENT <- pb_recurrent
    }
    
  }
  
  
  
  if (no_y) {
    LL <- lapply(LL, function(x) {
      x$Y <- NULL
      x$p <- NULL
      return (x)
    })
  }

  # mynames <- names(LL)
  # for (n in mynames) {
  #   if (is.null(LL[[n]]$X)) {
  #     LL[[n]] <- NULL
  #   }
  # }
  # 
  LL
    
}

filter_ds <- function(ds, keep, name) {
  ds$X <- ds$X[,keep, drop = FALSE]
  ds$Y <- ds$Y[keep,]
  if (!is.null(ds$p)) {
    ds$p <- ds$p[keep,]
  }
  ds$NAME <- name
  ds
}




# get gene from platform, if multiple probes/genes take one with
# highest mean
getGeneFromPlatform <- function(X, gene1, PL, gene2 = NULL) {
  if ('GPL570' == PL) {
    load('data/RData/GPL570.RData')
    PL <- GPL570
  } else if ('GPL96' == PL) {
    load('data/RData/GPL96.RData')
    PL <- GPL96
  } else if ('ENSEMBL' == PL) {
    PL <- readRDS('data/ensemble.rds')
  } else if ('GPL10558' == PL) {
    load('data/RData/GPL10558.RData')
    PL <- GPL10558
  } else {
    stop('platform ', PL, ' not found')
  }
  
  get_gene(X, PL, gene1)
}

sum_genes <- function(x1,x2) {
  if (is.null(x1) || is.null(x2)) {
    return(NULL)
  }
  #x1 + x2
  matrix(scale(as.double(x1)) + scale(as.double(x2)), 
         nrow = 1)
  
}

invalidParams <- function(gene1, gene2, gene3, sum) {
  missing_gene <- sum && (is.null(gene1) || is.null(gene2))
  extra_gene <- sum && !is.null(gene3)
  missing_gene | extra_gene
}
  
get_target <- function(gene1, gene2 = NULL, gene3 = NULL, sum = FALSE) {
  
  if (invalidParams(gene1, gene2, gene3, sum)) {
    stop("invalid params")
  }
  
  library(EnsDb.Hsapiens.v86)
  
  TARGET <- readRDS('data/target.rds')
  PL <- readRDS('data/ensemble.rds')

  x <- getGeneFromPlatform(TARGET$X, gene1, 'ENSEMBL')
  
  if (identical(gene2, NULL)) {
    TARGET$X <- x
  } else {
    x2 <- getGeneFromPlatform(TARGET$X, gene2, 'ENSEMBL')
    if (!is.null(gene3)) {
      x3 <- getGeneFromPlatform(TARGET$X, gene3, 'ENSEMBL')
      TARGET$X <- x - 0.5*(x2+x3)
    } else {
      TARGET$X <- x - x2
      if (sum) {
        TARGET$X <- sum_genes(x, x2)
      }
    }
    
  }
  
  return(TARGET)

}

subtract_genes <- function(x1,x2, x3, gene1, gene2, gene3) {
  
  gene_not_found <-function(x,gene) {
    if (!is.null(gene) && is.null(x)) {
      message('gene ', gene, '  not found')
      return(TRUE)
    }
    return (FALSE)
  }
  
  if (gene_not_found(x1,gene1) ||
      gene_not_found(x2,gene2) ||
      gene_not_found(x3,gene3) ) {
    return(NULL)
  }
  
  if (is.null(x1)) {    
    return(NULL)
  } else if (is.null(x2) && is.null(gene2)) {  # only x1 valid
    return(x1)
  } else if (is.null(x3) && is.null(gene3)) {  # only x1-x2 valid
    x <- x1 - x2
    rownames(x) <- paste0(gene1, '-', gene2)
  } else if (!is.null(x3)) { # x1 - x3 are valid
    x <- x1 - .5*(x2+x3)
    rownames(x) <- paste0(gene1, '-', gene2, ':', gene3)
  } else {
    stop('at least one gene is not found')
  }
  
  x
}

# returns row index of gene 'gene1' in DS$X
grep_gene <- function(DS, gene1, beg = '^', end = '$') {
  
  if (is.null(gene1)) {
    return(NULL)
  }
  
  qry <- paste0(beg, gene1,end)
  g <- grep(qry, rownames(DS$X))

  if (length(g) == 0) {
    message(gene1, ' not found in ', DS$NAME)
    return(NULL)
  } else if (length(g) > 1) {
    #stop(gene1, ' has multiple matches in', DS$NAME)
      w <- 1
      w <- which.max(rowMeans(DS$X[g,]))
      DS$X <- DS$X[w, , drop = FALSE]
      return(DS$X)
  }
  DS$X[g,,drop = FALSE]
}

get_tcga <- function(gene1, gene2 = NULL, gene3 = NULL, sum = FALSE) {

    if (invalidParams(gene1, gene2, gene3, sum)) {
      stop("invalid params")
    }
  
  TCGA <- readRDS('data/TCGA.rds')

  x1 <- grep_gene(TCGA, gene1, end = '\\|')
  x2 <- grep_gene(TCGA, gene2, end = '\\|')
  x3 <- grep_gene(TCGA, gene3, end = '\\|')  
  
  if (sum) {
    TCGA$X <- sum_genes(x1,x2) 
  } else {
    TCGA$X <- subtract_genes(x1,x2,x3, gene1, gene2,gene3) 
  }
  
  TCGA
}


get_beat_aml <- function(type, gene1, gene2 = NULL, gene3 = NULL, sum = FALSE) {
  
  stopifnot(type %in% c('BMA', 'PB'))
  
  if (invalidParams(gene1, gene2, gene3, sum)) {
    stop("invalid params")
  }
  
  
  if (type == 'BMA') {
    BEAT <- readRDS('data/BEAT_BMA.rds')  
  } else {
    BEAT <- readRDS('data/BEAT_PB.rds')  
  }
  
  x1 <- grep_gene(BEAT, gene1)
  x2 <- grep_gene(BEAT, gene2)
  x3 <- grep_gene(BEAT, gene3)

  if (sum) {
    BEAT$X <- sum_genes(x1,x2)
  } else {
    BEAT$X <- subtract_genes(x1,x2,x3, gene1, gene2, gene3)
  }
  BEAT
}


get_GSE37642 <- function(ds_num, gene1, gene2 = NULL, gene3 = NULL, sum = FALSE) {
  
  stopifnot(ds_num%in%1:2)
  
  if (invalidParams(gene1, gene2, gene3, sum)) {
    stop("invalid params")
  }
  
  
  if (ds_num == 1) {
    ds <- readRDS('data/GSE37642_1.rds')
    pl <- 'GPL570'
  } else {
    ds <- readRDS('data/GSE37642_2.rds')
    pl <- 'GPL96'
  }
  
  x1 <- getGeneFromPlatform(ds$X, gene1, pl)
  x2 <- NULL
  x3 <- NULL
  
  if (!is.null(gene2)) {
    x2 <- getGeneFromPlatform(ds$X, gene2, pl)
  }

    if (!is.null(gene3)) {
    x3 <- getGeneFromPlatform(ds$X, gene3, pl)
  }
  
  if (sum) {
    ds$X <- sum_genes(x1,x2)
  } else {
    ds$X <- subtract_genes(x1,x2,x3, gene1, gene2,gene3)
  }
  ds
}

## load('data/RData/GPL10558.RData', envir = e1)

get_GSE71014 <- function(gene1, gene2 = NULL, gene3 = NULL, sum = FALSE) {
  
  pl <- 'GPL10558'
  
  if (invalidParams(gene1, gene2, gene3, sum)) {
    stop("invalid params")
  }
  
  ds <- readRDS('data/GSE71014.rds')
  
  x1 <- getGeneFromPlatform(ds$X, gene1, pl)
  x2 <- NULL
  x3 <- NULL
  
  if (!is.null(gene2)) {
    x2 <- getGeneFromPlatform(ds$X, gene2, pl)
  }
  
  if (!is.null(gene3)) {
    x3 <- getGeneFromPlatform(ds$X, gene3, pl)
  }
  
  if (sum) {
    ds$X <- sum_genes(x1,x2)
  } else {
    ds$X <- subtract_genes(x1,x2, x3, gene1, gene2,gene3)
  }
  ds
}


get_GSE6891 <- function(ds_num, gene1, gene2 = NULL, gene3 = NULL, sum = FALSE) {
  
  stopifnot(ds_num%in%1:2)
  
  if (invalidParams(gene1, gene2, gene3, sum)) {
    stop("invalid params")
  }
  
  
  pl <- 'GPL570'
  
  if (ds_num == 1) {
    ds <- readRDS('data/GSE6891_1.rds')
  } else {
    ds <- readRDS('data/GSE6891_2.rds')
  }
  
  x1 <- getGeneFromPlatform(ds$X, gene1, pl)
  x2 <- NULL
  x3 <- NULL
  
  if (!is.null(gene2)) {
    x2 <- getGeneFromPlatform(ds$X, gene2, pl)
  }
  
  if (!is.null(gene3)) {
    x3 <- getGeneFromPlatform(ds$X, gene3, pl)
  }
  
  
  ds$X <- subtract_genes(x1,x2, x3, gene1, gene2, gene3)
  ds
  
}

get_genes_from_platform <- function(id, pl) {
  m <- match(id, pl$ID)
  genes <- pl[m,2]
  genes <- unlist(strsplit(genes, ' /// '))
  unique(genes)
}


get_all_genes <- function(ds) {
  if (ds == 'TCGA') {
    TCGA <- readRDS('data/TCGA.rds')
    r <- rownames(TCGA$X)
    r <- gsub('\\|.*$', '', r)
    return(r)
  } else if (ds == 'TARGET') {
    TARGET <- readRDS('data/target.rds')
    PL <- readRDS('data/ensemble.rds')
    genes <- get_genes_from_platform(rownames(TARGET$X), PL)
    return(genes)
  } else if (ds == 'BEAT_BMA') {
    BEAT <- readRDS('data/BEAT_BMA.rds')  
    return(rownames(BEAT$X))
  } else if (ds == 'BEAT_PB') {
    BEAT <- readRDS('data/BEAT_PB.rds')  
    return(rownames(BEAT$X))
  } else if (ds == 'GSE37642_1') {
    ds <- readRDS('data/GSE37642_1.rds')
    load('data/RData/GPL570.RData')
    genes <- get_genes_from_platform(rownames(ds$X), GPL570)
    return(genes)
  } else if (ds == 'GSE37642_2') {
    ds <- readRDS('data/GSE37642_2.rds')
    load('data/RData/GPL96.RData')
    genes <- get_genes_from_platform(rownames(ds$X), GPL96)
    return(genes)
  }else if (ds == 'GSE6891_1') {
    ds <- readRDS('data/GSE6891_1.rds')
    load('data/RData/GPL570.RData')
    genes <- get_genes_from_platform(rownames(ds$X), GPL570)
    return(genes)
  }else if (ds == 'GSE6891_2') {
    ds <- readRDS('data/GSE6891_2.rds')
    load('data/RData/GPL570.RData')
    genes <- get_genes_from_platform(rownames(ds$X), GPL570)
    return(genes)
  }else if (ds == 'GSE71014') {
    ds <- readRDS('data/GSE71014.rds')
    load('data/RData/GPL10558.RData')
    genes <- get_genes_from_platform(rownames(ds$X), GPL10558)
    return(genes)
  }else  {
    stop('ds not recognized: ', ds)
  }
}



get_data_lsc17 <- function(RNA_SEQ_only = FALSE,
                           SPLIT_DENOVO = FALSE,
                           NO_RECURRENT = FALSE,
                           wlab ='weights2') {

  ## LSC17 genes and weights from paper
  genes <- data.frame(gene = c('DNMT3B','ZBTB46', 'NYNRIN',
                               'ARHGAP22', 'LAPTM4B', 'MMRN1',
                               'DPYSL3','KIAA0125', 'CDK6', 'CPXM1',
                               'SOCS2','SMIM24','EMP1', 'NGFRAP1',
                               'CD34', 'AKR1C3', 'GPR56'),
                      weights = c(0.0874, -0.0347,0.00865, -0.0138,
                                  0.00582,0.0258, 0.0284, 0.0196,
                                  -0.0704, -0.0258, 0.0271, -0.0226,
                                  0.0146, 0.0465, 0.0338, -0.0402,
                                  0.0501)
  )
  
  genes$weights2 <- sign(genes$weights)
  
  # get genes
  lsc17 <- lapply(genes$gene, get_all_data, RNA_SEQ_only = RNA_SEQ_only, 
                  NO_RECURRENT = NO_RECURRENT, SPLIT_DENOVO = SPLIT_DENOVO)
  
  ds <- names(lsc17[[1]])
  
  extract_genes <- function(myds, lsc17) {
    lapply(lsc17, function(L) L[[myds]]$X)
  }
  
  calc_score_lsc17 <- function(ds, genes, lsc17, wlab) {
    target <- extract_genes(ds, lsc17)
    tt <- do.call(rbind, target)
    r <- rownames(tt)
    r <- gsub('\\|.*', '', r)
    m <- match(r, genes$gene)
    w <- genes[[wlab]][m]
    colMeans(tt*w)
  }
  
  # calculate lsc17 score for each dataset
  mm <- lapply(ds, calc_score_lsc17, genes, lsc17, wlab)
  names(mm) <- ds
  
  # create final list with clinical information
  final <- list()
  for (ds1 in ds) {
    final[[ds1]] <- lsc17[[1]][[ds1]]
    #final[[ds1]]$X <- t(as.matrix(mm[[ds1]]))
    #rownames[[ds1]]$X <- 'LSC17'
  }
  final
}

