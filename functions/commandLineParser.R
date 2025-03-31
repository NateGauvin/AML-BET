source('testDataFunctions.R')

parser <- ArgumentParser(
  description = paste0("Process data for mongo upload. 
                        Any data uploaded will overwrite existing data."))

parser$add_argument("--datasets", default = "all",
                    help = "comma-separated list of datasets to process.")
parser$add_argument("--expression", default = "yes",
                    help = "add expression data (yes/no)")
parser$add_argument("--clinical", default = "yes",
                    help = "add clinical data (yes/no)")
parser$add_argument("--db", default = "no",
                    help = "upload data to database (yes/no)")

args <- parser$parse_args()

for (p in c('expression', 'clinical', 'db')) {
  if (!args[[p]] %in% c('yes', 'no')) {
    stop('argument must be yes/no: ', p)
  }
}

all_datasets <- NULL
all_datasets <- Sys.glob('data/*.rds')
all_datasets <- sub('data/', '', all_datasets)
all_datasets <- sub('.rds', '', all_datasets)
if (args$datasets != "all") {
  datasets <- strsplit(args$datasets, ',')
  if (any(is.na(match(datasets[[1]], all_datasets)))) {
    stop('invalid datasets entered')
  }
  all_datasets <- datasets
}

if (args$expression == "yes") {
  apply(datasets,, expression_upload, db)
}

if (args$clinical == "yes") {
  apply(datasets,, upload_clinical_data, db)
}
