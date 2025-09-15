source('functions/processUploadFunctions.R')

#options(python_cmd = 'C:/Users/nateg/AppData/Local/Microsoft/WindowsApps/python.exe')


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
parser$add_argument("--genes", default = "yes",
                    help = "create collection for gene names (yes/no)")

args <- parser$parse_args()

for (p in c('expression', 'clinical', 'db')) {
  if (!args[[p]] %in% c('yes', 'no')) {
    stop('argument must be yes/no: ', p)
  }
}

if (args$db == "yes") {args$db <- TRUE}
if (args$db == "no") {args$db <- FALSE}

all_datasets <- NULL
all_datasets <- Sys.glob('data/*.rds')
all_datasets <- sub('data/', '', all_datasets)
all_datasets <- sub('.rds', '', all_datasets)

if (is.null(all_datasets)) {stop("No datasets found to process/upload.")}

if (args$datasets != "all") {
  datasets <- strsplit(args$datasets, ',')
  if (any(is.na(match(datasets[[1]], all_datasets)))) {
    stop('invalid datasets entered')
  }
  all_datasets <- datasets[[1]]
}

if (args$expression == "yes") {
  print(lapply(all_datasets, expression_upload, args$db))
}

if (args$clinical == "yes") {
  print(lapply(all_datasets, upload_clinical_data, args$db))
}

if (args$genes == "yes") {
  compile_genes()
}
