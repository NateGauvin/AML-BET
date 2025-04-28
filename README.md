[![all_datasets_no_db](https://github.com/NateGauvin/AML-BET/actions/workflows/process_upload.yml/badge.svg)](https://github.com/NateGauvin/AML-BET/actions/workflows/process_upload.yml)

# AML-BET (Acute Myeloid Leukemia Biomarker Evaluation Tool)

Prototype and thesis proposal can be found in AML-BET_Prototype.
Release versions of AML-BET will use source and modified code from leukemia_lib.


# leukemia_lib (Contributed by Dr. Garrett Dancik)

Functions and R data files for leukemia datasets used in the following publications:

- Dancik, G.M., Voutsas, I.F. & Vlahopoulos, S. Lower RNA expression of ALDH1A1 distinguishes the favorable risk group in acute myeloid leukemia. Mol Biol Rep 49, 3321–3331 (2022). https://doi.org/10.1007/s11033-021-07073-7 
- Dancik, G.M.; Varisli, L.; Tolan, V.; Vlahopoulos, S. Aldehyde Dehydrogenase Genes as Prospective Actionable Targets in Acute Myeloid Leukemia. Genes 2023, 14, 1807. https://doi.org/10.3390/genes14091807 


# Current Setup Instructions

At the moment, AML-BET allows you to populate a MongoDB instance with the complete expression and phenotype data from all currently supported datasets, as well as produce basic grouped boxplots. Instructions are as follows:

1. Download and install R from https://cloud.r-project.org/
2. Download and install RStudio from https://posit.co/download/rstudio-desktop/ (required to generate boxplot diagrams)
3. Download and install MongoDB from https://www.mongodb.com/docs/manual/installation/
4. Download and install MongoDB Compass from https://www.mongodb.com/try/download/compass (optional, allows you to view datasets)
5. Using RStudio or another text editor/IDE, navigate to functions/mongoConnection.R and add hashtags (#) to comment out the detailed connect_mongo() function.
6. Remove the respective hashtags that comment out the simpler connect_mongo() function.
7. Using the RStudio terminal or cmd, run the command "Rscript functions/commandLineParser --db yes" to populate your MongoDB instance with the expression and clinical data of all currently supported datasets.
8. Once populated, you can use the script functions/queryMongoGeneratePlots.R in RStudio to generate plots from the data.
