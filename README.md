[![all_datasets_no_db](https://github.com/NateGauvin/AML-BET/actions/workflows/process_upload.yml/badge.svg)](https://github.com/NateGauvin/AML-BET/actions/workflows/process_upload.yml)

# AML-BET (Acute Myeloid Leukemia Biomarker Evaluation Tool)

Prototype and thesis proposal can be found in AML-BET_Prototype.
Release versions of AML-BET will use source and modified code from leukemia_lib.


# leukemia_lib (Contributed by Dr. Garrett Dancik)

Functions and R data files for leukemia datasets used in the following publications:

- Dancik, G.M., Voutsas, I.F. & Vlahopoulos, S. Lower RNA expression of ALDH1A1 distinguishes the favorable risk group in acute myeloid leukemia. Mol Biol Rep 49, 3321–3331 (2022). https://doi.org/10.1007/s11033-021-07073-7 
- Dancik, G.M.; Varisli, L.; Tolan, V.; Vlahopoulos, S. Aldehyde Dehydrogenase Genes as Prospective Actionable Targets in Acute Myeloid Leukemia. Genes 2023, 14, 1807. https://doi.org/10.3390/genes14091807 


# Setup Instructions

Once R, RStudio, and MongoDB are installed on your machine and the repository is downloaded, navigate to file mongoConnection.R and comment out the connect_mongo() function and remove the comments for the connect_mongo() function immediately below it.
Then, from the Terminal, run "Rscript functions/commandLineParser --db yes" to process and upload all data to MongoDB. Finally, retrieval of the data can be done with Mongosh commands or within queryMongoGeneratePlots.R.
