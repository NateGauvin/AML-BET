[![all_datasets_no_db](https://github.com/NateGauvin/AML-BET/actions/workflows/process_upload.yml/badge.svg)](https://github.com/NateGauvin/AML-BET/actions/workflows/process_upload.yml)

# AML-BET (Acute Myeloid Leukemia Biomarker Evaluation Tool)

*AML-BET* is an online tool that allows for the evaluation of biomarkers within the RNA-Seq profiles of Acute Myeloid Leukemia patients. You can perform visual differential expression analysis within risk prognoses, generate Kaplan Meier curves between high and low expression groups, and view forest plot summarizations of survival data across all datasets.

The prototype and thesis proposal for *AML-BET* can be found in AML-BET_Prototype.

# leukemia_lib (Contributed by Dr. Garrett Dancik)

Functions and R data files for leukemia datasets used in the following publications:

- Dancik, G.M., Voutsas, I.F. & Vlahopoulos, S. Lower RNA expression of ALDH1A1 distinguishes the favorable risk group in acute myeloid leukemia. Mol Biol Rep 49, 3321–3331 (2022). https://doi.org/10.1007/s11033-021-07073-7 
- Dancik, G.M.; Varisli, L.; Tolan, V.; Vlahopoulos, S. Aldehyde Dehydrogenase Genes as Prospective Actionable Targets in Acute Myeloid Leukemia. Genes 2023, 14, 1807. https://doi.org/10.3390/genes14091807 


# Running AML-BET

The web version of AML-BET has yet to launch, but AML-BET can be run locally using Docker.

1. First, download and install docker from https://docs.docker.com/get-docker/.

2. Save the docker compose script, available within this repository as docker/docker-compose.yml.

3. From your terminal, run the following from the directory containing the *docker-compose.yml* file:

    ```
    docker compose up -d
    ```

    Initial startup

4. You can now access *AML-BET* by opening your browser and entering the following:

    ```
    http://localhost:3838
    ```

    It may take a while for the database to fully load.

5. To stop *AML-BET*, you can run the following docker command from the directory containing the *docker-compose.yml* file:

    ```
    docker compose down
    ```
