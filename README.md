# Characterizing the landscape of gene expression variance in humans

Repository contains scripts to reproduce results of the paper as below:

> *Characterizing the landscape of gene expression variance in humans.*
> Scott Wolf*, Diogo Melo*, Kristina Garske, Luisa Pallares, Julien Ayroles. 

# Workflow

The main pipeline for generating gene-level metrics is built using [Snakemake](https://snakemake.readthedocs.io/en/stable/). To run the pipeline, run the following command from the `snakemake` directory:

```{bash}
snakemake Snakefile
```


Many of the supporting functions can be found in [functions.R](functions.R).
## Snakemake pipeline

The following notes describe each piece of the Snakemake workflow in detail. The Snakemake file can be viewed at [snakemake/Snakefile](snakemake/Snakefile). The corresponding config file detailing the datasources is located at [snakemake/config.yaml](snakemake/config.yaml), and the config for running this on a Slurm cluster is located at [snakemake/config_slurm.yaml](snakemake/config_slurm.yaml).

### Data download

After reading each set of ids and appending the corresponding source, we download the data from the source. This is done by passing our study ids along with the corresponding source (Expression Atlas or recount3) and metadata for each. For sources found in Expression Atlas and recount3, we use the corresponding R packages (ExpressionAtlas and recount3) to pull the rawdata and corresponding metadata. For each of these, we output the data in `/snakemake/Rdatas/raw/{id}.rds`. The source for doing this can be found in: [download_datasets.R](snakemake/scripts/download_datasets.R).

### Preprocessing

After gathering the raw data, we need to preprocess the data. We initially passing the id and metadata. After loading the datasets, we preprocess the data by:

1. Filtering out manually curated columns of problematic metadata. The columns that were hand curated can be found in [snakemake/metadata/EA_metadata.csv](snakemake/metadata/EA_metadata.csv) and [snakemake/metadata/recount3_metadata.csv](snakemake/metadata/recount3_metadata.csv) accordingly.

2. Removing portions of studies where cell lines are used (i.e. Cultured Fibroblasts in the Skin tissue set from GTEx) and where multiple tissue types are pooled (i.e. Esophagus from GTEx).

3. Filtering to get only control samples from each study.

4. Summing technical replicates.

5. Drop gene with expression less than 1 cpm.

6. Drop genes with mean expression below 5 cpm

7. Drop mitochondrial genes

8. Remove extreme outliers from BLOOD and Stomach samples

9. Recalculate normalization factors on the remaining data

10. Remove redundant metadata

11. Remove factors with large number of levels



Following this, for each dataset, we output the filtered expression data and its corresponding filtered metadata at `/snakemake/Rdatas/preProcess/{id}.rds`. The source for doing this can be found in: [preprocessing.R](snakemake/scripts/preprocessing.R).

### Calculating residuals

Following preprocessing, we want to calculate residual expression after removing the confounding factors. We do this for each study by:

1. Making a design matrix using preprocessed data and ignoring manually identified problematic columns.

2. Calculate a variance stabilizing transform as implemented in [DESeq2](https://rdrr.io/bioc/DESeq2/man/varianceStabilizingTransformation.html).

3. Filter outlier individuals using a robust Principal Component Analysis approach as described by Chen et al. 2020.

4. Recalculate design matrix after filtering for outliers and stabilizing variance.



### Postprocessing 

#### CSV residuals

Generates a CSV file containing the residuals for each study.

#### Graph
Creates a a fully connected gene-by-gene graph in which each edge is weighted by the Spearman correlation between gene expression for each study.

#### Network stats

#### Gene metrics

#### Variance Scores

Calculates a standard, cross-study gene expression variation using the post correction expression residuals.