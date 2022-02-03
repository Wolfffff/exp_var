# Expression Variation

This repository contains code for exploring the landscape of variation across a multitude of studies. Specifically, we have a pipeline for pulling a multitude of RNA-Seq datasets, processing them, and we're starting to poke at a results.

## Pipeline and Data

For the exact pipeline we run, see the `snakemake` folder in the repository. The current set of studies from (Expression Atlas)[snakemake/metadata/EA_metadata.csv](snakemake/metadata/EA_metadata.csv) and [snakemake/metadata/recount3_metadata.csv](snakemake/metadata/recount3_metadata.csv). We currently have ~70 studies. The resulting data is placed in [snakemake/Rdatas/residuals/](snakemake/Rdatas/residuals/). After processing, we merge everything in [snakemake_processing.R](snakemake_processing.R).

## Plots
The plots are generated in [data/plots/](data/plots/). Check out the PCAs and mean-variance plots for the residuals that we've discussed in the past.

## Ranking

Implemented PCA scoring/ranking in [acrossStudiesVarScore.R](acrossStudiesVarScore.R) which yield a more as we discussed previously. The output is currently the best thing we have to run with and should be give us a decent first pass at gene lists for toying with.

## Evo

Here we've take a pass at getting a couple evolutionary metrics. We pulled $\omega$ and $\alpha$ from [PopHuman](https://academic.oup.com/nar/article/46/D1/D1003/4559406) with chimpanzees as the reference species. That can be found in [imkt_pull.R](imkt_pull.R) and a correlation correlation plot is generated in [evo_analysis.R](evo_analysis.R) to be found [data/plots/SpearmanCorrelations/corr_plot_with_pvals.png](data/plots/SpearmanCorrelations/corr_plot_with_pvals.png).


## Disease

Starting PTWAS analysis -- see [ptwas.R](ptwas.R). We have the table but need to deal with subsetting to disease phenotypes using their Table S1 from the paper.

## Other

It's currently a bit of a pain to pull new data as we have to find the right features to subset to (case in case control) and other features like multiple timepoints or technical replicates.


https://www.nature.com/articles/s41467-019-09583-2 has proved harder to use than expected and I'd appreciate guidance on the right way to use their data once you're back.

## Drivers of Gene Expression Variance difference

- Positive selection, continuously pushing advantageous alleles to fixation
    - dn/ds
- Purifying selection, removing deleterious variants
    - dn/ds, segregating variants, 
- Differential input of mutation, due to differential mutation rate or differences in the size of the mutational targets. 
    - Number of eQTLs, regulatory regions
- Canalization, modifiers that change the sensitivity of genes to genetic and environmental variation

Early-expressed genes should then be less variable, and late-expressed genes more variable due to the additive effect of several upstream regulators.

Is there any gene-level measure that has a more-upstream to less-upstream axis? (KEGG?)

### Canalization

Genes should have some expected level of phenotypic variation given their underlying genetic variation and eQTLs and all that. Genes that have less phenotypic variation than expected by this measure can be said to be more canalized, being robust to the input of genetic variation. What are the characteristics of these canalized genes?

Same thing for hyper-variable genes. Are they responding to the environment?

Pitch: We can test for substantial canalization by looking at the relation between genetic variation (broad sense: genetic diversity, segregating sites, eQTLs) and phenotypic variation in the expression data. If there is a relation between them, no evidence for canalization and phenotypic variation is proportional to genetic variation. Genes that deviate from this might be more or less canalized, and we can probe them for enrichments.
