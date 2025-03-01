# BayPassAcorn

A snakemake pipeline to run Genotype-Environment Association analysis (GEA) using the software BayPass (Gautier 2015) for Pool-Seq data. 

The pipeline:
- ensures close to linear multi-core performance for BayPass by splitting large-scale genomic data sets into user-defined smaller sets,
- can accomodate multiple environmental factors (each run as a separate univariate analysis),
- provides option to rank-transform environmental factors,
- optionally runs the window based Weighted-Z Analysis (WZA) (Booker et al. 2023) using the BayPass results,
- provides diagnostic plots (p-value distribution, QQplots, Manhattan plots) for each test employed via ggplot2.

All required software to run the pipeline is available via the [poolseq_tools](https://github.com/nikostourvas/poolseq_tools) Singularity/Apptainer container.



  
## Literature  
Booker TR, Yeaman S, Whiting JR, Whitlock MC (2024) The WZA : A window‐based method for characterizing genotype–environment associations. Molecular Ecology Resources 24:e13768. https://doi.org/10.1111/1755-0998.13768

Gautier M (2015) Genome-Wide Scan for Adaptive Divergence and Association with Population-Specific Covariates. Genetics 201:1555–1579. https://doi.org/10.1534/genetics.115.181453
