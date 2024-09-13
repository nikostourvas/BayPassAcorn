# Snakefile

configfile: "config.yaml"

#Define the software container
singularity: config["singularity"]

# Get list of sample names
SAMPLES = config["samples"]

# Define number of subs for splitting the baypass input file
N_SUBS=6
# Define number of replicates for the core model
N_CORE_REPLICATES=3
# Define the minimum haploid pool size
MIN_HAPLOID_POOL_SIZE=15
N_PILOT=100 # see model details below

rule all:
    input:
        expand("results/{sample}_baypassSplitOut_core/core_{i}_mat_omega.out", sample=SAMPLES,
            i=range(1, N_CORE_REPLICATES+1)),
        expand("results/{sample}_omega_comp.pdf", sample=SAMPLES),
        expand("results/{sample}_omega_comp.csv", sample=SAMPLES),
        expand("results/{sample}_baypassSplitOut_covariate/covariate_{i}_summary_betai_reg.out", 
            sample=SAMPLES, i=range(1, N_SUBS+1)),
        expand("results/{sample}_std_IS_model_diagnostics.pdf", sample=SAMPLES),
        expand("results/{sample}_xtxst_pvalue_dist.pdf", sample=SAMPLES),
        expand("results/{sample}_concatenated_res_covariate.csv", sample=SAMPLES),
    threads: 1
    resources:
        runtime=10,
        mem_mb=1000,
        slurm_partition="testing"

rule generate_complementary_inputs:
    input:
        envfactors="data/20240805_ACORN_dem_TOPO_by_POP_NT.csv",
        poolsizes="data/20240307_poolsizes_per_pop.csv",
    resources:
        runtime=5,
        mem_mb=1000,
        slurm_partition="batch"
    output:
        poolsizes="data/{sample}_poolsizes",
        efile="data/{sample}_efile",
        ecotype="data/{sample}_ecotype",
    script: "generate_complementary_inputs.R"

rule vcf2genobaypass:
    input:
        vcf="data/{sample}.vcf",
        poolsizes="data/{sample}_poolsizes",
        poolnames="data/{sample}_poolnames",
    params:
        prefix="{sample}",
        subs=N_SUBS,
    threads: 1
    resources:
        runtime=90,
        mem_mb=24000,
        slurm_partition="batch"
    output:
        genobaypass=["results/subsets/{{sample}}.genobaypass.sub{}".format(i) for i in range(1, N_SUBS+1)],
        snpdet=["results/subsets/{{sample}}.snpdet.sub{}".format(i) for i in range(1, N_SUBS+1)]
    script: "poolfstat_subsample.R"

#########################################
## Detecting outlier loci with baypass ##
#########################################

## Notes: 
## (1) This part of the code is run directly from the command line in the Terminal.
## (2)-d0yij = 1/5 of min haploid pool size; use npilot = 100 for final analysis (see the Manual for detailed parameter descriptions).
## (3) baypass is not scaling linearly with the no of threads; analyze subsets of data on single threads in parallel; 
## pooldata.subset() function in poolfstat can be used.
## (4) Option 3: Contrast analysis is applied on the subsetted ACORN data!

## Option 1. Scanning the genome for differentiation (without covariate) 
# running core model for scanning the genome for differentiation using the XtX statistics
rule run_baypass_core:
    input:
        sub="results/subsets/{sample}.genobaypass.sub{i}"
    params:
        threads=1,
        poolsizefile="data/{sample}_poolsizes",
        npop=lambda wildcards: config['sample_metadata'][wildcards.sample]['npop'],
        d0yij=MIN_HAPLOID_POOL_SIZE/5,
        npilot=N_PILOT,
    threads: 1
    resources:
        runtime=1440,
        mem_mb=1000,
        slurm_partition="batch"
    output:
        mat_omega = protected("results/{sample}_baypassSplitOut_core/core_{i}_mat_omega.out"),
        summary_lda_omega = protected("results/{sample}_baypassSplitOut_core/core_{i}_summary_lda_omega.out"),
        xtx_summary= protected("results/{sample}_baypassSplitOut_core/core_{i}_summary_pi_xtx.out"),
    shell:
        """
        g_baypass \
        -nthreads {params.threads} \
        -npop {params.npop} -gfile {input.sub} -poolsizefile {params.poolsizefile} \
        -d0yij {params.d0yij} -npilot {params.npilot} \
        -outprefix results/{wildcards.sample}_baypassSplitOut_core/core_{wildcards.i}
        """

rule compare_omega:
    input:
        omega1="results/{sample}_baypassSplitOut_core/core_1_mat_omega.out",
        omega2="results/{sample}_baypassSplitOut_core/core_2_mat_omega.out",
        omega3="results/{sample}_baypassSplitOut_core/core_3_mat_omega.out",
        xtx_summary="results/{sample}_baypassSplitOut_core/core_1_summary_pi_xtx.out"
    threads: 1
    resources:
        runtime=30,
        mem_mb=2000,
        slurm_partition="batch"
    output:
        omega_comp="results/{sample}_omega_comp.pdf",
        omega_comp_table="results/{sample}_omega_comp.csv",
        xtxst_pvalue_dist="results/{sample}_xtxst_pvalue_dist.pdf"
    script: "compare_omega.R"

## Option 2. Identifying SNPs associated with population covariate data
rule run_baypass_covariate:
    input:
        sub="results/subsets/{sample}.genobaypass.sub{i}",
        omegafile="results/{sample}_baypassSplitOut_core/core_1_mat_omega.out",
    params:
        threads=1,
        poolsizefile="data/{sample}_poolsizes",
        npop=lambda wildcards: config['sample_metadata'][wildcards.sample]['npop'],
        efile="data/{sample}_efile",
        d0yij=MIN_HAPLOID_POOL_SIZE/5,
        npilot=N_PILOT,
    threads: 1
    resources:
        runtime=12960,
        mem_mb=2000,
        slurm_partition="htc"
    output:
        protected("results/{sample}_baypassSplitOut_covariate/covariate_{i}_covariate.std"),
        protected("results/{sample}_baypassSplitOut_covariate/covariate_{i}_DIC.out"),
        protected("results/{sample}_baypassSplitOut_covariate/covariate_{i}_summary_beta_params.out"),
        protected("results/{sample}_baypassSplitOut_covariate/covariate_{i}_summary_betai_reg.out"),
        protected("results/{sample}_baypassSplitOut_covariate/covariate_{i}_summary_pi_xtx.out"),
        protected("results/{sample}_baypassSplitOut_covariate/covariate_{i}_summary_yij_pij.out")
    shell:    
        """
        g_baypass \
        -nthreads {params.threads} \
        -npop {params.npop} -gfile {input.sub} -poolsizefile {params.poolsizefile} \
        -d0yij {params.d0yij} -npilot {params.npilot} \
        -omegafile {input.omegafile} \
        -efile {params.efile} \
        -outprefix results/{wildcards.sample}_baypassSplitOut_covariate/covariate_{wildcards.i}
        """

rule baypass_covariate_diagnostics:
    input:
        summary_betai_reg = "results/{sample}_baypassSplitOut_covariate/covariate_1_summary_betai_reg.out",
    threads: 1
    resources:
        runtime=90,
        mem_mb=8000,
        slurm_partition="batch"
    output:
        diagnostics_is = "results/{sample}_std_IS_model_diagnostics.pdf",
    script: "model_diagnostics_covariate.R"

## Option 3. Running contrast analysis estimating C2 statistic: 
## population ecotype is a binary trait (dry = -1; moist = 1)
rule run_baypass_C2:
    input:
        sub="results/subsets/{sample}.genobaypass.sub{i}",
        omegafile="results/{sample}_baypassSplitOut_core/core_1_mat_omega.out",
    params:
        threads=1,
        poolsizefile="data/{sample}_poolsizes",
        npop=lambda wildcards: config['sample_metadata'][wildcards.sample]['npop'],
        ecotype="data/{sample}_ecotype",
        d0yij=MIN_HAPLOID_POOL_SIZE/5,
        npilot=N_PILOT,
    threads: 1
    resources:
        runtime=2880,
        mem_mb=1000,
        slurm_partition="htc"
    output:
        protected("results/{sample}_baypassSplitOut_contrast/contrast_{i}_covariate.std"),
        protected("results/{sample}_baypassSplitOut_contrast/contrast_{i}_DIC.out"),
        protected("results/{sample}_baypassSplitOut_contrast/contrast_{i}_summary_beta_params.out"),
        protected("results/{sample}_baypassSplitOut_contrast/contrast_{i}_summary_betai_reg.out"),
        protected("results/{sample}_baypassSplitOut_contrast/contrast_{i}_summary_contrast.out"),
        protected("results/{sample}_baypassSplitOut_contrast/contrast_{i}_summary_pi_xtx.out"),
        protected("results/{sample}_baypassSplitOut_contrast/contrast_{i}_summary_yij_pij.out")
    shell:
        """
        g_baypass \
        -nthreads {params.threads} \
        -npop {params.npop} -gfile {input.sub} -poolsizefile {params.poolsizefile} \
        -d0yij {params.d0yij} -npilot {params.npilot} \
        -omegafile {input.omegafile} \
        -contrastfile {params.ecotype} -efile {params.ecotype} \
        -outprefix results/{wildcards.sample}_baypassSplitOut_contrast/contrast_{wildcards.i}
        """

rule c2_diagnostics:
    input:
        summary_contrast = "results/{sample}_baypassSplitOut_contrast/contrast_1_summary_contrast.out",
    params:
        ecotype="data/ecotype",
    threads: 1
    resources:
        runtime=30,
        mem_mb=8000,
        slurm_partition="batch"
    output:
        diagnostics_c2 = "results/{sample}_C2_model_diagnostics.pdf",
    script: "model_diagnostics_c2.R"

rule concatenate_results_covariate:
    input:
        covariate = expand("results/{{sample}}_baypassSplitOut_covariate/covariate_{i}_summary_betai_reg.out", 
        i=range(1, N_SUBS+1)),
    params:
        prefixcovariate = "results/{sample}_baypassSplitOut_covariate/covariate",
        subs=N_SUBS,
        snpdetprefix = "results/subsets/{sample}.snpdet.sub",
        retrieve_c2= False
    threads: 1
    resources:
        runtime=30,
        mem_mb=8000,
        slurm_partition="batch"        
    output:
        covariateresults = "results/{sample}_concatenated_res_covariate.csv",
    script: "concatenate_res.R"

rule concatenate_results_c2:
    input:
        contrast = expand("results/{{sample}}_baypassSplitOut_contrast/contrast_{i}_summary_betai_reg.out", 
         i=range(1, N_SUBS+1))
    params:
        prefixcontrast = "results/{sample}_baypassSplitOut_contrast/contrast",
        subs=N_SUBS,
        snpdetprefix = "results/subsets/{sample}.snpdet.sub",
        retrieve_c2= True
    threads: 1
    resources:
        runtime=30,
        mem_mb=8000,
        slurm_partition="batch"        
    output:
        contrasttable = "results/{sample}_concatenated_res_contrast.csv"
    script: "concatenate_res.R"