# Snakefile

# Define configuration file that lists samples and metadata
# That is the file you should edit prior to analysis
configfile: "config.yaml"

# Get the software container from config file
singularity: config["singularity"]

# Get list of sample names from config file
SAMPLES = config["samples"]

# Get input file names from config file
ENVFACTORS = config["input_files"]["envfactors"]
POOLSIZES = config["input_files"]["poolsizes"]

# Get parameters from config file
N_SUBS = config["parameters"]["n_subs"]
N_CORE_REPLICATES = config["parameters"]["n_core_replicates"]
MIN_HAPLOID_POOL_SIZE = config["parameters"]["min_haploid_pool_size"]
N_PILOT = config["parameters"]["n_pilot"]

# Get resources from config file
RESOURCES = config["resources"]

# function to calculate the number of populations
import os

def count_words_in_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    return len(content.split())


rule all:
    input:
        expand("results/{sample}_baypassSplitOut_core/core_{i}_mat_omega.out", sample=SAMPLES,
            i=range(1, N_CORE_REPLICATES+1)),
        expand("results/{sample}_omega_comp.pdf", sample=SAMPLES),
        expand("results/{sample}_omega_comp.csv", sample=SAMPLES),
        expand("results/{sample}_baypassSplitOut_covariate/covariate_{i}_summary_betai_reg.out", 
            sample=SAMPLES, i=range(1, N_SUBS+1)),
        expand("results/{sample}_std_IS_model_BFis.png", sample=SAMPLES),
        expand("results/{sample}_C2_model_diagnostics.pdf", sample=SAMPLES),
        expand("results/{sample}_xtxst_pvalue_dist.pdf", sample=SAMPLES),
        expand("results/{sample}_concatenated_res_covariate.csv", sample=SAMPLES),
        expand("results/{sample}_concatenated_res_contrast.csv", sample=SAMPLES)
    resources:
        runtime=RESOURCES["all"]["runtime"],
        mem_mb=RESOURCES["all"]["mem_mb"],
        slurm_partition=RESOURCES["all"]["slurm_partition"]

rule generate_complementary_inputs:
    input:
        envfactors=ENVFACTORS,
        poolsizes=POOLSIZES,
    resources:
        runtime=RESOURCES["generate_complementary_inputs"]["runtime"],
        mem_mb=RESOURCES["generate_complementary_inputs"]["mem_mb"],
        slurm_partition=RESOURCES["generate_complementary_inputs"]["slurm_partition"]
    output:
        poolsizes="data/{sample}_poolsizes",
        efile="data/{sample}_efile",
        ecotype="data/{sample}_ecotype",
    script: "scripts/generate_complementary_inputs.R"

rule vcf2genobaypass:
    input:
        vcf="data/{sample}.vcf",
        poolsizes="data/{sample}_poolsizes",
        poolnames="data/{sample}_poolnames",
    params:
        prefix="{sample}",
        subs=N_SUBS,
    resources:
        runtime=RESOURCES["vcf2genobaypass"]["runtime"],
        mem_mb=RESOURCES["vcf2genobaypass"]["mem_mb"],
        slurm_partition=RESOURCES["vcf2genobaypass"]["slurm_partition"]
    output:
        genobaypass=["results/subsets/{{sample}}.genobaypass.sub{}".format(i) for i in range(1, N_SUBS+1)],
        snpdet=["results/subsets/{{sample}}.snpdet.sub{}".format(i) for i in range(1, N_SUBS+1)]
    script: "scripts/poolfstat_subsample.R"

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
        npop=lambda wildcards: count_words_in_file(f"data/{wildcards.sample}_poolnames"),
        d0yij=MIN_HAPLOID_POOL_SIZE/5,
        npilot=N_PILOT,
    resources:
        runtime=RESOURCES["run_baypass_core"]["runtime"],
        mem_mb=RESOURCES["run_baypass_core"]["mem_mb"],
        slurm_partition=RESOURCES["run_baypass_core"]["slurm_partition"]
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
    resources:
        runtime=RESOURCES["compare_omega"]["runtime"],
        mem_mb=RESOURCES["compare_omega"]["mem_mb"],
        slurm_partition=RESOURCES["compare_omega"]["slurm_partition"]
    output:
        omega_comp="results/{sample}_omega_comp.pdf",
        omega_comp_table="results/{sample}_omega_comp.csv",
        xtxst_pvalue_dist="results/{sample}_xtxst_pvalue_dist.pdf"
    script: "scripts/compare_omega.R"

## Option 2. Identifying SNPs associated with population covariate data
rule run_baypass_covariate:
    input:
        sub="results/subsets/{sample}.genobaypass.sub{i}",
        omegafile="results/{sample}_baypassSplitOut_core/core_1_mat_omega.out",
    params:
        threads=1,
        poolsizefile="data/{sample}_poolsizes",
        npop=lambda wildcards: count_words_in_file(f"data/{wildcards.sample}_poolnames"),
        efile="data/{sample}_efile",
        d0yij=MIN_HAPLOID_POOL_SIZE/5,
        npilot=N_PILOT,
    resources:
        runtime=RESOURCES["run_baypass_covariate"]["runtime"],
        mem_mb=RESOURCES["run_baypass_covariate"]["mem_mb"],
        slurm_partition=RESOURCES["run_baypass_covariate"]["slurm_partition"]
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
        envfactor_names = "data/{sample}_efile_envfactor_names"
    resources:
        runtime=RESOURCES["baypass_covariate_diagnostics"]["runtime"],
        mem_mb=RESOURCES["baypass_covariate_diagnostics"]["mem_mb"],
        slurm_partition=RESOURCES["baypass_covariate_diagnostics"]["slurm_partition"]
    output:
        BFis = "results/{sample}_std_IS_model_BFis.png",
        eBPis = "results/{sample}_std_IS_model_eBPis.png",
        Beta_is = "results/{sample}_std_IS_model_Beta_is.png",
    script: "scripts/model_diagnostics_covariate.R"

## Option 3. Running contrast analysis estimating C2 statistic: 
## population ecotype is a binary trait (dry = -1; moist = 1)
rule run_baypass_c2:
    input:
        sub="results/subsets/{sample}.genobaypass.sub{i}",
        omegafile="results/{sample}_baypassSplitOut_core/core_1_mat_omega.out",
    params:
        threads=1,
        poolsizefile="data/{sample}_poolsizes",
        npop=lambda wildcards: count_words_in_file(f"data/{wildcards.sample}_poolnames"),
        ecotype="data/{sample}_ecotype",
        d0yij=MIN_HAPLOID_POOL_SIZE/5,
        npilot=N_PILOT,
    resources:
        runtime=RESOURCES["run_baypass_c2"]["runtime"],
        mem_mb=RESOURCES["run_baypass_c2"]["mem_mb"],
        slurm_partition=RESOURCES["run_baypass_c2"]["slurm_partition"]
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
        summary_betai_reg = "results/{sample}_baypassSplitOut_contrast/contrast_1_summary_betai_reg.out",
        summary_contrast = "results/{sample}_baypassSplitOut_contrast/contrast_1_summary_contrast.out",
    params:
        ecotype="data/ecotype",
    resources:
        runtime=RESOURCES["c2_diagnostics"]["runtime"],
        mem_mb=RESOURCES["c2_diagnostics"]["mem_mb"],
        slurm_partition=RESOURCES["c2_diagnostics"]["slurm_partition"]
    output:
        diagnostics_c2 = "results/{sample}_C2_model_diagnostics.pdf",
    script: "scripts/model_diagnostics_c2.R"

rule concatenate_results_covariate:
    input:
        envfactor_names = "data/{sample}_efile_envfactor_names",
        covariate = expand("results/{{sample}}_baypassSplitOut_covariate/covariate_{i}_summary_betai_reg.out", 
        i=range(1, N_SUBS+1)),
    params:
        prefixcovariate = "results/{sample}_baypassSplitOut_covariate/covariate",
        subs=N_SUBS,
        snpdetprefix = "results/subsets/{sample}.snpdet.sub",
        retrieve_c2= False
    resources:
        runtime=RESOURCES["concatenate_results_covariate"]["runtime"],
        mem_mb=RESOURCES["concatenate_results_covariate"]["mem_mb"],
        slurm_partition=RESOURCES["concatenate_results_covariate"]["slurm_partition"]       
    output:
        covariateresults = "results/{sample}_concatenated_res_covariate.csv",
        manhattanplot = "results/{sample}_manhattanplot_covariate.png"
    script: "scripts/concatenate_res.R"

rule concatenate_results_c2:
    input:
        contrast = expand("results/{{sample}}_baypassSplitOut_contrast/contrast_{i}_summary_betai_reg.out", 
         i=range(1, N_SUBS+1))
    params:
        prefixcontrast = "results/{sample}_baypassSplitOut_contrast/contrast",
        subs=N_SUBS,
        snpdetprefix = "results/subsets/{sample}.snpdet.sub",
        retrieve_c2= True
    resources:
        runtime=RESOURCES["concatenate_results_c2"]["runtime"],
        mem_mb=RESOURCES["concatenate_results_c2"]["mem_mb"],
        slurm_partition=RESOURCES["concatenate_results_c2"]["slurm_partition"]       
    output:
        contrasttable = "results/{sample}_concatenated_res_contrast.csv"
    script: "scripts/concatenate_res.R"