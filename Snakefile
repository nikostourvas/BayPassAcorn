# Snakefile

configfile: "config.yaml"

#Define the software container
singularity: "~/other_apps/bin/poolseq_tools_0.2.11.sif" # this path is on my system, you need to change it to your path

# Get list of sample names
with open("data/sample_list.txt") as f:
     SAMPLES = f.read().splitlines()

# Define number of subs for splitting the baypass input file
N_SUBS=12

rule all:
    input:
        expand("results/{sample}_concatenated_res_contrast.csv", sample=SAMPLES),
        expand("results/{sample}_baypassSplitOut_core/core_{i}_mat_omega.out", sample=SAMPLES,
            i=range(1,4)),
        "results/{sample}_omega_comp.pdf".format(sample=SAMPLES[0]),
        "results/{sample}_omega_comp.csv".format(sample=SAMPLES[0]),
        "results/{sample}_std_IS_model_diagnostics.pdf".format(sample=SAMPLES[0])

rule vcf2genobaypass:
    input:
        vcf="data/{sample}.vcf"
    params:
        poolsizes="data/{sample}_poolsizes",
        poolnames="data/{sample}_poolnames",
        prefix="{sample}",
        subs=N_SUBS,
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
        poolsizefile="results/subsets/{sample}.poolsize",
        npop=lambda wildcards: config['samples'][wildcards.sample]['npop'],
        d0yij=3,
        npilot=20
    output:
        mat_omega = "results/{sample}_baypassSplitOut_core/core_{i}_mat_omega.out",
        summary_lda_omega = "results/{sample}_baypassSplitOut_core/core_{i}_summary_lda_omega.out",
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
        omega3="results/{sample}_baypassSplitOut_core/core_3_mat_omega.out"
    output:
        omega_comp="results/{sample}_omega_comp.pdf",
        omega_comp_table="results/{sample}_omega_comp.csv"
    script: "compare_omega.R"

## Option 3. Running contrast analysis estimating C2 statistic: 
## population ecotype is a binary trait (dry = -1; moist = 1)
rule run_baypass_C2:
    input:
        sub="results/subsets/{sample}.genobaypass.sub{i}",
        omegafile="results/{sample}_baypassSplitOut_core/core_1_mat_omega.out",
    params:
        threads=1,
        poolsizefile="results/subsets/{sample}.poolsize",
        npop=lambda wildcards: config['samples'][wildcards.sample]['npop'],
        ecotype="data/ecotype",
        d0yij=3,
        npilot=20,
    output:
        covariate = "results/{sample}_baypassSplitOut_contrast/contrast_{i}_covariate.std",
        dic = "results/{sample}_baypassSplitOut_contrast/contrast_{i}_DIC.out",
        summary_beta_params = "results/{sample}_baypassSplitOut_contrast/contrast_{i}_summary_beta_params.out",
        summary_betai_reg = "results/{sample}_baypassSplitOut_contrast/contrast_{i}_summary_betai_reg.out",
        summary_contrast = "results/{sample}_baypassSplitOut_contrast/contrast_{i}_summary_contrast.out",
        summary_pi_xtx = "results/{sample}_baypassSplitOut_contrast/contrast_{i}_summary_pi_xtx.out",
        summary_yij_pij = "results/{sample}_baypassSplitOut_contrast/contrast_{i}_summary_yij_pij.out"
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

rule model_diagnostics:
    input:
        summary_betai_reg = "results/{sample}_baypassSplitOut_contrast/contrast_1_summary_betai_reg.out",
        summary_contrast = "results/{sample}_baypassSplitOut_contrast/contrast_1_summary_contrast.out",
    output:
        diagnostics_is = "results/{sample}_std_IS_model_diagnostics.pdf",
        diagnostics_c2 = "results/{sample}_C2_model_diagnostics.pdf"
    script: "model_diagnostics.R"

rule concatenate_results:
    input:
        contrast = expand("results/{{sample}}_baypassSplitOut_contrast/contrast_{i}_summary_betai_reg.out", 
        i=range(1, N_SUBS+1))
    params:
        prefixcontrast = "results/{sample}_baypassSplitOut_contrast/contrast",
        subs=N_SUBS,
        snpdetprefix = "results/subsets/{sample}.snpdet.sub"
    output:
        contrasttable = "results/{sample}_concatenated_res_contrast.csv"
    script: "concatenate_res.R"