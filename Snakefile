# Snakefile

#Define the software container
singularity: "~/other_apps/bin/poolseq_tools_0.2.11.sif" # this path is on my system, you need to change it to your path

# Get list of sample names
with open("data/sample_list.txt") as f:
     SAMPLES = f.read().splitlines()

# Define number of chunks for splitting the baypass input file
N_CHUNKS=20

rule all:
    input:
        expand("results/{sample}_baypassSplitOut/chunk_{i}_DIC.out", sample=SAMPLES,
            i=range(1, N_CHUNKS))

rule vcf2counts:
    input:
        vcf="data/{sample}.vcf"
    output:
        counts=temp("results/{sample}_counts.txt")
    log:
        vcf2counts="logs/vcf2counts/{sample}_vcf2counts.log"
    shell:
        "bash vcf_to_counts.sh {input.vcf} > {output.counts} 2> {log.vcf2counts}"

# pool name fix - remove suffix from pool names
# If you do not to fix the pool names, you can skip this rule 
rule pool_name_fix:
    input:
        counts="results/{sample}_counts.txt"
    output:
        counts_fixed="results/{sample}_counts_fixed.txt"
    shell:
        "bash pool_name_fix.sh {input.counts} > {output.counts_fixed}"

# counts to baypass input file
# This is still a (computationally) costly trip to R, 
# although it is substantially more efficient than loading the VCF file
# to the poolfstat R package. 
# TODO: rewrite this with unix tools
rule counts2baypass:
    input:
        counts="results/{sample}_counts_fixed.txt"
    output:
        baypass="results/{sample}_geno_baypass.txt"
    script: "counts2baypass.R"

# generate subsets of the baypass input file
rule sample_baypass:
    input:
        "results/{sample}_geno_baypass.txt"
    output:
        ["results/{{sample}}_geno_baypass_{}.txt".format(i) for i in range(1, N_CHUNKS)]
    shell:
        "for i in {{1..{N_CHUNKS}}}; do awk -v i=$i 'NR % {N_CHUNKS} == i' {input} > results/{wildcards.sample}_geno_baypass_$i.txt; done"

#########################################
## Detecting outlier loci with baypass ##
#########################################

## Notes: 
## (1) This part of the code is run directly from the command line in the Terminal.
## (2)-d0yij = 1/5 of min haploid pool size; use npilot = 100 for final analysis (see the Manual for detailed parameter descriptions).
## (3) baypass is not scaling linearly with the no of threads; analyze subsets of data on single threads in parallel; pooldata.subset() function in poolfstat can be used.
## (4) Option 3: Contrast analysis is applied on the subsetted ACORN data!

## Option 1. Scanning the genome for differentiation (without covariate) 
# running core model for scanning the genome for differentiation using the XtX statistics
rule run_baypass_1:
    input:
        chunk="results/{sample}_geno_baypass_{i}.txt",
    params:
        threads=1,
        poolsizefile="data/poolsize",
        npop=40,
        d0yij=3,
        npilot=20
    output:
        dic = "results/{sample}_baypassSplitOut/chunk_{i}_DIC.out",
        mat_omega = "results/{sample}_baypassSplitOut/chunk_{i}_mat_omega.out",
        summary_beta_params = "results/{sample}_baypassSplitOut/chunk_{i}_summary_beta_params.out",
        summary_lda_omega = "results/{sample}_baypassSplitOut/chunk_{i}_summary_lda_omega.out",
        summary_pi_xtx = "results/{sample}_baypassSplitOut/chunk_{i}_summary_pi_xtx.out",
        summary_yij_pij = "results/{sample}_baypassSplitOut/chunk_{i}_summary_yij_pij.out"
    shell:
        """
        g_baypass \
        -nthreads {params.threads} \
        -npop {params.npop} -gfile {input.chunk} -poolsizefile {params.poolsizefile} \
        -d0yij {params.d0yij} -npilot {params.npilot} \
        -outprefix results/{wildcards.sample}_baypassSplitOut/chunk_{wildcards.i}
        """