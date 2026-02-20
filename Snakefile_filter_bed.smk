# Snakefile

# Define configuration file that lists samples and metadata
# That is the file you should edit prior to analysis
configfile: "config.yaml"

# Get the software container from config file
singularity: config["singularity"]

# Get list of sample names from config file
SAMPLES = config["samples"]

# Get list of environment factors from config file
ENVFACTOR_NAMES = config["envfactor_names"]

# Get input file names from config file
ENVFACTORS = config["input_files"]["envfactors"]
POOLSIZES = config["input_files"]["poolsizes"]

# Get parameters from config file
RANK_ANALYSIS_FLAG = config["parameters"]["rank_transform_covariates"]
N_SUBS = config["parameters"]["n_subs"]
N_CORE_REPLICATES = config["parameters"]["n_core_replicates"]
MIN_HAPLOID_POOL_SIZE = config["parameters"]["min_haploid_pool_size"]
N_PILOT = config["parameters"]["n_pilot"]
BF_THRESHOLD = config["parameters"]["bf_threshold"]
SPEARMAN_THRESHOLD = config["parameters"]["spearman_threshold"]
WZA_WINDOW_SIZE = config["parameters"]["WZA_window_size"]
WZA_MAF_THRESHOLD = config["parameters"]["WZA_maf_threshold"]
WZA_FDR = config["parameters"]["WZA_fdr"]

# Get resources from config file
RESOURCES = config["resources"]
localrules: all, generate_complementary_inputs, compare_omega

# function to calculate the number of populations
import os

def count_words_in_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    return len(content.split())

rule all:
    input:
        expand("results/{sample}_concatenated_res_covariate_filtered.csv", sample=SAMPLES),
        expand("results/{sample}_concatenated_res_covariate_spearman_filtered.csv", sample=SAMPLES),
        # expand("data/WZA/{sample}_WZA_input.csv", sample=SAMPLES),
        # expand("data/WZA/{sample}_WZA_input_filtered.csv", sample=SAMPLES),
        # expand("results/WZA_res/{sample}_{envfactor}_BF_WZA_output.csv", sample=SAMPLES, envfactor=ENVFACTOR_NAMES),
        # expand("results/WZA_res/{sample}_{envfactor}_spearman_WZA_output.csv", sample=SAMPLES, envfactor=ENVFACTOR_NAMES),
        # expand("results/WZA_res/{sample}_WZA_manhattan_plots_BF.png", sample=SAMPLES),
        # expand("results/WZA_res/{sample}_{envfactor}_WZA_output_fdr.csv", sample=SAMPLES, envfactor=ENVFACTOR_NAMES),



# The following 3 rules remove known Transposable Elements (TEs) from the BayPass results.
rule make_snps_bed:
    input:
        covariateresults = "results/{sample}_concatenated_res_covariate.csv",
        covariateresults_spearman = "results/{sample}_concatenated_res_covariate_spearman.csv",
    output:
        covariateresults_bed = "results/{sample}_concatenated_res_covariate.bed",
        covariateresults_spearman_bed = "results/{sample}_concatenated_res_covariate_spearman.bed",
    resources:
        runtime=RESOURCES["make_snps_bed"]["runtime"],
        mem_mb=RESOURCES["make_snps_bed"]["mem_mb"],
        slurm_partition=RESOURCES["make_snps_bed"]["slurm_partition"]
    shell:
        """
        awk 'BEGIN{{FS=","; OFS="\t"}} NR > 1 {{
            if ($1 ~ /^[0-9]+$/) {{
                chr=sprintf("Qrob_Chr%02d", $1); 
            }} else {{
                chr=$1;
                gsub(/\.Sc/, "_Sc", chr);
            }}
            print chr, $2-1, $2, $0
        }}' {input.covariateresults} > {output.covariateresults_bed}

        awk 'BEGIN{{FS=","; OFS="\t"}} NR > 1 {{
            if ($1 ~ /^[0-9]+$/) {{
                chr=sprintf("Qrob_Chr%02d", $1); 
            }} else {{
                chr=$1;
                gsub(/\.Sc/, "_Sc", chr);
            }}
            print chr, $2-1, $2, $0
        }}' {input.covariateresults_spearman} > {output.covariateresults_spearman_bed}
        """

rule intersect_snps_tes:
    input:
        covariateresults_bed = "results/{sample}_concatenated_res_covariate.bed",
        covariateresults_spearman_bed = "results/{sample}_concatenated_res_covariate_spearman.bed",
        te_bed = "data/tes.bed",
    output:
        filtered_bed = "results/{sample}_concatenated_res_covariate_filtered.bed",
        filtered_spearman_bed = "results/{sample}_concatenated_res_covariate_spearman_filtered.bed",
    resources:
        runtime=RESOURCES["intersect_snps_tes"]["runtime"],
        mem_mb=RESOURCES["intersect_snps_tes"]["mem_mb"],
        slurm_partition=RESOURCES["intersect_snps_tes"]["slurm_partition"]
    shell:
        """
        bedtools intersect -a {input.covariateresults_bed} -b {input.te_bed} -v > {output.filtered_bed}
        bedtools intersect -a {input.covariateresults_spearman_bed} -b {input.te_bed} -v > {output.filtered_spearman_bed}
        """

rule make_filtered_results_covariate:
    input:
        covariateresults = "results/{sample}_concatenated_res_covariate.csv",
        covariateresults_spearman = "results/{sample}_concatenated_res_covariate_spearman.csv",
        filtered_bed = "results/{sample}_concatenated_res_covariate_filtered.bed",
        filtered_spearman_bed = "results/{sample}_concatenated_res_covariate_spearman_filtered.bed",
    output:
        filtered_covariateresults = "results/{sample}_concatenated_res_covariate_filtered.csv",
        filtered_covariateresults_spearman = "results/{sample}_concatenated_res_covariate_spearman_filtered.csv",
    resources:
        runtime=RESOURCES["make_filtered_results_covariate"]["runtime"],
        mem_mb=RESOURCES["make_filtered_results_covariate"]["mem_mb"],
        slurm_partition=RESOURCES["make_filtered_results_covariate"]["slurm_partition"]
    shell:
        """
        head -n 1 {input.covariateresults} > {output.filtered_covariateresults}
        cut -f 4- {input.filtered_bed} >> {output.filtered_covariateresults}

        head -n 1 {input.covariateresults_spearman} > {output.filtered_covariateresults_spearman}
        cut -f 4- {input.filtered_spearman_bed} >> {output.filtered_covariateresults_spearman}
        """

# rule gea_scatter_plots:
#     input:
#         filtered_covariateresults = "results/{sample}_concatenated_res_covariate_filtered.csv",
#         filtered_covariateresults_spearman = "results/{sample}_concatenated_res_covariate_spearman_filtered.csv",
#         frequency_table = "data/{sample}_AlleleFrequencyTable_filter_TEs.txt",
#         envfactor_names = "data/{sample}_efile_envfactor_names",
#         efile="data/{sample}_efile",
#     resources:
#         runtime=RESOURCES["gea_scatter_plots"]["runtime"],
#         mem_mb=RESOURCES["gea_scatter_plots"]["mem_mb"],
#         slurm_partition=RESOURCES["gea_scatter_plots"]["slurm_partition"]
#     params:
#         BFthreshold = BF_THRESHOLD,
#         spearman_threshold = SPEARMAN_THRESHOLD,
#         pdfprefix = "results/{sample}_scatterplots/"
#     log:
#         "logs/gea_scatter_plots/{sample}.log"
#     output:
#         BFvsSpearman = "results/{sample}_BFvsSpearman.png",
#         significant_snps = "results/{sample}_significant_snps.csv",
#         significant_snps_bf_spearman2_5 = "results/{sample}_significant_snps_bf_spearman2_5.csv",
#         significant_snps_bf_spearman1 = "results/{sample}_significant_snps_bf_spearman1.csv",
#         significant_snps_bf_spearman0_1 = "results/{sample}_significant_snps_bf_spearman0_1.csv",
#         significant_snps_0_1 = "results/{sample}_significant_snps_top0_1_percent.csv",
#         significant_snps_1 = "results/{sample}_significant_snps_top1_percent.csv",
#         significant_snps_5 = "results/{sample}_significant_snps_top5_percent.csv",
#         significant_snps_10 = "results/{sample}_significant_snps_top10_percent.csv",
#         scatterplots = expand("results/{{sample}}_scatterplots/{envfactor}_scatterplots.pdf", envfactor=ENVFACTOR_NAMES)
#     script: "scripts/scatter_plots.R"

# rule concatenate_results_c2:
#     input:
#         contrast = expand("results/{{sample}}_baypassSplitOut_contrast/contrast_{i}_summary_betai_reg.out", 
#          i=range(1, N_SUBS+1))
#     params:
#         prefixcontrast = "results/{sample}_baypassSplitOut_contrast/contrast",
#         subs=N_SUBS,
#         snpdetprefix = "results/subsets/{sample}.snpdet.sub",
#         retrieve_c2= True
#     resources:
#         runtime=RESOURCES["concatenate_results_c2"]["runtime"],
#         mem_mb=RESOURCES["concatenate_results_c2"]["mem_mb"],
#         slurm_partition=RESOURCES["concatenate_results_c2"]["slurm_partition"]       
#     output:
#         contrasttable = "results/{sample}_concatenated_res_contrast.csv"
#     script: "scripts/concatenate_res.R"

# rule prepare_WZA:
#     input:
#         filtered_covariateresults = "results/{sample}_concatenated_res_covariate_filtered.csv",
#         filtered_covariateresults_spearman = "results/{sample}_concatenated_res_covariate_spearman_filtered.csv",
#         frequency_table = "data/{sample}_AlleleFrequencyTable_filter_TEs.txt",
#     params:
#         window_size = WZA_WINDOW_SIZE,
#     resources:
#         runtime=RESOURCES["prepare_WZA"]["runtime"],
#         mem_mb=RESOURCES["prepare_WZA"]["mem_mb"],
#         slurm_partition=RESOURCES["prepare_WZA"]["slurm_partition"]
#     output:
#         WZA_input = "data/WZA/{sample}_WZA_input.csv"
#     shell:
#         """
#         paste -d',' {input.filtered_covariateresults_spearman} <(cut -d ',' -f 6- {input.filtered_covariateresults}) \
#             | awk -f scripts/Make_Genomic_Windows.sh -v window_size={params.window_size} > tmp1_{wildcards.sample}

#         awk -f scripts/Calculate_Mean_MAF.sh {input.frequency_table} | cut -d',' -f3 > tmpMAF_{wildcards.sample}

#         paste -d',' tmp1_{wildcards.sample} tmpMAF_{wildcards.sample} > {output}
#         rm tmp1_{wildcards.sample} tmpMAF_{wildcards.sample}
#         """

# rule filter_WZA_input:
#     input:
#         WZA_input = "data/WZA/{sample}_WZA_input.csv",
#     params:
#         wza_maf_threshold = WZA_MAF_THRESHOLD,
#     output:
#         WZA_input_filtered = "data/WZA/{sample}_WZA_input_filtered.csv"
#     script:
#         "scripts/filter_WZA_input.R"

# rule WZA:
#     input:
#         WZA_input_filtered = "data/WZA/{sample}_WZA_input_filtered.csv",
#     resources:
#         runtime=RESOURCES["WZA"]["runtime"],
#         mem_mb=RESOURCES["WZA"]["mem_mb"],
#         slurm_partition=RESOURCES["WZA"]["slurm_partition"]
#     output:
#         WZA_output_BF = "results/WZA_res/{sample}_{envfactor}_BF_WZA_output.csv",
#     shell:
#         """
#         python3 scripts/general_WZA_script.py \
#             --correlations {input.WZA_input_filtered} \
#             --summary_stat {wildcards.envfactor}_BF \
#             --window window_id \
#             --MAF MAF \
#             --sep "," \
#             --large_i_small_p \
#             --min_snps 3 \
#             --retain POS \
#             --output {output.WZA_output_BF}
#         """

# rule WZA_spearman:
#     input:
#         WZA_input_filtered = "data/WZA/{sample}_WZA_input_filtered.csv",
#     resources:
#         runtime=RESOURCES["WZA"]["runtime"],
#         mem_mb=RESOURCES["WZA"]["mem_mb"],
#         slurm_partition=RESOURCES["WZA"]["slurm_partition"]
#     output:
#         WZA_output_spearman = "results/WZA_res/{sample}_{envfactor}_spearman_WZA_output.csv",
#     shell:
#         """
#         python3 scripts/general_WZA_script.py \
#             --correlations {input.WZA_input_filtered} \
#             --summary_stat {wildcards.envfactor}_spearman \
#             --window window_id \
#             --MAF MAF \
#             --sep "," \
#             --large_i_small_p \
#             --min_snps 3 \
#             --retain POS \
#             --output {output.WZA_output_spearman}
#         """

# rule WZA_diagnostics:
#     input:
#         envfactor_names = "data/{sample}_efile_envfactor_names",
#         WZA_output_BF = expand("results/WZA_res/{{sample}}_{envfactor}_BF_WZA_output.csv", envfactor=ENVFACTOR_NAMES),
#         WZA_output_spearman = expand("results/WZA_res/{{sample}}_{envfactor}_spearman_WZA_output.csv", envfactor=ENVFACTOR_NAMES),
#     params:
#         prefix = "results/WZA_res/{sample}_",
#         FDR_level = WZA_FDR,
#     resources:
#         runtime=RESOURCES["WZA_diagnostics"]["runtime"],
#         mem_mb=RESOURCES["WZA_diagnostics"]["mem_mb"],
#         slurm_partition=RESOURCES["WZA_diagnostics"]["slurm_partition"]
#     log:
#         "logs/WZA_diagnostics/{sample}.log"
#     output:
#         WZA_manhattan_plots_BF = "results/WZA_res/{sample}_WZA_manhattan_plots_BF.png",
#         WZA_manhattan_plots_BF_wo_GIF = "results/WZA_res/{sample}_WZA_manhattan_plots_BF_wo_GIF.png",
#         WZA_manhattan_plots_spearman = "results/WZA_res/{sample}_WZA_manhattan_plots_spearman.png",
#         WZA_manhattan_plots_spearman_wo_GIF = "results/WZA_res/{sample}_WZA_manhattan_plots_spearman_wo_GIF.png",
#         WZA_q0_1 = "results/WZA_res/{sample}_significant_windows_q0.1.csv",
#         WZA_q0_05 = "results/WZA_res/{sample}_significant_windows_q0.05.csv",
#         WZA_q0_01 = "results/WZA_res/{sample}_significant_windows_q0.01.csv",
#         WZA_q0_001 = "results/WZA_res/{sample}_significant_windows_q0.001.csv",
#         WZA_q0_1_noGIF = "results/WZA_res/{sample}_significant_windows_q0.1_noGIF.csv",
#         WZA_q0_05_noGIF = "results/WZA_res/{sample}_significant_windows_q0.05_noGIF.csv",
#         WZA_q0_01_noGIF = "results/WZA_res/{sample}_significant_windows_q0.01_noGIF.csv",
#         WZA_q0_001_noGIF = "results/WZA_res/{sample}_significant_windows_q0.001_noGIF.csv",
#         WZA_q0_1_rho_top2_5 = "results/WZA_res/{sample}_significant_windows_q0.1_rho_top2_5.csv",
#         WZA_q0_05_rho_top2_5 = "results/WZA_res/{sample}_significant_windows_q0.05_rho_top2_5.csv",
#         WZA_q0_01_rho_top2_5 = "results/WZA_res/{sample}_significant_windows_q0.01_rho_top2_5.csv",
#         WZA_q0_001_rho_top2_5 = "results/WZA_res/{sample}_significant_windows_q0.001_rho_top2_5.csv",
#         WZA_q0_1_noGIF_rho_top2_5 = "results/WZA_res/{sample}_significant_windows_q0.1_noGIF_rho_top2_5.csv",
#         WZA_q0_05_noGIF_rho_top2_5 = "results/WZA_res/{sample}_significant_windows_q0.05_noGIF_rho_top2_5.csv",
#         WZA_q0_01_noGIF_rho_top2_5 = "results/WZA_res/{sample}_significant_windows_q0.01_noGIF_rho_top2_5.csv",
#         WZA_q0_001_noGIF_rho_top2_5 = "results/WZA_res/{sample}_significant_windows_q0.001_noGIF_rho_top2_5.csv",
#         WZA_top0_1 = "results/WZA_res/{sample}_significant_windows_top0_1.csv",
#         WZA_top0_05 = "results/WZA_res/{sample}_significant_windows_top0_05.csv",
#         WZA_top0_01 = "results/WZA_res/{sample}_significant_windows_top0_01.csv",
#         WZA_top0_001 = "results/WZA_res/{sample}_significant_windows_top0_001.csv",
#         WZA_top0_01_rho_top2_5 = "results/WZA_res/{sample}_significant_windows_top0_01_rho_top2_5.csv",
#         WZA_top0_001_rho_top2_5 = "results/WZA_res/{sample}_significant_windows_top0_001_rho_top2_5.csv",
#         WZA_output_fdr = expand("results/WZA_res/{{sample}}_{envfactor}_WZA_output_fdr.csv", envfactor=ENVFACTOR_NAMES)
#     script: "scripts/WZA_diagnostics.R"