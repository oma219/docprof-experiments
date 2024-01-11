##################################################
# Name: exp_15.smk
# Description: Contains the workflow and methods
#              needed for experiment 15.
#
#              This experiment is intended to 
#              compare docprofiles to Kraken
#              on the problem of 16S rRNA 
#              classification.
#
# Date: Nov 21st, 2023
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_paths_to_input_fastq_exp15(wildcards):
    """ Return paths to input paired-end fastq files """
    file_list = ["", ""]
    for data_file in os.listdir(f"exp15_data/"):
        if data_file.endswith(".fastq") and "_mate1" in data_file:
            file_list[0] = f"exp15_data/" + data_file
        elif data_file.endswith(".fastq") and "_mate2" in data_file:
            file_list[1] = f"exp15_data/" + data_file
    assert file_list[0] != "" and file_list[1] != ""
    print(file_list)
    return file_list

def get_input_files_for_analysis_exp15(wildcards):
    """ Returns the names of files that we want to use as input """
####################################################
# Section 2: Rules needed for this experiment type
####################################################

rule remove_ambiguous_reads_exp15:
    input:
        get_paths_to_input_fastq_exp15
    output:
        "exp15_input_reads/{sample_name_exp15}_mate1.fastq",
        "exp15_input_reads/{sample_name_exp15}_mate2.fastq"
    shell:
        """
        python3 {repo_dir}/src/remove_ambiguous_reads.py \
                --mate1 {input[0]} \
                --mate2 {input[1]} \
                --output-dir exp15_input_reads/ \
                --read-taxrank {path_to_read_tax_exp15} \
                --silva-taxrank {silva_tax_ranks_exp15} \
                --sample {sample_name_exp15}
        """

rule run_kraken2_and_bracken_exp15:
    input:
        "exp15_input_reads/{sample_name_exp15}_mate1.fastq",
        "exp15_input_reads/{sample_name_exp15}_mate2.fastq"
    output:
        "exp15_results/kraken2_minimizers/{sample_name_exp15}.kreport2",
        "exp15_results/kraken2_minimizers/{sample_name_exp15}.kraken2",
        "exp15_results/kraken2_minimizers/{sample_name_exp15}_genus.bracken",
        "exp15_results/kraken2_minimizers/{sample_name_exp15}_kraken2.time",
        "exp15_results/kraken2_minimizers/{sample_name_exp15}_bracken.time"
    shell:
        """
        {time_prog} {time_format} --output={output[3]} kraken2 --db silva \
                                                        --threads 1 \
                                                        --report {output[0]} \
                                                        --paired {input[0]} {input[1]} \
                                                        > {output[1]}
        {time_prog} {time_format} --output={output[4]} bracken -d {kraken_db_path_exp15} \
                                                        -r 250 \
                                                        -l G \
                                                        -i {output[0]} \
                                                        -o {output[2]} 
        """

rule run_kraken2_and_bracken_with_no_minimizers_exp15:
    input:
        "exp15_input_reads/{sample_name_exp15}_mate1.fastq",
        "exp15_input_reads/{sample_name_exp15}_mate2.fastq"
    output:
        "exp15_results/kraken2_no_minimizers/{sample_name_exp15}.kreport2",
        "exp15_results/kraken2_no_minimizers/{sample_name_exp15}.kraken2",
        "exp15_results/kraken2_no_minimizers/{sample_name_exp15}_genus.bracken",
        "exp15_results/kraken2_no_minimizers/{sample_name_exp15}_kraken2.time",
        "exp15_results/kraken2_no_minimizers/{sample_name_exp15}_bracken.time"
    shell:
        """
        {time_prog} {time_format} --output={output[3]} kraken2 --db silva_nominimizer \
                                                        --threads 1 \
                                                        --report {output[0]} \
                                                        --paired {input[0]} {input[1]} \
                                                        > {output[1]}
        {time_prog} {time_format} --output={output[4]} bracken -d {kraken_nominimizer_db_path_exp15} \
                                                        -r 250 \
                                                        -l G \
                                                        -i {output[0]} \
                                                        -o {output[2]} 
        """

rule run_pfpdoc_no_minimizers_exp15:
    input:
        "exp15_input_reads/{sample_name_exp15}_mate1.fastq",
        "exp15_input_reads/{sample_name_exp15}_mate2.fastq"
    output:
        "exp15_results/docprof_no_minimizers/{sample_name_exp15}_mate1.fastq",
        "exp15_results/docprof_no_minimizers/{sample_name_exp15}_mate2.fastq",
        "exp15_results/docprof_no_minimizers/{sample_name_exp15}_mate1.fastq.listings",
        "exp15_results/docprof_no_minimizers/{sample_name_exp15}_mate2.fastq.listings",
        "exp15_results/docprof_no_minimizers/{sample_name_exp15}_mate1.log",
        "exp15_results/docprof_no_minimizers/{sample_name_exp15}_mate2.log",
        "exp15_results/docprof_no_minimizers/{sample_name_exp15}_mate1.time",
        "exp15_results/docprof_no_minimizers/{sample_name_exp15}_mate2.time"
    shell:
        """
        cp {input[0]} {output[0]}
        pfp_doc64 run \
				  --ref {pfpdoc_db_exp15} \
				  --pattern {output[0]} \
				  --taxcomp  \
				  --num-col 7 2> {output[4]}
        grep "processing" {output[4]} | awk '{{print substr($7, 2)}}' >> {output[6]}

        cp {input[1]} {output[1]}
        pfp_doc64 run \
				  --ref {pfpdoc_db_exp15} \
				  --pattern {output[1]} \
				  --taxcomp \
				  --num-col 7 2> {output[5]}
        grep "processing" {output[5]} | awk '{{print substr($7, 2)}}' >> {output[7]}
        """


rule analyze_results_from_docprof_approaches_exp15:
    input:
        "exp15_results/docprof_no_minimizers/{sample_name_exp15}_mate1.fastq.listings",
        "exp15_results/docprof_no_minimizers/{sample_name_exp15}_mate2.fastq.listings"
    output:
        "exp15_results/docprof_no_minimizers/{sample_name_exp15}.classification_results.csv",
        "exp15_results/docprof_no_minimizers/{sample_name_exp15}.abundance_results.csv"
        # "exp15_results/docprof_no_minimizers/{sample_name_exp15}.debug.txt"
    shell:
        """
        touch {output[1]}
        python3 {repo_dir}/src/classify_silva_readset.py \
                        --mate1-listings {input[0]} \
                        --mate2-listings {input[1]} \
                        --doc-id-to-traversal {doc_to_trav_exp15} \
                        --silva-tax-ranks {silva_tax_ranks_exp15} \
                        --readset-truthset {path_to_read_tax_exp15} \
                        --output-dir exp15_results/docprof_no_minimizers/ \
                        --trav-to-length {trav_to_length_exp15} \
                        --sample-name {sample_name_exp15} \
                        --classify-docprof 

        """

rule analyze_results_from_kraken2_approaches_exp15:
    input:
        "exp15_results/{kraken_approach}/{sample_name_exp15}.kreport2",
        "exp15_results/{kraken_approach}/{sample_name_exp15}.kraken2",
        "exp15_results/{kraken_approach}/{sample_name_exp15}_genus.bracken"
    output:
        "exp15_results/{kraken_approach}/{sample_name_exp15}.classification_results.csv",
        "exp15_results/{kraken_approach}/{sample_name_exp15}.abundance_results.csv"
    shell:
        """
        python3 {repo_dir}/src/classify_silva_readset.py \
                        --doc-id-to-traversal {doc_to_trav_exp15} \
                        --silva-tax-ranks {silva_tax_ranks_exp15} \
                        --readset-truthset {path_to_read_tax_exp15} \
                        --output-dir exp15_results/{wildcards.kraken_approach}/ \
 				        --kraken2-read-class {input[1]} \
 				        --bracken-output {input[2]} \
                        --trav-to-length {trav_to_length_exp15} \
                        --sample-name {sample_name_exp15} \
                        --classify-kraken
        """

rule run_exp15:
    params:
        sample_name_exp15=config["SAMPLE_NAME_EXP15"]
    input:
        expand("exp15_results/{method}/" + sample_name_exp15 + ".classification_results.csv", method=["kraken2_minimizers", "docprof_no_minimizers", "kraken2_no_minimizers"]),
        expand("exp15_results/{method}/" + sample_name_exp15 + ".abundance_results.csv", method=["kraken2_minimizers", "docprof_no_minimizers", "kraken2_no_minimizers"])

