##################################################
# Name: exp_16.smk
# Description: Contains the workflow and methods
#              needed for experiment 16.
#
#              This experiment is intended to 
#              compare docprofiles query speed
#              when using the original mode,
#              optimized mode, and optimized
#              mode with ftab querying.
#
# Date: Jan 10th, 2024
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_paths_to_input_fasta_exp16(wildcards):
    """ Return paths to input paired-end fastq files """
    file_list = []
    for data_file in os.listdir(f"exp16_data/"):
        if data_file.endswith(".fna") or data_file.endswith(".fa"):
            file_list.append(f"exp16_data/" + data_file)
    return file_list

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Generate the input files for pfp_doc64

rule build_complete_reference_from_all_documents_exp16:
    input:
        get_paths_to_input_fasta_exp16
    output:
        "exp16_full_reference/output.fa"
    shell:
        """
        cat {input} > {output}
        """

rule build_filelist_for_index_exp16:
    input:
        get_paths_to_input_fasta_exp16
    output:
        "exp16_filelist/filelist.txt"
    run:
        with open(output[0], "w") as out_fd:
            for i, path in enumerate(input):
                out_fd.write(f"{path} {i+1}\n")

rule simulate_input_reads_at_certain_length_and_error_rate_exp16:
    input:
        "exp16_full_reference/output.fa"
    output:
        "exp16_input_reads/{percent}_percent_error_{length}_length.fa"
    shell:
        """
        python3 {repo_dir}/src/gen_reads.py --input {input} \
                                            --error {wildcards.percent} \
                                            --num 100000 \
                                            --length {wildcards.length} > {output}
        """

rule build_pfpdoc_index_for_document_listing_exp16:
    input:
        "exp16_filelist/filelist.txt"
    output:
        "exp16_index/output.fna",
        "exp16_index/output.fna.sdap",
        "exp16_index/output.fna.edap"
    shell:
        """
        pfp_doc64 build --filelist {input} \
                        --output exp16_index/output \
                        --revcomp \
                        --two-pass exp16_index/temp \
                        --tmp-size 4GB
        """

# Section 2.2: Run the different modes of pfp_doc64

rule run_pfpdoc_in_normal_mode_exp16:
    input:
        "exp16_index/output.fna",
        "exp16_input_reads/{percent}_percent_error_{length}_length.fa"
    output:
        "exp16_results/{percent}_percent_error_{length}_length/normal.listings",
        "exp16_results/{percent}_percent_error_{length}_length/normal.log"
    shell:
        """
        pfp_doc64 run --ref exp16_index/output \
                      --pattern {input[1]} \
                      --output exp16_results/{wildcards.percent}_percent_error_{wildcards.length}_length/normal 2> {output[1]}
        """

rule run_pfpdoc_in_optimized_mode_exp16:
    input:
        "exp16_index/output.fna",
        "exp16_input_reads/{percent}_percent_error_{length}_length.fa"
    output:
        "exp16_results/{percent}_percent_error_{length}_length/optimized.listings",
        "exp16_results/{percent}_percent_error_{length}_length/optimized.log"
    shell:
        """
        pfp_doc64 run --ref exp16_index/output \
                      --pattern {input[1]} \
                      --output exp16_results/{wildcards.percent}_percent_error_{wildcards.length}_length/optimized \
                      --optimized 2> {output[1]}
        """

rule run_pfpdoc_in_ftab_mode_exp16:
    input:
        "exp16_index/output.fna",
        "exp16_input_reads/{percent}_percent_error_{length}_length.fa"
    output:
        "exp16_results/{percent}_percent_error_{length}_length/ftab.listings",
        "exp16_results/{percent}_percent_error_{length}_length/ftab.log"
    shell:
        """
        pfp_doc64 run --ref exp16_index/output \
                      --pattern {input[1]} \
                      --output exp16_results/{wildcards.percent}_percent_error_{wildcards.length}_length/ftab \
                      --ftab 2> {output[1]}
        """

rule generate_exp16_results:
    input:
        expand("exp16_results/{percent}_percent_error_{length}_length/{query_type}.log", percent=[0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1], 
                                                                                   length=[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000], 
                                                                                   query_type=["normal", "optimized", "ftab"])
    output:
        "exp16_final_csv/output.csv"
    shell:
        """
        # generate the ftab vs. normal lines
        for err in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1;
        do
            for len in 100 200 300 400 500 600 700 800 900 1000;
            do
                ftab_file="exp16_results/${{err}}_percent_error_${{len}}_length/ftab.log"
                normal_file="exp16_results/${{err}}_percent_error_${{len}}_length/normal.log"

                ftab_time=$(grep "querying the patterns" $ftab_file | awk  '{{print substr($(NF-1), 2)}}')
                normal_time=$(grep "processing the patterns" $normal_file | awk  '{{print substr($(NF-1), 2)}}')

                printf "ftab,%s,%s,%s,%s\n" $err $len $ftab_time $normal_time >> {output}
            done
        done

        # generate the optimized vs. normal lines
        for err in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1;
        do
            for len in 100 200 300 400 500 600 700 800 900 1000;
            do
                optimized_file="exp16_results/${{err}}_percent_error_${{len}}_length/optimized.log"
                normal_file="exp16_results/${{err}}_percent_error_${{len}}_length/normal.log"

                optimized_time=$(grep "processing the patterns" $optimized_file | awk  '{{print substr($(NF-1), 2)}}')
                normal_time=$(grep "processing the patterns" $normal_file | awk  '{{print substr($(NF-1), 2)}}')

                printf "optimized,%s,%s,%s,%s\n" $err $len $optimized_time $normal_time >> {output} 
            done
        done
        """
    

