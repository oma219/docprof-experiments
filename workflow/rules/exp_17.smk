##################################################
# Name: exp_17.smk
# Description: Contains the workflow and methods
#              needed for experiment 17.
#
#              This experiment is intended to 
#              compare docprofiles query speed
#              when using the original mode,
#              optimized mode, and optimized
#              mode with ftab querying.
#
# Date: Feb 22nd, 2024
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_paths_to_input_fasta_exp17(wildcards):
    """ Return paths to input paired-end fastq files """
    file_list = []
    for data_file in os.listdir(f"exp17_data/"):
        if data_file.endswith(".fna") or data_file.endswith(".fa"):
            file_list.append(f"exp17_data/" + data_file)
    return file_list

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Generate the input files for pfp_doc64

rule build_complete_reference_from_all_documents_exp17:
    input:
        get_paths_to_input_fasta_exp17
    output:
        "exp17_full_reference/output.fa"
    shell:
        """
        cat {input} > {output}
        """

rule build_filelist_for_index_exp17:
    input:
        get_paths_to_input_fasta_exp17
    output:
        "exp17_filelist/filelist.txt"
    run:
        with open(output[0], "w") as out_fd:
            for i, path in enumerate(input):
                out_fd.write(f"{path} {i+1}\n")

rule simulate_input_reads_at_certain_length_and_error_rate_exp17:
    input:
        "exp17_full_reference/output.fa"
    output:
        "exp17_input_reads/{percent}_percent_error_{length}_length.fa"
    shell:
        """
        python3 {repo_dir}/src/gen_reads.py --input {input} \
                                            --error {wildcards.percent} \
                                            --num 250000 \
                                            --length {wildcards.length} > {output}
        """

# Section 2.2: Build the different indexes needed

rule build_pfpdoc_index_for_document_listing_exp17:
    input:
        "exp17_filelist/filelist.txt"
    output:
        "exp17_index/output.fna",
        "exp17_index/output.fna.sdap",
        "exp17_index/output.fna.edap"
    shell:
        """
        pfp_doc64 build --filelist {input} \
                        --output exp17_index/output \
                        --revcomp \
                        --two-pass exp17_index/temp \
                        --tmp-size 4GB
        """

rule build_pfpdoc_minimizer_index_for_document_listing_exp17:
    input:
        "exp17_filelist/filelist.txt"
    output:
        "exp17_minimizer_index/output.fna",
        "exp17_minimizer_index/output.fna.sdap",
        "exp17_minimizer_index/output.fna.edap"
    shell:
        """
        pfp_doc64 build --filelist {input} \
                        --output exp17_minimizer_index/output \
                        --revcomp \
                        --two-pass exp17_minimizer_index/temp \
                        --tmp-size 4GB \
                        --minimizers \
                        --small-window 4 \
                        --large-window 11
        """

# Section 2.3: Run the different modes of pfp_doc64 (optimized, optimized+ftab, optimized+minimizer, optimized+minimizer+ftab)

rule run_pfpdoc_in_normal_mode_exp17:
    input:
        "exp17_index/output.fna",
        "exp17_input_reads/{percent}_percent_error_{length}_length.fa"
    output:
        "exp17_results/{percent}_percent_error_{length}_length/normal.listings",
        "exp17_results/{percent}_percent_error_{length}_length/normal.log"
    shell:
        """
        pfp_doc64 run --ref exp17_index/output \
                      --pattern {input[1]} \
                      --output exp17_results/{wildcards.percent}_percent_error_{wildcards.length}_length/normal \
                      2> {output[1]}
        """

rule run_pfpdoc_in_optimized_mode_exp17:
    input:
        "exp17_index/output.fna",
        "exp17_input_reads/{percent}_percent_error_{length}_length.fa"
    output:
        "exp17_results/{percent}_percent_error_{length}_length/optimized.listings",
        "exp17_results/{percent}_percent_error_{length}_length/optimized.log"
    shell:
        """
        pfp_doc64 run --ref exp17_index/output \
                      --pattern {input[1]} \
                      --output exp17_results/{wildcards.percent}_percent_error_{wildcards.length}_length/optimized \
                      --optimized 2> {output[1]}
        """

rule run_pfpdoc_in_ftab_mode_exp17:
    input:
        "exp17_index/output.fna",
        "exp17_input_reads/{percent}_percent_error_{length}_length.fa"
    output:
        "exp17_results/{percent}_percent_error_{length}_length/ftab.listings",
        "exp17_results/{percent}_percent_error_{length}_length/ftab.log"
    shell:
        """
        pfp_doc64 run --ref exp17_index/output \
                      --pattern {input[1]} \
                      --output exp17_results/{wildcards.percent}_percent_error_{wildcards.length}_length/ftab \
                      --ftab 2> {output[1]}
        """

rule run_pfpdoc_in_minimizer_mode_exp17:
    input:
        "exp17_minimizer_index/output.fna",
        "exp17_input_reads/{percent}_percent_error_{length}_length.fa"
    output:
        "exp17_results/{percent}_percent_error_{length}_length/minimizer.listings",
        "exp17_results/{percent}_percent_error_{length}_length/minimizer.log"
    shell:
        """
        pfp_doc64 run --ref exp17_minimizer_index/output \
                      --pattern {input[1]} \
                      --output exp17_results/{wildcards.percent}_percent_error_{wildcards.length}_length/minimizer \
                      --minimizers \
                      --small-window 4 \
                      --large-window 11 \
                      --optimized 2> {output[1]}
        """

rule run_pfpdoc_in_minimizer_ftab_mode_exp17:
    input:
        "exp17_minimizer_index/output.fna",
        "exp17_input_reads/{percent}_percent_error_{length}_length.fa"
    output:
        "exp17_results/{percent}_percent_error_{length}_length/minimizer_ftab.listings",
        "exp17_results/{percent}_percent_error_{length}_length/minimizer_ftab.log"
    shell:
        """
        pfp_doc64 run --ref exp17_minimizer_index/output \
                      --pattern {input[1]} \
                      --output exp17_results/{wildcards.percent}_percent_error_{wildcards.length}_length/minimizer_ftab \
                      --minimizers \
                      --small-window 4 \
                      --large-window 11 \
                      --ftab  2> {output[1]}
        """

# Section 2.4: Gather all the results

rule generate_exp17_results:
    input:
        expand("exp17_results/{percent}_percent_error_{length}_length/{query_type}.log", percent=[0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1], 
                                                                                   length=[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 4000], 
                                                                                   query_type=["normal", "optimized", "minimizer", "ftab", "minimizer_ftab"])
    output:
        "exp17_final_csv/output.csv"
    shell:
        """
        echo helloworld1
        # generate the normal vs. optimized lines
        for err in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1;
        do
            for len in 100 200 300 400 500 600 700 800 900 1000 1500 2000 4000;
            do
                printf "%s,%s\n" $err $len
                optimized_file="exp17_results/${{err}}_percent_error_${{len}}_length/optimized.log"
                normal_file="exp17_results/${{err}}_percent_error_${{len}}_length/normal.log"

                optimized_time=$(grep "processing the patterns" $optimized_file | awk  '{{print substr($(NF-1), 2)}}')
                normal_time=$(grep "processing the patterns" $normal_file | awk  '{{print substr($(NF-1), 2)}}')

                printf "optimized,%s,%s,%s,%s\n" $err $len $optimized_time $normal_time >> {output}
            done
        done
        echo helloworld1
        # generate the normal vs. minimizer lines
        for err in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1;
        do
            for len in 100 200 300 400 500 600 700 800 900 1000 1500 2000 4000;
            do
                minimizer_file="exp17_results/${{err}}_percent_error_${{len}}_length/minimizer.log"
                normal_file="exp17_results/${{err}}_percent_error_${{len}}_length/normal.log"

                minimizer_time=$(grep "processing the patterns" $minimizer_file | awk  '{{print substr($(NF-1), 2)}}')
                normal_time=$(grep "processing the patterns" $normal_file | awk  '{{print substr($(NF-1), 2)}}')

                printf "minimizer,%s,%s,%s,%s\n" $err $len $minimizer_time $normal_time >> {output}
            done
        done

        # generate the normal vs. ftab lines
        for err in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1;
        do
            for len in 100 200 300 400 500 600 700 800 900 1000 1500 2000 4000;
            do
                ftab_file="exp17_results/${{err}}_percent_error_${{len}}_length/ftab.log"
                normal_file="exp17_results/${{err}}_percent_error_${{len}}_length/normal.log"

                ftab_time=$(grep "querying the patterns" $ftab_file | awk  '{{print substr($(NF-1), 2)}}')
                normal_time=$(grep "processing the patterns" $normal_file | awk  '{{print substr($(NF-1), 2)}}')

                printf "ftab,%s,%s,%s,%s\n" $err $len $ftab_time $normal_time >> {output} 
            done
        done

        # generate the normal vs. minimizer_ftab lines
        for err in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1;
        do
            for len in 100 200 300 400 500 600 700 800 900 1000 1500 2000 4000;
            do
                minimizer_ftab_file="exp17_results/${{err}}_percent_error_${{len}}_length/minimizer_ftab.log"
                normal_file="exp17_results/${{err}}_percent_error_${{len}}_length/normal.log"

                minimizer_ftab_time=$(grep "querying the patterns" $minimizer_ftab_file | awk  '{{print substr($(NF-1), 2)}}')
                normal_time=$(grep "processing the patterns" $normal_file | awk  '{{print substr($(NF-1), 2)}}')

                printf "minimizer+ftab,%s,%s,%s,%s\n" $err $len $minimizer_ftab_time $normal_time >> {output} 
            done
        done
        """
    

