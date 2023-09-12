##################################################
# Name: exp_10.smk
# Description: Contains the workflow and methods
#              needed for experiment 10.
#
#              This experiment verifies that the 
#              output of the taxonomic document
#              array is correct by double checking
#              it with the r-index.
#
# Date: Sept 8th, 2023
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_input_fasta_file_exp10(wildcards):
    """ Returns path to SILVA database with taxonomic path in each header """
    file_list = []
    for data_file in os.listdir(f"exp10_data/"):
        if data_file.endswith(".fasta"):
            file_list.append(f"exp10_data/" + data_file)
    assert len(file_list) == 1
    return file_list[0]

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Generate the n genera files

rule convert_to_uppercase_and_samelines_exp10:
    input:
        get_input_fasta_file_exp10
    output:
        "exp10_data/silva_database.fa"
    shell:
        "seqtk seq -U {input} > {output}"

rule generate_separate_fasta_files_exp10:
    input:
        "exp10_data/silva_database.fa"
    output:
        "exp10_input_files/filelist.txt"
    shell:
        """
        python3 {repo_dir}/src/separate_silva_classes.py \
        -i {input} \
        -o exp10_input_files/ \
        --use-order \
        --tree exp10_data/tax_slv_ssu_138.1.tre \
        --tree-map exp10_data/tax_slv_ssu_138.1.map \
        -n {num_genera_exp10} \
        """


# Section 2.2: Build SPUMONI and pfp_doc index and r-index 

rule build_spumoni_index_exp10:
    input:
        "exp10_input_files/filelist.txt"
    output:
        "exp10_indexes/spumoni_index/output.fa"
    shell:
        """
        spumoni build --filelist {input} \
                      --no-rev-comp \
                      --MS --PML \
                      --no-digest \
                      --prefix exp10_indexes/spumoni_index/output 
        """

rule build_pfpdoc_tax_index_exp10:
    input:
        "exp10_input_files/filelist.txt"
    output:
        "exp10_indexes/doc_index/filelist.txt",
        "exp10_indexes/doc_index/output.fna.taxcomp.sdap",
        "exp10_indexes/doc_index/output.fna.taxcomp.edap",
        "exp10_indexes/doc_index/output.fna.taxcomp.of.sdap",
        "exp10_indexes/doc_index/output.fna.taxcomp.of.edap"
    shell:
        """
        cp {input[0]} {output[0]}
        pfp_doc64 build -f {output} \
                        -o exp10_indexes/doc_index/output \
                        --taxcomp \
                        --num-col 7 -n
        """

rule build_rindex_exp10:
    input:
        "exp10_indexes/spumoni_index/output.fa"
    output:
        "exp10_indexes/r_index/output.fa",
        "exp10_indexes/r_index/output.fa.ri"
    shell:
        """
        cp {input[0]} {output[0]} 
        cd {r_dir} 
        ri-buildfasta {base_dir}/{output[0]}
        cd {base_dir}
        """


# Section 2.3: Simulate reads from the full reference, 
# and compute MS with respect to it.

rule simulate_short_reads_exp10:
    input:
        "exp10_indexes/spumoni_index/output.fa"
    output:
        "exp10_reads/reads.fa"
    shell:
        """
        python3 {repo_dir}/src/gen_reads.py -i {input} -e 0.05 -n 1000000 --RNA > {output}
        """

rule run_spumoni_to_generate_ms_exp10:
    input:
        "exp10_indexes/spumoni_index/output.fa",
        "exp10_reads/reads.fa"
    output:
        "exp10_reads/reads.fa.lengths"
    shell:
        """
        # Don't make this multi-threaded, it will mess up the next step
        spumoni run -r exp10_indexes/spumoni_index/output -p {input[1]} -M -n 
        """

#   Section 2.4: Extract MEMs using the MS lengths
#   and the pivot reads itself. 

rule extract_mems_or_halfmems_based_on_ms_exp10:
    input:
        "exp10_reads/reads.fa.lengths",
        "exp10_reads/reads.fa"
    output:
        "exp10_reads/half_mems.fastq"
    shell:
        """
        python3 {repo_dir}/src/extract_mems.py \
        -l {input[0]} \
        -p {input[1]} \
        -c {num_mems_exp10} \
        --mems \
        -t 10 \
        -o {output} 
        """

# Section 2.5: Compute the leftmost, rightmost occurrences using 
# the document array profiles and use r-index to find all occurrences

rule determined_doc_array_listings_for_mems_exp10:
    input:
        "exp10_indexes/doc_index/output.fna.taxcomp.sdap",
        "exp10_indexes/doc_index/output.fna.taxcomp.edap",
        "exp10_indexes/doc_index/output.fna.taxcomp.of.sdap",
        "exp10_indexes/doc_index/output.fna.taxcomp.of.edap",
        "exp10_reads/half_mems.fastq"
    output:
        "exp10_results/doc_array/half_mems.fastq",
        "exp10_results/doc_array/half_mems.fastq.listings"
    shell:
        """
        cp {input[4]} {output[0]}
        pfp_doc64 run --ref exp10_indexes/doc_index/output \
                      --pattern {output[0]} \
                      --num-col 7 \
                      --taxcomp
        """

rule align_against_rindex_of_database_exp10:
    input:
        "exp10_indexes/r_index/output.fa",
        "exp10_indexes/r_index/output.fa.ri",
        "exp10_reads/half_mems.fastq"
    output:
        "exp10_results/r_index/align_to_database.sam"
    shell:
        """
        cd {r_dir} 
        ri-align locate {base_dir}/{input[0]} {base_dir}/{input[2]} \
        > {base_dir}/{output}
        """

# Section 2.6: Compare the listings

rule compare_doc_listings_with_rindex_exp10:
    input:
        "exp10_input_files/filelist.txt",
        "exp10_results/r_index/align_to_database.sam",
        "exp10_results/doc_array/half_mems.fastq.listings"
    output:
        "exp10_final_results/output.txt"
    shell:
        """
        python3 {repo_dir}/src/compare_lca_queries.py \
        --filelist {input[0]} \
        --sam {input[1]} \
        --listings {input[2]} \
        --output {output[0]}
        """


