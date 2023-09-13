##################################################
# Name: exp_11.smk
# Description: Contains the workflow and methods
#              needed for experiment 11.
#
#              This experiment focuses on a small
#              handmade taxonomy with a sub-selection
#              of species and want to visualize
#              the LCA hits.
#
# Date: Sept 12th, 2023
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_input_fasta_files_exp11(wildcards):
    """ Returns path to all genomes files that will be used in database """
    file_list = []
    for dataset_num in range(1, num_docs_exp11+1):
        for data_file in os.listdir(f"exp11_data/dataset_{dataset_num}"):
            if data_file.endswith(".fna"):
                file_list.append(f"exp11_data/dataset_{dataset_num}/" + data_file)
    return file_list

def get_specific_fasta_file_exp11(wildcards):
    """ Returns a specific genome associated with a certains dataset number """
    file_list = []
    for data_file in os.listdir(f"exp11_data/dataset_{wildcards.num}"):
        if data_file.endswith(".fna"):
            file_list.append(f"exp11_data/dataset_{wildcards.num}/" +  data_file)
    assert len(file_list) == 1
    return file_list

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Create a filelist to be used for index building

rule create_filelist_for_indexes_exp11:
    input:
        get_input_fasta_files_exp11
    output:
        "exp11_indexes/filelist.txt"
    run:
        with open(output[0], "w") as out_fd:
            for i, path in enumerate(input):
                out_fd.write(f"{path} {i+1}\n")

# Section 2.2: Simulate reads from each leaf node

rule simulate_short_reads_from_each_database_exp11:
    input:
        get_specific_fasta_file_exp11
    output:
        "exp11_reads/dataset_{num}_reads.fa"
    shell:
        """
        mason_simulator -ir {input} -n {num_reads_exp11} -v -o {output[0]} --illumina-read-length 150
        """

# Section 2.3: Build SPUMONI index and DocProfiles index

rule build_spumoni_index_exp11:
    input:
        "exp11_indexes/filelist.txt"
    output:
        "exp11_indexes/spumoni_index/output.fa"
    shell:
        """
        spumoni build --filelist {input} \
                      --MS --PML \
                      --no-digest \
                      --prefix exp11_indexes/spumoni_index/output 
        """

rule build_pfpdoc_tax_index_exp11:
    input:
        "exp11_indexes/filelist.txt"
    output:
        "exp11_indexes/doc_index/filelist.txt",
        "exp11_indexes/doc_index/output.fna.taxcomp.sdap",
        "exp11_indexes/doc_index/output.fna.taxcomp.edap",
        "exp11_indexes/doc_index/output.fna.taxcomp.of.sdap",
        "exp11_indexes/doc_index/output.fna.taxcomp.of.edap"
    shell:
        """
        cp {input[0]} {output[0]}
        pfp_doc64 build -f {output[0]} \
                        -o exp11_indexes/doc_index/output \
                        --taxcomp \
                        --num-col 7 -n \
                        --revcomp
        """

# Section 2.4: Generate MS with respect to a database of all genomes

rule run_spumoni_to_generate_ms_exp11:
    input:
        "exp11_indexes/spumoni_index/output.fa",
        "exp11_reads/dataset_{num}_reads.fa"
    output:
        "exp11_results/ms_lengths/dataset_{num}.fa",
        "exp11_results/ms_lengths/dataset_{num}.fa.lengths"
    shell:
        """
        cp {input[1]} {output[0]}

        # Don't make this multi-threaded, it will mess up the next step
        spumoni run -r exp11_indexes/spumoni_index/output -p {output[0]} -M -n 
        """

#   Section 2.5: Extract MEMs using the MS lengths
#   and the pivot reads itself. 

rule extract_mems_or_halfmems_based_on_ms_exp11:
    input:
        "exp11_results/ms_lengths/dataset_{num}.fa",
        "exp11_results/ms_lengths/dataset_{num}.fa.lengths"
    output:
        "exp11_results/extracted_mems/half_mems_dataset_{num}.fastq"
    shell:
        """
        python3 {repo_dir}/src/extract_mems.py \
        -l {input[1]} \
        -p {input[0]} \
        -c {num_reads_exp11} \
        --mems \
        -t 10 \
        -o {output[0]} 
        """

# Section 2.6: Compute the leftmost, rightmost occurrences using 
# the document array profiles

rule determined_doc_array_listings_for_mems_exp11:
    input:
        "exp11_indexes/doc_index/output.fna.taxcomp.sdap",
        "exp11_indexes/doc_index/output.fna.taxcomp.edap",
        "exp11_indexes/doc_index/output.fna.taxcomp.of.sdap",
        "exp11_indexes/doc_index/output.fna.taxcomp.of.edap",
        "exp11_results/extracted_mems/half_mems_dataset_{num}.fastq"
    output:
        "exp11_results/listings/half_mems_dataset_{num}.fastq",
        "exp11_results/listings/half_mems_dataset_{num}.fastq.listings"
    shell:
        """
        cp {input[4]} {output[0]}
        pfp_doc64 run --ref exp11_indexes/doc_index/output \
                      --pattern {output[0]} \
                      --num-col 7 \
                      --taxcomp
        """

# Section 2.7: Take the listings and plot the LCA on the
# taxonomy on a tree diagram.

rule plot_lca_on_taxonomy_exp11:
    input:
        expand("exp11_results/listings/half_mems_dataset_{num}.fastq.listings", num=range(1, num_docs_exp11+1))
    output:
        expand("exp11_results/plots/doc_{num}_results.txt", num=range(1, num_docs_exp11+1)),
        expand("exp11_results/plots/doc_{num}_plots.png", num=range(1, num_docs_exp11+1))
    shell:
        """
        echo helloworld > {output[0]}
        python3 {repo_dir}/src/plot_lca_on_tree_exp11.py \
                --taxonomy exp11_data/exp11_taxonomy.txt \
                --input {input} \
                --num-of-docs {num_docs_exp11} \
                --doc-to-name exp11_data/exp11_doc_names.txt \
                --output-dir exp11_results/plots/
        """

