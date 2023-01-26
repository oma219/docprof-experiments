##################################################
# Name: exp_8.smk
# Description: Contains the workflow and methods
#              needed for experiment 8.
#
#              Analyze the monotonic increases
#              of document array profiles from
#              left to right, and visa-versa.
#
# Date: December 4th, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_input_fasta_file_exp8(wildcards):
    """ Returns path to SILVA database with taxonomic path in each header """
    file_list = []
    for data_file in os.listdir(f"data/"):
        if data_file.endswith(".fasta"):
            file_list.append(f"data/" + data_file)
    assert len(file_list) == 1
    return file_list[0]


####################################################
# Section 2: Rules needed for this experiment type
####################################################

rule convert_to_uppercase_and_samelines_exp8:
    input:
        get_input_fasta_file_exp8
    output:
        "data/silva_database.fa"
    shell:
        "seqtk seq -U {input} > {output}"

rule generate_separate_fasta_files_exp8:
    input:
        "data/silva_database.fa"
    output:
        "exp8_input_files/filelist.txt"
    shell:
        """
        python3 {repo_dir}/src/separate_silva_classes.py \
        -i {input} \
        -o exp8_input_files/ \
        --use-order \
        --tree data/tax_slv_ssu_138.1.tre \
        --tree-map data/tax_slv_ssu_138.1.map \
        -n {num_genera_exp8} \
        """

rule build_docprofs_for_silva_exp8:
    input:
        "exp8_input_files/filelist.txt"
    output:
        "exp8_index/output.fna.edap",
        "exp8_index/output.fna.sdap"
    shell:
        """
        pfp_doc build -f {input} -o exp8_index/output
        """

rule print_out_profiles_exp8:
    input:
        "exp8_index/output.fna.edap",
        "exp8_index/output.fna.sdap"
    output:
        "exp8_profiles/output.sdap.csv",
        "exp8_profiles/output.edap.csv"
    shell:
        """
        pfp_doc info -r exp8_index/output -o exp8_profiles/output -n 1000000
        """

rule analyze_the_profiles_exp8:
    input:
        "exp8_profiles/output.sdap.csv",
        "exp8_profiles/output.edap.csv"
    output:
        "exp8_profiles/output.monotonic_increases_plot.png",
        "exp8_profiles/output.monotonic_profiles_plot.png"
    shell:
        """
        python3 {repo_dir}/src/analyze_docprofs.py -i exp8_profiles/output
        """



