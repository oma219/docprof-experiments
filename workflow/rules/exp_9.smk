##################################################
# Name: exp_9.smk
# Description: Contains the workflow and methods
#              needed for experiment 9.
#
#              Analyze the size of the taxonomic
#              compressed data-structure at different
#              number of columnms.
#
# Date: April 4th, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_input_fasta_file_exp9(wildcards):
    """ Returns path to SILVA database with taxonomic path in each header """
    file_list = []
    for data_file in os.listdir(f"exp9_data/"):
        if data_file.endswith(".fasta"):
            file_list.append(f"exp9_data/" + data_file)
    assert len(file_list) == 1
    return file_list[0]


####################################################
# Section 2: Rules needed for this experiment type
####################################################

rule convert_to_uppercase_and_samelines_exp9:
    input:
        get_input_fasta_file_exp9
    output:
        "exp9_data/silva_database.fa"
    shell:
        "seqtk seq -U {input} > {output}"

rule generate_separate_fasta_files_exp9:
    input:
        "exp9_data/silva_database.fa"
    output:
        "exp9_input_files/filelist.txt"
    shell:
        """
        python3 {repo_dir}/src/separate_silva_classes.py \
        -i {input} \
        -o exp9_input_files/ \
        --use-order \
        --tree exp9_data/tax_slv_ssu_138.1.tre \
        --tree-map exp9_data/tax_slv_ssu_138.1.map \
        -n {num_genera_exp9} \
        """

rule build_docprofs_with_taxonomic_comp_exp9:
    input:
        "exp9_input_files/filelist.txt"
    output:
        "exp9_indexes/num_cols_{col}/filelist.txt"
    shell:
        """
        cp {input[0]} {output}
        pfp_doc64 build -f {output} \
                        -o exp9_indexes/num_cols_{wildcards.col}/output \
                        --taxcomp \
                        --num-col {wildcards.col}
        """

rule build_docprofs_with_topk_comp_exp9:
    input:
        "exp9_input_files/filelist.txt"
    output:
        "exp9_indexes/topk_num_cols_{col}/filelist.txt"
    shell:
        """
        cp {input[0]} {output}
        pfp_doc64 build -f {output} \
                        -o exp9_indexes/topk_num_cols_{wildcards.col}/output \
                        --top-k \
                        --num-col {wildcards.col}
        """

rule build_docprofs_indexes_for_full_profile_exp9:
    input:
        "exp9_input_files/filelist.txt"
    output:
        "exp9_indexes/full_structure/filelist.txt"
    shell:
        """
        cp {input[0]} {output}
        pfp_doc64 build -f {output} \
                        -o exp9_indexes/full_structure/output \
        """

rule analyze_docprofs_indexes_with_taxonomic_comp_exp9:
    input:
        "exp9_indexes/num_cols_{col}/filelist.txt"
    output:
        "exp9_indexes/num_cols_{col}/index_stats.txt"
    shell:
        """
        sdap=$(ls -l exp9_indexes/num_cols_{wildcards.col}/output.fna.taxcomp.sdap | awk '{{print $5}}')
        edap=$(ls -l exp9_indexes/num_cols_{wildcards.col}/output.fna.taxcomp.edap | awk '{{print $5}}')
        sdap_of=$(ls -l exp9_indexes/num_cols_{wildcards.col}/output.fna.taxcomp.of.sdap | awk '{{print $5}}')
        edap_of=$(ls -l exp9_indexes/num_cols_{wildcards.col}/output.fna.taxcomp.of.edap | awk '{{print $5}}')

        table_sum=$(($sdap + $edap))
        of_sum=$(($sdap_of + $edap_of))

        printf "%s,%d,%d\n" {wildcards.col} $table_sum $of_sum > {output}
        """

rule analyze_docprofs_indexes_with_topk_comp_exp9:
    input:
        "exp9_indexes/topk_num_cols_{col}/filelist.txt"
    output:
        "exp9_indexes/topk_num_cols_{col}/index_stats.txt"
    shell:
        """
        sdap=$(ls -l exp9_indexes/topk_num_cols_{wildcards.col}/output.fna.topk.sdap | awk '{{print $5}}')
        edap=$(ls -l exp9_indexes/topk_num_cols_{wildcards.col}/output.fna.topk.edap | awk '{{print $5}}')

        table_sum=$(($sdap + $edap))

        printf "%d\n" $table_sum > {output}
        """

rule analyze_docprofs_indexes_for_full_profile_exp9:
    input:
        "exp9_indexes/full_structure/filelist.txt"
    output:
        "exp9_indexes/full_structure/index_stats.txt"
    shell:
        """
        sdap=$(ls -l exp9_indexes/full_structure/output.fna.sdap | awk '{{print $5}}')
        edap=$(ls -l exp9_indexes/full_structure/output.fna.edap | awk '{{print $5}}')

        total_sum=$(($sdap + $edap))

        printf "%d\n" $total_sum > {output}
        """

rule generate_exp9_files:
    input:
        expand("exp9_indexes/num_cols_{col}/index_stats.txt", col=range(2,16)),
        expand("exp9_indexes/topk_num_cols_{col}/index_stats.txt", col=range(2,16)),
        "exp9_indexes/full_structure/index_stats.txt"
    output:
        "exp9_output/analysis.csv"
    shell:
        """
        for col in $(seq 2 15); do 
            taxfile="exp9_indexes/num_cols_${{col}}/index_stats.txt"
            kfile="exp9_indexes/topk_num_cols_${{col}}/index_stats.txt"
            
            table=$(cat $taxfile | awk -F, '{{print $2}}')
            overflow=$(cat $taxfile | awk -F, '{{print $3}}')
            ktable=$(head -n1 $kfile | awk '{{print $1}}')
            
            printf "%d,%d,%d,%d\n" $col $table $overflow $ktable >> {output}
        done
        cat exp9_indexes/full_structure/index_stats.txt >> {output}
        """
