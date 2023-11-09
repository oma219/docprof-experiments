##################################################
# Name: exp_14.smk
# Description: Contains the workflow and methods
#              needed for experiment 14.
#
#              This experiment uses a newer
#              script to separate out the genera
#              in the SILVA database and then 
#              build document array profiles over it.
#
# Date: Nov 9th, 2023
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_input_fasta_file_exp14(wildcards):
    """ Returns path to SILVA database with taxonomic path in each header """
    file_list = []
    for data_file in os.listdir(f"exp14_data/"):
        if data_file.endswith(".fasta") and data_file.startswith("SILVA"):
            file_list.append(f"exp14_data/" + data_file)
    assert len(file_list) == 1
    return file_list[0]

####################################################
# Section 2: Rules needed for this experiment type
####################################################

rule convert_to_uppercase_and_samelines_exp14:
    input:
        get_input_fasta_file_exp14
    output:
        "exp14_silva_database/silva_database.fa"
    shell:
        "seqtk seq -U {input} > {output}"

rule generate_separate_fasta_files_exp14:
    input:
        "exp14_silva_database/silva_database.fa"
    output:
        "exp14_input_files/filelist.txt"
    shell:
        """
        python3 {repo_dir}/src/separate_silva_genera_v2.py \
        -i {input} \
        -o exp14_input_files/ \
        --tree exp14_data/tax_slv_ssu_138.1.tre \
        --tree-map exp14_data/tax_slv_ssu_138.1.map \
        --tree-rank exp14_data/tax_slv_ssu_138.1.txt \
        -n {num_genera_exp14}
        """

rule generate_dna_version_of_fasta_files_exp14:
    input:
        "exp14_input_files/filelist.txt"
    output:
        "exp14_input_files_DNA/filelist.txt"
    shell:
        """
        work_dir=$(pwd)
        files=(exp14_input_files/*.fa)
        num_docs=${{#files[@]}}
        
        touch {output}
        for i in $(seq 1 $num_docs); do
            input_file="${{work_dir}}/exp14_input_files/doc_${{i}}_seq.fa"
            output_file="${{work_dir}}/exp14_input_files_DNA/doc_${{i}}_seq.fa"

            seqtk seq -r $input_file > $output_file
            printf "%s %d\n" $output_file $i >> {output}
        done
        """

rule build_index_with_two_pass_exp14:
    input:
        "exp14_input_files_DNA/filelist.txt"
    output:
        "exp14_index/filelist.txt",
        "exp14_index/build.log",
        "exp14_index/output.fna.taxcomp.sdap",
        "exp14_index/output.fna.taxcomp.edap",
        "exp14_index/output.fna.taxcomp.of.sdap",
        "exp14_index/output.fna.taxcomp.of.edap",
        "exp14_index/time_and_mem.log",
        "exp14_index/read_and_write.log"
    shell:
        """
        cp {input} {output[0]}
        {time_prog} {time_format} --output={output[6]} \
        {strace_prog} {strace_args} -o {output[7]} \
        pfp_doc64 build --filelist {output[0]} \
                        --revcomp \
                        --taxcomp \
                        --num-col 7 \
                        --two-pass exp14_index/temp  \
                        --tmp-size {tmp_mem_used_exp14} \
                        --output exp14_index/output 2> {output[1]}
        """