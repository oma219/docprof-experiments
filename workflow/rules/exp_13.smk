##################################################
# Name: exp_13.smk
# Description: Contains the workflow and methods
#              needed for experiment 10.
#
#              This experiment builds progressively
#              larger document array profiles 
#              over a group of FASTAS provided
#              and profiles the time/memory.
#
# Date: Oct 14th, 2023
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_all_genomes_for_exp13(wildcards):
    """ Returns all the genomes in the data/ folder """
    file_list = []
    for data_file in os.listdir(f"exp13_data/"):
        if data_file.endswith(".fna") or data_file.endswith(".fa"):
            file_list.append(f"exp13_data/" + data_file)
    return file_list

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Create a filelist of different sizes

rule create_filelists_for_different_sizes_exp13:
    input:
        get_all_genomes_for_exp13
    output:
        "exp13_filelist/filelist_{num}_genomes.txt"
    run:
        with open(output[0], "w") as out_fd:
            for i, path in enumerate(input):
                out_fd.write(f"{path} {(i+1)}\n")
                if (i+1) == int(wildcards.num):
                    break

# Section 2.2: Build the index in different ways

rule build_index_with_two_pass_exp13:
    input:
        "exp13_filelist/filelist_{num}_genomes.txt"
    output:
        "exp13_index/{num}_genomes/filelist.txt",
        "exp13_index/{num}_genomes/build.log",
        "exp13_index/{num}_genomes/output.fna.sdap",
        "exp13_index/{num}_genomes/output.fna.edap",
        "exp13_index/{num}_genomes/time_and_mem.log",
        "exp13_index/{num}_genomes/read_and_write.log"
    shell:
        """
        cp {input} {output[0]}
        {time_prog} {time_format} --output={output[4]} \
        {strace_prog} {strace_args} -o {output[5]} \
        pfp_doc64 build --filelist {output[0]} \
                        --two-pass exp13_index/{wildcards.num}_genomes/temp  \
                        --tmp-size 400GB \
                        --output exp13_index/{wildcards.num}_genomes/output 2> {output[1]}
        """