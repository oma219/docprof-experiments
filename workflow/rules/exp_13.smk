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

# Section 2.3: Gather all the time and memory results

rule gather_data_of_two_pass_exp13:
    input:
        expand("exp13_index/{num}_genomes/time_and_mem.log", num=range(2,11))
    output:
        "exp13_results/output.csv"
    shell:
        """
        printf "num,n,time,pass1,pass2,maxrss,tmpmem\n" > {output}
        for path in {input}; do
            num=$(echo $path | grep -oP '[0-9]*_genomes' | grep -oP '[0-9]*')
            time=$(cat $path | awk '{{print $6}}')
            raw_ram=$(cat $path | awk '{{print $10}}')
            ram=$(printf "%.2f" $(echo $raw_ram/1048576 | bc -l))

            firstpass=$(cat "exp13_index/${{num}}_genomes/build.log" | grep "1st-pass" | grep -oP '\([0-9]*.[0-9]* ')
            secondpass=$(cat "exp13_index/${{num}}_genomes/build.log" | grep "2nd pass" | grep -oP '\([0-9]*.[0-9]* ')
            firstpass=${{firstpass:1}}
            secondpass=${{secondpass:1}}

            tmp_mem=$(cat "exp13_index/${{num}}_genomes/build.log" | grep "stats:" | awk '{{print $11}}')
            tmp_mem=$(printf "%.2f" $(echo $tmp_mem/1073741824 | bc -l))

            # Get n, remove command at end, convert to GB
            n=$(cat "exp13_index/${{num}}_genomes/build.log" | grep "stats:" | awk '{{print $5}}')
            n=${{n::-1}}
            n=$(printf "%.2f" $(echo "$n/1073741824" | bc -l))

            printf "%d,%.2f,%.2f,%.3f,%.3f,%.2f,%.2f\n" $num $n $time $firstpass $secondpass $ram $tmp_mem >> {output}
        done
        """