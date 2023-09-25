##################################################
# Name: exp_12.smk
# Description: Contains the workflow and methods
#              needed for experiment 10.
#
#              This experiment compares the construction
#              time for document array profiles when 
#              using heuristic/no-heuristic/2-pass 
#              algorithm.
#
# Date: Sept 24th, 2023
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_all_genomes_for_exp12(wildcards):
    """ Returns all the genomes in the data/ folder """
    file_list = []
    for data_file in os.listdir(f"exp12_data/"):
        if data_file.endswith(".fna"):
            file_list.append(f"exp12_data/" + data_file)
    return file_list

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Create a filelist of different sizes

rule create_filelists_for_different_sizes_exp12:
    input:
        get_all_genomes_for_exp12
    output:
        "exp12_filelist/filelist_{num}_genomes.txt"
    run:
        with open(output[0], "w") as out_fd:
            for i, path in enumerate(input):
                out_fd.write(f"{path} {(i+1)}\n")
                if (i+1) == int(wildcards.num):
                    break

# Section 2.2: Build the index in different ways

rule build_index_with_heuristics_exp12:
    input:
        "exp12_filelist/filelist_{num}_genomes.txt"
    output:
        "exp12_index/heuristic/{num}_genomes/filelist.txt",
        "exp12_index/heuristic/{num}_genomes/build.log",
        "exp12_index/heuristic/{num}_genomes/output.fna.sdap",
        "exp12_index/heuristic/{num}_genomes/output.fna.edap"
    shell:
        """
        cp {input} {output[0]}
        pfp_doc64 build --filelist {output[0]} \
                        --output exp12_index/heuristic/{wildcards.num}_genomes/output 2> {output[1]}
        """

rule build_index_with_no_heuristics_exp12:
    input:
        "exp12_filelist/filelist_{num}_genomes.txt"
    output:
        "exp12_index/no_heuristic/{num}_genomes/filelist.txt",
        "exp12_index/no_heuristic/{num}_genomes/build.log",
        "exp12_index/no_heuristic/{num}_genomes/output.fna.sdap",
        "exp12_index/no_heuristic/{num}_genomes/output.fna.edap"
    shell:
        """
        cp {input} {output[0]}
        pfp_doc64 build --filelist {output[0]} \
                        --no-heuristic \
                        --output exp12_index/no_heuristic/{wildcards.num}_genomes/output 2> {output[1]}
        """

rule build_index_with_two_pass_exp12:
    input:
        "exp12_filelist/filelist_{num}_genomes.txt"
    output:
        "exp12_index/two_pass/{num}_genomes/filelist.txt",
        "exp12_index/two_pass/{num}_genomes/build.log",
        "exp12_index/two_pass/{num}_genomes/output.fna.sdap",
        "exp12_index/two_pass/{num}_genomes/output.fna.edap"
    shell:
        """
        cp {input} {output[0]}
        pfp_doc64 build --filelist {output[0]} \
                        --two-pass ./temp \
                        --tmp-size 10GB \
                        --output exp12_index/two_pass/{wildcards.num}_genomes/output 2> {output[1]}
        """

# Section 2.3: Combine all the results into one file

rule compile_results_from_different_methods_exp12:
    input:
        expand("exp12_index/{type}/{num}_genomes/build.log", 
               type=["heuristic", "no_heuristic","two_pass"],
               num=[2,3,4,5,10,15,20,25,30])
    output:
        "exp12_results/output.csv"
    shell:
        """
        printf "num,heuristic,noheuristic,twopasstotal,firstpass,secondpass\n" > {output}
        for num in 2 3 4 5 10 15 20 25 30; do
            heur=$(grep 'finished' "exp12_index/heuristic/${{num}}_genomes/build.log" | awk '{{print $NF}}')
            no_heur=$(grep 'finished' "exp12_index/no_heuristic/${{num}}_genomes/build.log" | awk '{{print $NF}}')
            tpass=$(grep 'finished' "exp12_index/two_pass/${{num}}_genomes/build.log" | awk '{{print $NF}}')
            onepass=$(grep "1st-pass" "exp12_index/two_pass/${{num}}_genomes/build.log" | awk '{{print substr($(NF-1), 2)}}')
            twopass=$(grep "2nd pass" "exp12_index/two_pass/${{num}}_genomes/build.log" | awk '{{print substr($(NF-1), 2)}}')
            printf "%d,%.2f,%.2f,%.2f,%.2f,%.2f\n" $num $heur $no_heur $tpass $onepass $twopass >> {output}
        done
        """



            
