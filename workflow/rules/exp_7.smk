##################################################
# Name: exp_7.smk
# Description: Contains the workflow and methods
#              needed for experiment 7.
#
#              Generates ANI measures between all
#              pairs of genomes in dataset
#
# Date: October 13, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_group_lists_exp7(wildcards):
    """ Returns a list of genomes given a dataset number """
    file_list = []
    for i in range(1, num_classes_exp7+1):
        for data_file in os.listdir(f"data/dataset_{i}"):
            if data_file.endswith(".fna"):
                file_list.append(f"data/dataset_{i}/" + data_file)
    return file_list

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Build a file-list to feed into fastANI

rule build_filelist_exp7:
    input:
        get_group_lists_exp7
    output:
        "exp7_filelist/input_list.txt"
    run:
        with open(output[0], "w") as out_fd:
            for path in input:
                out_fd.write(f"{path}\n")


# Section 2.2: Generate the ANI values for this dataset

rule compute_ANI_with_fastANI_exp7:
    input:
        "exp7_filelist/input_list.txt"
    output:
        "exp7_results/fastANI_output.txt"
    shell:
        """
        fastANI -t 16 --ql {input} --rl {input} -o {output}
        """

# Section 2.3: Analyze the output from fastANI, only consider
# comparisons that are between different datasets and write
# that to a csv file

rule analyze_ANI_from_fastANI_exp7:
    input:
        "exp7_results/fastANI_output.txt"
    output:
        "exp7_results/output.csv"
    run:
        with open(output[0], "w") as out_fd:
            out_fd.write("name,ANI\n")

            # Go through fastANI output and only consider relevant
            # lines
            with open(input[0], "r") as input_fd:
                for line in input_fd:
                    line_split = line.split()
                    dataset_1 = line_split[0].split("/")[1]
                    dataset_2 = line_split[1].split("/")[1]

                    # Across classes ...
                    if dataset_1 != dataset_2:
                        out_fd.write(f"{dataset_name_exp7},{line_split[2]}\n")
