##################################################
# Name: exp_6.smk
# Description: Contains the workflow and methods
#              needed for experiment 5.
#
#              Compares read classification of 
#              SPUMONI 2 to document array profiles
#              on real mock community reads.
#
# Date: October 13, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_group_lists_exp6(wildcards):
    """ Returns a list of genomes given a dataset number """
    file_list = []
    # Add the other strains ...
    for data_file in os.listdir(f"data/dataset_{wildcards.num}"):
        if data_file.endswith(".fna"):
            file_list.append(f"data/dataset_{wildcards.num}/" + data_file)
    # Add the true strain ...
    for data_file in os.listdir(f"data/class_{wildcards.num}/true_strain"):
        if data_file.endswith(".fasta"):
            file_list.append(f"data/class_{wildcards.num}/true_strain/" + data_file)
    return file_list

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Build a file-list to be used for SPUMONI 2 and pfp_doc
# to build their respective indexes. One index for each
# each class.

rule build_filelist_for_each_species_exp6:
    input:
        get_group_lists_exp6
    output:
        "exp6_filelists/class_{num}_filelist.txt"
    run:
        with open(output[0], "w") as out_fd:
            for i, strain_genome in enumerate(input):
                out_fd.write(f"{strain_genome} {i+1}\n")

# Section 2.2: Copy reads to working and directory and subset them
# during testing.

rule copy_read_files_for_each_species_exp6:
    input:
        "data/class_{num}/reads/class_{num}_reads.fa"
    output:
        "exp6_read_sets/class_{num}_reads.fa"
    shell:
        """
        cp {input} {output}
        """

# Section 2.3: Build index for both SPUMONI 2 and the document 
# array profiles for each species

rule build_spumoni_index_for_exp6:
    input:
        "exp6_filelists/class_{num}_filelist.txt"
    output:
        "exp6_indexes/class_{num}/spumoni_index/spumoni_full_ref.fa.thrbv.ms"
    shell:
        """
        spumoni build -i {input} -b exp6_indexes/class_{wildcards.num}/spumoni_index/ -n -M -d
        """

rule build_doc_array_profiles_for_listing_exp6:
    input:
        "exp6_filelists/class_{num}_filelist.txt"
    output:
        "exp6_indexes/class_{num}/docprofiles_index/full_ref.fna.sdap",
        "exp6_indexes/class_{num}/docprofiles_index/full_ref.fna.edap"
    shell:
        """
        pfp_doc build -f {input} -o exp6_indexes/class_{wildcards.num}/docprofiles_index/full_ref -r
        """

# Section 2.4: Run SPUMONI on all the reads in order to extract MEMs
# and compute read classifications for SPUMONI

rule run_query_with_spumoni_for_doc_labels_exp6:
    input:
        "exp6_indexes/class_{num}/spumoni_index/spumoni_full_ref.fa.thrbv.ms",
        "exp6_read_sets/class_{num}_reads.fa"
    output:
        "exp6_results/class_{num}/spumoni/ms_results/class_{num}_reads.fna",
        "exp6_results/class_{num}/spumoni/ms_results/class_{num}_reads.fna.lengths",
        "exp6_results/class_{num}/spumoni/ms_results/class_{num}_reads.fna.doc_numbers"
    shell:
        """
        cp {input[1]} {output[0]}
        spumoni run -r exp6_indexes/class_{wildcards.num}/spumoni_index/spumoni_full_ref.fa -p {output[0]} -M -d -n
        """

# Section 2.5: Extract MEMs from each class of reads using the matching statistics
# computed by SPUMONI to feed into document array profile, and query each MEM

rule extract_mems_based_on_ms_exp6:
    input:
        "exp6_results/class_{num}/spumoni/ms_results/class_{num}_reads.fna",
        "exp6_results/class_{num}/spumoni/ms_results/class_{num}_reads.fna.lengths",
    output:
        "exp6_results/class_{num}/doclist/ms_results/class_{num}_reads.fna",
        "exp6_results/class_{num}/doclist/ms_results/class_{num}_reads.fna.lengths",
        "exp6_results/class_{num}/doclist/mem_file/class_{num}_mems.fastq"
    shell:
        """
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        python3 {repo_dir}/src/extract_mems_exp3.py \
        -l {output[1]} \
        -p {output[0]} \
        --mems \
        -t 0 \
        -o {output[2]} 
        """

rule determine_doc_array_listings_for_mems_exp6:
    input:
        "exp6_indexes/class_{num}/docprofiles_index/full_ref.fna.sdap",
        "exp6_results/class_{num}/doclist/mem_file/class_{num}_mems.fastq"
    output:
        "exp6_results/class_{num}/doclist/listings/class_{num}_reads.fastq",
        "exp6_results/class_{num}/doclist/listings/class_{num}_reads.fastq.listings"
    shell:
        """
        cp {input[1]} {output[0]}
        pfp_doc run -r exp6_indexes/class_{wildcards.num}/docprofiles_index/full_ref -p {output[0]}
        """

# Section 2.6: Generate confusion matix for SPUMONI and document array
# profile separately by analyzing the reads from each class.

rule classify_spumoni_reads_exp6:
    input:
        expand("exp6_results/class_{num}/spumoni/ms_results/class_{num}_reads.fna.doc_numbers", num=range(1, num_classes_exp6+1)),
        expand("exp6_results/class_{num}/spumoni/ms_results/class_{num}_reads.fna.lengths", num=range(1, num_classes_exp6+1))
    output:
        "exp6_final_results/spumoni.csv"
    run:
        def get_class_counts_for_read(doc_list, length_list):
            """ Takes in document labels for read and return counts """
            assert len(doc_list) == len(length_list), "mismatch in list lengths!"
            counts = [0 for i in range(num_strains_exp6)]

            prev_length = -1
            for doc_num, length_num in zip(doc_list, length_list):
                if length_num >= 15 and length_num >= prev_length:
                    counts[doc_num] += length_num - 14
                prev_length = length_num
            return counts

        confusion_matrix = [[0 for x in range(num_strains_exp6)] for i in range(num_classes_exp6)]
        
        # Goes through each *.doc_number files
        for class_num, input_file in enumerate(input[:num_classes_exp6]):
            # Open both the *.doc_number and *.lengths file for certain class
            with open(input_file, "r") as input_fd, open(input[class_num+num_classes_exp6], "r") as lengths_fd:
                for line, length_line in zip(input_fd, lengths_fd):
                    if '>' not in line:
                        counts = get_class_counts_for_read([int(x) for x in line.split()], [int(x) for x in length_line.split()])
                        max_count = max(counts)
                        #pred_class = counts.index(max_count)

                        # This section just handles case where two classes have same weight
                        pred_class = -1
                        pred_class_list = [i for i in range(len(counts)) if counts[i] == max_count]
                        if len(pred_class_list) > 1:
                            pred_class = random.choice(pred_class_list)
                        else:
                            pred_class = pred_class_list[0]

                        confusion_matrix[class_num][pred_class] += 1
        
        # Write out the csv data file ...
        with open(output[0], "w") as out_fd:
            out_fd.write("tool,dataset,tp,fn\n")
            for i, curr_row in enumerate(confusion_matrix):
                tp = curr_row[-1]
                fn = sum(curr_row) - tp
                results = [tp, fn]
                out_fd.write(",".join(["spumoni", str(i)] + [str(x) for x in results]) + "\n")


rule classify_docarray_profile_reads_exp6:
    input:
        expand("exp6_results/class_{num}/doclist/listings/class_{num}_reads.fastq.listings", num=range(1, num_classes_exp6+1))
    output:
        "exp6_final_results/doclist.csv"
    run:        
        confusion_matrix = [[0 for x in range(num_strains_exp6)] for i in range(num_classes_exp6)]

        for class_num, input_file in enumerate(input):
            with open(input_file, "r") as input_fd:
                curr_read_num = 0
                curr_mem_length = 0
                curr_read_stats = [0 for i in range(num_strains_exp6)]
                for line in input_fd:
                    if '>' in line:
                        line_split = line.strip().split("_")
                        curr_mem_length = int(line_split[5])
                        
                        # Need to classify the read now ...
                        if int(line_split[1]) != curr_read_num:
                            pred_class_weight = max(curr_read_stats)
                            #pred_class = curr_read_stats.index(pred_class_weight)

                            # This section just handles case where two classes have same weight
                            pred_class = -1
                            pred_class_list = [i for i in range(len(curr_read_stats)) if curr_read_stats[i] == pred_class_weight]
                            if len(pred_class_list) > 1:
                                pred_class = random.choice(pred_class_list)
                            else:
                                pred_class = pred_class_list[0]
                            confusion_matrix[class_num][pred_class] += 1
                            
                            curr_read_num += 1
                            curr_read_stats = [0 for i in range(num_strains_exp6)]
                    else:
                        line_split = line.strip().split()
                        classes = [int(x) for x in line_split[1][1:-1].split(",")]
                        assert len(classes) > 0, "no classes found for MEMs"

                        # Weighting classification by MEM length
                        if curr_mem_length >= 15:
                            for class_id in classes:
                                curr_read_stats[class_id] += (curr_mem_length - 14) * (1.0/len(classes))
                
                # Classify the last read ...
                pred_class_weight = max(curr_read_stats)
                pred_class = curr_read_stats.index(pred_class_weight)
                confusion_matrix[class_num][pred_class] += 1
                
                curr_read_num += 1
                curr_read_stats = [0 for i in range(num_strains_exp6)]
        
        # Write out the csv data file ...
        with open(output[0], "w") as out_fd:
            out_fd.write("tool,dataset,tp,fn\n")
            for i, curr_row in enumerate(confusion_matrix):
                tp = curr_row[-1]
                fn = sum(curr_row) - tp
                results = [tp, fn]
                out_fd.write(",".join(["docprofile", str(i)] + [str(x) for x in results]) + "\n")