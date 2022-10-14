##################################################
# Name: exp_5.smk
# Description: Contains the workflow and methods
#              needed for experiment 5.
#
#              Compares read classification of 
#              SPUMONI 2 to document array profiles
#              on simulated reads from different datasets.
#
# Date: October 2, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_group_lists_exp5(wildcards):
    """ Returns a list of genomes given a dataset number """
    file_list = []
    for data_file in os.listdir(f"data/dataset_{wildcards.num}"):
        if data_file.endswith(".fna"):
            file_list.append(f"data/dataset_{wildcards.num}/" + data_file)
    return file_list

def list_of_genomes_to_index_exp5(wildcards):
    """
    Returns all the paths to files to be included
    in the index.
    """
    input_files = []
    for num in range(1, num_classes_exp5+1):
        for data_file in os.listdir(f"data/dataset_{num}"):
            if data_file.endswith(".fna"):
                input_files.append(f"{base_dir}/data/dataset_{num}/{data_file}")
    return input_files

####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Build a file-list to be used for SPUMONI 2 and pfp_doc
# to build their respective indexes.

rule build_filelist_for_tools_exp5:
    input:
        list_of_genomes_to_index_exp5
    output:
        "exp5_filelists/index_filelist.txt"
    run:
        with open(output[0], "w") as out_fd:
            for num in range(1, num_classes_exp5+1):
                for data_file in os.listdir(f"data/dataset_{num}"):
                    if data_file.endswith(".fna"):
                        out_fd.write(f"{base_dir}/data/dataset_{num}/{data_file} {num}\n")


# Section 2.2: Simulate long reads from a random genome from each dataset, and 
# subset the reads to match the number of reads required per class.

rule generate_raw_positive_long_reads_exp5:
    input:
        get_group_lists_exp5
    output:
        "exp5_raw_reads/dataset_{num}/dataset_{num}_reads.fastq"
    shell:
        """
        set +o pipefail;
        positive_genome=$(ls data/dataset_{wildcards.num}/*.fna | shuf | head -n1)

        pbsim --depth 95.0 --prefix exp5_raw_reads/dataset_{wildcards.num}/dataset_{wildcards.num}_reads \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 4800 --length-max 5200 $positive_genome

        cat 'exp5_raw_reads/dataset_{wildcards.num}/dataset_{wildcards.num}_reads'*.fastq > {output}
        ls  'exp5_raw_reads/dataset_{wildcards.num}/dataset_{wildcards.num}_reads_'*.fastq | xargs rm
        ls  'exp5_raw_reads/dataset_{wildcards.num}/dataset_{wildcards.num}_reads'*.ref | xargs rm
        ls  'exp5_raw_reads/dataset_{wildcards.num}/dataset_{wildcards.num}_reads'*.maf | xargs rm
        """

rule convert_long_reads_to_fasta_and_subset_exp5:
    input:
        "exp5_raw_reads/dataset_{num}/dataset_{num}_reads.fastq"
    output:
        "exp5_reads/dataset_{num}_reads.fa"
    shell:
        """
        num_lines=$(({num_reads_per_class_exp5} * 4))
        head -n $num_lines {input} > exp5_raw_reads/dataset_{wildcards.num}/dataset_{wildcards.num}_subset_reads.fastq

        seqtk seq -a exp5_raw_reads/dataset_{wildcards.num}/dataset_{wildcards.num}_subset_reads.fastq > {output}

        if [ $(grep -c '>' {output}) != {num_reads_per_class_exp5} ]; then echo 'number of reads assertion failed.'; exit 1; fi
        rm exp5_raw_reads/dataset_{wildcards.num}/dataset_{wildcards.num}_subset_reads.fastq
        """

# Section 2.5: Build index for both SPUMONI 2 and the document 
# array profiles

rule build_spumoni_index_for_exp5:
    input:
        "exp5_filelists/index_filelist.txt"
    output:
        "exp5_indexes/spumoni_index/spumoni_full_ref.fa.thrbv.ms"
    shell:
        """
        spumoni build -i {input} -b exp5_indexes/spumoni_index/ -n -M -d
        """

rule build_doc_array_profiles_for_listing_exp5:
    input:
        "exp5_filelists/index_filelist.txt"
    output:
        "exp5_indexes/docprofiles_index/full_ref.fna.sdap",
        "exp5_indexes/docprofiles_index/full_ref.fna.edap"
    shell:
        """
        pfp_doc build -f {input} -o exp5_indexes/docprofiles_index/full_ref -r
        """

# Section 2.6: Run SPUMONI on all the reads in order to extract MEMs
# and compute read classifications for SPUMONI

rule run_query_with_spumoni_for_doc_labels_exp5:
    input:
        "exp5_indexes/spumoni_index/spumoni_full_ref.fa.thrbv.ms",
        "exp5_reads/dataset_{num}_reads.fa"
    output:
        "exp5_results/spumoni/ms_results/class_{num}_reads.fna",
        "exp5_results/spumoni/ms_results/class_{num}_reads.fna.lengths",
        "exp5_results/spumoni/ms_results/class_{num}_reads.fna.doc_numbers"
    shell:
        """
        cp {input[1]} {output[0]}
        spumoni run -r exp5_indexes/spumoni_index/spumoni_full_ref.fa -p {output[0]} -M -d -n
        """

# Section 2.7: Extract MEMs from each class of reads using the matching statistics
# computed by SPUMONI to feed into document array profile, and query each MEM

rule extract_mems_based_on_ms_exp5:
    input:
        "exp5_results/spumoni/ms_results/class_{num}_reads.fna",
        "exp5_results/spumoni/ms_results/class_{num}_reads.fna.lengths"
    output:
        "exp5_results/doclist/ms_results/class_{num}_reads.fna",
        "exp5_results/doclist/ms_results/class_{num}_reads.fna.lengths",
        "exp5_results/doclist/mem_file/class_{num}_mems.fastq"
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

rule determine_doc_array_listings_for_mems_exp5:
    input:
        "exp5_indexes/docprofiles_index/full_ref.fna.sdap",
        "exp5_results/doclist/mem_file/class_{num}_mems.fastq"
    output:
        "exp5_results/doclist/listings/class_{num}_reads.fastq",
        "exp5_results/doclist/listings/class_{num}_reads.fastq.listings"
    shell:
        """
        cp {input[1]} {output[0]}
        pfp_doc run -r exp5_indexes/docprofiles_index/full_ref -p {output[0]}
        """


# Section 2.8: Generate confusion matix for SPUMONI and document array
# profile separately by analyzing the reads from each class.

rule classify_spumoni_reads_exp5:
    input:
        expand("exp5_results/spumoni/ms_results/class_{num}_reads.fna.doc_numbers", num=range(1, num_classes_exp5+1)),
        expand("exp5_results/spumoni/ms_results/class_{num}_reads.fna.lengths", num=range(1, num_classes_exp5+1))
    output:
        "exp5_results/spumoni/csv_files/classification_results.csv"
    run:
        def calculate_accuracy_values(confusion_matrix, num_datasets):
            """ Calculates accuracy score given confusion matrix """
            accuracies = []
            for pivot in range(num_datasets):
                tp = confusion_matrix[pivot][pivot]
                fp = fn = tn = 0
                for row in range(num_datasets):
                    for column in range(num_datasets):
                        curr = confusion_matrix[row][column]
                        if(column == pivot and row != pivot):
                            fp += curr
                        elif(row == pivot and column != pivot):
                            fn += curr
                        elif(row != pivot):
                            tn += curr
                accuracies.append([pivot,tp,tn,fp,fn])
            return accuracies

        def get_class_counts_for_read(doc_list, length_list):
            """ Takes in document labels for read and return counts """
            assert len(doc_list) == len(length_list), "mismatch in list lengths!"
            counts = [0 for i in range(num_classes_exp5)]

            prev_length = -1
            for doc_num, length_num in zip(doc_list, length_list):
                if length_num >= 15 and length_num >= prev_length:
                    counts[doc_num] += length_num - 14
                prev_length = length_num
            return counts

        confusion_matrix = [[0 for x in range(num_classes_exp5)] for i in range(num_classes_exp5)]
        
        # Goes through each *.doc_number files
        for class_num, input_file in enumerate(input[:num_classes_exp5]):
            # Open both the *.doc_number and *.lengths file for certain class
            with open(input_file, "r") as input_fd, open(input[class_num+num_classes_exp5], "r") as lengths_fd:
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
            out_fd.write("tool,dataset,tp,tn,fp,fn\n")
            for results in calculate_accuracy_values(confusion_matrix, num_classes_exp5):
                out_fd.write(",".join(["spumoni"] + [str(x) for x in results]) + "\n")


rule classify_docarray_profile_reads_exp5:
    input:
        expand("exp5_results/doclist/listings/class_{num}_reads.fastq.listings", num=range(1, num_classes_exp5+1))
    output:
        "exp5_results/doclist/csv_files/classification_results.csv"
    run:
        def calculate_accuracy_values(confusion_matrix, num_datasets):
            """ Calculates accuracy score given confusion matrix """
            accuracies = []
            for pivot in range(num_datasets):
                tp = confusion_matrix[pivot][pivot]
                fp = fn = tn = 0
                for row in range(num_datasets):
                    for column in range(num_datasets):
                        curr = confusion_matrix[row][column]
                        if(column == pivot and row != pivot):
                            fp += curr
                        elif(row == pivot and column != pivot):
                            fn += curr
                        elif(row != pivot):
                            tn += curr
                accuracies.append([pivot,tp,tn,fp,fn])
            return accuracies
        
        confusion_matrix = [[0 for x in range(num_classes_exp5)] for i in range(num_classes_exp5)]

        for class_num, input_file in enumerate(input):
            with open(input_file, "r") as input_fd:
                curr_read_num = 0
                curr_mem_length = 0
                curr_read_stats = [0 for i in range(num_classes_exp5)]
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
                            curr_read_stats = [0 for i in range(num_classes_exp5)]
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
                curr_read_stats = [0 for i in range(num_classes_exp5)]
        
        # Write out the csv data file ...
        with open(output[0], "w") as out_fd:
            out_fd.write("tool,dataset,tp,tn,fp,fn\n")
            for results in calculate_accuracy_values(confusion_matrix, num_classes_exp5):
                out_fd.write(",".join(["docprofile"] + [str(x) for x in results]) + "\n")