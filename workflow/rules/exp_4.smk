##################################################
# Name: exp_4.smk
# Description: Contains the workflow and methods
#              needed for experiment 4.
#
#              Compares read classification of 
#              SPUMONI 2 to document array profiles
#              on real mock community reads.
#
# Date: October 11, 2022
##################################################

####################################################
# Section 1: Helper functions needed for these
#            experiment rules.
####################################################

def get_readset_for_exp4(wildcards):
    """
    Return a path to the read-set with real ONT reads
    from a Mock Community. This experiment assumes we
    only have 1 FASTQ file to process. 
    """
    file_list = []
    for data_file in os.listdir(f"data/full_read_set"):
        if data_file.endswith(".fastq"):
            file_list.append(f"{base_dir}/data/full_read_set/" + data_file)
    
    if len(file_list) > 1:
        print("Error: Expected only one input Fastq file.")
        exit(-1)

    return file_list

def get_all_zymo_genomes_exp4(wildcards):
    """
    Return paths to all eight of the genomes that were truly
    sequenced in the Zymo Mock Community experiment.
    """
    file_list = []
    for data_file in os.listdir(f"data/true_set"):
        if data_file.endswith(".fasta"):
            file_list.append(f"{base_dir}/data/true_set/" + data_file)

    if len(file_list) != 8:
        print("Error: Expected incorrect number of genomes found.")
        exit(-1)
    
    # Move the yeast to be last in list since it important for later rules
    yeast_file = f"{base_dir}/data/true_set/Saccharomyces_cerevisiae_complete_genome.fasta"
    assert yeast_file in file_list

    file_list.remove(yeast_file)
    file_list.append(yeast_file)

    return file_list

def list_of_genomes_to_index_exp4(wildcards):
    """
    Returns all the paths to files to be included
    in the index.
    """
    input_files = []
    for num in range(1, num_classes_exp4+1):
        for data_file in os.listdir(f"data/dataset_{num}"):
            if data_file.endswith(".fna"):
                input_files.append(f"{base_dir}/data/dataset_{num}/{data_file}")
    return input_files


####################################################
# Section 2: Rules needed for this experiment type
####################################################

# Section 2.1: Sub-sample ONT reads that are long enough to classify confidently (>5000 bp) and
#              trim the beginning of the ONT to exclude the bases that will be used for experiment

rule sample_from_full_ONT_reads_exp4:
    input:
        get_readset_for_exp4
    output:
        "exp4_intermediate/step_1/sampled_reads.fastq"
    shell:
        """
        seqtk seq -L 5000 {input[0]} > {output}
        """

# Section 2.2: Build a multi-FASTA file for the 8 genomes in the Zymo Mock Community and build
#              a minimap2 index for it.

rule build_ZymoMC_ref_exp4:
    input:
        get_all_zymo_genomes_exp4
    output:
        "exp4_intermediate/step_3/zymo_mc.fa"
    shell:
        """
        cat {input} > {output}
        """

rule build_minimap2_index_over_zymo_mc_exp4:
    input:
        "exp4_intermediate/step_3/zymo_mc.fa"
    output:
        "exp4_intermediate/step_4/zymo_mc.mmi"
    shell:
        """
        minimap2 -x map-ont -d {output} {input}
        """

# Section 2.3: Align all of the reads to develop the gold standard classifications, 
#              and filter them based on MAPQ. Generate files that contain the 
#              reference headers to separate out the microbial reads versus yeast 
#              reads.

rule align_all_reads_to_zymo_mc_exp4:
    input:
        "exp4_intermediate/step_1/sampled_reads.fastq",
        "exp4_intermediate/step_4/zymo_mc.mmi"
    output:
        "exp4_intermediate/step_5/zymo_mc_read_aln.sam"
    shell:
        """
        minimap2 -p 0.6 -a {input[1]} {input[0]} > {output}
        """

rule filter_all_reads_by_mapq_exp4:
    input:
        "exp4_intermediate/step_5/zymo_mc_read_aln.sam"
    output:
        "exp4_intermediate/step_6/zymo_mc_read_aln.filtered.sam",
        "exp4_intermediate/step_6/zymo_mc_read_aln.filtered.sorted.sam"
    shell:
        """
        samtools view -u -q 30 {input} > {output[0]}
        samtools sort {output[0]} -o {output[1]}
        """

rule generate_reference_seq_name_list_exp4:
    input:
        get_all_zymo_genomes_exp4
    output:
        expand("exp4_ref_name_lists/class_{n}.txt", n=range(1, 9)),
        "exp4_ref_name_lists/README"
    shell:
        """
        # IMPORTANT: this produces the class lists in the same order as the
        # database genomes to make sure the class numbers are aligned.
        i=1
        for genus in $(tail -n +2 data/README_dataset_summary.txt | awk '{{print $2}}')
        do
            genome_path=$(realpath data/true_set/*.fasta | grep "$genus")
            grep '^>' $genome_path | awk '{{print substr($1,2)}}' > "exp4_ref_name_lists/class_$i.txt"
            echo $genome_path >> "exp4_ref_name_lists/README"
            i=$((i+1))
        done
        """

rule generate_separate_read_files_for_each_class_exp4:
    input:
        expand("exp4_ref_name_lists/class_{n}.txt", n=range(1, 9)),
        "exp4_intermediate/step_6/zymo_mc_read_aln.filtered.sorted.sam"
    output:
        expand("exp4_read_sets/class_{n}_reads.fa", n=range(1, 9))
    shell:
        """
        python3 {repo_dir}/src/classify_reads_sam.py \
        -i exp4_intermediate/step_6/zymo_mc_read_aln.filtered.sorted.sam \
        -c exp4_ref_name_lists/class_{{1..8}}.txt \
        -o exp4_read_sets/ \
        -r exp4_intermediate/step_1/sampled_reads.fastq
        """

# Section 2.4: Build a file-list to be used for SPUMONI 2 and pfp_doc
# to build their respective indexes.

rule build_filelist_for_tools_exp4:
    input:
        list_of_genomes_to_index_exp4
    output:
        "exp4_filelists/index_filelist.txt"
    run:
        with open(output[0], "w") as out_fd:
            for num in range(1, num_classes_exp4+1):
                for data_file in os.listdir(f"data/dataset_{num}"):
                    if data_file.endswith(".fna"):
                        out_fd.write(f"{base_dir}/data/dataset_{num}/{data_file} {num}\n")



# Section 2.5: Build index for both SPUMONI 2 and the document 
# array profiles

rule build_spumoni_index_for_exp4:
    input:
        "exp4_filelists/index_filelist.txt"
    output:
        "exp4_indexes/spumoni_index/spumoni_full_ref.fa.thrbv.ms"
    shell:
        """
        spumoni build -i {input} -b exp4_indexes/spumoni_index/ -n -M -d
        """

rule build_doc_array_profiles_for_listing_exp4:
    input:
        "exp4_filelists/index_filelist.txt"
    output:
        "exp4_indexes/docprofiles_index/full_ref.fna.sdap",
        "exp4_indexes/docprofiles_index/full_ref.fna.edap"
    shell:
        """
        pfp_doc build -f {input} -o exp4_indexes/docprofiles_index/full_ref -r
        """

# Section 2.6: Run SPUMONI on all the reads in order to extract MEMs
# and compute read classifications for SPUMONI

rule run_query_with_spumoni_for_doc_labels_exp4:
    input:
        "exp4_indexes/spumoni_index/spumoni_full_ref.fa.thrbv.ms",
        "exp4_read_sets/class_{num}_reads.fa"
    output:
        "exp4_results/spumoni/ms_results/class_{num}_reads.fna",
        "exp4_results/spumoni/ms_results/class_{num}_reads.fna.lengths",
        "exp4_results/spumoni/ms_results/class_{num}_reads.fna.doc_numbers"
    shell:
        """
        cp {input[1]} {output[0]}
        spumoni run -r exp4_indexes/spumoni_index/spumoni_full_ref.fa -p {output[0]} -M -d -n
        """

# Section 2.7: Extract MEMs from each class of reads using the matching statistics
# computed by SPUMONI to feed into document array profile, and query each MEM

rule extract_mems_based_on_ms_exp4:
    input:
        "exp4_results/spumoni/ms_results/class_{num}_reads.fna",
        "exp4_results/spumoni/ms_results/class_{num}_reads.fna.lengths"
    output:
        "exp4_results/doclist/ms_results/class_{num}_reads.fna",
        "exp4_results/doclist/ms_results/class_{num}_reads.fna.lengths",
        "exp4_results/doclist/mem_file/class_{num}_mems.fastq"
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

rule determine_doc_array_listings_for_mems_exp4:
    input:
        "exp4_indexes/docprofiles_index/full_ref.fna.sdap",
        "exp4_results/doclist/mem_file/class_{num}_mems.fastq"
    output:
        "exp4_results/doclist/listings/class_{num}_reads.fastq",
        "exp4_results/doclist/listings/class_{num}_reads.fastq.listings"
    shell:
        """
        cp {input[1]} {output[0]}
        pfp_doc run -r exp4_indexes/docprofiles_index/full_ref -p {output[0]}
        """

#   Section 2.8: Generate confusion matix for SPUMONI and document array
#   profile separately by analyzing the reads from each class.

rule classify_docarray_profile_reads_exp4:
    input:
        expand("exp4_results/doclist/listings/class_{num}_reads.fastq.listings", num=range(1, num_classes_exp4+1))
    output:
        "exp4_results/doclist/csv_files/classification_results.csv"
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
        
        confusion_matrix = [[0 for x in range(num_classes_exp4)] for i in range(num_classes_exp4)]

        for class_num, input_file in enumerate(input):
            with open(input_file, "r") as input_fd:
                curr_read_num = 0
                curr_mem_length = 0
                curr_read_stats = [0 for i in range(num_classes_exp4)]
                for line in input_fd:
                    if '>' in line:
                        line_split = line.strip().split("_")
                        curr_mem_length = int(line_split[5])
                        
                        # Need to classify the read now ...
                        if int(line_split[1]) != curr_read_num:
                            pred_class_weight = max(curr_read_stats)
                            pred_class = curr_read_stats.index(pred_class_weight)
                            confusion_matrix[class_num][pred_class] += 1
                            
                            curr_read_num += 1
                            curr_read_stats = [0 for i in range(num_classes_exp4)]
                    else:
                        line_split = line.strip().split()
                        classes = [int(x) for x in line_split[1][1:-1].split(",")]
                        assert len(classes) > 0, "no classes found for MEMs"

                        # Weighting classificaiton by MEM length
                        if curr_mem_length >= 15:
                            for class_id in classes:
                                curr_read_stats[class_id] += curr_mem_length
                
                # Classify the last read ...
                pred_class_weight = max(curr_read_stats)
                pred_class = curr_read_stats.index(pred_class_weight)
                confusion_matrix[class_num][pred_class] += 1
                
                curr_read_num += 1
                curr_read_stats = [0 for i in range(num_classes_exp4)]
        
        # Write out the csv data file ...
        with open(output[0], "w") as out_fd:
            out_fd.write("tool,dataset,tp,tn,fp,fn\n")
            for results in calculate_accuracy_values(confusion_matrix, num_classes_exp4):
                out_fd.write(",".join(["docprofile"] + [str(x) for x in results]) + "\n")


rule classify_spumoni_reads_exp4:
    input:
        expand("exp4_results/spumoni/ms_results/class_{num}_reads.fna.doc_numbers", num=range(1, num_classes_exp4+1)),
        expand("exp4_results/spumoni/ms_results/class_{num}_reads.fna.lengths", num=range(1, num_classes_exp4+1))
    output:
        "exp4_results/spumoni/csv_files/classification_results.csv"
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
            counts = [0 for i in range(num_classes_exp4)]

            prev_length = -1
            for doc_num, length_num in zip(doc_list, length_list):
                if length_num >= 15 and length_num >= prev_length:
                    counts[doc_num] += length_num
                prev_length = length_num
            return counts

        confusion_matrix = [[0 for x in range(num_classes_exp4)] for i in range(num_classes_exp4)]
        
        # Goes through each *.doc_number files
        for class_num, input_file in enumerate(input[:num_classes_exp4]):
            # Open both the *.doc_number and *.lengths file for certain class
            with open(input_file, "r") as input_fd, open(input[class_num+num_classes_exp4], "r") as lengths_fd:
                for line, length_line in zip(input_fd, lengths_fd):
                    if '>' not in line:
                        counts = get_class_counts_for_read([int(x) for x in line.split()], [int(x) for x in length_line.split()])
                        max_count = max(counts)
                        pred_class = counts.index(max_count)
                        confusion_matrix[class_num][pred_class] += 1
        
        # Write out the csv data file ...
        with open(output[0], "w") as out_fd:
            out_fd.write("tool,dataset,tp,tn,fp,fn\n")
            for results in calculate_accuracy_values(confusion_matrix, num_classes_exp4):
                out_fd.write(",".join(["spumoni"] + [str(x) for x in results]) + "\n")