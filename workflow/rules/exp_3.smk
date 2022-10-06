###############################################################
# Name: exp_3.smk
# Description: Contains functions and rules for
#              the type of experiment type 3:
# 
#              Simulates reads that contain a mixture of E. coli
#              and AMR genes from different classes. Then tries to
#              classify them using this document array listing
#              approach and SPUMONI 2 usind sampled document array.
#
# Date: October 5, 2022
###############################################################

####################################################
# Section 1: Code that is always run, when this 
#            experiment is called. It genereates the
#            needed files.
####################################################

####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_complete_gene_set_exp3(wildcards):
    """ Grab the file the data/genes folder """
    input_files = []
    for data_file in os.listdir(f"data/genes/"):
        if data_file.endswith(".fna"):
            input_files.append("data/genes/" + data_file)
    return input_files[0]

def get_bacteria_genome_exp3(wildcards):
    """ Grab the bacteria (E. coli) genome file in folder """
    input_files = []
    for data_file in os.listdir(f"data/dataset_1/"):
        if data_file.endswith(".fna"):
            input_files.append("data/dataset_1/" + data_file)
    return input_files[0]

####################################################
# Section 3: Rules needed for this experiment type
####################################################

#   Section 3.1: Separate out the genes from each sub-class
#   into separate FASTA files

rule generate_separate_gene_subclasses_exp3:
    input:
        get_complete_gene_set_exp3
    output:
        expand("exp3_gene_data/raw_seq/class_{num}.fna", num=range(1, num_gene_classes_exp3+1))
    shell:
        """
        python3 {repo_dir}/src/separate_gene_classes.py -n {num_gene_classes_exp3} -i {input} -o exp3_gene_data/raw_seq/
        """

#   Section 3.2: Simulate nanopore reads from the E. coli genome, and subset it
#   make sure we have enough reads for each class.

rule generate_long_reads_from_bacteria_exp3:
    input:
        get_bacteria_genome_exp3
    output:
        "exp3_bacteria_reads/raw_reads/all_classes.fna"
    shell:
        """
        pbsim --depth 100.0 --prefix exp3_bacteria_reads/raw_reads/all_classes \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 900 --length-max 1100 {input}
        
        cat 'exp3_bacteria_reads/raw_reads/all_classes'*.fastq > exp3_bacteria_reads/raw_reads/all_classes.fastq
        seqtk seq -a exp3_bacteria_reads/raw_reads/all_classes.fastq > {output}
        rm 'exp3_bacteria_reads/raw_reads/all_classes_'*.fastq
        ls  'exp3_bacteria_reads/raw_reads/all_classes'*.ref | xargs rm
        ls  'exp3_bacteria_reads/raw_reads/all_classes'*.maf | xargs rm
        """

rule subset_long_reads_from_bacteria_exp3:
    input:
        "exp3_bacteria_reads/raw_reads/all_classes.fna"
    output:
        "exp3_bacteria_reads/raw_reads/all_classes_subset.fna"
    shell:
        """
        num_reads_present=$(grep -c '>' {input})
        needed_reads=$(({reads_per_class_exp3} * {num_gene_classes_exp3}))
        if (( $num_reads_present < $needed_reads )); then echo 'not enough reads for all classes!'; exit 1; fi

        num_lines=$(({reads_per_class_exp3} * {num_gene_classes_exp3} * 2))
        head -n $num_lines {input} > {output}
        """

#   Section 3.3: Simulate reads from each class of genes in order
#   mix them into the reads from the Bacteria

rule generate_long_reads_from_genes_exp3:
    input:
        "exp3_gene_data/raw_seq/class_{num}.fna"
    output:
        "exp3_gene_data/sim_reads/class_{num}_reads.fna"
    shell:
        """
        pbsim --depth 10.0 --prefix exp3_gene_data/sim_reads/class_{wildcards.num}_reads \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 900 --length-max 1100 {input}
        
        cat 'exp3_gene_data/sim_reads/class_{wildcards.num}_reads'*.fastq > exp3_gene_data/sim_reads/class_{wildcards.num}_reads.fastq
        seqtk seq -a exp3_gene_data/sim_reads/class_{wildcards.num}_reads.fastq > {output}
        rm 'exp3_gene_data/sim_reads/class_{wildcards.num}_reads_'*.fastq
        ls  'exp3_gene_data/sim_reads/class_{wildcards.num}_reads'*.ref | xargs rm
        ls  'exp3_gene_data/sim_reads/class_{wildcards.num}_reads'*.maf | xargs rm
        rm exp3_gene_data/sim_reads/class_{wildcards.num}_reads.fastq 
        """  

#   Section 3.4: Generate a mix of reads for each class of genes, basically
#   I would take a bacteria and randomly insert an AMR gene into it. They were 
#   simulated to be roughly the same size so the read will be 50:50 bacteria and 
#   and AMR gene.

rule generate_mixed_reads_exp3:
    input:
        "exp3_bacteria_reads/raw_reads/all_classes_subset.fna",
        expand("exp3_gene_data/sim_reads/class_{num}_reads.fna", num=range(1, num_gene_classes_exp3+1))
    output:
        expand("exp3_mixed_reads/class_{num}_reads.fna", num=range(1, num_gene_classes_exp3+1))
    shell:
        """
        python3 {repo_dir}/src/generate_mixed_reads.py -n {num_gene_classes_exp3} \
                                                       -r {input[0]} \
                                                       -i exp3_gene_data/sim_reads/ \
                                                       -o exp3_mixed_reads/ \
                                                       --num-reads {reads_per_class_exp3}
        """

#   Section 3.5: Build a file-list for database, build SPUMONI index over
#   AMR genes, and generate document labels for each read to use for 
#   classification

rule generate_gene_file_list_exp3:
    input:
        expand("exp3_gene_data/raw_seq/class_{num}.fna", num=range(1, num_gene_classes_exp3+1))
    output:
        "exp3_filelist/gene_filelist.txt"
    run:
        with open(output[0], "w") as out_fd:
            pos = 1
            for in_file in input:
                out_fd.write(f"{in_file} {pos}\n")
                pos += 1

rule build_spumoni_index_for_exp3:
    input:
        "exp3_filelist/gene_filelist.txt"
    output:
        "exp3_indexes/spumoni_index/spumoni_full_ref.fa.thrbv.ms"
    shell:
        """
        spumoni build -i {input} -b exp3_indexes/spumoni_index/ -n -M -d
        """

rule run_query_with_spumoni_for_doc_labels_exp3:
    input:
        "exp3_indexes/spumoni_index/spumoni_full_ref.fa.thrbv.ms",
        "exp3_mixed_reads/class_{num}_reads.fna"
    output:
        "exp3_results/spumoni/ms_results/class_{num}_reads.fna",
        "exp3_results/spumoni/ms_results/class_{num}_reads.fna.lengths",
        "exp3_results/spumoni/ms_results/class_{num}_reads.fna.doc_numbers"
    shell:
        """
        cp {input[1]} {output[0]}
        spumoni run -r exp3_indexes/spumoni_index/spumoni_full_ref.fa -p {output[0]} -M -d -n
        """

#   Section 3.6: Extract MEMs from each read, and then list the documents
#   for each MEM in order to "classify" each read

rule run_query_with_spumoni_for_MEMs_exp3:
    input:
        "exp3_indexes/spumoni_index/spumoni_full_ref.fa.thrbv.ms",
        "exp3_mixed_reads/class_{num}_reads.fna"
    output:
        "exp3_results/doclist/ms_results/class_{num}_reads.fna",
        "exp3_results/doclist/ms_results/class_{num}_reads.fna.lengths"
    shell:
        """
        cp {input[1]} {output[0]}
        spumoni run -r exp3_indexes/spumoni_index/spumoni_full_ref.fa -p {output[0]} -M -n
        """

rule extract_mems_based_on_ms_exp3:
    input:
        "exp3_results/doclist/ms_results/class_{num}_reads.fna.lengths",
        "exp3_results/doclist/ms_results/class_{num}_reads.fna"
    output:
        "exp3_results/doclist/mem_file/class_{num}_mems.fastq"
    shell:
        """
        python3 {repo_dir}/src/extract_mems_exp3.py \
        -l {input[0]} \
        -p {input[1]} \
        --mems \
        -t 0 \
        -o {output} 
        """

rule build_doc_array_profiles_for_listing_exp3:
    input:
        "exp3_filelist/gene_filelist.txt"
    output:
        "exp3_indexes/docprofiles_index/full_ref.fna.sdap",
        "exp3_indexes/docprofiles_index/full_ref.fna.edap"
    shell:
        """
        pfp_doc build -f {input} -o exp3_indexes/docprofiles_index/full_ref -r
        """

rule determine_doc_array_listings_for_mems_exp3:
    input:
        "exp3_indexes/docprofiles_index/full_ref.fna.sdap",
        "exp3_results/doclist/mem_file/class_{num}_mems.fastq"
    output:
        "exp3_results/doclist/listings/class_{num}_reads.fastq",
        "exp3_results/doclist/listings/class_{num}_reads.fastq.listings"
    shell:
        """
        cp {input[1]} {output[0]}
        pfp_doc run -r exp3_indexes/docprofiles_index/full_ref -p {output[0]}
        """

#   Section 3.7: Generate confusion matix for SPUMONI and document array
#   profile separately by analyzing the reads from each class.

rule classify_spumoni_reads_exp3:
    input:
        expand("exp3_results/spumoni/ms_results/class_{num}_reads.fna.doc_numbers", num=range(1, num_gene_classes_exp3+1))
    output:
        "exp3_results/spumoni/csv_files/classification_results.csv"
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

        def get_class_counts_for_read(doc_list):
            """ Takes in document labels for read and return counts """
            counts = [0 for i in range(num_gene_classes_exp3)]
            for doc_num in doc_list:
                counts[doc_num] += 1
            return counts
        
        confusion_matrix = [[0 for x in range(num_gene_classes_exp3)] for i in range(num_gene_classes_exp3)]
        
        for class_num, input_file in enumerate(input):
            with open(input_file, "r") as input_fd:
                for line in input_fd:
                    if '>' not in line:
                        counts = get_class_counts_for_read([int(x) for x in line.split()])
                        max_count = max(counts)
                        pred_class = counts.index(max_count)
                        confusion_matrix[class_num][pred_class] += 1
        
        # Write out the csv data file ...
        with open(output[0], "w") as out_fd:
            out_fd.write("tool,dataset,tp,tn,fp,fn\n")
            for results in calculate_accuracy_values(confusion_matrix, num_gene_classes_exp3):
                out_fd.write(",".join(["spumoni"] + [str(x) for x in results]) + "\n")


rule classify_docarray_profile_reads_exp3:
    input:
        expand("exp3_results/doclist/listings/class_{num}_reads.fastq.listings", num=range(1, num_gene_classes_exp3+1))
    output:
        "exp3_results/doclist/csv_files/classification_results.csv"
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
        
        confusion_matrix = [[0 for x in range(num_gene_classes_exp3)] for i in range(num_gene_classes_exp3)]

        for class_num, input_file in enumerate(input):
            with open(input_file, "r") as input_fd:
                curr_read_num = 0
                curr_mem_length = 0
                curr_read_stats = [0 for i in range(num_gene_classes_exp3)]
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
                            curr_read_stats = [0 for i in range(num_gene_classes_exp3)]
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
                curr_read_stats = [0 for i in range(num_gene_classes_exp3)]
        
        # Write out the csv data file ...
        with open(output[0], "w") as out_fd:
            out_fd.write("tool,dataset,tp,tn,fp,fn\n")
            for results in calculate_accuracy_values(confusion_matrix, num_gene_classes_exp3):
                out_fd.write(",".join(["docprofile"] + [str(x) for x in results]) + "\n")

