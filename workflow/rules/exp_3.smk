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
            print(data_file)
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
        pbsim --depth 10.0 --prefix exp3_bacteria_reads/raw_reads/all_classes \
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
        echo "hello" > exp3_mixed_reads/class_1_reads.fna
        echo "hello" > exp3_mixed_reads/class_2_reads.fna
        echo "hello" > exp3_mixed_reads/class_3_reads.fna
        echo "hello" > exp3_mixed_reads/class_4_reads.fna

        python3 {repo_dir}/src/generate_mixed_reads.py -n {num_gene_classes_exp3} \
                                                       -r {input[0]} \
                                                       -i exp3_gene_data/sim_reads/ \
                                                       -o exp3_mixed_reads/ \
                                                       --num-reads {reads_per_class_exp3}
        """


