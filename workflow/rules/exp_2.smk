###############################################################
# Name: exp_2.smk
# Description: Contains functions and rules for
#              the type of experiment type 2:
# 
#              Simulates reads from different metagenomic 
#              classes, and then extracts MEMs, then
#              identifies which classes each MEM occurs
#              using document array profiles, and the
#              r-index. This is DIFFERENT from exp_0 
#              because I will put all the classes into one
#              r-index and time it in comparison to doc-profiles.
#
#              This workflow is heavily based on the 
#              exp7 workflow in the oma219/khoice repo.
#
# Date: September 29, 2022
###############################################################

####################################################
# Section 1: Code that is always run, when this 
#            experiment is called. It genereates the
#            needed files.
####################################################

if exp_type == 2:
    if not os.path.isdir("exp2_non_pivot_data/"):
        for i in range(1, num_datasets_exp2 + 1):
            dir_prefix = f"data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"exp2_non_pivot_data/dataset_{i}"):
                os.makedirs(f"exp2_non_pivot_data/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"exp2_non_pivot_data/dataset_{i}")
            
            # Extract the pivot and move to new directory
            dir_prefix = f"exp2_non_pivot_data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            pivot = random.choice(list_of_files)

            if not os.path.isdir(f"exp2_pivot_data/dataset_{i}/"):
                os.makedirs(f"exp2_pivot_data/dataset_{i}/")
            
            # Remove pivot from other set of genomes
            shutil.copy(pivot, f"exp2_pivot_data/dataset_{i}/pivot_{i}.fna")
            os.remove(pivot)
    
    # Build the file-list needed by pfp_doc
    if not os.path.isdir("exp2_file_lists/"):
        os.makedirs(f"exp2_file_lists")

        with open("exp2_file_lists/input_list.txt", "w") as input_fd:
            for i in range(1, num_datasets_exp2 + 1):
                dir_prefix = f"exp2_non_pivot_data/dataset_{i}/"
                list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix) if x.endswith(".fna")]
                
                for path in list_of_files:
                    input_fd.write("{} {}\n".format(path, i))

####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_dataset_non_pivot_genomes_exp2(wildcards):
    """ Returns list of non-pivot genomes for a given dataset """
    input_files = []
    for data_file in os.listdir(f"exp2_non_pivot_data/dataset_{wildcards.num}/"):
        if data_file.endswith(".fna"):
            file_name = data_file.split(".fna")[0]
            input_files.append(f"exp2_non_pivot_data/dataset_{wildcards.num}/{file_name}.fna")
    return input_files

def get_combined_ref_dataset_exp2(wildcards):
    """ Returns list of forward combined refs for each dataset """
    input_files = []
    for i in range(1, num_datasets_exp2 + 1):
        input_files.append(f"exp2_combined_refs/dataset_{i}/combined_ref_forward.fna")
    return input_files

def all_relevant_log_files_exp2(wildcards):
    input_files = []

    # Add files for doclist approach
    for i in range(1, num_datasets_exp2+1):
        input_files.append("exp2_results/doc_lists/{}/{}/pivot_{}.fastq.listings.log".format(wildcards.mem_type, wildcards.read_type, i))
    input_files.append("exp2_results/doc_lists/{}/{}/pivot_{}.docprofiles.bwt".format(wildcards.mem_type, wildcards.read_type, 1))
    input_files.append("exp2_results/doc_lists/{}/{}/pivot_{}.docprofiles".format(wildcards.mem_type, wildcards.read_type, 1))

    # Add files for doclist_optimize approach
    for i in range(1, num_datasets_exp2+1):
        input_files.append("exp2_results/doc_lists/{}/{}/pivot_{}_optimize.fastq.listings.log".format(wildcards.mem_type, wildcards.read_type, i))
    input_files.append("exp2_results/doc_lists/{}/{}/pivot_{}_optimize.docprofiles".format(wildcards.mem_type, wildcards.read_type, 1))
    input_files.append("exp2_results/doc_lists/{}/{}/pivot_{}_optimize.docprofiles.bwt".format(wildcards.mem_type, wildcards.read_type, 1))

    # Add files for rindex approach
    for i in range(1, num_datasets_exp2+1):
        input_files.append("exp2_results/r_index/{}/{}/pivot_{}.sam.log".format(wildcards.mem_type, wildcards.read_type, i))
    input_files.append("exp2_indexes/r_index/combined_ref_all.fna.ri")
    input_files.append("exp2_indexes/r_index/combined_ref_all.fna.1.ri")

    return input_files

####################################################
# Section 3: Rules needed for this experiment type
####################################################

#   Section 3.1: Generating long reads from each 
#   out-pivot genome for each dataset, focused on 
#   long reads since I expect a wider range of MEM
#   lengths, so a more variety of class memberships.

rule generate_raw_positive_long_reads_exp2:
    input:
        "exp2_pivot_data/dataset_{num}/pivot_{num}.fna"
    output:
        "exp2_pivot_data/dataset_{num}/ont/pivot_{num}_reads.fna"
    shell:
        """
        pbsim --depth 5.0 --prefix exp2_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 200 exp2_pivot_data/dataset_{wildcards.num}/pivot_{wildcards.num}.fna
        
        cat 'exp2_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads'*.fastq > exp2_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads.fastq
        seqtk seq -a exp2_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads.fastq > {output}
        rm 'exp2_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads_'*.fastq
        ls  'exp2_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads'*.ref | xargs rm
        ls  'exp2_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads'*.maf | xargs rm
        """

#   Section 3.2: Contains rules for building a concatenated reference of all
#   the genomes in each database. Then, building reverse complements, and 
#   concatenating those as well. SPUMONI just uses the forward sequences, but
#   we want the forward + rev. comp. for the r-index

rule build_forward_ref_dataset_exp2:
    input:
        get_dataset_non_pivot_genomes_exp2
    output:
        "exp2_combined_refs/dataset_{num}/combined_ref_forward.fna"
    shell:
        "cat {input} > {output}"

rule concat_only_forward_refs_for_all_datasets_exp2:
    input:
        get_combined_ref_dataset_exp2
    output:
        "exp2_combined_refs/combined_ref_forward_only.fna"
    shell:
        "cat {input} > {output}"

rule build_revcomp_for_all_datasets_exp2:
    input:
        "exp2_combined_refs/combined_ref_forward_only.fna"
    output:
        "exp2_combined_refs/combined_ref_revcomp_only.fna"
    shell:
        "seqtk seq -r {input} | sed 's/^>/>revcomp_/' > {output}"

rule combine_forward_and_revcomp_for_all_datasets_exp2:
    input:
        "exp2_combined_refs/combined_ref_forward_only.fna",
        "exp2_combined_refs/combined_ref_revcomp_only.fna"
    output:
        "exp2_combined_refs/combined_ref_all.fna"
    shell:
        "cat {input} > {output}"


#   Section 3.3: Build the SPUMONI index and then use it
#   to compute MS wrt to the database to find MEMs

rule build_spumoni_index_exp2:
    input:
        "exp2_combined_refs/combined_ref_forward_only.fna"
    output:
        "exp2_indexes/spumoni_index/spumoni_full_ref.fa",
        "exp2_indexes/spumoni_index/spumoni_full_ref.fa.thrbv.ms"
    shell:
        # Only use forward references for spumoni
        """
        cp {input} exp2_indexes/spumoni_index/combined_ref_forward_only.fna
        spumoni build -r exp2_indexes/spumoni_index/combined_ref_forward_only.fna -M -n
        """

rule run_spumoni_to_generate_ms_exp2:
    input:
        "exp2_indexes/spumoni_index/spumoni_full_ref.fa",
        "exp2_pivot_data/dataset_{num}/{read_type}/pivot_{num}_reads.fna"
    output:
        "exp2_pivot_data/dataset_{num}/{read_type}/pivot_{num}_reads.fna.lengths"
    shell:
        """
        # Don't make this multi-threaded, it will mess up the next step
        spumoni run -r {input[0]} -p {input[1]} -M -n 
        """

#   Section 3.4: Extract MEMs using the pivot's MS lengths
#   and the pivot reads itself.

rule extract_mems_or_halfmems_based_on_ms_exp2:
    input:
        "exp2_pivot_data/dataset_{pivot}/{read_type}/pivot_{pivot}_reads.fna.lengths",
        "exp2_pivot_data/dataset_{pivot}/{read_type}/pivot_{pivot}_reads.fna"
    output:
        "exp2_{mem_type}_data/{read_type}/pivot_{pivot}.fastq"
    shell:
        """
        python3 {repo_dir}/src/extract_mems.py \
        -l {input[0]} \
        -p {input[1]} \
        -c {mem_count_exp2} \
        --{wildcards.mem_type} \
        -t {thresh_exp2} \
        -o {output} 
        """

#   Section 3.5: Build r-index for each dataset individually over
#   the reference containing both the forward and reverse complement.

rule build_rindex_over_entire_database_exp2:
    input:
        "exp2_combined_refs/combined_ref_all.fna"
    output:
        "exp2_indexes/r_index/combined_ref_all.fna.ri",
        "exp2_indexes/r_index/combined_ref_all.fna.1.ri"
    shell:
        """
        cp {input} exp2_indexes/r_index/combined_ref_all.fna

        cd {r_dir} 
        ri-buildfasta {base_dir}/exp2_indexes/r_index/combined_ref_all.fna
        cd {base_dir}
        """

# Section 3.6: Build a pfp_doc index in order to prepare for computing
# document listings for each provided MEM, and then query the index with
# MEMs

rule build_doc_array_profiles_for_full_database_exp2:
    input:
        "exp2_file_lists/input_list.txt"
    output:
        "exp2_indexes/docprofiles_index/full_ref.fna.sdap",
        "exp2_indexes/docprofiles_index/full_ref.fna.edap"
    shell:
        """
        pfp_doc build -f {input} -o exp2_indexes/docprofiles_index/full_ref -r
        """

# Section 3.7: Run queries through both the r-index and the 
# document array profiles and time them.

rule determine_doc_array_listings_for_mems_exp2:
    input:
        "exp2_indexes/docprofiles_index/full_ref.fna.sdap",
        "exp2_{mem_type}_data/{read_type}/pivot_{pivot}.fastq"
    output:
        "exp2_results/doc_lists/{mem_type}/{read_type}/pivot_{pivot}.fastq",
        "exp2_results/doc_lists/{mem_type}/{read_type}/pivot_{pivot}_optimize.fastq",
        "exp2_results/doc_lists/{mem_type}/{read_type}/pivot_{pivot}.fastq.listings",
        "exp2_results/doc_lists/{mem_type}/{read_type}/pivot_{pivot}.fastq.listings.resources",
        "exp2_results/doc_lists/{mem_type}/{read_type}/pivot_{pivot}_optimize.fastq.listings",
        "exp2_results/doc_lists/{mem_type}/{read_type}/pivot_{pivot}_optimize.fastq.listings.resources",
        "exp2_results/doc_lists/{mem_type}/{read_type}/pivot_{pivot}.fastq.listings.log",
        "exp2_results/doc_lists/{mem_type}/{read_type}/pivot_{pivot}_optimize.fastq.listings.log",
        "exp2_results/doc_lists/{mem_type}/{read_type}/pivot_{pivot}.docprofiles.bwt",
        "exp2_results/doc_lists/{mem_type}/{read_type}/pivot_{pivot}.docprofiles",
        "exp2_results/doc_lists/{mem_type}/{read_type}/pivot_{pivot}_optimize.docprofiles.bwt",
        "exp2_results/doc_lists/{mem_type}/{read_type}/pivot_{pivot}_optimize.docprofiles"
    shell:
        """
        cp {input[1]} {output[0]}
        {time_prog} {time_format} --output={output[3]} pfp_doc run -r exp2_indexes/docprofiles_index/full_ref -p {output[0]} -s 2> {output[2]}.log
        cp exp2_indexes/docprofiles_index/full_ref.fna.docprofiles.bwt {output[8]}
        cp exp2_indexes/docprofiles_index/full_ref.fna.docprofiles {output[9]}

        cp {input[1]} {output[1]}
        {time_prog} {time_format} --output={output[5]} pfp_doc run -r exp2_indexes/docprofiles_index/full_ref -p {output[1]} -s -l 150 2> {output[4]}.log
        cp exp2_indexes/docprofiles_index/full_ref.fna.docprofiles.bwt {output[10]}
        cp exp2_indexes/docprofiles_index/full_ref.fna.docprofiles {output[11]}
        """

rule align_pivot_against_rindex_of_database_exp2:
    input:
        "exp2_indexes/r_index/combined_ref_all.fna.ri",
        "exp2_{mem_type}_data/{read_type}/pivot_{pivot}.fastq"
    output:
        "exp2_results/r_index/{mem_type}/{read_type}/pivot_{pivot}.sam",
        "exp2_results/r_index/{mem_type}/{read_type}/pivot_{pivot}.sam.resources",
        "exp2_results/r_index/{mem_type}/{read_type}/pivot_{pivot}.sam.log"
    shell:
        """
        cd {r_dir} 
        {time_prog} {time_format} --output="{base_dir}/{output[1]}" ri-align locate {base_dir}/exp2_indexes/r_index/combined_ref_all.fna {base_dir}/{input[1]} \
        > {base_dir}/{output[0]} 2> {base_dir}/{output[2]}
        """

# Section 3.8: Build a csv-file for a single pivot results
# from both the r-index and doc-profiles

rule build_csv_results_for_single_pivot_exp2:
    input:
        all_relevant_log_files_exp2
    output:
        "exp2_final_output/{mem_type}/{read_type}/output.csv"
    shell:
        """        
        # Build a bash array of files
        FILES=()        
        for file in {input}
        do
            FILES+=("$file")
        done

        # Add header to csv file
        echo "approach,numclasses,n,r,time,indexsize" >> {output[0]}

        # Analyze the doclisting approaches
        pos=0; n=0; r=0
        for approach in "doclist" "doclist_optimize" 
        do
            total_time=0.0
            for i in $(seq 1 {num_datasets_exp2})
            do
                curr_file=${{FILES[$pos]}}

                curr_time=$(grep "run_main" $curr_file | awk '{{print substr($(NF-1), 2)}}')
                total_time=$(echo "$total_time + $curr_time" | bc)

                n=$(grep "statistics" $curr_file | awk '{{print substr($6, 1, length($6)-1)}}')
                r=$(grep "statistics" $curr_file | awk '{{print $9}}')
                
                pos=$((pos+1))
            done

            index_size_1=0; index_size_2=0

            # Get size of index files
            index_size_1=$(ls -l ${{FILES[$pos]}} | awk '{{print $5}}')
            pos=$((pos+1))

            index_size_2=$(ls -l ${{FILES[$pos]}} | awk '{{print $5}}')
            pos=$((pos+1))
            total_index=$((index_size_1+index_size_2))

            echo "$approach,{num_datasets_exp2},$n,$r,$total_time,$total_index" >> {output[0]}
        done

        # Analyze the rindex approach
        for approach in "rindex"
        do
            total_time=0.0
            for i in $(seq 1 {num_datasets_exp2})
            do
                curr_file=${{FILES[$pos]}}

                curr_time=$(grep "time" $curr_file | awk '{{print substr($3, 1, length($3)-1)}}')
                total_time=$(echo "$total_time + $curr_time" | bc)
                
                pos=$((pos+1))
            done

            index_size_1=0; index_size_2=0

            # Get size of index files
            index_size_1=$(ls -l ${{FILES[$pos]}} | awk '{{print $5}}')
            pos=$((pos+1))

            index_size_2=$(ls -l ${{FILES[$pos]}} | awk '{{print $5}}')
            pos=$((pos+1))
            total_index=$((index_size_1+index_size_2))

            echo "$approach,{num_datasets_exp2},$n,$r,$total_time,$total_index" >> {output[0]}
        done
        """

