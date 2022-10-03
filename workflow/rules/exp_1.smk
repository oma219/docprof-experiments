###############################################################
# Name: exp_1.smk
# Description: Contains functions and rules for
#              the type of experiment type 1:
# 
#              Simulates reads from different metagenomic 
#              classes, and then extracts MEMs, then
#              identifies which classes each MEM occurs
#              using document array profiles, and the
#              r-index.
#
#              This workflow is heavily based on the 
#              exp7 workflow in the oma219/khoice repo.
#
# Date: September 12, 2022
###############################################################

####################################################
# Section 1: Code that is always run, when this 
#            experiment is called. It genereates the
#            needed files.
####################################################

if exp_type == 1:
    if not os.path.isdir("exp1_non_pivot_data/"):
        for i in range(1, num_datasets_exp1 + 1):
            dir_prefix = f"data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"exp1_non_pivot_data/dataset_{i}"):
                os.makedirs(f"exp1_non_pivot_data/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"exp1_non_pivot_data/dataset_{i}")
            
            # Extract the pivot and move to new directory
            dir_prefix = f"exp1_non_pivot_data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            pivot = random.choice(list_of_files)

            if not os.path.isdir(f"exp1_pivot_data/dataset_{i}/"):
                os.makedirs(f"exp1_pivot_data/dataset_{i}/")
            
            # Remove pivot from other set of genomes
            shutil.copy(pivot, f"exp1_pivot_data/dataset_{i}/pivot_{i}.fna")
            os.remove(pivot)
    
    # Build the file-list needed by pfp_doc
    if not os.path.isdir("exp1_file_lists/"):
        os.makedirs(f"exp1_file_lists")

        with open("exp1_file_lists/input_list.txt", "w") as input_fd:
            for i in range(1, num_datasets_exp1 + 1):
                dir_prefix = f"exp1_non_pivot_data/dataset_{i}/"
                list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix) if x.endswith(".fna")]
                
                for path in list_of_files:
                    input_fd.write("{} {}\n".format(path, i))


####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_dataset_non_pivot_genomes_exp1(wildcards):
    """ Returns list of non-pivot genomes for a given dataset """
    input_files = []
    for data_file in os.listdir(f"exp1_non_pivot_data/dataset_{wildcards.num}/"):
        if data_file.endswith(".fna"):
            file_name = data_file.split(".fna")[0]
            input_files.append(f"exp1_non_pivot_data/dataset_{wildcards.num}/{file_name}.fna")
    return input_files

def get_combined_ref_dataset_exp1(wildcards):
    """ Returns list of forward combined refs for each dataset """
    input_files = []
    for i in range(1, num_datasets_exp1 + 1):
        input_files.append(f"exp1_combined_refs/dataset_{i}/combined_ref_forward.fna")
    return input_files

def get_python_input_exp1(wildcards):
    """ Returns a list of all pivot SAM files for a read type """
    input_files = []
    curr_pivot = wildcards.num
    for i in range(1,num_datasets_exp1+1):
        input_files.append(f"exp1_sam_files/{wildcards.mem_type}/{wildcards.read_type}/pivot_{curr_pivot}_align_dataset_{i}.sam")
    return input_files

####################################################
# Section 3: Rules needed for this experiment type
####################################################

#   Section 3.1: Generating long reads from each 
#   out-pivot genome for each dataset, focused on 
#   long reads since I expect a wider range of MEM
#   lengths, so a more variety of class memberships.

rule generate_raw_positive_long_reads_exp1:
    input:
        "exp1_pivot_data/dataset_{num}/pivot_{num}.fna"
    output:
        "exp1_pivot_data/dataset_{num}/ont/pivot_{num}_reads.fna"
    shell:
        """
        pbsim --depth 20.0 --prefix exp1_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 200 exp1_pivot_data/dataset_{wildcards.num}/pivot_{wildcards.num}.fna
        
        cat 'exp1_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads'*.fastq > exp1_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads.fastq
        seqtk seq -a exp1_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads.fastq > {output}
        rm 'exp1_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads_'*.fastq
        ls  'exp1_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads'*.ref | xargs rm
        ls  'exp1_pivot_data/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_reads'*.maf | xargs rm
        """

#   Section 3.2: Contains rules for building a concatenated reference of all
#   the genomes in each database. Then, building reverse complements, and 
#   concatenating those as well. SPUMONI just uses the forward sequences, but
#   we want the forward + rev. comp. for the r-index

rule build_forward_ref_dataset_exp1:
    input:
        get_dataset_non_pivot_genomes_exp1
    output:
        "exp1_combined_refs/dataset_{num}/combined_ref_forward.fna"
    shell:
        "cat {input} > {output}"

rule build_rev_comp_ref_dataset_exp1:
    input:
        "exp1_combined_refs/dataset_{num}/combined_ref_forward.fna"
    output:
        "exp1_combined_refs/dataset_{num}/combined_ref_rcomp.fna"
    shell:
        "seqtk seq -r {input} | sed 's/^>/>revcomp_/' > {output}"

rule concat_forward_and_rev_comp_dataset_exp1:
    input:
        "exp1_combined_refs/dataset_{num}/combined_ref_forward.fna",
        "exp1_combined_refs/dataset_{num}/combined_ref_rcomp.fna"
    output:
        "exp1_combined_refs/dataset_{num}/combined_ref_all.fna"
    shell:
        "cat {input} > {output}"

rule concat_only_forward_refs_for_all_datasets_exp1:
    input:
        get_combined_ref_dataset_exp1
    output:
        "exp1_combined_refs/combined_ref_forward_only.fna"
    shell:
        "cat {input} > {output}"

#   Section 3.3: Build the SPUMONI index and then use it
#   to compute MS wrt to the database to find MEMs

rule build_spumoni_index_exp1:
    input:
        "exp1_combined_refs/combined_ref_forward_only.fna"
    output:
        "exp1_indexes/spumoni_index/spumoni_full_ref.fa",
        "exp1_indexes/spumoni_index/spumoni_full_ref.fa.thrbv.ms"
    shell:
        # Only use forward references for spumoni
        """
        cp {input} exp1_indexes/spumoni_index/combined_ref_forward_only.fna
        spumoni build -r exp1_indexes/spumoni_index/combined_ref_forward_only.fna -M -n
        """

rule run_spumoni_to_generate_ms_exp1:
    input:
        "exp1_indexes/spumoni_index/spumoni_full_ref.fa",
        "exp1_pivot_data/dataset_{num}/{read_type}/pivot_{num}_reads.fna"
    output:
        "exp1_pivot_data/dataset_{num}/{read_type}/pivot_{pivot}_reads.fna.lengths"
    shell:
        """
        # Don't make this multi-threaded, it will mess up the next step
        spumoni run -r {input[0]} -p {input[1]} -M -n 
        """

#   Section 3.4: Extract MEMs using the pivot's MS lengths
#   and the pivot reads itself. 

rule extract_mems_or_halfmems_based_on_ms_exp1:
    input:
        "exp1_pivot_data/dataset_{pivot}/{read_type}/pivot_{pivot}_reads.fna.lengths",
        "exp1_pivot_data/dataset_{pivot}/{read_type}/pivot_{pivot}_reads.fna"
    output:
        "exp1_{mem_type}_data/{read_type}/pivot_{pivot}.fastq"
    shell:
        """
        python3 {repo_dir}/src/extract_mems.py \
        -l {input[0]} \
        -p {input[1]} \
        -c {mem_count_exp1} \
        --{wildcards.mem_type} \
        -t {thresh_exp1} \
        -o {output} 
        """

#   Section 3.5: Build r-index for each dataset individually over
#   the reference containing both the forward and reverse complement.

rule build_rindex_over_individual_databases_exp1:
    input:
        "exp1_combined_refs/dataset_{num}/combined_ref_all.fna"
    output:
        "exp1_combined_refs/dataset_{num}/combined_ref_all.fna.ri"
    shell:
        """
        cd {r_dir} 
        ri-buildfasta {base_dir}/{input}
        cd {base_dir}
        """

#   Section 3.6: Create SAM file using r-index to locate all the 
#   half-MEMs or MEMs in each of the individual databases.

rule align_pivot_against_rindex_of_database_exp1:
    input:
        "exp1_combined_refs/dataset_{num}/combined_ref_all.fna",
        "exp1_combined_refs/dataset_{num}/combined_ref_all.fna.ri",
        "exp1_{mem_type}_data/{read_type}/pivot_{pivot}.fastq"
    output:
        "exp1_sam_files/{mem_type}/{read_type}/pivot_{pivot}_align_dataset_{num}.sam"
    shell:
        """
        cd {r_dir} 
        ri-align -m 1 locate {base_dir}/{input[0]} {base_dir}/{input[2]} \
        > {base_dir}/{output}
        """

# Section 3.7: Build a pfp_doc index in order to prepare for computing
# document listings for each provided MEM, and then query the index with
# MEMs

rule build_doc_array_profiles_for_full_database_exp1:
    input:
        "exp1_file_lists/input_list.txt"
    output:
        "exp1_indexes/docprofiles_index/full_ref.fna.sdap",
        "exp1_indexes/docprofiles_index/full_ref.fna.edap"
    shell:
        """
        pfp_doc build  -f {input} -o exp1_indexes/docprofiles_index/full_ref -r
        """

rule determined_doc_array_listings_for_mems_exp1:
    input:
        "exp1_indexes/docprofiles_index/full_ref.fna.sdap",
        "exp1_{mem_type}_data/{read_type}/pivot_{pivot}.fastq"
    output:
        "exp1_doclistings/{mem_type}/{read_type}/pivot_{pivot}.fastq",
        "exp1_doclistings/{mem_type}/{read_type}/pivot_{pivot}.fastq.listings"
    shell:
        """
        cp {input[1]} {output[0]}
        pfp_doc run -r exp1_indexes/docprofiles_index/full_ref -p {output[0]} -l 1000
        """


#   Section 3.7: Analyze the SAM files with respect to each database
#   to determine how we should split each 

rule analyze_sam_for_each_type_exp1:
    input:
        get_python_input_exp1,
        "exp1_doclistings/{mem_type}/{read_type}/pivot_{num}.fastq.listings"
    output:
        "exp1_output_data/{mem_type}/{read_type}/pivot_{num}/output.txt"
    shell:
        """
        python3 {repo_dir}/src/analyze_sam_exp1.py \
        -n {num_datasets_exp1} \
        -c {wildcards.num} \
        -s {base_dir}/exp1_sam_files/{wildcards.mem_type}/{wildcards.read_type}/ \
        -o {base_dir}/exp1_output_data/{wildcards.mem_type}/{wildcards.read_type}/pivot_{wildcards.num}/ \
        --doclist exp1_doclistings/{wildcards.mem_type}/{wildcards.read_type}/pivot_{wildcards.num}.fastq.listings
        """