"""
Purpose: Script to remove reads in 16S rRNA
         read dataset whose taxonomic path is 
         ambiguous and cannot be clearly mapped 
         onto SILVA databaes.
"""

import os
import sys
import argparse

###########################################
# Section 1: Parse SILVA taxonomy and 
# find the all nodes at each level of the tree
############################################

def get_traversal_from_full_line(line):
    # Get the traversal part of line, the presence of spaces in
    # in traversal string makes it slightly less than trivial
    line_split = line.split()
    last_pos_of_traversal = max([i for i in range(len(line_split)) if ";" in line_split[i]])
    traversal = ' '.join(line_split[0:last_pos_of_traversal+1])
    
    assert (len(line_split) - last_pos_of_traversal) <= 4   
    return last_pos_of_traversal, traversal

def get_all_nodes_at_this_level_in_silva(input_path, requested_rank):
    # Format of outlist (node #, full traversal, name)
    out_list = []
    with open(input_path, "r") as in_fd:
        for line in in_fd:
            line_split = line.split()
            last_pos_of_traversal, traversal = get_traversal_from_full_line(line)
                        
            if line_split[last_pos_of_traversal+2] == requested_rank:
                name = traversal.split(";")[-2]
                node_id = int(line_split[last_pos_of_traversal+1])
                out_list.append((node_id, traversal, name))
    return out_list

def get_silva_dict_for_level_to_list(silva_rank_to_id_path):
    """
    Creates a dictionary thats takes a taxonomic rank like
    'domain', 'phylum' and returns a list of the names in the
    SILVA taxonomy. 
    """
    tree_levels = ["domain", "phylum", "class", "order", "family", "genus"]
    silva_ranks_dict = {}

    for level in tree_levels:
        silva_level_tuples = get_all_nodes_at_this_level_in_silva(silva_rank_to_id_path, level)
        silva_level_names = [tup[2] for tup in silva_level_tuples]
        silva_ranks_dict[level] = silva_level_names
    return silva_ranks_dict

###########################################
# Section 2: Extract the taxonomic
# traversal information for read set
############################################

def get_dict_mapping_read_name_to_traversal(input_file, sample_name=""):
    seqtax = {}
    with open(input_file, "r") as in_fd:
        for line in in_fd:
            line_split = line.split()
            assert len(line_split) == 2
            
            if sample_name in line_split[0]:
                seqtax[line_split[0]] = line_split[1]
    return seqtax


############################################
# Section 3: Verify that a read set traversal
# is unambiguous with respect to SILVA
############################################

def extract_specific_level_from_read_set_traversal(traversal, level):
    level_to_info = {"domain": ("sk__", 0), 
                     "phylum": ("p__", 2),
                     "class": ("c__", 3),
                     "order": ("o__", 4),
                     "family": ("f__", 5),
                     "genus": ("g__", 6)}
    assert level in level_to_info
    prefix, relevant_pos = level_to_info[level]

    trav_split = traversal.split(";")
    assert prefix in trav_split[relevant_pos]

    return trav_split[relevant_pos][len(prefix):]

def traversal_is_unambiguous(traversal, silva_ranks):
    levels = ["domain", "phylum", "class", "order", "family", "genus"]
    for level in levels:
        curr_name = extract_specific_level_from_read_set_traversal(traversal, level)
        # num_matches = sum([curr_name in x for x in silva_ranks[level]])
        num_matches = sum([curr_name in x for x in silva_ranks[level]])

        #if curr_name not in silva_ranks[level]:
        if num_matches == 0 or num_matches > 1:
            return False
        elif curr_name not in silva_ranks[level]:
            return False
        else:
            assert num_matches == 1
            assert curr_name in silva_ranks[level]
        
        # Special case (class = Clostridia)
        if level == "class" and curr_name == "Clostridia":
            return False
        # if level == "class" and curr_name == "Bacilli":
        #     return False
        # if level == "order" and curr_name == "Bacillales":
        #     return False
    return True

############################################
# Main method
############################################

def go_through_reads_and_remove_ambiguous_reads(args):
    # Generate output file paths
    mate1_file_out = args.output_dir + os.path.basename(args.mate1)
    mate2_file_out = args.output_dir + os.path.basename(args.mate2)

    # Get a dictionary to a list of names for each rank
    silva_ranks_dict = get_silva_dict_for_level_to_list(args.silva_tax_file)

    # Get a dictionary mapping read names to traversals
    read_to_traversal_dict = get_dict_mapping_read_name_to_traversal(args.read_tax_file, args.sample_name)

    # Load in all the lines for both files
    with open(args.mate1, "r") as mate1_fd, open(args.mate2, "r") as mate2_fd:
        mate1_lines = [x.strip() for x in  mate1_fd.readlines()]
        mate2_lines = [x.strip() for x in  mate2_fd.readlines()]
    
    # Open the output files, and print out only read sequences that are
    # unambiguous in the mate1 file
    with open(mate1_file_out, "w") as mate1_out, open(mate2_file_out, "w") as mate2_out:
        def print_read_in_fastq(header, seq, qual, out_fd):
            out_fd.write(f"{header}\n{seq}\n+\n{qual}\n")

        header_m1 = ""; seq_m1 = ""; qual_m1 = ""; pos = 0
        header_m2 = ""; seq_m2 = ""; qual_m2 = ""; pos = 0
        seq_list = set()

        for mate1_line, mate2_line in zip(mate1_lines, mate2_lines):
            if mate1_line.startswith(f"@{args.sample_name}"):
                assert mate2_line.startswith(f"@{args.sample_name}")
                header_m1 = mate1_line; header_m2 = mate2_line 
                pos += 1
            elif pos == 1:
                seq_m1 = mate1_line; seq_m2 = mate2_line
                pos += 1
            elif pos == 2:
                pos += 1
            elif pos == 3:
                qual_m1 = mate1_line; qual_m2 = mate2_line
                pos = 0
                
                # key statement: only keep sequence
                # if traversal is unambiguous so each 
                # member in traversal exists in SILVA
                mate1_read_group = header_m1.split("-")[0][1:]
                mate2_read_group = header_m2.split("-")[0][1:]
                assert mate1_read_group == mate2_read_group
                assert mate1_read_group in read_to_traversal_dict

                if traversal_is_unambiguous(
                        read_to_traversal_dict[mate1_read_group],
                        silva_ranks_dict):
                    print_read_in_fastq(header_m1, seq_m1, qual_m1, mate1_out)
                    print_read_in_fastq(header_m2, seq_m2, qual_m2, mate2_out)


def parse_arguments():
    parser = argparse.ArgumentParser(description="remove ambiguous reads from 16S rRNA files ...")
    parser.add_argument("--mate1", dest="mate1", help="path to mate1 file", required=True)
    parser.add_argument("--mate2", dest="mate2", help="path to mate2 file", required=True)
    parser.add_argument("--output-dir", dest="output_dir", help="path to output_dir", required=True)
    parser.add_argument("--sample", dest="sample_name", help="name of sample (e.g. A500_V12)", required=True)
    parser.add_argument("--silva-taxrank", dest="silva_tax_file", help="path to silva taxrank file (e.g. *txt)", required=True)
    parser.add_argument("--read-taxrank", dest="read_tax_file", help="path to read taxrank file (e.g. seqtax_*.tab)", required=True)
    # parser.add_argument("--human", dest="human_dataset", action="store_true", default=False, help="readset is human-gut")
    # parser.add_argument("--soil", dest="soil_dataset", action="store_true", default=False, help="readset is soil")
    # parser.add_argument("--ocean", dest="ocean_dataset", action="store_true", default=False, help="readset is ocean")
    args = parser.parse_args()
    return args

def check_args(args):
    if not os.path.isfile(args.mate1):
        print("Error: input file 1 provided is not valid.")
        exit(1)
    if not os.path.isfile(args.mate2):
        print("Error: input file 2 provided is not valid.")
        exit(1)
    
    if not os.path.isdir(args.output_dir):
        print("Error: output directory is not valid.")
        exit(1)
    elif args.output_dir[-1] != "/":
        args.output_dir += "/"
    if not os.path.isfile(args.silva_tax_file):
        print("Error: silva taxonomy rank file provided is not valid.")
        exit(1)
    if not os.path.isfile(args.read_tax_file):
        print("Error: readset taxonomy rank file provided is not valid.")
        exit(1)
    
    # dataset_count = sum([args.human_dataset, args.soil_dataset, args.ocean_dataset])
    # if dataset_count != 1:
    #     print("Error: must specify exactly one readset type")
    #     exit(1)

if __name__ == "__main__":
    args = parse_arguments()
    check_args(args)
    go_through_reads_and_remove_ambiguous_reads(args)