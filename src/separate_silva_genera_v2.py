import os
import sys
import argparse
import multiprocess as mp

########################################
# Section 1: SILVA Tree-related methods
########################################

class TreeNode(object):   
    def __init__(self, data, traversal="", rank=""):
        self.data = data
        self.children = []
        self.traversal = traversal
        self.rank = rank

    def add_child(self, obj):
        self.children.append(obj)
    
    def add_traversal(self, traversal):
        self.traversal = traversal
    
    def add_rank(self, rank):
        self.rank = rank

def strip_parentheses(tree_str):
    return tree_str[1:-1]

def get_silva_newick_tree(input_filepath):
    """ Loads the input in Newick Tree format """
    in_fd = open(input_filepath, "r")
    lines = [x.strip() for x in in_fd.readlines()]
    assert len(lines) == 1
    
    # Remove two extra characters at end
    silva_tree_str = lines[0][:-2]
    return silva_tree_str

def validate_subtree(tree_str):
    assert tree_str.count("(") == tree_str.count(")")

def node_split(input_str):
    """ 
    Split str based on comma but go around interior nodes so 
    ((One,Two),(Three,Four)),Five would return a list of 
    [((One,Two),(Three,Four)), Five]
    """
    out_list = []; pos = 0
    while pos < len(input_str):
        curr_str = ""
        # Case 1: interior node 
        if input_str[pos] == "(":
            curr_str += input_str[pos]
            status = 1
            pos += 1

            while status != 0:
                if input_str[pos] == "(":
                    status += 1
                elif input_str[pos] == ")":
                    status -= 1
                curr_str += input_str[pos]
                pos += 1
            out_list.append(curr_str)
        # Case 2: comma
        elif input_str[pos] == ",":
              pos += 1
        # Case 3: tip
        else:
            while pos < len(input_str) and input_str[pos] != ",":
                curr_str += input_str[pos]
                pos += 1
            out_list.append(curr_str)
    return out_list

def build_newick_tree(curr_node, tree_str, node_map):
    """ builds the tree based on newick format """
    validate_subtree(tree_str)
    if "(" not in tree_str and ")" not in tree_str:
        children = tree_str.split(",")
        for child in children:
            assert int(child) in node_map
            curr_traversal, curr_rank = node_map[int(child)]
            curr_node.add_child(TreeNode(child, curr_traversal, curr_rank))
    else:
        interior_children = node_split(tree_str)
        for child in interior_children:
            if "(" not in child:
                assert int(child) in node_map
                curr_traversal, curr_rank = node_map[int(child)]
                curr_node.add_child(TreeNode(child, curr_traversal, curr_rank))
            else:
                curr_node.add_child(build_newick_tree(TreeNode("Interior"), strip_parentheses(child), node_map))
    return curr_node

def get_genus_from_traversal(traversal):
    return traversal.split(";")[-2]

def get_domain_from_traversal(traversal):
    return traversal.split(";")[0]

def get_traversal_from_full_line(line):
    # Get the traversal part of line, the presence of spaces in
    # in traversal string makes it slightly less than trivial
    line_split = line.split()
    last_pos_of_traversal = max([i for i in range(len(line_split)) if ";" in line_split[i]])
    traversal = ' '.join(line_split[0:last_pos_of_traversal+1])
    
    assert (len(line_split) - last_pos_of_traversal) <= 4   
    return last_pos_of_traversal, traversal

def get_genera_list(input_path):
    # Each genus is a tuple (node #, full_traversal, domain, genus name)
    genera_list = []
    with open(input_path, "r") as in_fd:
        for line in in_fd:
            line_split = line.split()
            last_pos_of_traversal, traversal = get_traversal_from_full_line(line)
    
            # If this a genus line, store the information
            if line_split[last_pos_of_traversal+2] == "genus":
                domain = get_domain_from_traversal(traversal)
                genus = get_genus_from_traversal(traversal)
                node_id = int(line_split[last_pos_of_traversal+1])
                genera_list.append((node_id, traversal, domain, genus))
        return genera_list 

def generate_node_id_to_traversal_dict(traversal_to_level_path):
    """ create a dictionary mapping node id to taxonomy traversal """
    node_map = {}
    with open(traversal_to_level_path, "r") as in_fd:
        for line in in_fd:
            line_split = line.split()
            last_pos_of_traversal, traversal = get_traversal_from_full_line(line)
            
            node_id = int(line_split[last_pos_of_traversal+1])
            tax_rank = line_split[last_pos_of_traversal+2]
            
            node_map[node_id] = (traversal, tax_rank)
    return node_map

def extract_traversal_only_up_to_genus(full_traversal):
    line_split = full_traversal.split(";")
    new_traversal = ";".join(line_split[0:len(line_split)-1]) + ";"
    return new_traversal

def get_ref_seq_headers_in_silva(input_path):
    seq_head_list = []
    with open(input_path, "r") as in_fd:
        for line in in_fd:
            if ">" in line:
                new_trav = extract_traversal_only_up_to_genus(" ".join(line.split()[1:]))
                seq_head_list.append(new_trav)
    return seq_head_list

def traversal_newick_tree(curr_node, tree_traversal):
    """ performs a pre-order traversal of tree """
    tree_traversal.append((curr_node.traversal, curr_node.rank))
    if len(curr_node.children) == 0:
        return tree_traversal
    else:
        for child in curr_node.children:
            tree_traversal = traversal_newick_tree(child, tree_traversal)
    return tree_traversal

def generate_genus_to_raw_sequences_dict(genera_order, input_file_path):
    """ create dictionary going from genus name to list of all raw sequences for it """
    genus_to_raw_seq_list = {}
    for curr_genera_tup in genera_order:
        genus_to_raw_seq_list[curr_genera_tup[0]] = []
    
    def add_raw_seq_to_dict(trav, header, seq, seq_dict):
        if len(trav) > 0 and len(header) > 0 and len(seq) > 0:
            if trav in seq_dict:
                output_str = f"{header}\n{seq}\n"
                seq_dict[trav].append(output_str)
        return seq_dict

    with open(input_file_path, "r") as silva_db:
        curr_trav = ""; curr_header = ""; curr_seq = ""
        for line in silva_db:
            if ">" in line:
                genus_to_raw_seq_list = add_raw_seq_to_dict(curr_trav, 
                                                            curr_header, 
                                                            curr_seq,
                                                            genus_to_raw_seq_list)
                curr_trav = ""; curr_header = ""; curr_seq = ""

                curr_trav = extract_traversal_only_up_to_genus(" ".join(line.split()[1:]))
                curr_header = line.strip()
            else:
                curr_seq += line.strip()
        genus_to_raw_seq_list = add_raw_seq_to_dict(curr_trav, 
                                                    curr_header, 
                                                    curr_seq,
                                                    genus_to_raw_seq_list)
    return genus_to_raw_seq_list

def write_individual_genus_to_file(input_tup):
    genus_num, traversal, genus_to_seq_dict, output_dir = input_tup
    with open(output_dir+f"doc_{genus_num+1}_seq.fa", "w") as out_fd:
        assert traversal in genus_to_seq_dict
        if len(genus_to_seq_dict[traversal]) == 0:
            return 0
        for seq in genus_to_seq_dict[traversal]:
            out_fd.write(seq)
    return 1

def create_filelist(output_dir, num_genera_written):
    with open(output_dir+"filelist.txt", "w") as out_fd:
        for i in range(1, num_genera_written+1):
            out_fd.write(f"{output_dir}doc_{i}_seq.fa {i}\n")

########################################
# Section 2: Main method, general helper
########################################

def main(args):
    # Step 1: Load in SILVA tree structure, determine a list of the genera
    #         ordered based on the tree structure.
    newick_tree = get_silva_newick_tree(args.tree_path)
    node_map = generate_node_id_to_traversal_dict(args.tree_rank_path)
    root = build_newick_tree(TreeNode("Root"), strip_parentheses(newick_tree), node_map)
    print("\n[log] build the SILVA taxonomy")

    # Step 2: Get the pre-order traversal, take only genus leaf nodes
    traversal = traversal_newick_tree(root, [])
    genera_in_order = [tup for tup in traversal if tup[1] == "genus"]
    print(f"[log] computed the list of {len(genera_in_order)} genera in tree-order")

    # Step 3: Count number of sequences per genus
    genus_to_seq_count = {}
    for curr_genus in genera_in_order:
        genus_to_seq_count[curr_genus[0]] = 0
    
    seq_headers = get_ref_seq_headers_in_silva(args.input_file)
    for curr_seq_header in seq_headers:
        if curr_seq_header in genus_to_seq_count:
            genus_to_seq_count[curr_seq_header] += 1
    no_seq_count = sum([1 for key in genus_to_seq_count.keys() if genus_to_seq_count[key] == 0])
    print(f"[log] counted number of sequence per genus. Found {no_seq_count} with no sequences.")

    # Step 4: Generate a dictionary mapping each genus to list of sequences
    genus_to_sequences = generate_genus_to_raw_sequences_dict(genera_in_order, args.input_file)
    total_number_of_seqs_at_genera_level = sum([len(genus_to_sequences[key]) for key in genus_to_sequences.keys()])
    print(f"[log] out of a total of {len(seq_headers)}, {total_number_of_seqs_at_genera_level} are specified to genus level.")

    # Step 5: Write out n genera to files
    genera_to_write = args.num_genera if args.num_genera < len(genera_in_order) else len(genera_in_order)
    genera_write_list = [(i, genera_in_order[i][0], genus_to_sequences, args.output_dir) for i in range(genera_to_write)]
    
    with mp.Pool(8) as pool:
        res = pool.map(write_individual_genus_to_file, genera_write_list)
    print(f"[log] finished writing {sum(res)} genera to files.")

    create_filelist(args.output_dir, sum(res))
    print(f"[log] finished writing a filelist.")
    

def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script allows users to separate a SILVA database at genus-level.")
    parser.add_argument("-i", "--input", dest="input_file", required=True, help="input SILVA database")
    parser.add_argument("-o", "--output-dir", dest="output_dir", required=True, help="output directory for files")
    parser.add_argument("--tree", dest="tree_path", default="", required=True, help="path to *.tre file for the SILVA database")
    parser.add_argument("--tree-map", dest="tree_map_path", default="", required=True, help="path to *.map file for the SILVA database")
    parser.add_argument("--tree-rank", dest="tree_rank_path", default="", required=True, help="path to *.txt file for the SILVA database")
    parser.add_argument("-n", "--num", dest="num_genera", default=10, type=int, help="number of genera to write out to separate files")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Validate the command-line arguments """
    if not os.path.isfile(args.input_file):
        print("Error: the path for the provided input file is not valid.")
        exit(1)
    if not os.path.isdir(args.output_dir):
        print("error: the output directory path is not valid.")
        exit(1)
    elif args.output_dir[-1] != '/':
        args.output_dir += "/"
    
    if not os.path.isfile(args.tree_path):
        print("Error: path to *.tre file is not valid")
        exit(1)
    if not os.path.isfile(args.tree_map_path):
        print("Error: path to *.map file is not valid")
        exit(1)
    if not os.path.isfile(args.tree_rank_path):
        print("Error: path to *.txt file is not valid")
        exit(1)
    
    if ".tre" not in args.tree_path:
        print("Error: the tree file does not have valid extension (*.tre)")
        exit(1)
    if ".map" not in args.tree_map_path:
        print("Error: the tree map file does not have valid extension (*.map)")
        exit(1)
    if ".txt" not in args.tree_rank_path:
        print("Error: the tree map file does not have valid extension (*.map)")
        exit(1)

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)